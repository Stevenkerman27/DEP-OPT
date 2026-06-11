import time
import threading
import tkinter as tk

# 确保使用 pyOptSparse 2.11.3
from pyoptsparse import Optimization, SLSQP


# 全局变量，用于 objfunc 中更新 GUI
root = None
labels = {}  # 键为变量名，值为对应的 Label 控件


def update_gui(current_dvs):
    """
    current_dvs: dict，比如 {"x0": 1.2345, "x1": -0.5678}
    把每个变量的当前数值更新到 Tkinter 窗口上的对应 Label 中。
    """
    for name, val in current_dvs.items():
        labels[name].config(text=f"{name}: {val:.6f}")
    # 立即刷新界面
    root.update_idletasks()


def objfunc(xdict):
    """
    pyOptSparse 要求的回调接口：返回 (funcs, fail)
    - xdict: {'x0': 当前 x0, 'x1': 当前 x1}，值都是标量。
    - funcs: 一个字典，包含目标和所有约束（本例无约束，只返回目标）
          funcs = {'f': f_val}
    - fail: 布尔，False 表示评估成功，True/1 表示评估失败
    """
    # 读取当前尝试的 x0, x1
    x0 = xdict["x0"]
    x1 = xdict["x1"]

    # 计算目标函数值
    f_val = (x0 - 3.0) ** 2 + (x1 + 1.0) ** 2

    # 在 GUI 上实时更新当前变量值，并加延时
    update_gui({"x0": x0, "x1": x1})
    time.sleep(0.3)  # 延时 0.3 秒，便于观察

    # 构造 funcs 字典（本例无约束）
    funcs = {"f": f_val}
    fail = False
    return funcs, fail


def run_optimization():
    """
    在后台线程中启动优化，指定 SLSQP 用有限差分 (sens="FD")，
    差分步长 sensStep=0.005。避免阻塞 Tkinter 主线程。
    """
    # 1. 构造优化问题
    opt_prob = Optimization("demo_problem", objfunc)

    # 添加两个连续型设计变量 x0、x1，范围 [-10, 10]，初始值都设为 0.0
    opt_prob.addVar("x0", "c", lower=-10.0, upper=10.0, value=0.0)
    opt_prob.addVar("x1", "c", lower=-10.0, upper=10.0, value=0.0)

    # 添加目标函数 f，对应 objfunc 返回的 funcs["f"]
    opt_prob.addObj("f")

    # 2. 选择 SLSQP 优化器
    optimizer = SLSQP()
    solution = optimizer(opt_prob,
                         sens="FD",
                         sensStep=0.005)

    # 优化结束后在控制台输出最终结果
    print(solution)
    # 此时 Tkinter 窗口依然保留，方便查看最后一次输出的 x0, x1


def main():
    global root, labels

    # 1. 创建 Tkinter 主窗口
    root = tk.Tk()
    root.title("优化变量实时监控")

    # 2. 在窗口里创建两个 Label，分别显示 x0、x1 的当前值
    lbl_x0 = tk.Label(root, text="x0: --", font=("Arial", 14))
    lbl_x0.pack(padx=10, pady=5)
    lbl_x1 = tk.Label(root, text="x1: --", font=("Arial", 14))
    lbl_x1.pack(padx=10, pady=5)

    labels["x0"] = lbl_x0
    labels["x1"] = lbl_x1

    # 3. 启动后台线程跑优化，避免阻塞主线程的 GUI 事件循环
    opt_thread = threading.Thread(target=run_optimization, daemon=True)
    opt_thread.start()

    # 4. 启动 Tkinter 主循环
    root.mainloop()

main()
