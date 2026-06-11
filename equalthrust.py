import openvsp as vsp
import os
import numpy as np
import threading
import tkinter as tk
import prop
import infrastructure as opb
import csv
cfg = {
    "wing_S": 0, 
    "bref":   0,
    "cref":   0
}

airfoil_cfg = {
    "Camber" : 0.04,
    "CamberLoc": 0.4,
    "ThickChord": 0.12
}

sub_cfg = {
    "Length_Start" : 0.06,
    "Length_End" : 0.06,
    "length" : 0.8
}

# 全局变量，用于 objfunc 中更新 GUI
root = None
labels = {}  # 键为变量名，值为对应的 Label 控件

# 优化历史
ld_hst = []
LP_hst =[]
thrust_hst=[]
landing_thrust_hst=[]
power_hst = []
RPMC_hst = []
RPML_hst = []
liftL_hst = []
propC_hst = []
propL_hst = []

# 参数
mass = 36
opb.mass = mass
cruise_spd = 15
min_speed = 7
max_AOA = 10
flap_angle = 30
density = 1.225
opb.density = density
tess_interval = 0.006
iter = 0
sen_step = 0.02
drag_maxit = 5

Nprops = 4
mode = 0 # 0 for VLM, 1 for panel
span = 1.015
prop_choice = prop.prop_choice
liftprop_Dia = 8
tipprop_Dia = 13

# 获取当前脚本所在目录
script_dir = os.path.dirname(os.path.abspath(__file__))
# 创建 outputs 子文件夹（如果不存在）
output_dir = os.path.join(script_dir, "outputs")
os.makedirs(output_dir, exist_ok=True)
os.chdir(script_dir)
prop_data = prop.read_apce_grouped('data')
# 切换当前工作目录
os.chdir(output_dir)

def update_gui(current_dvs):
    for name, val in current_dvs.items():
        labels[name].config(text=f"{name}: {val}")
    # 立即刷新界面
    root.update_idletasks()

def objfunc():
    global iter
    # chord and taper
    Mean_chord = 0.24
    taper = 1
    root = Mean_chord * 2  / (1 + taper )
    tip = Mean_chord * 2 - root
    #twist
    twist = 0
    #prop ele
    prop_ele = 0
    # RPS ratio
    thrust_ratio = 0.8
    thrust_ratio_landing = 0.4
    #定义机翼
    spanlist, chordlist, twistlist, wing_S, Cl_target = opb.think_trapwing(root, tip, span, twist, cruise_spd)
    #定义螺旋桨
    prop_pos, prop_D_inch = opb.place_prop(span, Nprops, liftprop_Dia, tipprop_Dia)
    prop_D = np.array(prop_D_inch) * opb.inch_in_m
    opb.create_wing(spanlist, chordlist, twistlist, prop_pos, prop_D, prop_ele, tess_interval, airfoil_cfg, sub_cfg)
    # 气动分析
    cfg["wing_S"] = wing_S
    cfg["bref"] = 2 * span
    cfg["cref"] = Mean_chord
    drag, alpha, lift, netdrag0, power, _ = opb.runaero(mode, 0, max_AOA, 2, cruise_spd, cfg, Cl_target)
    RPM, Ct, Cp = prop.equal_thrust(prop_data, netdrag0, cruise_spd, prop_D_inch, thrust_ratio) #大螺旋桨转速
    #修正升力
    total_lift = lift + netdrag0 * np.sin(alpha/360 * 2 * np.pi)
    Cl_target = mass / total_lift * Cl_target
    #第二次气动
    drag, alpha, lift, netdrag1, power, _ = opb.runaero(mode, max(0, alpha-1), min(10,alpha+1), 2, cruise_spd, cfg, Cl_target, 0, prop_D, RPM, Ct, Cp)
    thrust1 = netdrag0
    thrust0 = 0 
    drag_tol = drag * 0.01
    for _ in range(drag_maxit):
        if (abs(netdrag1) < drag_tol):
            print("Converged!")
            break
        thrust2 = (netdrag0 * thrust1 - netdrag1 * thrust0)/(netdrag0-netdrag1)
        RPM, Ct, Cp = prop.equal_thrust(prop_data, thrust2, cruise_spd, prop_D_inch, thrust_ratio) #大螺旋桨转速
        #修正升力
        total_lift = lift + thrust2 * np.sin(alpha/360 * 2 * np.pi)
        Cl_target = mass / total_lift * Cl_target
        #气动
        drag, alpha, lift, netdrag2, power, _ = opb.runaero(mode, max(0, alpha-1), min(10,alpha+1), 2, cruise_spd, cfg, Cl_target, 0, prop_D, RPM, Ct, Cp)
        # 更新下次迭代点
        netdrag0 = netdrag1
        netdrag1 = netdrag2
        thrust0 = thrust1
        thrust1 = thrust2

    drag_res = netdrag1 #巡航阻力残差
    thrust_hst.append(thrust1)
    
    LD = mass / drag
    eff = mass / power
    power_hst.append(power)
    AR = (span*2)**2/wing_S
    propC_hst.append([Ct, Cp])

    # 升力约束
    drag, _, lift_L, netdrag0, power_L, _ = opb.runaero(mode, max_AOA, max_AOA, 1, min_speed, cfg, _, flap_angle)
    RPM_L, Ct, Cp = prop.equal_thrust(prop_data, netdrag0, min_speed, prop_D_inch, thrust_ratio_landing) #大螺旋桨转速
    #第二次气动
    drag, _, lift_L, netdrag1, power_L, Cl_target_landing = opb.runaero(mode, max_AOA, max_AOA, 1, min_speed, cfg, _, flap_angle, prop_D, RPM_L, Ct, Cp)
    total_lift = lift_L + netdrag0 * np.sin(max_AOA/360 * 2 * np.pi)
    #RPM循环
    thrust1 = netdrag0
    thrust0 = 0 
    drag_tol = drag * 0.02
    for _ in range(drag_maxit):
        if (abs(netdrag1) < drag_tol):
            break
        thrust2 = (netdrag0 * thrust1 - netdrag1 * thrust0)/(netdrag0-netdrag1)
        RPM_L, Ct, Cp = prop.equal_thrust(prop_data, thrust2, min_speed, prop_D_inch, thrust_ratio_landing) #大螺旋桨转速
        #气动
        drag, _, lift_L, netdrag2, power_, Cl_target_landing= opb.runaero(mode, max_AOA, max_AOA, 1, min_speed, cfg, _, flap_angle, prop_D, RPM_L, Ct, Cp)
        #修正升力
        total_lift = lift_L + thrust2 * np.sin(max_AOA/360 * 2 * np.pi)
        # 更新下次迭代点
        netdrag0 = netdrag1
        netdrag1 = netdrag2
        thrust0 = thrust1
        thrust1 = thrust2

    #添加历史
    propL_hst.append([Ct, Cp])
    ld_hst.append(LD)
    LP_hst.append(eff)
    landing_thrust_hst.append(thrust1)
    RPMC_hst.append(str(RPM[0])+"-"+str(RPM[-1]))
    RPML_hst.append(str(RPM_L[0])+"-"+str(RPM_L[-1]))
    liftL_hst.append(total_lift)

    # GUI
    update_gui({"root": root, "tip": tip, "LD":LD, "L/P": eff, "alpha":alpha, "minlift":total_lift,"thrust ratio":thrust_ratio, "thrust_ratio_landing": thrust_ratio_landing,
                "Cl_target": Cl_target, "RPM_C": RPMC_hst[-1], "RPM_L": RPML_hst[-1], "Cl_landing":  Cl_target_landing,
                "根梢比": taper, "平均弦长": Mean_chord, "展弦比": AR, "巡航阻力残差": drag_res , "降落阻力残差": netdrag1, "power": power })


def main():
    global root, labels
    root = tk.Tk()
    root.title("优化变量实时监控")
    lbl_x15 = tk.Label(root, text="--", font=("Arial", 14))#平均弦长
    lbl_x15.pack(padx=10, pady=5)
    lbl_x14 = tk.Label(root, text="--", font=("Arial", 14))#根梢比
    lbl_x14.pack(padx=10, pady=5)
    lbl_x0 = tk.Label(root, text="--", font=("Arial", 14)) #Root
    lbl_x0.pack(padx=10, pady=5)
    lbl_x1 = tk.Label(root, text="--", font=("Arial", 14)) #tip
    lbl_x1.pack(padx=10, pady=5)
    lbl_x2 = tk.Label(root, text="--", font=("Arial", 14))#LD
    lbl_x2.pack(padx=10, pady=5)
    lbl_x3 = tk.Label(root, text="--", font=("Arial", 14))#迎角
    lbl_x3.pack(padx=10, pady=5)
    lbl_x9 = tk.Label(root, text="--", font=("Arial", 14))#Cl
    lbl_x9.pack(padx=10, pady=5)
    lbl_x21 = tk.Label(root, text="--", font=("Arial", 14))#总效率
    lbl_x21.pack(padx=10, pady=5)
    lbl_x8 = tk.Label(root, text="--", font=("Arial", 14))#RPM_C
    lbl_x8.pack(padx=10, pady=5)
    lbl_x18 = tk.Label(root, text="--", font=("Arial", 14))#RPM_L
    lbl_x18.pack(padx=10, pady=5)
    lbl_x13 = tk.Label(root, text="--", font=("Arial", 14))#螺旋桨推力比
    lbl_x13.pack(padx=10, pady=5)
    lbl_x26 = tk.Label(root, text="--", font=("Arial", 14))#螺旋桨推力比降落
    lbl_x26.pack(padx=10, pady=5)
    lbl_x23 = tk.Label(root, text="--", font=("Arial", 14))#阻力残差
    lbl_x23.pack(padx=10, pady=5)
    lbl_x24 = tk.Label(root, text="--", font=("Arial", 14))#降落阻力残差
    lbl_x24.pack(padx=10, pady=5)
    lbl_x7 = tk.Label(root, text="--", font=("Arial", 14))#升力约束
    lbl_x7.pack(padx=10, pady=5)
    lbl_x16 = tk.Label(root, text="--", font=("Arial", 14))#展弦比
    lbl_x16.pack(padx=10, pady=5)
    lbl_x19 = tk.Label(root, text="--", font=("Arial", 14))#总效率
    lbl_x19.pack(padx=10, pady=5)
    lbl_x27 = tk.Label(root, text="--", font=("Arial", 14))#power
    lbl_x27.pack(padx=10, pady=5)
    labels["root"] = lbl_x0
    labels["tip"] = lbl_x1
    labels["LD"] = lbl_x2
    labels["alpha"] = lbl_x3
    labels["minlift"] = lbl_x7
    labels["RPM_C"] = lbl_x8
    labels["Cl_target"] = lbl_x9
    labels["thrust ratio"] = lbl_x13
    labels["根梢比"] = lbl_x14
    labels["平均弦长"] = lbl_x15
    labels["展弦比"] = lbl_x16
    labels["RPM_L"] = lbl_x18
    labels["L/P"] = lbl_x19
    labels["Cl_landing"] = lbl_x21
    labels["巡航阻力残差"] = lbl_x23
    labels["降落阻力残差"] = lbl_x24
    labels["thrust_ratio_landing"] = lbl_x26
    labels["power"] = lbl_x27
    # 启动后台线程跑优化，避免阻塞主线程的 GUI 事件循环
    opt_thread = threading.Thread(target=objfunc, daemon=True)
    opt_thread.start()
    # 启动 Tkinter 主循环
    root.mainloop()

main()
vsp.WriteVSPFile("equal.vsp3", vsp.SET_ALL)
print(f"模型已保存")
with open('history.csv', 'w', encoding='utf-8', newline='') as f:
    writer = csv.writer(f, delimiter=',')
    writer.writerow(['LD_hst'] + ld_hst)
    writer.writerow(['LP_hst'] + LP_hst)
    writer.writerow(['RPMC'] + RPMC_hst)
    writer.writerow(['RPML'] + RPML_hst)
    writer.writerow(['Thrust'] + thrust_hst)
    writer.writerow(['Thrust_L'] + landing_thrust_hst)
    writer.writerow(['Lift_L'] + liftL_hst)
    writer.writerow(['POWER'] + power_hst)
    writer.writerow(['propC'] + propC_hst)
    writer.writerow(['propL'] + propL_hst)