import matplotlib.pyplot as plt
import numpy as np

def plot_wing_planforms(
    initial_params,
    optimized_params,
    initial_color='cornflowerblue',
    optimized_color='lightcoral'
):
    """
    绘制并比较优化前后的机翼半翼平面图。

    参数:
    initial_params (dict): 初始机翼的参数字典，包含 'span', 'mean_chord', 'taper_ratio'。
    optimized_params (dict): 优化后机翼的参数字典，与 initial_params 结构相同。
    initial_color (str): 初始机翼的绘图颜色。
    optimized_color (str): 优化后机翼的绘图颜色。
    """

    def get_wing_coords(span, mean_chord, taper_ratio):
        """根据给定的空气动力学参数计算机翼四个角的坐标。"""
        half_span = span / 2.0

        # 根据梢根比(λ)和平均弦长(c_avg)计算翼根和翼尖弦长
        # c_avg = (c_root + c_tip) / 2
        # taper_ratio = c_tip / c_root
        # 解方程组可得：
        if (1 + taper_ratio) == 0: # 避免除以零
            raise ValueError("梢根比不能为 -1")
        
        c_root = (2 * mean_chord) / (1 + taper_ratio)
        c_tip = c_root * taper_ratio

        # 定义半翼平面的四个顶点坐标 (假设前缘无后掠)
        # 顺序: 翼根前缘 -> 翼尖前缘 -> 翼尖后缘 -> 翼根后缘
        x_coords = np.array([0, half_span, half_span, 0, 0])
        y_coords = np.array([0, 0, c_tip, c_root, 0])
        
        return x_coords, y_coords

    # --- 设置绘图样式以匹配参考图 ---
    #plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(10, 7))

    # --- 计算并绘制初始机翼 ---
    x_initial, y_initial = get_wing_coords(
        initial_params['span'],
        initial_params['mean_chord'],
        initial_params['taper_ratio']
    )
    ax.plot(x_initial, y_initial, color=initial_color, linewidth=2, label='Initial wing')

    # --- 计算并绘制优化后机翼 ---
    x_optimized, y_optimized = get_wing_coords(
        optimized_params['span'],
        optimized_params['mean_chord'],
        optimized_params['taper_ratio']
    )
    ax.plot(x_optimized, y_optimized, color=optimized_color, linewidth=2, label='Optimized wing')

    # --- 格式化图表 ---
    ax.set_xlabel('Span [m]', fontsize=14)
    ax.set_ylabel('Chord [m]', fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, which='both', linestyle='-', linewidth=0.5, color='gray', alpha=0.4)
    
    # 将Y轴反转，使机翼向下延伸，与参考图一致
    ax.invert_yaxis()
    
    # 设置坐标轴比例为1:1，确保几何形状不失真
    ax.set_aspect('equal', adjustable='box')
    
    # 调整坐标轴范围，使其更美观
    all_x = np.concatenate((x_initial, x_optimized))
    all_y = np.concatenate((y_initial, y_optimized))
    ax.set_xlim(0, all_x.max() * 1.1)
    ax.set_ylim(all_y.max() * 1.1, -0.1)


    plt.show()


# --- 主程序入口 ---
if __name__ == '__main__':
    # --- 在这里定义你的机翼参数 ---
    # 初始机翼 (长方形，梢根比为1.0)
    initial_wing = {
        'span': 1.8,         # 总翼展 (m)
        'mean_chord': 0.24,    # 平均弦长 (m)
        'taper_ratio': 1.0    # 梢根比 (c_tip / c_root)
    }

    # 优化后机翼 (直角梯形)
    # 根据参考图，优化后的翼根弦长变大，翼尖弦长变小，但翼展几乎不变
    optimized_wing = {
        'span': 1.86,         # 总翼展 (m)
        'mean_chord': 0.197,    # 平均弦长 (m)
        'taper_ratio': 0.59   # 梢根比 (c_tip / c_root)
    }

    # 调用函数绘图
    plot_wing_planforms(initial_wing, optimized_wing)