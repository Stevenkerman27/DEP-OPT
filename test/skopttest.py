import numpy as np
import matplotlib.pyplot as plt

from skopt import gp_minimize
from skopt.plots import (
    plot_convergence,
    plot_objective_2D,
    plot_evaluations,
    plot_regret
)
from skopt.space import Real

# 定义典型的二维多峰函数：Himmelblau 函数
def himmelblau(x):
    x1, x2 = x
    return (x1**2 + x2 - 11)**2 + (x1 + x2**2 - 7)**2

# 搜索空间：x1, x2 都在 [-5, 5] 区间
bounds = [Real(-5.0, 5.0, name='x1'),
          Real(-5.0, 5.0, name='x2')]

# 调用 gp_minimize
res = gp_minimize(
    func=himmelblau,
    initial_point_generator="lhs",
    dimensions=bounds,
    acq_func="LCB",          # 选择 Expected Improvement
    n_calls=40,             # 总评估次数
    random_state=0,
)

# 输出最优结果
print(f"最佳目标值: {res.fun:.6f}")
print(f"对应位置: x1 = {res.x[0]:.6f}, x2 = {res.x[1]:.6f}")

# —— 绘图 —— #
# 创建一个 1x3 的 Figure 和 Axes
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# 子图 1：收敛曲线
plot_convergence(res, ax=axes[0])
axes[0].set_title("Convergence Plot")

# 子图 2：二维等高面和采样点
plot_objective_2D(
    result=res,
    dimension_identifier1='x1',  # 或者 0
    dimension_identifier2='x2',  # 或者 1
    levels=30,
    ax=axes[1]
)
axes[1].set_title("2D Objective Surface")

# 子图 3：采样分布
plot_evaluations(res, ax=axes[2])
axes[2].set_title("Sampling Locations")

# 调整布局并一起显示
plt.tight_layout()
plt.show()
