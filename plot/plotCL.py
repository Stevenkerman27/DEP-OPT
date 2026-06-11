import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
import io
import os
from scipy.signal import savgol_filter   # ✨ 用于平滑曲线

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# 设置全局字体为 Times New Roman
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams["mathtext.fontset"] = "stix"

# ---------------------- 处理 VSPAero 数据 ----------------------
def process_wing_data(file_path):
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        return None, None

    # --- 提取 Wing 段数据 ---
    wing_header_index = -1
    for i, line in enumerate(lines):
        if line.strip().startswith("Wing"):
            wing_header_index = i
            break
    if wing_header_index == -1:
        return None, None

    wing_data_lines = []
    for i in range(wing_header_index + 1, len(lines)):
        line_parts = lines[i].strip().split()
        if not line_parts or not line_parts[0].isdigit():
            break
        wing_data_lines.append(lines[i])

    wing_data_string = lines[wing_header_index] + "".join(wing_data_lines)
    df = pd.read_csv(io.StringIO(wing_data_string), sep=r"\s+")

    # --- 仅保留 Wing 1 数据 ---
    df["Wing"] = pd.to_numeric(df["Wing"], errors="coerce")
    wing1_data = df[df["Wing"] == 1].copy()

    required_cols = ["Cl", "Yavg", "V/Vref", "Chord"]
    if not all(col in wing1_data.columns for col in required_cols):
        print(f"错误: 文件 {file_path} 缺少以下列: "
              f"{[c for c in required_cols if c not in wing1_data.columns]}")
        return None, None

    for col in required_cols:
        wing1_data[col] = pd.to_numeric(wing1_data[col], errors="coerce")
    wing1_data.dropna(subset=required_cols, inplace=True)
    if wing1_data.empty:
        return None, None

    # 按展向排序
    wing1_data = wing1_data.sort_values("Yavg").reset_index(drop=True)

    # --- 提取整体 CL ---
    comp_header_index = -1
    CL_overall = None
    for i, line in enumerate(lines):
        if line.strip().startswith("Comp"):
            comp_header_index = i
            break

    if comp_header_index != -1:
        comp_string = "".join(lines[comp_header_index:])
        comp_df = pd.read_csv(io.StringIO(comp_string), sep=r"\s+")
        if "Comp" in comp_df.columns and "CL" in comp_df.columns:
            comp_df["Comp"] = pd.to_numeric(comp_df["Comp"], errors="coerce")
            target_row = comp_df[comp_df["Comp"] == 1]
            if not target_row.empty:
                CL_overall = target_row["CL"].iloc[0]

    if CL_overall is None:
        print(f"警告: 文件 {file_path} 未找到 Overall CL，使用 1.0 替代。")
        CL_overall = 1.0

    # --- 计算展向升力分布 ---
    wing1_data["y_value_final"] = (
        #wing1_data["Cl"] * wing1_data["Chord"] * CL_overall * (wing1_data["V/Vref"] ** 2)
        wing1_data["Cl"] * wing1_data["Chord"]
    )
    print("avg speed: " + str(np.sum(wing1_data["V/Vref"]* wing1_data["Chord"]) / np.sum(wing1_data["Chord"])))

    # --- 归一化 X ---
    max_yavg = wing1_data["Yavg"].max()
    wing1_data["x_value_final"] = (
        wing1_data["Yavg"] / max_yavg if max_yavg > 0 else 0
    )

    return wing1_data[["x_value_final", "y_value_final"]], CL_overall


# ---------------------- 读取 CFD 数据 ----------------------
def read_cfd_load(file_path,target_area=1.0):
    """
    读取 CFD 载荷分布文件 load.csv。
    文件前两列分别是 位置(y) 和 力(lift)。
    """
    try:
        df = pd.read_csv(file_path)
    except FileNotFoundError:
        print(f"未找到 CFD 文件: {file_path}")
        return None

    if df.shape[1] < 2:
        print(f"CFD 文件 {file_path} 格式错误（至少需要两列）。")
        return None

    # 假定第1列是位置 y，第2列是升力
    y = df.iloc[:, 0].values
    L = df.iloc[:, 1].values

    # --- 归一化 ---
    y_norm = (y - y.min()) / (y.max() - y.min())

    current_area = np.trapz(L, y_norm)
    L_scaled = L * (target_area / current_area)

    return pd.DataFrame({"x_value_final": y_norm, "y_value_final": L_scaled})


# ---------------------- 主绘图函数 ----------------------
def plot_multiple_spanwise_distributions(files_to_plot,
                                         include_cfd=False,       #  是否绘制CFD曲线
                                         smooth_vlm=True,        #  是否平滑VLM曲线
                                         smooth_window=7,        #  平滑窗口长度（奇数）
                                         smooth_poly=3):  #  Savitzky-Golay滤波多项式阶数        
    """
    绘制多个归一化升力分布，可选CFD曲线和平滑处理。
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown']

    vlm_area_for_cfd = 1.0  # 默认面积
    first_vlm_processed = False

    print("-" * 30)

    # --- 绘制 VSPAero (VLM) 结果 ---
    for i, (filename, label) in enumerate(files_to_plot.items()):
        try:
            wing_data, cl_overall = process_wing_data(filename)
            if cl_overall is not None:
                print(f"文件 '{filename}': 读取到 Overall CL = {cl_overall:.4f}")

            if wing_data is not None and not wing_data.empty:
                y_val = wing_data["y_value_final"].copy()
                x_val = wing_data["x_value_final"].copy()

                # 平滑处理
                if smooth_vlm and len(y_val) >= smooth_window:
                    y_val = savgol_filter(y_val, smooth_window, smooth_poly)

                # 计算该曲线的面积 (integral of c*Cl over 2y/b)
                current_vlm_area = np.trapz(y_val, x_val)
                
                # 记录第一个 VLM 文件的面积作为 CFD 的基准
                if not first_vlm_processed:
                    vlm_area_for_cfd = current_vlm_area
                    first_vlm_processed = True
                    print(f"CFD 将参考 '{label}' 的面积进行归一化: {vlm_area_for_cfd:.4f}")

                ax.plot(
                    x_val,
                    y_val,
                    linestyle='-',
                    color=colors[i % len(colors)],
                    linewidth=3,
                    #label=f"{label} (VSPAero)"
                    label=f"{label} "
                )
            else:
                print(f"错误: 无法处理或绘制文件 '{filename}' 的数据。")
        except Exception as e:
            print(f"处理文件 '{filename}' 时发生未知错误: {e}")

    # --- 可选：绘制 CFD 结果 ---
    if include_cfd:
        cfd_file = os.path.join(script_dir, "load.csv")
        cfd_data = read_cfd_load(cfd_file, target_area=vlm_area_for_cfd)
        if cfd_data is not None:
            ax.plot(
                cfd_data["x_value_final"],
                cfd_data["y_value_final"],
                linestyle='--',
                color='black',
                linewidth=3,
                label="RANS CFD method"
            )
            print(f"CFD 文件 '{cfd_file}' 已成功读取并绘制。")
        else:
            print("未找到或无法读取 CFD 数据。")
    else:
        print("已关闭 CFD 曲线绘制。")

    print("-" * 30)

    # --- 图表标签 ---
    ax.set_xlabel("2y/b [-]", fontsize=30)
    #ax.set_ylabel("Lift [N]", fontsize=22)
    ax.set_ylabel(r"$c c_l$ [m]", fontsize=30)
    ax.grid(True, linestyle='--', linewidth=0.5, color='gray', alpha=0.6)
    ax.tick_params(axis='both', which='major', labelsize=28)
    ax.legend(fontsize=30,loc='lower left')
    ax.set_xlim(0, 1)
    plt.tight_layout()
    plt.show()


# ---------------------- 主程序入口 ----------------------
if __name__ == "__main__":
    files = {
        "optL.lod": "Optimized",
        #"optC.lod": "VLM method",
        "baseL.lod": "Baseline"
    }

    # 调整下面两个布尔开关即可：
    plot_multiple_spanwise_distributions(
        files_to_plot=files,
        include_cfd=False,    # ← 设置为 False 可关闭 CFD 曲线
        smooth_vlm=True,     # ← 设置为 False 不平滑 VLM 结果
        smooth_window=7,     # ← 平滑窗口，建议 5~11
        smooth_poly=3        # ← 拟合阶数，通常 2~3 即可
    )
