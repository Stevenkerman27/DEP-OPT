import openvsp as vsp
import os
import numpy as np
import pandas as pd
import infrastructure as opb
import matplotlib.pyplot as plt
cfg = {"wing_S": 0, "bref": 0,"cref": 0}

wing_pos = {"name": "Mainwing", "x":0, "y":0, "z":0, "yr": 1}

airfoil_cfg = {"filename": None, "Camber" : 0.04,"CamberLoc": 0.4,"ThickChord": 0.12}

flap_cfg = {"name": "Flaperon", "c" :1, "Length_Start" : 0.3,"Length_End" : 0.25,"length" : 0.8}

def_cfg = {"Flaperon": 30}

weight = {"S_density": 1.75, "CF_Strength": 400e6, "CF_rho": 2000, "fuse_weight": 1.2, "prop-DW": 0.1, "payload": 0.5} # weight in kg, strength in Pa

#setting
include_weight = 1
include_TL = 1
reuse_pre = 1
mass = 1.2

root = None
labels = {}  # 键为变量名，值为对应的 Label 控件

# 优化历史
ld_hst = []
LP_hst =[]
power_hst = []

# 参数
case_name = "convergence"
file_name = case_name + ".vsp3"

opb.case_name = case_name
opb.file_name = file_name

cruise_spd = 15

G = 6
SF = 2
opb.SF = SF
opb.G = G

tess_interval = 0.016
prop_Dia = 9

tess_choice = [0.04, 0.03, 0.02, 0.013, 0.01, 0.0083, 0.0075, 0.0067, 0.0062]
conve_hst = []

# 获取当前脚本所在目录
script_dir = os.path.dirname(os.path.abspath(__file__))
# 创建 outputs 子文件夹（如果不存在）
output_dir = os.path.join(script_dir, "outputs")
os.chdir(output_dir)


def converge(tess):
    opb.ini_geom()
    global iter
    global mass
    opb.mass = mass
    # span, chord and taper
    span = 0.6
    Mean_chord = 0.16
    CG = Mean_chord / 3 
    taper = 0.4
    #root and tip
    root = 2 * Mean_chord / (1 + taper)
    tip = root * taper
    #twist
    twist = 0
    #prop ele
    prop_ele = 0.01
    #wing angle
    wing_angle = 0
    wing_pos["yr"] = wing_angle
    #定义机翼
    spanlist, chordlist, twistlist, wing_S, Cl_target = opb.think_trapwing(root, tip, span, 0.03, twist, cruise_spd)
    #定义螺旋桨
    prop_D= opb.place_single_prop(prop_Dia, 0.2, tess, prop_ele)
    #创建主翼
    Npanel = opb.create_wing(wing_pos, spanlist, chordlist, twistlist, 0, tess, airfoil_cfg)
    vsp.Update()
    vsp.WriteVSPFile(file_name)

    # 气动分析
    cfg["wing_S"] = wing_S
    cfg["bref"] = 2 * span
    cfg["cref"] = Mean_chord
    ele_def0 = 0
    def_cfg["Flaperon"] = 0
    def_cfg["elevator"] = ele_def0

    RPM = [4000]
    Ct = [0.02]
    Cp = [0.02]
    drag, alpha, lift, netdrag1, power, _, CMy1 = opb.runaero(CG, 2, 2, 1, cruise_spd, cfg, Cl_target, opb.solver_config1,def_cfg, prop_D, RPM, Ct, Cp)
    return(Npanel, drag, lift, CMy1)
    

def main():
    for i in range(0,  len(tess_choice)):
        tess = tess_choice[i]
        Npanel, drag, lift, CMy1 = converge(tess)
        conve_hst.append([Npanel, drag, lift, CMy1])

main()
vsp.WriteVSPFile(file_name)
print(f"模型已保存")

res = np.asarray(conve_hst, dtype=float)  # shape = [N, 4] -> [Npanel, drag, lift, CMy1]
df = pd.DataFrame(res, columns=["Npanel", "drag", "lift", "CMy1"])

# === 保存 CSV ===
csv_path = os.path.join(output_dir, "convergence_results.csv")
df.to_csv(csv_path, index=False, float_format="%.8g")
print(f"[Convergence] 结果已写入: {csv_path}")

# ==== 方式1：四张单独图 ====
for col in ["drag", "lift", "CMy1"]:
    plt.figure(figsize=(5, 4))
    plt.plot(df["Npanel"].to_numpy(), df[col].to_numpy(), marker="o")
    plt.xlabel("Npanel")
    plt.ylabel(col)
    plt.title(f"{col} vs Npanel")
    plt.grid(True, linestyle="--", alpha=0.4)

plt.show()

