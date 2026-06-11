import openvsp as vsp
import os
import numpy as np
import pandas as pd
import infrastructure as opb
import prop
import matplotlib.pyplot as plt

# Global Config (Mirroring opt)
wing_pos = {"name": "Mainwing", "x":0, "y":0, "z":0, "yr": 2} # 2 deg for opt
tail_pos = {"name": "HT", "x":0.615, "y":0, "z":-0.03, "yr": 0}
airfoil_cfg = {"filename": None, "Camber" : 0.04,"CamberLoc": 0.4,"ThickChord": 0.12}
airfoiltail_cfg = {"filename": None, "Camber" : 0, "CamberLoc": 0, "ThickChord": 0.12}
tail_cfg = {"root": 0.168,"tip": 0.095,"span": 0.24}

# Params
cruise_spd = 15
fuse_w = 0.07
liftprop_Dia = 8
tipprop_Dia = 13
prop_ele = 0.01
CG = 0.072748
G = 5 # Match evaluate.py
SF = 1.5
opb.SF = SF
opb.G = G
opb.g = 9.8
opb.density = 1.225
mass_prop = {"S_density": 1.6, "CF_Strength": 450e6, "CF_rho": 2000, "fuse_mass": 1.2, "prop": [0.05,0.05,0.05,0.12], "payload": 0.3}
opb.mass_prop = mass_prop

# settings (unused or legacy removed if redundant)
case_name = "convergence"
file_name = case_name + ".vsp3"

opb.case_name = case_name
opb.file_name = file_name

tess_choice = [0.04, 0.03, 0.02, 0.013, 0.01, 0.0083, 0.0075, 0.0067]
conve_hst = []

# 获取当前脚本所在目录
script_dir = os.path.dirname(os.path.abspath(__file__))
# 创建 outputs 子文件夹（如果不存在）
output_dir = os.path.join(script_dir, "outputs")
os.makedirs(output_dir, exist_ok=True)
os.chdir(output_dir)


def converge(tess):
    vsp.ClearVSPModel()
    opb.ini_geom()
    
    # Aircraft Parameters (opt)
    Mean_chord = 0.197
    taper = 0.59
    span = 1.0
    root = 2 * Mean_chord / (1 + taper)
    tip = root * taper
    twist = 0

    # 1. Define Geometry
    spanlist, chordlist, twistlist, wing_S, Cl_target = opb.think_trapwing(root, tip, span, fuse_w, twist, cruise_spd)
    
    # 2. Place Propellers (4 props)
    prop_pos, prop_D_inch, Nprops = opb.place_prop(span, liftprop_Dia, tipprop_Dia, prop_ele, tess, fuse_w)
    prop_D = np.array(prop_D_inch) * opb.inch_in_m
    
    # 3. Create Wing and HT
    Npanel_wing = opb.create_wing(wing_pos, spanlist, chordlist, twistlist, 0, tess, airfoil_cfg)
    
    tail_span = [tail_cfg["span"]]
    tail_chord = [tail_cfg["root"], tail_cfg["tip"]]
    tail_twist = [0]
    Npanel_tail = opb.create_wing(tail_pos, tail_span, tail_chord, tail_twist, opb.max_sweeploc, tess, airfoiltail_cfg)
    
    Npanel = Npanel_wing + Npanel_tail
    vsp.Update()
    vsp.WriteVSPFile(file_name)

    # 4. Mass and Aero Setup
    # Use mass_sim_iter to get a consistent mass for target Cl
    y_cp = span * 0.4 # approximation for initial mass
    mass, _ = opb.mass_sim_iter(span, wing_S, y_cp, prop_pos, spanlist, chordlist)
    opb.mass = mass
    
    cfg_aero = {"wing_S": wing_S, "bref": 2 * span, "cref": Mean_chord}
    Cl_req = mass * 9.8 / (0.5 * 1.225 * cruise_spd**2 * wing_S)
    
    # Calculate Prop Parameters (fixed point for convergence)
    RPM = [4300,4400,4400,4400]
    Ct = [0.015,0.015,0.015,0.015]
    Cp = [0.01,0.01,0.01,0.01] 
    
    # 5. Run Aero
    drag, alpha, lift, netdrag, power, _, CMy = opb.runaero(CG, 2, 2, 1, cruise_spd, cfg_aero, Cl_req, opb.solver_config1, {}, prop_D, RPM, Ct, Cp)
    
    return(Npanel, drag, lift, CMy)
    

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
    plt.savefig(os.path.join(output_dir, f"convergence_{col}.png"))
    print(f"Saved plot: convergence_{col}.png")

# plt.show()

