import openvsp as vsp
import os
from pyoptsparse import SLSQP, Optimization
import numpy as np
import threading
import tkinter as tk
import prop
import infrastructure as opb
import csv
import warnings
# This is for vsp3.41!!!
cfg = {"wing_S": 0, "bref": 0,"cref": 0}

wing_pos = {"name": "Mainwing", "x":0, "y":0, "z":0, "yr": 0}
tail_pos = {"name": "Horizontal_stab", "x":0.615, "y":0, "z":-0.03, "yr": 0}

airfoil_cfg = {"filename": None, "Camber" : 0.04,"CamberLoc": 0.4,"ThickChord": 0.12}

airfoiltail_cfg = {"filename": None, "Camber" : 0,"CamberLoc": 0,"ThickChord": 0.12}

flap_cfg = {"name": "Flaperon", "c" :1, "Length_Start" : 0.3,"Length_End" : 0.25,"EtaStart" : 0.8, "EtaEnd": 0.07}

ELE_cfg = {"name": "elevator", "c" :0, "Length_Start" : 0.05,"Length_End" : 0.05,"EtaStart" : 0.727000, "EtaEnd": 0}

def_cfg = {"elevator": -5, "Flaperon": -20}

tail_cfg = {"root": 0.168,"tip": 0.095,"span": 0.24}

mass_prop = {"S_density": 1.6, "CF_Strength": 450e6, "CF_rho": 2000, "fuse_mass": 1.2, "prop": [0.05,0.05,0.05,0.12], "payload": 0.3}
opb.mass_prop = mass_prop

#setting
include_weight = 1
include_TL = 1

mass = 1

root = None
labels = {}  # 键为变量名，值为对应的 Label 控件

# 优化历史
ld_hst = []
LP_hst =[]
ratio_hst = []

power_hst = []
target_hst = []
RPMC_hst = []
RPML_hst = []
liftL_hst = []
taper_hst =[]
chord_hst = []
span_hst = []
wingang_hst = []
ele_C = []
ele_L = []
ele_hst = []
alpha_hst = []

# 参数
case_name = "auto"
file_name = case_name + ".vsp3"

opb.case_name = case_name
opb.file_name = file_name
ACC = 0.1
g = 9.8
opb.g = g
cruise_spd = 15
min_speed = 7
max_AOA = 10

density = 1.225
opb.density = density

#flight setting
G = 5
SF = 1.5
n_ult = G * SF
opb.SF = SF
opb.G = G
g = 9.8
opb.g = g
cD0_C = 0.08
cD0_L = 0.1
cD0_S = 0.06
d0_C = opb.D0(cD0_C, cruise_spd, cD0_S) 
d0_L = opb.D0(cD0_L, min_speed, cD0_S) 

#analysis parameter
tess_interval = 0.01
iter = 0
sen_step = 0.02
drag_maxit = 9
far_1 = 3
far_2 = 10
wakeN_1 = 16
wakeN_2 = 32
omega_lift = 0.7
omega_drag = 0.7
max_step = 0.3  # deg

#Normalization
thrust_mul = 5
max_thr_ratio = 0.82
min_thr_ratio = 0.3
taper_mul = 5
chord_div = 3
span_mul = 2.5
angle_mul = 10

prop_choice = prop.prop_choice
liftprop_Dia = 8
tipprop_Dia = 13
fuse_w = 0.07

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

def reuse(value, seq, tar_seq):
    if not seq:
        return None
    idx = min(range(len(seq)), key=lambda i: abs(seq[i] - value))
    return tar_seq[idx]

def objfunc(x_dict):
    opb.ini_geom()
    global iter
    global mass
    # span, chord and taper
    span = float(x_dict["span"]) * span_mul
    Mean_chord = float(x_dict["Mean_chord"]) / chord_div
    CG = Mean_chord / 3 
    taper = float(x_dict["taper"])
    taper = taper * taper_mul
    #root and tip
    root = 2 * Mean_chord / (1 + taper)
    tip = root * taper
    #twist
    twist = 0
    #prop ele
    prop_ele = 0.01
    #wing angle
    wing_angle = float(x_dict["wing_angle"]) * angle_mul
    wing_pos["yr"] = wing_angle
    
    # RPS ratio
    thrust_ratio = float(x_dict["thrust_ratio"]) * thrust_mul
    thrust_ratio_landing = float(x_dict["thrust_ratio_landing"]) * thrust_mul
    # Record
    taper_hst.append(taper)
    span_hst.append(span)
    chord_hst.append(Mean_chord)
    ratio_hst.append(thrust_ratio)
    wingang_hst.append(wing_angle)
    # 定义问题
    funcs = {}
    #定义机翼
    spanlist, chordlist, twistlist, wing_S, Cl_target = opb.think_trapwing(root, tip, span, fuse_w, twist, cruise_spd)
    #定义螺旋桨
    prop_pos, prop_D_inch, Nprops = opb.place_prop(span, liftprop_Dia, tipprop_Dia, prop_ele, tess_interval, fuse_w)
    prop_D = np.array(prop_D_inch) * opb.inch_in_m
    #创建主翼
    Npanel = opb.create_wing(wing_pos, spanlist, chordlist, twistlist, 0, tess_interval, airfoil_cfg, flap_cfg)
    #创建尾翼
    tail_span = [tail_cfg["span"]]
    tail_chord = [tail_cfg["root"],tail_cfg["tip"]]
    tail_twist = [0]
    opb.create_wing(tail_pos, tail_span, tail_chord, tail_twist, opb.max_sweeploc, tess_interval, airfoiltail_cfg, ELE_cfg)
    vsp.Update()
    vsp.WriteVSPFile(file_name)
    #计算重量
    y_cp = (4/np.pi + (n_ult-1) * (1 + 2 * taper)/(1 + taper))/ (3 * n_ult)
    mass, wing_mass = opb.mass_sim_iter(span, wing_S, y_cp, prop_pos, spanlist, chordlist)
    opb.mass = mass

    flap_angle = def_cfg["Flaperon"]
    def_cfg["Flaperon"] = 0

    # GUI
    update_gui({"root": root, "tip": tip, "span": span, "根梢比": taper, "平均弦长": Mean_chord, "Nprops": Nprops, "wing_angle": wing_angle, "当前循环": iter, 
                "总质量": mass, "机翼质量": wing_mass, "thrust ratio":thrust_ratio, "thrust_ratio_landing": thrust_ratio_landing})

    flight_condition = {"speed": cruise_spd, "max_AOA":max_AOA, "Cl_target":Cl_target, "TR": thrust_ratio, "d0":d0_C }
    geo_info = {"spanlist":spanlist, "chordlist":chordlist, "span":span, "wing_S":wing_S, "bref":Mean_chord, "cref":2*span,
                "CG": CG, "def_cfg": def_cfg, "prop_D":prop_D,"prop_D_inch":prop_D_inch, "prop_pos":prop_pos}
    config={"max_it":drag_maxit, "tol":[0.01, 0.01, 0.01], "relax":[omega_drag, omega_lift], "max_alpha_step":max_step,"propdata":prop_data}

    lift, drag, power, alpha, RPM, thrust, mass_result, ele_def = opb.single_point(flight_condition, geo_info, config)

    mass = mass_result["mass"]
    wing_mass = mass_result["wing_mass"]
    drag_resC = thrust - drag #巡航阻力残差
    lift_res = lift - mass*g
    
    LD = lift / (drag + d0_C)
    eff = mass / power
    power_hst.append(power)
    AR = (span*2)**2/wing_S
    ele_C.append(ele_def)
    alpha_hst.append(alpha)
    ld_hst.append(LD)
    LP_hst.append(eff)

    if include_TL:
        # 升力约束
        max_L_AOA = max_AOA - wing_angle 
        def_cfg["Flaperon"] = flap_angle

        flight_condition = {"speed": min_speed, "max_AOA":max_L_AOA, "Cl_target":-1, "TR": thrust_ratio_landing, "d0":d0_L }
        geo_info["def_cfg"] = def_cfg

        lift_L, drag, power_L, alpha, RPM_L, thrust, mass_L, ele_def = opb.single_point(flight_condition, geo_info, config)

        drag_resL = thrust - drag #最小速度阻力残差

        # 返回目标值
        penalty = lift_L - mass * g #升力不足时penalty负
        if penalty < 0:
            penalty = penalty * - 2
        else: 
            penalty = penalty * 0.5
        target = power + penalty # 升力不足target更大
    else: #不考虑起降
        RPM_L = ["-"]
        target = power
    funcs["obj"] = target
    iter = iter + 1

    #添加历史
    ele_L.append(ele_def)
    ele_hst.append(str(ele_C[-1]) + " - " + str(ele_L[-1]))
    RPMC_hst.append(str(RPM[0])+"-"+str(RPM[-1]))
    RPML_hst.append(str(RPM_L[0])+"-"+str(RPM_L[-1])) 
    
    liftL_hst.append(lift_L)
    target_hst.append(target)
    
    # GUI
    update_gui({"LD":LD, "L/P": eff, "alpha":alpha_hst[-1], "minlift":lift_L, "Cl_target": Cl_target, "RPM_C": RPMC_hst[-1], "RPM_L": RPML_hst[-1],
                 "ele_angle": ele_hst[-1], "升阻比提升": ld_hst[-1] - ld_hst[0], "功率变化":power_hst[-1] - power_hst[0], 
                "阻力残差": str(drag_resC) + " - " + str(drag_resL), "功率-目标": str(power) + " - " + str(target)})
    return funcs, False

def run_optimization():
    # 初始化问题
    optProb = Optimization("Auto", objfunc)
    # 添加设计变量
    optProb.addVar("span", "c", lower = 0.75 / span_mul, upper = 0.93 / span_mul, value = 0.9 / span_mul)
    optProb.addVar("Mean_chord", "c", lower = 0.12 * chord_div, upper = 0.22 * chord_div, value = 0.2 * chord_div)
    optProb.addVar("taper", "c", lower=0.5 / taper_mul, upper = 0.9 /taper_mul, value = 0.6 / taper_mul) #tip/root
    optProb.addVar("wing_angle", "c", lower= 0 / angle_mul, upper = 2 / angle_mul, value = 1/  angle_mul)
    optProb.addVar("thrust_ratio", "c", lower = min_thr_ratio / thrust_mul, upper = max_thr_ratio / thrust_mul, value = 0.7 /thrust_mul)
    optProb.addVar("thrust_ratio_landing", "c", lower = min_thr_ratio / thrust_mul, upper = max_thr_ratio / thrust_mul, value = 0.4 /thrust_mul)
    # rst begin addObj
    optProb.addObj("obj")

    # Check optimization problem
    print(optProb)
    optProb.printSparsity()

    # 配置SLSQP参数（更大胆的收敛策略）
    optOptions = {"ACC": ACC, "MAXIT": 15, "IPRINT": 2}
    opt = SLSQP(options=optOptions)
    sol = opt(optProb, sens="FD", sensStep = sen_step, storeHistory="opt_hist.hst")
    print(sol)

def main():
    global root, labels
    root = tk.Tk()
    root.title("优化变量实时监控")
    lbl_x24 = tk.Label(root, text="--", font=("Arial", 13))#span
    lbl_x24.pack(padx=10, pady=5)
    lbl_x15 = tk.Label(root, text="--", font=("Arial", 13))#平均弦长
    lbl_x15.pack(padx=10, pady=5)
    lbl_x14 = tk.Label(root, text="--", font=("Arial", 13))#根梢比
    lbl_x14.pack(padx=10, pady=5)
    lbl_x0 = tk.Label(root, text="--", font=("Arial", 13)) #Root
    lbl_x0.pack(padx=10, pady=5)
    lbl_x1 = tk.Label(root, text="--", font=("Arial", 13)) #tip
    lbl_x1.pack(padx=10, pady=5)
    lbl_x4 = tk.Label(root, text="--", font=("Arial", 13)) #总质量
    lbl_x4.pack(padx=10, pady=5)
    lbl_x5 = tk.Label(root, text="--", font=("Arial", 13)) #机翼质量
    lbl_x5.pack(padx=10, pady=5)
    lbl_x27 = tk.Label(root, text="--", font=("Arial", 13))#Nprops
    lbl_x27.pack(padx=10, pady=5)
    lbl_x29 = tk.Label(root, text="--", font=("Arial", 13)) #安装角
    lbl_x29.pack(padx=10, pady=5)
    lbl_x13 = tk.Label(root, text="--", font=("Arial", 13))#螺旋桨推力比
    lbl_x13.pack(padx=10, pady=5)
    lbl_x26 = tk.Label(root, text="--", font=("Arial", 13))#螺旋桨推力比降落
    lbl_x26.pack(padx=10, pady=5)
    lbl_x6 = tk.Label(root, text="0", font=("Arial", 13))#迭代数
    lbl_x6.pack(padx=10, pady=5)
    lbl_x2 = tk.Label(root, text="--", font=("Arial", 13))#LD
    lbl_x2.pack(padx=10, pady=5)
    lbl_x3 = tk.Label(root, text="--", font=("Arial", 13))#迎角
    lbl_x3.pack(padx=10, pady=5)
    lbl_x9 = tk.Label(root, text="--", font=("Arial", 13))#Cl
    lbl_x9.pack(padx=10, pady=5)
    lbl_x8 = tk.Label(root, text="--", font=("Arial", 13))#RPM_C
    lbl_x8.pack(padx=10, pady=5)
    lbl_x18 = tk.Label(root, text="--", font=("Arial", 13))#RPM_L
    lbl_x18.pack(padx=10, pady=5)
    lbl_x23 = tk.Label(root, text="--", font=("Arial", 13))#阻力残差
    lbl_x23.pack(padx=10, pady=5)
    lbl_x7 = tk.Label(root, text="--", font=("Arial", 13))#升力约束
    lbl_x7.pack(padx=10, pady=5)
    lbl_x28 = tk.Label(root, text="--", font=("Arial", 13))#升降舵偏转
    lbl_x28.pack(padx=10, pady=5)
    lbl_x17 = tk.Label(root, text="--", font=("Arial", 13))#升阻比提升
    lbl_x17.pack(padx=10, pady=5)
    lbl_x19 = tk.Label(root, text="--", font=("Arial", 13))#总效率
    lbl_x19.pack(padx=10, pady=5)
    lbl_x20 = tk.Label(root, text="--", font=("Arial", 13))#效率提升
    lbl_x20.pack(padx=10, pady=5)
    lbl_x25 = tk.Label(root, text="--", font=("Arial", 13))#优化目标
    lbl_x25.pack(padx=10, pady=5)
    labels["root"] = lbl_x0
    labels["tip"] = lbl_x1
    labels["LD"] = lbl_x2
    labels["alpha"] = lbl_x3
    labels["总质量"] = lbl_x4
    labels["机翼质量"] = lbl_x5
    labels["当前循环"] = lbl_x6
    labels["minlift"] = lbl_x7
    labels["RPM_C"] = lbl_x8
    labels["Cl_target"] = lbl_x9
    labels["thrust ratio"] = lbl_x13
    labels["根梢比"] = lbl_x14
    labels["平均弦长"] = lbl_x15
    labels["升阻比提升"] = lbl_x17
    labels["RPM_L"] = lbl_x18
    labels["L/P"] = lbl_x19
    labels["功率变化"] = lbl_x20
    labels["阻力残差"] = lbl_x23
    labels["span"] = lbl_x24
    labels["功率-目标"] = lbl_x25
    labels["thrust_ratio_landing"] = lbl_x26
    labels["Nprops"] = lbl_x27
    labels["ele_angle"] = lbl_x28
    labels["wing_angle"] = lbl_x29
    # 启动后台线程跑优化，避免阻塞主线程的 GUI 事件循环
    opt_thread = threading.Thread(target=run_optimization, daemon=True)
    opt_thread.start()
    # 启动 Tkinter 主循环
    root.mainloop()

main()
vsp.WriteVSPFile(file_name, vsp.SET_ALL)
print(f"模型已保存")
with open('history.csv', 'w', encoding='utf-8', newline='') as f:
    writer = csv.writer(f, delimiter=',')
    writer.writerow(['span'] + span_hst)
    writer.writerow(['taper'] + taper_hst)
    writer.writerow(['chord'] + chord_hst)
    writer.writerow(['angle'] + wingang_hst)
    writer.writerow(['LD_hst'] + ld_hst)
    writer.writerow(['LP_hst'] + LP_hst)
    writer.writerow(['RPMC'] + RPMC_hst)
    writer.writerow(['thrust_ratio'] + ratio_hst)
    writer.writerow(['Lift_L'] + liftL_hst)
    writer.writerow(['Target'] + target_hst)
    writer.writerow(['POWER'] + power_hst)