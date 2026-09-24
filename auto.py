import openvsp as vsp
import os
from pyoptsparse import SLSQP, Optimization
import numpy as np
import threading
import tkinter as tk
import prop
import infrastructure as opb
import config as project_config
import csv
import warnings
# This is for vsp3.41!!!
cfg = {"wing_S": 0, "bref": 0,"cref": 0}

wing_pos = dict(project_config.AUTO_WING_POS)
tail_pos = dict(project_config.AUTO_TAIL_POS)
airfoil_cfg = dict(project_config.AIRFOIL_CFG)
airfoiltail_cfg = dict(project_config.TAIL_AIRFOIL_CFG)
flap_cfg = dict(project_config.FLAP_CFG)
ELE_cfg = dict(project_config.ELEVATOR_CFG)
def_cfg = dict(project_config.AUTO_DEF_CFG)
tail_cfg = dict(project_config.TAIL_CFG)
mass_prop = dict(project_config.MASS_PROP)
opb.mass_prop = mass_prop

#setting
include_weight = project_config.INCLUDE_WEIGHT
include_TL = project_config.INCLUDE_TAKEOFF_LANDING

mass = 1

root = None
labels = {}  # 键为变量名，值为对应的 Label 控件

# 优化历史
ld_hst = []
LP_hst =[]
ratio_hst = []
landing_ratio_hst = []

power_hst = []
target_hst = []
target_delta_hst = []
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
ACC = project_config.SLSQP_ACC
g = project_config.G
opb.g = g
cruise_spd = project_config.CRUISE_SPEED
min_speed = project_config.LANDING_SPEED
max_AOA = project_config.MAX_AOA

density = project_config.DENSITY
opb.density = density

#flight setting
G = project_config.ULTIMATE_LOAD_FACTOR
SF = project_config.SAFETY_FACTOR
n_ult = G * SF
opb.SF = SF
opb.G = G
cD0_C = project_config.CRUISE_CD0
cD0_L = project_config.LANDING_CD0
cD0_S = project_config.REFERENCE_CD0_S
d0_C = opb.D0(cD0_C, cruise_spd, cD0_S) 
d0_L = opb.D0(cD0_L, min_speed, cD0_S) 

#analysis parameter
tess_interval = project_config.MESH_INTERVAL
iter = 0
sen_step = project_config.SLSQP_SENS_STEP
sens_scale = project_config.SLSQP_SENSITIVITY

span_min, span_max, span_value = project_config.DESIGN_BOUNDS["span"]
chord_min, chord_max, chord_value = project_config.DESIGN_BOUNDS["Mean_chord"]
taper_min, taper_max, taper_value = project_config.DESIGN_BOUNDS["taper"]
angle_min, angle_max, angle_value = project_config.DESIGN_BOUNDS["wing_angle"]
thrust_min, thrust_max, thrust_value = project_config.DESIGN_BOUNDS["thrust_ratio"]
_, _, landing_thrust_value = project_config.DESIGN_BOUNDS["thrust_ratio_landing"]

if not span_min <= span_value <= span_max:
    raise ValueError(f"span initial value {span_value} is outside [{span_min}, {span_max}]")
if not chord_min <= chord_value <= chord_max:
    raise ValueError(f"Mean_chord initial value {chord_value} is outside [{chord_min}, {chord_max}]")
if not taper_min <= taper_value <= taper_max:
    raise ValueError(f"taper initial value {taper_value} is outside [{taper_min}, {taper_max}]")
if not angle_min <= angle_value <= angle_max:
    raise ValueError(f"wing_angle initial value {angle_value} is outside [{angle_min}, {angle_max}]")
if not thrust_min <= thrust_value <= thrust_max:
    raise ValueError(f"thrust_ratio initial value {thrust_value} is outside [{thrust_min}, {thrust_max}]")
landing_thrust_min, landing_thrust_max, _ = project_config.DESIGN_BOUNDS["thrust_ratio_landing"]
if not landing_thrust_min <= landing_thrust_value <= landing_thrust_max:
    raise ValueError(f"thrust_ratio_landing initial value {landing_thrust_value} is outside [{landing_thrust_min}, {landing_thrust_max}]")

prop_choice = prop.prop_choice
liftprop_Dia = project_config.PROP_DIAMETERS_INCH["lift"]
tipprop_Dia = project_config.PROP_DIAMETERS_INCH["tip"]
fuse_w = project_config.FUSELAGE_HALF_WIDTH

# 获取当前脚本所在目录
script_dir = os.path.dirname(os.path.abspath(__file__))
# 创建 outputs 子文件夹（如果不存在）
output_dir = os.path.join(script_dir, "outputs")
os.makedirs(output_dir, exist_ok=True)
os.chdir(script_dir)
prop_data = prop.read_apce_grouped('data')
# 切换当前工作目录
os.chdir(output_dir)
slsqp_ifile = os.path.join(output_dir, f"SLSQP_{os.getpid()}.out")
slsqp_history = os.path.join(output_dir, f"opt_hist_{os.getpid()}.hst")

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
    span = span_min + float(x_dict["span"]) * (span_max - span_min)
    Mean_chord = chord_min + float(x_dict["Mean_chord"]) * (chord_max - chord_min)
    CG = Mean_chord / 3 
    taper = taper_min + float(x_dict["taper"]) * (taper_max - taper_min)
    #root and tip
    root = 2 * Mean_chord / (1 + taper)
    tip = root * taper
    #twist
    twist = 0
    #prop ele
    prop_ele = 0.01
    #wing angle
    wing_angle = angle_min + float(x_dict["wing_angle"]) * (angle_max - angle_min)
    wing_pos["yr"] = wing_angle
    
    # RPS ratio
    thrust_ratio = thrust_min + float(x_dict["thrust_ratio"]) * (thrust_max - thrust_min)
    thrust_ratio_landing = landing_thrust_min + float(x_dict["thrust_ratio_landing"]) * (landing_thrust_max - landing_thrust_min)
    # Record
    taper_hst.append(taper)
    span_hst.append(span)
    chord_hst.append(Mean_chord)
    ratio_hst.append(thrust_ratio)
    landing_ratio_hst.append(thrust_ratio_landing)
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

    max_AOA_case = max_AOA - wing_angle
    flight_condition = {"speed": cruise_spd, "max_AOA":max_AOA_case, "Cl_target":Cl_target, "TR": thrust_ratio, "d0_others":d0_C }
    geo_info = {"spanlist":spanlist, "chordlist":chordlist, "span":span, "wing_S":wing_S, "bref":2*span, "cref":Mean_chord,
                "CG": CG, "def_cfg": def_cfg, "prop_D":prop_D,"prop_D_inch":prop_D_inch, "prop_pos":prop_pos,
                "t_over_c": airfoil_cfg["ThickChord"]}
    broyden_config = opb.make_broyden_config(prop_data, max_AOA_case, mass * g)
    lift, drag, power, alpha, RPM, thrust, mass_result, ele_def = opb.single_point(
        flight_condition, geo_info, broyden_config
    )
    RPM, Ct, Cp = prop.equal_thrust(
        prop_data, thrust, cruise_spd, prop_D_inch, thrust_ratio
    )
    cruise_stall_limit = opb.run_stall_limit_sweep(
        CG,
        cruise_spd,
        geo_info,
        airfoil_cfg,
        angle=def_cfg.copy(),
        prop_D=prop_D,
        RPM=RPM,
        Ct=Ct,
        Cp=Cp,
        sol_config=opb.solver_config0,
        limit_config=dict(project_config.STALL_LIMIT_CONFIG),
        flap_cfg=flap_cfg,
        flap_deflection=0.0,
    )
    if not cruise_stall_limit["usable"]:
        raise RuntimeError(
            f"Cruise stall limit is unusable: {cruise_stall_limit['reason']}"
        )
    cruise_alpha_limit = cruise_stall_limit["alpha_limit"]
    flight_condition["max_AOA"] = cruise_alpha_limit
    broyden_config = opb.make_broyden_config(
        prop_data, cruise_alpha_limit, mass * g
    )
    lift, drag, power, alpha, RPM, thrust, mass_result, ele_def = opb.single_point(
        flight_condition, geo_info, broyden_config
    )

    mass = mass_result["mass"]
    wing_mass = mass_result["wing_mass"]
    d0_wing_C = opb.calculate_parasite_drag_from_lod(cruise_spd, airfoil_cfg["ThickChord"])
    drag_resC = thrust - drag - d0_C - d0_wing_C #巡航阻力残差
    lift_res = lift - mass*g
    
    LD = lift / (drag + d0_C + d0_wing_C)
    eff = mass / power
    power_hst.append(power)
    AR = (span*2)**2/wing_S
    ele_C.append(ele_def)
    alpha_hst.append(alpha)
    ld_hst.append(LD)
    LP_hst.append(eff)

    if include_TL:
        # 升力约束
        max_L_AOA = max_AOA_case
        def_cfg["Flaperon"] = flap_angle

        flight_condition = {"speed": min_speed, "max_AOA":max_L_AOA, "Cl_target":-1, "TR": thrust_ratio_landing, "d0_others":d0_L }
        geo_info["def_cfg"] = def_cfg

        broyden_config = opb.make_broyden_config(
            prop_data, max_L_AOA, mass * g
        )
        lift_L, drag, power_L, alpha, RPM_L, thrust, mass_L, ele_def = opb.single_point(
            flight_condition, geo_info, broyden_config
        )
        RPM_L, Ct_L, Cp_L = prop.equal_thrust(
            prop_data, thrust, min_speed, prop_D_inch, thrust_ratio_landing
        )
        landing_stall_limit = opb.run_stall_limit_sweep(
            CG,
            min_speed,
            geo_info,
            airfoil_cfg,
            angle=def_cfg.copy(),
            prop_D=prop_D,
            RPM=RPM_L,
            Ct=Ct_L,
            Cp=Cp_L,
            sol_config=opb.solver_config0,
            limit_config=dict(project_config.STALL_LIMIT_CONFIG),
            flap_cfg=flap_cfg,
            flap_deflection=flap_angle,
        )
        if not landing_stall_limit["usable"]:
            raise RuntimeError(
                f"Landing stall limit is unusable: {landing_stall_limit['reason']}"
            )
        max_L_AOA = landing_stall_limit["alpha_limit"]
        flight_condition["max_AOA"] = max_L_AOA
        broyden_config = opb.make_broyden_config(
            prop_data, max_L_AOA, mass * g
        )
        lift_L, drag, power_L, alpha, RPM_L, thrust, mass_L, ele_def = opb.single_point(
            flight_condition, geo_info, broyden_config
        )

        d0_wing_L = opb.calculate_parasite_drag_from_lod(min_speed, airfoil_cfg["ThickChord"])
        drag_resL = thrust - drag - d0_L - d0_wing_L #最小速度阻力残差

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

    if target_hst:
        target_delta_hst.append(target - target_hst[-1])
    else:
        target_delta_hst.append("")

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

def objfunc_sens(x_dict, funcs):
    dv_names = [
        "span",
        "Mean_chord",
        "taper",
        "wing_angle",
        "thrust_ratio",
        "thrust_ratio_landing",
    ]
    funcs_sens = {"obj": {}}
    fail = False
    base_obj = float(funcs["obj"])

    for name in dv_names:
        step = sen_step * sens_scale[name]
        x_perturbed = dict(x_dict)
        x_perturbed[name] = float(x_dict[name]) + step
        funcs_perturbed, fail_perturbed = objfunc(x_perturbed)
        funcs_sens["obj"][name] = (float(funcs_perturbed["obj"]) - base_obj) / step
        fail = fail or fail_perturbed

    return funcs_sens, fail

def run_optimization():
    # 初始化问题
    optProb = Optimization("Auto", objfunc)
    # 添加设计变量
    optProb.addVar("span", "c", lower = 0.0, upper = 1.0, value = (span_value - span_min) / (span_max - span_min))
    optProb.addVar("Mean_chord", "c", lower = 0.0, upper = 1.0, value = (chord_value - chord_min) / (chord_max - chord_min))
    optProb.addVar("taper", "c", lower = 0.0, upper = 1.0, value = (taper_value - taper_min) / (taper_max - taper_min))
    optProb.addVar("wing_angle", "c", lower = 0.0, upper = 1.0, value = (angle_value - angle_min) / (angle_max - angle_min))
    optProb.addVar("thrust_ratio", "c", lower = 0.0, upper = 1.0, value = (thrust_value - thrust_min) / (thrust_max - thrust_min))
    optProb.addVar("thrust_ratio_landing", "c", lower = 0.0, upper = 1.0, value = (landing_thrust_value - landing_thrust_min) / (landing_thrust_max - landing_thrust_min))
    # rst begin addObj
    optProb.addObj("obj")

    # Check optimization problem
    print(optProb)
    optProb.printSparsity()

    # 配置SLSQP参数（更大胆的收敛策略）
    optOptions = {
        "ACC": ACC,
        "MAXIT": project_config.SLSQP_MAXIT,
        "IPRINT": 2,
        "IFILE": slsqp_ifile,
    }
    opt = SLSQP(options=optOptions)
    sol = opt(optProb, sens=objfunc_sens, storeHistory=slsqp_history)
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
    writer.writerow(['thrust_ratio_landing'] + landing_ratio_hst)
    writer.writerow(['Lift_L'] + liftL_hst)
    writer.writerow(['Target'] + target_hst)
    writer.writerow(['Target_delta'] + target_delta_hst)
    writer.writerow(['POWER'] + power_hst)
