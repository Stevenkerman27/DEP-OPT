import openvsp as vsp
import os
import numpy as np
import prop
import infrastructure as opb
import csv
cfg = {"wing_S": 0, "bref": 0,"cref": 0}

wing_pos = {"name": "wing", "x":0, "y":0, "z":0, "yr": 2} # 2 degre for opt
tail_pos = {"name": "HT", "x":0.615, "y":0, "z":-0.03, "yr": 0}

airfoil_cfg = {"filename": None, "Camber" : 0.04, "CamberLoc": 0.4, "ThickChord": 0.12}

airfoiltail_cfg = {"filename": None, "Camber" : 0, "CamberLoc": 0, "ThickChord": 0.12}

flap_cfg = {"name": "Flaperon", "c" :1, "Length_Start" : 0.3,"Length_End" : 0.25,"EtaStart" : 0.8, "EtaEnd": 0.07}

ELE_cfg = {"name": "elevator", "c" :0, "Length_Start" : 0.05,"Length_End" : 0.05,"EtaStart" : 0.727000, "EtaEnd": 0}

def_cfg = {"elevator": -5, "Flaperon": 20}

tail_cfg = {"root": 0.168,"tip": 0.095,"span": 0.24}

mass_prop = {"S_density": 1.6, "CF_Strength": 450e6, "CF_rho": 2000, "fuse_mass": 1.2, "prop": [0.05,0.05,0.05,0.12], "payload": 0.3}
opb.mass_prop = mass_prop
# 历史
ld_hst = []
LP_hst =[]
alpha_hst = []
thrust_hst=[]
lift_hst = []
clean_drag_hst = []
landing_drag_hst = []
landing_thrust_hst=[]
power_hst = []
target_hst = []
RPM_hst = []

ele_C = []
ele_L = []
ele_hst = []

case_name = "evaluate"
file_name = case_name + ".vsp3"

opb.case_name = case_name
opb.file_name = file_name

# 参数
typ_speed = 12
max_AOA = 10
density = 1.225
opb.density = density
cD0 = 0.12
cD0_S = 0.06
G = 5
SF = 1.5
n_ult = G * SF
opb.SF = SF
opb.G = G
g = 9.8
opb.g = g

#analysis parameter
tess_interval = 0.0093
drag_maxit = 9
omega_lift = 0.7
omega_drag = 0.7
max_step = 0.4  # deg

#定义飞机
Kn = 0.1
fuse_w = 0.07 # half of width of fuselage
Nprops = 4
prop_choice = prop.prop_choice
liftprop_Dia = 8
tipprop_Dia = 13

current_config = "opt"
if (current_config == "opt"):
    Mean_chord = 0.197
    taper = 0.59
    span = 1
    wing_pos["yr"] = 2
    CG = 0.072748 #opt
elif(current_config == "base"):
    Mean_chord = 0.24
    taper = 1
    span = 0.97
    wing_pos["yr"] = 0
    CG = 0.072192137 # baseline

root = 2 * Mean_chord / (1 + taper)
tip = root * taper
#twist
twist = 0
#prop ele
prop_ele = 0.01
#[15, 0, 0.73], [7.5, -20, 0.5]]
condition = [[7.5, -20, 0.5]] # speed, flap, thrust_ratio
# 获取当前脚本所在目录
script_dir = os.path.dirname(os.path.abspath(__file__))
# 创建 outputs 子文件夹（如果不存在）
output_dir = os.path.join(script_dir, "outputs")
os.makedirs(output_dir, exist_ok=True)
os.chdir(script_dir)
prop_data = prop.read_apce_grouped('data')
# 切换当前工作目录
os.chdir(output_dir)
vsp.ClearVSPModel()
opb.ini_geom()
#定义机翼
spanlist, chordlist, twistlist, wing_S, Cl_target = opb.think_trapwing(root, tip, span, fuse_w, twist, typ_speed)
flap_cfg["start_l"] = fuse_w
#定义螺旋桨
prop_pos, prop_D_inch, Nprops = opb.place_prop(span, liftprop_Dia, tipprop_Dia, prop_ele, tess_interval, fuse_w)
prop_D = np.array(prop_D_inch) * opb.inch_in_m
#创建主翼
opb.create_wing(wing_pos, spanlist, chordlist, twistlist, 0, tess_interval, airfoil_cfg, flap_cfg)
#创建尾翼
tail_span = [tail_cfg["span"]]
tail_chord = [tail_cfg["root"],tail_cfg["tip"]]
tail_twist = [0] 
opb.create_wing(tail_pos, tail_span, tail_chord, tail_twist, opb.max_sweeploc, tess_interval, airfoiltail_cfg, ELE_cfg)

#计算重量
y_cp = (4/np.pi + (n_ult-1) * (1 + 2 * taper)/(1 + taper))/ (3 * n_ult) * span
mass, wing_mass = opb.mass_sim_iter(span, wing_S, y_cp, prop_pos, spanlist, chordlist)
opb.mass = mass

# find NP

cfg["wing_S"] = wing_S
cfg["bref"] = 2 * span
cfg["cref"] = Mean_chord
vsp.WriteVSPFile("evaluate.vsp3")
print(f"模型已保存")
if (CG):
    print("CG known, skip")
else:
    CG = opb.cal_cg(Kn, cfg, max_AOA, typ_speed)
print("CG: " + str(CG))


# evaluate performance
for n in condition:
    max_AOA = max_AOA - wing_pos["yr"]
    speed = n[0]
    d0 = opb.D0(cD0, speed, cD0_S) 
    thrust_ratio = n[2]
    Cl_target =  mass * g / (0.5 * density * speed**2 * wing_S)
    def_cfg["Flaperon"] = n[1]

    flight_condition = {"speed": speed, "max_AOA":max_AOA, "Cl_target":Cl_target, "TR": thrust_ratio, "d0":d0 }
    geo_info = {"spanlist":spanlist, "chordlist":chordlist, "span":span, "wing_S":wing_S, "bref":Mean_chord, "cref":2*span,
                "CG": CG, "def_cfg": def_cfg, "prop_D":prop_D,"prop_D_inch":prop_D_inch, "prop_pos":prop_pos}
    config={"max_it":drag_maxit, "tol":[0.01, 0.01, 0.01], "relax":[omega_drag, omega_lift], "max_alpha_step":max_step,"propdata":prop_data} 
    # 2 options of accuracy for wake, tol and fcators in sequence of drag, lift, moment, relaxtion factor in sequence of drag and lift
    lift, drag, power, alpha, RPM, thrust, mass_result, _ = opb.single_point(flight_condition, geo_info, config)

    mass = mass_result["mass"]
    wing_mass = mass_result["wing_mass"]
    #修正升力
    lift_hst.append(lift)
    LD = lift / (drag + d0)
    ld_hst.append(LD)
    eff = mass / power
    power_hst.append(power)
    LP_hst.append(eff)
    RPM_hst.append(str(RPM[0])+"-"+str(RPM[-1]))
    alpha_hst.append(alpha)

if os.path.exists('history.csv'):
    os.remove('history.csv')
    
with open('history.csv', 'w', encoding='utf-8', newline='') as f:
    writer = csv.writer(f, delimiter=',')
    writer.writerow(['LD_hst'] + ld_hst)
    writer.writerow(['LP_hst'] + LP_hst)
    writer.writerow(['Lift_hst'] + lift_hst)
    writer.writerow(['RPM'] + RPM_hst)
    writer.writerow(['POWER'] + power_hst)
    writer.writerow(['AOA'] + alpha_hst)
    writer.writerow(['CG'] + [CG])
    writer.writerow(['Mass'] + [mass, wing_mass])
    writer.writerow(['Area'] + [wing_S])