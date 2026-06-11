import openvsp as vsp
import os
import math
import numpy as np
import infrastructure as opb
tess_int = 0.006

inch_in_m = 0.0254
prop_D = inch_in_m * 10
prop_pos = 0.16
prop_ele = 0.01
root = 0.2
tip = 0.05
span = 0.7
nSecs = 10
mass = 20
cruise_spd =15

wing_pos = {"name": "Mainwing", "x":0, "y":0, "z":0, "yr": 0}
flap_cfg = {
        "name" : "flap",
        "c": 1,
        "Length_Start" : 0.3,
        "Length_End" : 0.25,
        "length" : 0.8
    }

airfoil_cfg = {"Camber" : 0.04,"CamberLoc": 0.4,"ThickChord": 0.12}

# 获取当前脚本所在目录
script_dir = os.path.dirname(os.path.abspath(__file__))
# 创建 outputs 子文件夹（如果不存在）
output_dir = os.path.join(script_dir, "outputs")
os.makedirs(output_dir, exist_ok=True)
# 切换当前工作目录
os.chdir(output_dir)

spanlist, chordlist, twistlist, area, Cl = opb.generate_elliptical_wing(root, tip, span, nSecs, cruise_spd)
opb.create_wing(wing_pos, spanlist, chordlist, twistlist, 0, tess_int, airfoil_cfg, flap_cfg)
print("wing area is: " + str(area))
print("Coefficient of lift is: " + str(Cl))

# 找到 VSPAERO 设置的容器 ID
cs_group_container_id = vsp.FindContainer("VSPAEROSettings", 0)
defl_parm = vsp.FindParm(cs_group_container_id,
                            "DeflectionAngle",
                            "ControlSurfaceGroup_0")
vsp.SetParmVal(defl_parm, 30)  # 单位：度
vsp.Update()

vsp.WriteVSPFile("Ellip.vsp3", vsp.SET_ALL)
print(f"模型已保存")