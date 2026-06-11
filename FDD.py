import numpy as np
import os
import infrastructure as opb
import matplotlib.pyplot as plt
import openvsp as vsp
span = 0.9
Nprops = 4
liftprop_Dia = 8
tipprop_Dia = 13

# 获取当前脚本所在目录
script_dir = os.path.dirname(os.path.abspath(__file__))
# 创建 outputs 子文件夹（如果不存在）
output_dir = os.path.join(script_dir, "outputs")
os.makedirs(output_dir, exist_ok=True)
os.chdir(output_dir)

wing_pos = {"name": "Mainwing", "x":0, "y":0, "z":0, "yr": 0}
airfoil_cfg = {
    "filename":None,
    "Camber" : 0.04,
    "CamberLoc": 0.4,
    "ThickChord": 0.12
}

sub_cfg = {
    "name": "flap",
    "c": 0,
    "Length_Start" : 0.07,
    "Length_End" : 0.04,
    "length" : 0.8
}

tess_int = 0.012

prop_pos, prop_D_inch, Nprops = opb.place_prop(span, liftprop_Dia, tipprop_Dia, 0, tess_int)
prop_D = np.array(prop_D_inch) * opb.inch_in_m

spanlist, chordlist, twistlist, wing_S, Cl_target = opb.FFD_wing(
    root=0.25, tip=0.15, span=span, twist=None,
    prop_pos=prop_pos, prop_D=prop_D,
    ctrl=[0.2, 0.2, 0.2, 0.1], ctrl_dist=0.6, air_spd = 15)

_ = opb.create_wing(wing_pos, spanlist, chordlist, twistlist, 0.01, tess_int, airfoil_cfg, sub_cfg)

vsp.WriteVSPFile("FDD.vsp3", vsp.SET_ALL)
print(f"模型已保存")