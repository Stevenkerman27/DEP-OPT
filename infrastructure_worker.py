import json
import os
import sys
import ctypes

import numpy as np
import openvsp as vsp

import infrastructure as opb


if os.name == "nt":
    error_mode = 0x0001 | 0x0002 | 0x8000
    ctypes.windll.kernel32.SetErrorMode(error_mode)


request_path = sys.argv[1]
response_path = sys.argv[2]

with open(request_path, "r", encoding="utf-8") as request_file:
    request = json.load(request_file)

os.chdir(request["workdir"])
opb.case_name = request["case_name"]
opb.file_name = request["file_name"]
opb.density = request["density"]

vsp.ClearVSPModel()
vsp.ReadVSPFile(request["model_file"])

drag, alpha, lift, netdrag, power, Cl_list, CMy = opb.runaero(
    request["CG"],
    request["AlphaStart_input"],
    request["AlphaEnd_input"],
    request["AlphaNpts_input"],
    request["air_spd"],
    request["wing_cfg"],
    request["Cl_target"],
    request["sol_config"],
    request["angle"],
    np.asarray(request["prop_D"], dtype=float),
    np.asarray(request["RPM"], dtype=float),
    np.asarray(request["Ct"], dtype=float),
    np.asarray(request["Cp"], dtype=float),
    reynolds=request["reynolds"],
    wake_num_iter=request["wake_num_iter"],
)

response = {
    "drag": drag,
    "alpha": alpha,
    "lift": lift,
    "netdrag": netdrag,
    "power": power,
    "Cl_list": Cl_list,
    "CMy": CMy,
}
with open(response_path, "w", encoding="utf-8") as response_file:
    json.dump(response, response_file, ensure_ascii=False)
