import openvsp as vsp
import os
import numpy as np
import prop
import infrastructure as opb
import config as project_config
import csv
import json
import shutil
from pathlib import Path
import stall_limit
from spanwise_lift import (
    plot_spanwise_lift_curves,
    smooth_vlm_cl_by_sheet,
    xfoil_clmax_by_strip,
)
cfg = {"wing_S": 0, "bref": 0,"cref": 0}

wing_pos = dict(project_config.EVALUATION_WING_POS)
tail_pos = dict(project_config.EVALUATION_TAIL_POS)
airfoil_cfg = dict(project_config.AIRFOIL_CFG)
airfoiltail_cfg = dict(project_config.TAIL_AIRFOIL_CFG)
flap_cfg = dict(project_config.FLAP_CFG)
ELE_cfg = dict(project_config.ELEVATOR_CFG)
def_cfg = dict(project_config.EVALUATION_DEF_CFG)
tail_cfg = dict(project_config.TAIL_CFG)
mass_prop = dict(project_config.MASS_PROP)
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
lift_at_vlm_stall_limit_hst = []
lift_margin_hst = []

current_config = os.environ.get("DEP_CONFIG", "base")
case_speed = float(os.environ.get("DEP_SPEED", "8"))
case_flap = float(os.environ.get("DEP_FLAP", "-15"))
case_thrust_ratio = float(os.environ.get("DEP_THRUST_RATIO", "0.5"))

case_name = os.environ.get("DEP_CASE_NAME", f"evaluate_{current_config}_{case_speed:g}")
file_name = case_name + ".vsp3"

opb.case_name = case_name
opb.file_name = file_name

# 参数
typ_speed = project_config.EVALUATION_TYPICAL_SPEED
max_AOA = project_config.MAX_AOA
density = project_config.DENSITY
opb.density = density
cD0 = project_config.EVALUATION_CD0
cD0_S = project_config.REFERENCE_CD0_S
G = project_config.ULTIMATE_LOAD_FACTOR
SF = project_config.SAFETY_FACTOR
n_ult = G * SF
opb.SF = SF
opb.G = G
g = project_config.G
opb.g = g

#analysis parameter
tess_interval = project_config.EVALUATION_MESH_INTERVAL

#定义飞机
Kn = project_config.EVALUATION_KN
fuse_w = project_config.FUSELAGE_HALF_WIDTH # half of width of fuselage
Nprops = 4
prop_choice = prop.prop_choice
liftprop_Dia = project_config.PROP_DIAMETERS_INCH["lift"]
tipprop_Dia = project_config.PROP_DIAMETERS_INCH["tip"]

evaluation_config = project_config.EVALUATION_CONFIGURATIONS[current_config]
Mean_chord = evaluation_config["Mean_chord"]
taper = evaluation_config["taper"]
span = evaluation_config["span"]
wing_pos["yr"] = evaluation_config["wing_angle"]
CG = evaluation_config["CG"]

root = 2 * Mean_chord / (1 + taper)
tip = root * taper
#twist
twist = 0
#prop ele
prop_ele = 0.01
#[15, 0, 0.73], [7.5, -20, 0.5]]
condition = [[case_speed, case_flap, case_thrust_ratio]] # speed, flap, thrust_ratio
# 获取当前脚本所在目录
script_dir = os.path.dirname(os.path.abspath(__file__))
# 创建 outputs 子文件夹（如果不存在）
output_dir = os.path.join(script_dir, "outputs", case_name)
os.makedirs(output_dir, exist_ok=True)
os.chdir(script_dir)
prop_data_dir = os.environ.get("DEP_PROP_DATA_DIR", "data")
if not os.path.isdir(prop_data_dir):
    raise FileNotFoundError(f"Propeller data directory does not exist: {os.path.abspath(prop_data_dir)}")
prop_data = prop.read_apce_grouped(prop_data_dir)
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

    flight_condition = {"speed": speed, "max_AOA":max_AOA, "Cl_target":Cl_target, "TR": thrust_ratio, "d0_others":d0 }
    geo_info = {"spanlist":spanlist, "chordlist":chordlist, "span":span, "wing_S":wing_S, "bref":2*span, "cref":Mean_chord,
                "CG": CG, "def_cfg": def_cfg, "prop_D":prop_D,"prop_D_inch":prop_D_inch, "prop_pos":prop_pos,
                "t_over_c": airfoil_cfg["ThickChord"]}
    config = opb.make_broyden_config(prop_data, max_AOA, mass * g)
    # 2 options of accuracy for wake, tol and fcators in sequence of drag, lift, moment, relaxtion factor in sequence of drag and lift
    lift, drag, power, alpha, RPM, thrust, mass_result, balanced_elevator = opb.single_point(flight_condition, geo_info, config)

    mass = mass_result["mass"]
    wing_mass = mass_result["wing_mass"]
    def_cfg["elevator"] = balanced_elevator

    lod_path = Path(opb._find_result_file((".lod", "_DegenGeom.lod")))
    lod_cases = stall_limit.read_lod_cases(lod_path)
    if len(lod_cases) != 1:
        raise RuntimeError(f"Low-speed evaluation LOD must contain one case: {lod_path}")
    lift_case = lod_cases[0]
    wing_rows = stall_limit.lod_case_wing_rows(lift_case)
    vortex_sheets = sorted({int(row["VortexSheet"]) for row in wing_rows})
    main_wing_sheet = max(
        vortex_sheets,
        key=lambda sheet_id: max(
            row["Yavg"] for row in wing_rows
            if int(row["VortexSheet"]) == sheet_id
        ),
    )
    lift_rows = [
        row for row in wing_rows
        if int(row["VortexSheet"]) == main_wing_sheet
    ]
    lift_span = np.asarray([row["Yavg"] for row in lift_rows], dtype=float)
    lift_v_ratio = np.asarray([row["V/Vref"] for row in lift_rows], dtype=float)
    lift_cl_vlm = np.asarray([
        row["Cl"] / row["V/Vref"] ** 2
        for row in lift_rows
    ], dtype=float)
    lift_cl_vlm_smoothed = smooth_vlm_cl_by_sheet(
        lift_rows,
        lift_cl_vlm,
        project_config.STALL_LIMIT_CONFIG["vlm_cl_smooth_window"],
        project_config.STALL_LIMIT_CONFIG["vlm_cl_smooth_polyorder"],
    )
    lift_reynolds = np.asarray([
        density * speed * row["V/Vref"] * row["Chord"] / project_config.MU
        for row in lift_rows
    ], dtype=float)

    flap_angles = np.zeros(len(lift_rows), dtype=float)
    hinge_points = np.ones(len(lift_rows), dtype=float)
    for index, row in enumerate(lift_rows):
        eta = row["Yavg"] / span
        eta_start = flap_cfg["EtaStart"]
        eta_end = flap_cfg["EtaEnd"]
        if min(eta_start, eta_end) <= eta <= max(eta_start, eta_end):
            span_fraction = (eta - eta_start) / (eta_end - eta_start)
            flap_length = flap_cfg["Length_Start"] + span_fraction * (
                flap_cfg["Length_End"] - flap_cfg["Length_Start"]
            )
            hinge_points[index] = 1.0 - flap_length
            flap_angles[index] = -case_flap

    xfoil_result = xfoil_clmax_by_strip(
        lift_reynolds,
        flap_angles,
        hinge_points,
        airfoil_cfg,
        project_config.STALL_LIMIT_CONFIG,
        project_config.XFOIL_CONFIG,
        output_dir,
    )
    lift_clmax_xfoil = xfoil_result["clmax"]
    if not np.isfinite(lift_clmax_xfoil).all():
        raise RuntimeError(
            "XFOIL returned no lift capacity for one or more main-wing strips"
        )
    strip_area = np.asarray([row["dArea"] for row in lift_rows], dtype=float)
    balanced_lod_path = Path(output_dir) / f"{case_name}_balanced.lod"
    shutil.copy2(lod_path, balanced_lod_path)
    stall_scan_case_name = f"{case_name}_stall_scan"
    stall_scan_model_path = Path(output_dir) / f"{stall_scan_case_name}.vsp3"
    vsp.WriteVSPFile(str(stall_scan_model_path), vsp.SET_ALL)

    stall_scan_rpm, stall_scan_ct, stall_scan_cp = prop.equal_thrust(
        prop_data,
        thrust,
        speed,
        prop_D_inch,
        thrust_ratio,
    )
    opb.case_name = stall_scan_case_name
    opb.file_name = str(stall_scan_model_path)
    (
        _,
        _,
        _,
        _,
        _,
        stall_scan_polar_cl,
        _,
    ) = opb._runaero_isolated(
        CG,
        project_config.STALL_LIMIT_CONFIG["alpha_min"],
        project_config.STALL_LIMIT_CONFIG["alpha_max"],
        project_config.STALL_LIMIT_CONFIG["alpha_points"],
        speed,
        cfg,
        0.0,
        opb.solver_config0,
        def_cfg.copy(),
        prop_D,
        stall_scan_rpm,
        stall_scan_ct,
        stall_scan_cp,
    )
    opb.case_name = case_name
    opb.file_name = file_name

    stall_scan_lod_path = Path(output_dir) / f"{stall_scan_case_name}.lod"
    stall_scan_cases = sorted(
        stall_limit.read_lod_cases(stall_scan_lod_path),
        key=stall_limit.lod_case_alpha,
    )
    if len(stall_scan_cases) < 2:
        raise RuntimeError("VLM stall scan returned fewer than two alpha cases")

    stall_scan_strip_keys = [
        (int(row["VortexSheet"]), int(row["TrailVort"]))
        for row in lift_rows
    ]
    stall_scan_cl = []
    for stall_scan_case in stall_scan_cases:
        stall_scan_rows = [
            row for row in stall_limit.lod_case_wing_rows(stall_scan_case)
            if int(row["VortexSheet"]) == main_wing_sheet
        ]
        current_strip_keys = [
            (int(row["VortexSheet"]), int(row["TrailVort"]))
            for row in stall_scan_rows
        ]
        if current_strip_keys != stall_scan_strip_keys:
            raise ValueError("VLM stall scan strip ordering does not match the balanced LOD")
        local_cl = np.asarray([
            row["Cl"] / row["V/Vref"] ** 2
            for row in stall_scan_rows
        ], dtype=float)
        local_cl = smooth_vlm_cl_by_sheet(
            stall_scan_rows,
            local_cl,
            project_config.STALL_LIMIT_CONFIG["vlm_cl_smooth_window"],
            project_config.STALL_LIMIT_CONFIG["vlm_cl_smooth_polyorder"],
        )
        stall_scan_cl.append(local_cl)

    stall_scan_alpha = np.asarray([
        stall_limit.lod_case_alpha(case)
        for case in stall_scan_cases
    ], dtype=float)
    stall_scan_cl = np.asarray(stall_scan_cl, dtype=float)
    stall_scan_crossings = []
    for case_index in range(len(stall_scan_cases) - 1):
        left_residual = stall_scan_cl[case_index] - lift_clmax_xfoil
        right_residual = stall_scan_cl[case_index + 1] - lift_clmax_xfoil
        for strip_index in range(len(lift_rows)):
            if left_residual[strip_index] < 0.0 <= right_residual[strip_index]:
                crossing_fraction = (
                    -left_residual[strip_index]
                    / (right_residual[strip_index] - left_residual[strip_index])
                )
                crossing_alpha = stall_scan_alpha[case_index] + (
                    stall_scan_alpha[case_index + 1] - stall_scan_alpha[case_index]
                ) * crossing_fraction
                stall_scan_crossings.append(
                    (crossing_alpha, case_index, strip_index, crossing_fraction)
                )
    if not stall_scan_crossings:
        raise RuntimeError(
            "VLM alpha sweep did not cross any main-wing XFOIL Clmax before "
            f"{project_config.STALL_LIMIT_CONFIG['alpha_max']:g} deg"
        )

    (
        stall_limit_alpha,
        stall_limit_case_index,
        critical_strip_index,
        stall_limit_fraction,
    ) = min(stall_scan_crossings, key=lambda crossing: crossing[0])
    stall_limit_cl_by_strip = (
        stall_scan_cl[stall_limit_case_index]
        + stall_limit_fraction
        * (
            stall_scan_cl[stall_limit_case_index + 1]
            - stall_scan_cl[stall_limit_case_index]
        )
    )
    critical_vlm_cl = stall_limit_cl_by_strip[critical_strip_index]
    critical_xfoil_cl = lift_clmax_xfoil[critical_strip_index]
    stall_scan_polar_cl = np.asarray(stall_scan_polar_cl, dtype=float)
    if len(stall_scan_polar_cl) != len(stall_scan_alpha):
        raise ValueError("VSPAERO polar count does not match the alpha scan cases")
    cl_at_stall_limit = float(
        np.interp(stall_limit_alpha, stall_scan_alpha, stall_scan_polar_cl)
    )
    lift_at_stall_limit = (
        0.5 * density * speed**2 * wing_S * cl_at_stall_limit
        + thrust * np.sin(np.deg2rad(stall_limit_alpha))
    )
    lift_margin = lift_at_stall_limit - lift
    lift_at_vlm_stall_limit_hst.append(lift_at_stall_limit)
    lift_margin_hst.append(lift_margin)

    spanwise_png = Path(output_dir) / f"spanwise_lift_V{speed:.1f}.png"
    spanwise_csv = Path(output_dir) / f"spanwise_lift_V{speed:.1f}.csv"
    lift_margin_json = Path(output_dir) / f"lift_margin_V{speed:.1f}.json"
    plot_spanwise_lift_curves(
        lift_span,
        lift_cl_vlm_smoothed,
        lift_clmax_xfoil,
        spanwise_png,
        (
            f"Low-speed lift at {speed:g} m/s\n"
            f"Balanced lift = {lift:.2f} N | "
            f"VLM lift at {stall_limit_alpha:.2f} deg = "
            f"{lift_at_stall_limit:.2f} N | "
            f"margin = {lift_margin:.2f} N"
        ),
        span,
        xfoil_transition=(
            project_config.XFOIL_CONFIG["xtr_upper"],
            project_config.XFOIL_CONFIG["xtr_lower"],
        ),
        xfoil_alpha_limit=project_config.XFOIL_CONFIG["alpha_max_deg"],
    )
    with spanwise_csv.open("w", newline="", encoding="utf-8-sig") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow([
            "VortexSheet",
            "TrailVort",
            "Yavg_m",
            "V_over_Vref",
            "Cl_VLM_local_raw",
            "Cl_VLM_local_smoothed",
            "Cl_capacity_XFOIL_peak_or_alpha_limit",
            "Reynolds",
            "Reynolds_XFOIL",
            "XFOIL_alpha_used_deg",
            "XFOIL_peak_captured",
            "XFOIL_alpha_limit_fallback",
            "Cl_VLM_at_stall_limit",
            "Strip_area_m2",
            "Flap_deflection_deg",
            "Hinge_point_chord_fraction",
        ])
        for index, row in enumerate(lift_rows):
            writer.writerow([
                int(row["VortexSheet"]),
                int(row["TrailVort"]),
                lift_span[index],
                lift_v_ratio[index],
                lift_cl_vlm[index],
                lift_cl_vlm_smoothed[index],
                lift_clmax_xfoil[index],
                lift_reynolds[index],
                xfoil_result["reynolds"][index],
                xfoil_result["alpha_peak"][index],
                bool(xfoil_result["peak_captured"][index]),
                bool(xfoil_result["alpha_limit_fallback"][index]),
                stall_limit_cl_by_strip[index],
                strip_area[index],
                flap_angles[index],
                hinge_points[index],
            ])

    lift_margin_result = {
        "speed_mps": float(speed),
        "balanced_lift_n": float(lift),
        "lift_at_vlm_stall_limit_n": float(lift_at_stall_limit),
        "lift_margin_n": float(lift_margin),
        "lift_margin_percent_of_balanced_lift": float(lift_margin / lift * 100.0),
        "stall_limit_alpha_deg": float(stall_limit_alpha),
        "stall_limit_critical_vortex_sheet": int(main_wing_sheet),
        "stall_limit_critical_trail_vortex": int(
            lift_rows[critical_strip_index]["TrailVort"]
        ),
        "stall_limit_critical_span_m": float(lift_span[critical_strip_index]),
        "stall_limit_critical_vlm_cl": float(critical_vlm_cl),
        "stall_limit_critical_xfoil_clmax": float(critical_xfoil_cl),
        "lift_margin_model": (
            "VSPAERO total lift at the first alpha where a main-wing VLM strip "
            "reaches XFOIL Clmax, minus balanced lift"
        ),
        "definition": (
            "Find the first main-wing strip crossing between VLM alpha-sweep cases "
            "using smoothed local-q-corrected VLM Cl and its XFOIL Clmax, linearly "
            "interpolate the crossing alpha, then evaluate total VLM lift at that "
            "alpha; subtract the balanced lift"
        ),
        "stall_scan_alpha_bracket_deg": [
            float(stall_scan_alpha[stall_limit_case_index]),
            float(stall_scan_alpha[stall_limit_case_index + 1]),
        ],
        "stall_scan_lod_path": str(stall_scan_lod_path.resolve()),
        "balanced_lod_path": str(balanced_lod_path.resolve()),
        "xfoil_peak_captured_strip_count": int(
            np.count_nonzero(xfoil_result["peak_captured"])
        ),
        "xfoil_alpha_limit_fallback_strip_count": int(
            np.count_nonzero(xfoil_result["alpha_limit_fallback"])
        ),
        "main_wing_vortex_sheet": main_wing_sheet,
        "main_wing_strip_count": len(lift_rows),
        "lod_path": str(lod_path.resolve()),
        "spanwise_csv_path": str(spanwise_csv.resolve()),
        "spanwise_plot_path": str(spanwise_png.resolve()),
    }
    lift_margin_json.write_text(
        json.dumps(lift_margin_result, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    #修正升力
    lift_hst.append(lift)
    d0_wing = opb.calculate_parasite_drag_from_lod(speed, airfoil_cfg["ThickChord"])
    LD = lift / (drag + d0 + d0_wing)
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
    writer.writerow(['LiftAtVLMStallLimit_N'] + lift_at_vlm_stall_limit_hst)
    writer.writerow(['LiftMargin_N'] + lift_margin_hst)
