import csv
import json
import os
from pathlib import Path
import shutil
import sys

PROJECT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_DIR))
os.environ["PATH"] = str(Path(sys.prefix) / "Library" / "bin") + os.pathsep + os.environ["PATH"]

import numpy as np
import openvsp as vsp

import config
import infrastructure as opb
import prop
import stall_limit
from spanwise_lift import (
    plot_spanwise_lift_curves,
    smooth_vlm_cl_by_sheet,
    xfoil_clmax_by_strip,
)


# Geometry values follow config.py; regenerate the LOD after changing the design.
# Edit DEFAULT_SPEED_MPS, or pass speed in m/s as the first command-line argument.
DEFAULT_SPEED_MPS = 9.0
SPEED_MPS = float(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_SPEED_MPS
WING_GEOMETRY = {
    "semispan_m": 0.9,
    "mean_chord_m": 0.2,
    "taper_ratio": 0.8,
    "incidence_deg": 0.0,
}
FLAP_DEFLECTION_DEG = -15.0
FLAP_GEOMETRY = dict(config.FLAP_CFG)
AIRFOIL_GEOMETRY = dict(config.AIRFOIL_CFG)
TAIL_AIRFOIL_GEOMETRY = dict(config.TAIL_AIRFOIL_CFG)
NEURALFOIL_XTR_UPPER = 1.0
NEURALFOIL_XTR_LOWER = 1.0
XFOIL_XTR_UPPER = 1.0
XFOIL_XTR_LOWER = 1.0
PLOT_XFOIL = True
LANDING_THRUST_RATIO = config.DESIGN_BOUNDS["thrust_ratio_landing"][2]
BOOTSTRAP_ALPHA_DEG = config.MAX_AOA - WING_GEOMETRY["incidence_deg"]


SPEED_TAG = f"{SPEED_MPS:g}mps"
OUTPUT_DIR = PROJECT_DIR / "outputs" / "dep_stall_limit_probe" / f"speed_{SPEED_TAG}"
CASE_STEM = f"landing_15deg_dep_balanced_{SPEED_TAG}_{os.getpid()}"
LOD_PATH = OUTPUT_DIR / f"{CASE_STEM}.lod"
SCAN_LOD_PATH = OUTPUT_DIR / f"{CASE_STEM}_stall_scan.lod"
PNG_PATH = OUTPUT_DIR / f"landing_spanwise_load_V{SPEED_MPS:g}mps.png"
CSV_PATH = OUTPUT_DIR / f"landing_spanwise_load_V{SPEED_MPS:g}mps.csv"
REPORT_PATH = OUTPUT_DIR / f"spanwise_load_analysis_V{SPEED_MPS:.1f}.json"
SEMISPAN_M = WING_GEOMETRY["semispan_m"]
MEAN_CHORD_M = WING_GEOMETRY["mean_chord_m"]

def short_chord_mask_by_sheet(rows):
    short_chord = np.zeros(len(rows), dtype=bool)
    for sheet_id in sorted({int(row["VortexSheet"]) for row in rows}):
        sheet_indices = np.asarray([
            index for index, row in enumerate(rows)
            if int(row["VortexSheet"]) == sheet_id
        ], dtype=int)
        sheet_chord = np.asarray(
            [rows[index]["Chord"] for index in sheet_indices],
            dtype=float,
        )
        short_chord[sheet_indices] = sheet_chord < 0.5 * np.median(sheet_chord)
    return short_chord


limit_config = dict(config.STALL_LIMIT_CONFIG)
limit_config["xtr_upper"] = NEURALFOIL_XTR_UPPER
limit_config["xtr_lower"] = NEURALFOIL_XTR_LOWER
xfoil_config = dict(config.XFOIL_CONFIG)
xfoil_config["xtr_upper"] = XFOIL_XTR_UPPER
xfoil_config["xtr_lower"] = XFOIL_XTR_LOWER
cache_step = limit_config["reynolds_cache_step"]

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
os.chdir(PROJECT_DIR)
prop_data = prop.read_apce_grouped(str(PROJECT_DIR / "data"))
os.chdir(OUTPUT_DIR)

opb.case_name = CASE_STEM
opb.file_name = CASE_STEM + ".vsp3"
opb.density = config.DENSITY
opb.g = config.G
opb.G = config.ULTIMATE_LOAD_FACTOR
opb.SF = config.SAFETY_FACTOR
opb.mass_prop = dict(config.MASS_PROP)
opb.ini_geom()

wing_position = dict(config.AUTO_WING_POS)
wing_position["yr"] = WING_GEOMETRY["incidence_deg"]
tail_position = dict(config.AUTO_TAIL_POS)
tail_geometry = dict(config.TAIL_CFG)
fuselage_half_width = config.FUSELAGE_HALF_WIDTH
root_chord = 2.0 * MEAN_CHORD_M / (1.0 + WING_GEOMETRY["taper_ratio"])
tip_chord = root_chord * WING_GEOMETRY["taper_ratio"]
twist = 0.0
prop_elevation = 0.01
wing_span_sections, wing_chord_sections, wing_twist_sections, wing_area, _ = (
    opb.think_trapwing(
        root_chord,
        tip_chord,
        SEMISPAN_M,
        fuselage_half_width,
        twist,
        SPEED_MPS,
    )
)
EXPECTED_WING_AREA_M2 = float(wing_area)
prop_position, prop_diameter_inch, prop_count = opb.place_prop(
    SEMISPAN_M,
    config.PROP_DIAMETERS_INCH["lift"],
    config.PROP_DIAMETERS_INCH["tip"],
    prop_elevation,
    config.MESH_INTERVAL,
    fuselage_half_width,
)
prop_position = prop_position.tolist()
prop_diameter = (
    np.asarray(prop_diameter_inch, dtype=float) * config.INCH_IN_M
).tolist()
opb.create_wing(
    wing_position,
    wing_span_sections,
    wing_chord_sections,
    wing_twist_sections,
    0.0,
    config.MESH_INTERVAL,
    AIRFOIL_GEOMETRY,
    FLAP_GEOMETRY,
)
opb.create_wing(
    tail_position,
    [tail_geometry["span"]],
    [tail_geometry["root"], tail_geometry["tip"]],
    [0.0],
    opb.max_sweeploc,
    config.MESH_INTERVAL,
    TAIL_AIRFOIL_GEOMETRY,
    dict(config.ELEVATOR_CFG),
)
vsp.Update()

ultimate_load_factor = config.ULTIMATE_LOAD_FACTOR * config.SAFETY_FACTOR
y_cp = (
    4.0 / np.pi
    + (ultimate_load_factor - 1.0)
    * (1.0 + 2.0 * WING_GEOMETRY["taper_ratio"])
    / (1.0 + WING_GEOMETRY["taper_ratio"])
) / (3.0 * ultimate_load_factor)
aircraft_mass_kg, wing_mass_kg = opb.mass_sim_iter(
    SEMISPAN_M,
    wing_area,
    y_cp,
    prop_position,
    wing_span_sections,
    wing_chord_sections,
)
opb.mass = aircraft_mass_kg
aircraft_weight_n = aircraft_mass_kg * config.G

deflection_config = dict(config.AUTO_DEF_CFG)
deflection_config["Flaperon"] = FLAP_DEFLECTION_DEG
wing_configuration = {
    "spanlist": wing_span_sections,
    "chordlist": wing_chord_sections,
    "span": SEMISPAN_M,
    "wing_S": wing_area,
    "bref": 2.0 * SEMISPAN_M,
    "cref": MEAN_CHORD_M,
    "CG": MEAN_CHORD_M / 3.0,
    "def_cfg": deflection_config,
    "prop_D": prop_diameter,
    "prop_D_inch": prop_diameter_inch,
    "prop_pos": prop_position,
    "t_over_c": AIRFOIL_GEOMETRY["ThickChord"],
}
landing_condition = {
    "speed": SPEED_MPS,
    "max_AOA": BOOTSTRAP_ALPHA_DEG,
    "Cl_target": -1,
    "TR": LANDING_THRUST_RATIO,
    "d0_others": opb.D0(
        config.LANDING_CD0,
        SPEED_MPS,
        config.REFERENCE_CD0_S,
    ),
}
broyden_config = opb.make_broyden_config(
    prop_data,
    BOOTSTRAP_ALPHA_DEG,
    aircraft_weight_n,
)
(
    balanced_lift,
    balanced_drag,
    balanced_power,
    balanced_alpha,
    balanced_rpm,
    balanced_thrust,
    mass_result,
    balanced_elevator,
) = opb.single_point(landing_condition, wing_configuration, broyden_config)
aircraft_mass_kg = float(mass_result["mass"])
wing_mass_kg = float(mass_result["wing_mass"])
aircraft_weight_n = aircraft_mass_kg * config.G
balanced_rpm, prop_ct, prop_cp = prop.equal_thrust(
    prop_data,
    balanced_thrust,
    SPEED_MPS,
    prop_diameter_inch,
    LANDING_THRUST_RATIO,
)

opb._runaero_isolated(
    wing_configuration["CG"],
    limit_config["alpha_min"],
    limit_config["alpha_max"],
    limit_config["alpha_points"],
    SPEED_MPS,
    wing_configuration,
    0.0,
    opb.solver_config0,
    deflection_config.copy(),
    prop_diameter,
    balanced_rpm,
    prop_ct,
    prop_cp,
)
LOD_PATH = Path(opb._find_result_file((".lod", "_DegenGeom.lod")))
shutil.copy2(LOD_PATH, SCAN_LOD_PATH)

scan_cases = stall_limit.read_lod_cases(SCAN_LOD_PATH)
for case in scan_cases:
    condition = case["condition"]
    if not np.isclose(condition["Vinf_"], SPEED_MPS, rtol=0.0, atol=1.0e-5):
        raise ValueError(f"LOD freestream speed does not match SPEED_MPS: {condition['Vinf_']}")
    if not np.isclose(condition["Bref_"], 2.0 * SEMISPAN_M, rtol=1.0e-3):
        raise ValueError(f"LOD span does not match WING_GEOMETRY: {condition['Bref_']}")
    if not np.isclose(condition["Cref_"], MEAN_CHORD_M, rtol=1.0e-3):
        raise ValueError(f"LOD mean chord does not match WING_GEOMETRY: {condition['Cref_']}")
    if not np.isclose(condition["Sref_"], EXPECTED_WING_AREA_M2, rtol=1.0e-3):
        raise ValueError(f"LOD reference area does not match WING_GEOMETRY: {condition['Sref_']}")

scan_cases.sort(key=stall_limit.lod_case_alpha)
capacity_cache = {}
capacity_by_case = {}
local_q_cases = []
local_q_capacity_by_case = {}

for case in scan_cases:
    rows = stall_limit.lod_case_main_wing_rows(case)
    reynolds = stall_limit.lod_case_reynolds(case, SPEED_MPS)
    flap_angles = np.zeros(len(rows), dtype=float)
    hinge_points = np.ones(len(rows), dtype=float)

    for index, row in enumerate(rows):
        eta = row["Yavg"] / SEMISPAN_M
        eta_start = FLAP_GEOMETRY["EtaStart"]
        eta_end = FLAP_GEOMETRY["EtaEnd"]
        if min(eta_start, eta_end) <= eta <= max(eta_start, eta_end):
            span_fraction = (eta - eta_start) / (eta_end - eta_start)
            flap_length = FLAP_GEOMETRY["Length_Start"] + span_fraction * (
                FLAP_GEOMETRY["Length_End"] - FLAP_GEOMETRY["Length_Start"]
            )
            hinge_points[index] = 1.0 - flap_length
            flap_angles[index] = -FLAP_DEFLECTION_DEG

    capacities = np.empty(len(rows), dtype=float)
    capacity_usable = np.empty(len(rows), dtype=bool)
    for index, (reynolds_value, deflection, hinge) in enumerate(
        zip(reynolds, flap_angles, hinge_points, strict=True)
    ):
        cache_key = (
            float(np.round(reynolds_value / cache_step) * cache_step),
            float(deflection),
            float(hinge),
        )
        if cache_key not in capacity_cache:
            capacity_cache[cache_key] = stall_limit.predict_section_clmax(
                AIRFOIL_GEOMETRY,
                cache_key[0],
                limit_config,
                flap_deflection=cache_key[1],
                hinge_point=cache_key[2],
            )
        capacities[index] = capacity_cache[cache_key]["cl_max"]
        capacity_usable[index] = capacity_cache[cache_key]["usable"]

    capacity_by_case[id(case)] = {
        "valid": True,
        "cl_max": capacities,
        "capacity_usable": capacity_usable,
    }

    local_q_rows = [row.copy() for row in case["rows"]]
    main_sheets = sorted({
        row["VortexSheet"] for row in local_q_rows if row["IsARotor"] == 0.0
    })[:config.MAIN_WING_SHEET_COUNT]
    for row in local_q_rows:
        if row["IsARotor"] == 0.0 and row["VortexSheet"] in main_sheets:
            row["Cl"] /= row["V/Vref"] ** 2
    main_wing_row_indices = np.asarray([
        index for index, row in enumerate(local_q_rows)
        if row["IsARotor"] == 0.0 and row["VortexSheet"] in main_sheets
    ], dtype=int)
    main_wing_rows = [local_q_rows[index] for index in main_wing_row_indices]
    smoothed_cl = smooth_vlm_cl_by_sheet(
        main_wing_rows,
        [row["Cl"] for row in main_wing_rows],
        limit_config["vlm_cl_smooth_window"],
        limit_config["vlm_cl_smooth_polyorder"],
    )
    for index, smoothed_value in zip(main_wing_row_indices, smoothed_cl, strict=True):
        local_q_rows[index]["Cl"] = float(smoothed_value)
    local_q_case = {"condition": case["condition"], "rows": local_q_rows}
    local_q_cases.append(local_q_case)
    local_q_capacity_by_case[id(local_q_case)] = capacity_by_case[id(case)]

strip_crossings = []
for case in local_q_cases:
    case_rows = stall_limit.lod_case_main_wing_rows(case)
    short_chord = short_chord_mask_by_sheet(case_rows)
    capacity_usable = local_q_capacity_by_case[id(case)]["capacity_usable"]
    if not capacity_usable[~short_chord].all():
        raise RuntimeError(
            "Natural-transition Clmax is unusable for a non-short-chord strip"
        )

for case_index in range(len(local_q_cases) - 1):
    left_case = local_q_cases[case_index]
    right_case = local_q_cases[case_index + 1]
    left_rows = stall_limit.lod_case_main_wing_rows(left_case)
    right_rows = stall_limit.lod_case_main_wing_rows(right_case)
    left_capacity = local_q_capacity_by_case[id(left_case)]["cl_max"]
    right_capacity = local_q_capacity_by_case[id(right_case)]["cl_max"]
    left_cl = np.asarray([row["Cl"] for row in left_rows], dtype=float)
    right_cl = np.asarray([row["Cl"] for row in right_rows], dtype=float)
    left_residual = left_cl - left_capacity
    right_residual = right_cl - right_capacity
    left_short_chord = short_chord_mask_by_sheet(left_rows)
    right_short_chord = short_chord_mask_by_sheet(right_rows)
    if len(left_rows) != len(right_rows):
        raise ValueError("Alpha-scan cases have different main-wing strip counts")
    for strip_index, (left_row, right_row) in enumerate(
        zip(left_rows, right_rows, strict=True)
    ):
        left_key = (left_row["VortexSheet"], left_row["TrailVort"])
        right_key = (right_row["VortexSheet"], right_row["TrailVort"])
        if left_key != right_key:
            raise ValueError("Alpha-scan main-wing strip ordering changed between cases")
        short_chord = (
            left_short_chord[strip_index] or right_short_chord[strip_index]
        )
        if (
            not short_chord
            and left_residual[strip_index] < 0.0 <= right_residual[strip_index]
        ):
            left_alpha = stall_limit.lod_case_alpha(left_case)
            right_alpha = stall_limit.lod_case_alpha(right_case)
            alpha_crossing = left_alpha + (right_alpha - left_alpha) * (
                -left_residual[strip_index]
                / (right_residual[strip_index] - left_residual[strip_index])
            )
            strip_crossings.append(
                (alpha_crossing, case_index, strip_index, left_key)
            )
if not strip_crossings:
    raise RuntimeError("No per-strip smoothed VLM Cl / Clmax crossing was found")
strip_crossings.sort(key=lambda crossing: crossing[0])
critical_alpha, crossing_case_index, crossing_strip_index, crossing_strip_key = (
    strip_crossings[0]
)
below_alpha = stall_limit.lod_case_alpha(local_q_cases[crossing_case_index])
above_alpha = stall_limit.lod_case_alpha(local_q_cases[crossing_case_index + 1])
critical_alpha = float(critical_alpha)
deflection_config["elevator"] = balanced_elevator
landing_condition["max_AOA"] = critical_alpha
broyden_config = opb.make_broyden_config(
    prop_data,
    critical_alpha,
    aircraft_weight_n,
)
(
    balanced_lift,
    balanced_drag,
    balanced_power,
    balanced_alpha,
    balanced_rpm,
    balanced_thrust,
    mass_result,
    balanced_elevator,
) = opb.single_point(landing_condition, wing_configuration, broyden_config)
aircraft_mass_kg = float(mass_result["mass"])
wing_mass_kg = float(mass_result["wing_mass"])
aircraft_weight_n = aircraft_mass_kg * config.G
LOD_PATH = Path(opb._find_result_file((".lod", "_DegenGeom.lod")))
cases = stall_limit.read_lod_cases(LOD_PATH)
if len(cases) != 1:
    raise RuntimeError(f"Final balanced LOD must contain exactly one case: {LOD_PATH}")

final_case = cases[0]
final_condition = final_case["condition"]
if not np.isclose(final_condition["Vinf_"], SPEED_MPS, rtol=0.0, atol=1.0e-5):
    raise ValueError(f"Final LOD freestream speed does not match SPEED_MPS: {final_condition['Vinf_']}")
if not np.isclose(final_condition["Bref_"], 2.0 * SEMISPAN_M, rtol=1.0e-3):
    raise ValueError(f"Final LOD span does not match WING_GEOMETRY: {final_condition['Bref_']}")
if not np.isclose(final_condition["Cref_"], MEAN_CHORD_M, rtol=1.0e-3):
    raise ValueError(f"Final LOD mean chord does not match WING_GEOMETRY: {final_condition['Cref_']}")
if not np.isclose(final_condition["Sref_"], EXPECTED_WING_AREA_M2, rtol=1.0e-3):
    raise ValueError(f"Final LOD reference area does not match generated wing: {final_condition['Sref_']}")

rows = stall_limit.lod_case_main_wing_rows(final_case)
span = np.asarray([row["Yavg"] for row in rows], dtype=float)
chord = np.asarray([row["Chord"] for row in rows], dtype=float)
v_ratio = np.asarray([row["V/Vref"] for row in rows], dtype=float)
cl_local_q = np.asarray([
    row["Cl"] / row["V/Vref"] ** 2
    for row in rows
], dtype=float)
cl_local_q_raw = cl_local_q.copy()
cl_local_q = smooth_vlm_cl_by_sheet(
    rows,
    cl_local_q,
    limit_config["vlm_cl_smooth_window"],
    limit_config["vlm_cl_smooth_polyorder"],
)
reynolds = config.DENSITY * SPEED_MPS * v_ratio * chord / config.MU

flap_angles = np.zeros(len(rows), dtype=float)
hinge_points = np.ones(len(rows), dtype=float)
for index, row in enumerate(rows):
    eta = row["Yavg"] / SEMISPAN_M
    eta_start = FLAP_GEOMETRY["EtaStart"]
    eta_end = FLAP_GEOMETRY["EtaEnd"]
    if min(eta_start, eta_end) <= eta <= max(eta_start, eta_end):
        span_fraction = (eta - eta_start) / (eta_end - eta_start)
        flap_length = FLAP_GEOMETRY["Length_Start"] + span_fraction * (
            FLAP_GEOMETRY["Length_End"] - FLAP_GEOMETRY["Length_Start"]
        )
        hinge_points[index] = 1.0 - flap_length
        flap_angles[index] = -FLAP_DEFLECTION_DEG

clmax = np.empty(len(rows), dtype=float)
for index, (reynolds_value, deflection, hinge) in enumerate(
    zip(reynolds, flap_angles, hinge_points, strict=True)
):
    cache_key = (
        float(np.round(reynolds_value / cache_step) * cache_step),
        float(deflection),
        float(hinge),
    )
    if cache_key not in capacity_cache:
        capacity_cache[cache_key] = stall_limit.predict_section_clmax(
            AIRFOIL_GEOMETRY,
            cache_key[0],
            limit_config,
            flap_deflection=cache_key[1],
            hinge_point=cache_key[2],
        )
    clmax[index] = capacity_cache[cache_key]["cl_max"]

short_chord_all = short_chord_mask_by_sheet(rows)
eligible_indices = np.flatnonzero(~short_chord_all)
if len(eligible_indices) == 0:
    raise RuntimeError("No non-short-chord strips are available for stall assessment")
critical_index_all = int(
    eligible_indices[np.argmax((cl_local_q - clmax)[eligible_indices])]
)
critical_sheet = int(rows[critical_index_all]["VortexSheet"])
sheet_mask = np.asarray([
    int(row["VortexSheet"]) == critical_sheet for row in rows
], dtype=bool)
rows = [row for index, row in enumerate(rows) if sheet_mask[index]]
span = span[sheet_mask]
chord = chord[sheet_mask]
v_ratio = v_ratio[sheet_mask]
cl_local_q = cl_local_q[sheet_mask]
cl_local_q_raw = cl_local_q_raw[sheet_mask]
reynolds = reynolds[sheet_mask]
flap_angles = flap_angles[sheet_mask]
hinge_points = hinge_points[sheet_mask]
clmax = clmax[sheet_mask]

clmax_xfoil = np.full(len(rows), np.nan, dtype=float)
xfoil_reynolds = np.full(len(rows), np.nan, dtype=float)
xfoil_alpha_peak = np.full(len(rows), np.nan, dtype=float)
xfoil_peak_captured = np.zeros(len(rows), dtype=bool)

if PLOT_XFOIL:
    xfoil_capacity = xfoil_clmax_by_strip(
        reynolds,
        flap_angles,
        hinge_points,
        AIRFOIL_GEOMETRY,
        limit_config,
        xfoil_config,
        OUTPUT_DIR,
    )
    clmax_xfoil = xfoil_capacity["clmax"]
    xfoil_reynolds = xfoil_capacity["reynolds"]
    xfoil_alpha_peak = xfoil_capacity["alpha_peak"]
    xfoil_peak_captured = xfoil_capacity["peak_captured"]

plot_order = np.argsort(span)
span = span[plot_order]
cl_local_q = cl_local_q[plot_order]
cl_local_q_raw = cl_local_q_raw[plot_order]
clmax = clmax[plot_order]
clmax_xfoil = clmax_xfoil[plot_order]
xfoil_reynolds = xfoil_reynolds[plot_order]
xfoil_alpha_peak = xfoil_alpha_peak[plot_order]
xfoil_peak_captured = xfoil_peak_captured[plot_order]
chord = chord[plot_order]
reynolds = reynolds[plot_order]
flap_angles = flap_angles[plot_order]
hinge_points = hinge_points[plot_order]
rows = [rows[index] for index in plot_order]
margin = clmax - cl_local_q

sref = float(final_condition["Sref_"])
wing_cl = 2.0 * np.sum(
    cl_local_q_raw * np.asarray([row["dArea"] for row in rows], dtype=float)
) / sref
wing_lift_n = 0.5 * config.DENSITY * SPEED_MPS**2 * sref * wing_cl
lift_deficit_n = aircraft_weight_n - wing_lift_n
short_chord = short_chord_mask_by_sheet(rows)
cl_local_q_line = cl_local_q.copy()
clmax_line = clmax.copy()
clmax_xfoil_line = clmax_xfoil.copy()
cl_local_q_line[short_chord] = np.nan
clmax_line[short_chord] = np.nan
clmax_xfoil_line[short_chord] = np.nan

eligible_indices = np.flatnonzero(~short_chord)
critical_index = int(
    eligible_indices[np.argmax((cl_local_q - clmax)[eligible_indices])]
)
plot_spanwise_lift_curves(
    span,
    cl_local_q_line,
    clmax_xfoil_line if PLOT_XFOIL else np.full(len(span), np.nan),
    PNG_PATH,
    (
        f"Landing at {SPEED_MPS:g} m/s: VSPAERO Cl vs section Clmax "
        f"at alpha = {critical_alpha:.2f} deg\n"
        f"VSPAERO integrated wing CL = {wing_cl:.3f} | "
        f"wing lift deficit = {lift_deficit_n:.2f} N"
    ),
    SEMISPAN_M,
    reference_clmax=clmax_line,
    xfoil_transition=(
        XFOIL_XTR_UPPER,
        XFOIL_XTR_LOWER,
    ) if PLOT_XFOIL else None,
    xfoil_alpha_limit=xfoil_config["alpha_max_deg"] if PLOT_XFOIL else None,
)

with CSV_PATH.open("w", newline="", encoding="utf-8-sig") as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow([
        "VortexSheet",
        "TrailVort",
        "Yavg_m",
        "alpha_deg",
        "Cl_VSPAERO_local_q_corrected",
        "Cl_VSPAERO_raw_local_q_corrected",
        "Clmax_NeuralFoil_natural_transition",
        "Cl_capacity_XFOIL_peak_or_alpha_limit",
        "XFOIL_alpha_used_deg",
        "XFOIL_peak_captured",
        "XFOIL_alpha_limit_fallback",
        "Clmax_minus_Cl",
        "Reynolds",
        "Reynolds_XFOIL",
        "Chord_m",
        "short_chord_lod_strip",
        "flap_deflection_deg",
        "hinge_point_chord_fraction",
    ])
    for index, row in enumerate(rows):
        writer.writerow([
            int(row["VortexSheet"]),
            int(row["TrailVort"]),
            span[index],
            critical_alpha,
            cl_local_q[index],
            cl_local_q_raw[index],
            clmax[index],
            clmax_xfoil[index],
            xfoil_alpha_peak[index],
            bool(xfoil_peak_captured[index]),
            bool(xfoil_capacity["alpha_limit_fallback"][plot_order[index]])
            if PLOT_XFOIL else False,
            margin[index],
            reynolds[index],
            xfoil_reynolds[index],
            chord[index],
            bool(short_chord[index]),
            flap_angles[index],
            hinge_points[index],
        ])

analysis_report = {
    "speed_mps": SPEED_MPS,
    "lod_path": str(LOD_PATH),
    "wing_geometry": WING_GEOMETRY,
    "flap_deflection_deg": FLAP_DEFLECTION_DEG,
    "neuralfoil_transition": {
        "xtr_upper": NEURALFOIL_XTR_UPPER,
        "xtr_lower": NEURALFOIL_XTR_LOWER,
    },
    "xfoil_transition": {
        "xtr_upper": XFOIL_XTR_UPPER,
        "xtr_lower": XFOIL_XTR_LOWER,
    },
    "alpha_case_count": len(scan_cases),
    "alpha_limit_deg": critical_alpha,
    "alpha_bracket_deg": [below_alpha, above_alpha],
    "critical_crossing_strip": {
        "vortex_sheet": int(crossing_strip_key[0]),
        "trail_vortex": int(crossing_strip_key[1]),
        "lower_alpha_deg": float(below_alpha),
        "upper_alpha_deg": float(above_alpha),
    },
    "vlm_cl_smoothing": {
        "method": "Savitzky-Golay, independently by VortexSheet",
        "window_length": limit_config["vlm_cl_smooth_window"],
        "polynomial_order": limit_config["vlm_cl_smooth_polyorder"],
    },
    "stall_scan_lod_path": str(SCAN_LOD_PATH),
    "final_balanced_lod_path": str(LOD_PATH),
    "critical_sheet": critical_sheet,
    "critical_trail_vortex": int(rows[critical_index]["TrailVort"]),
    "integrated_wing_cl": float(wing_cl),
    "landing_balance": {
        "bootstrap_alpha_deg": float(BOOTSTRAP_ALPHA_DEG),
        "trim_alpha_deg": float(balanced_alpha),
        "elevator_deflection_deg": float(balanced_elevator),
        "thrust_n": float(balanced_thrust),
        "thrust_ratio": float(LANDING_THRUST_RATIO),
        "rpm": [float(value) for value in balanced_rpm],
        "mass_kg": aircraft_mass_kg,
        "wing_mass_kg": wing_mass_kg,
    },
    "lift_budget": {
        "aircraft_mass_kg": aircraft_mass_kg,
        "aircraft_weight_n": float(aircraft_weight_n),
        "wing_lift_n": float(wing_lift_n),
        "wing_lift_deficit_n": float(lift_deficit_n),
        "definition": "aircraft weight minus main-wing aerodynamic lift; positive means wing lift is short",
        "excludes": "horizontal-tail lift and propeller vertical force",
    },
    "xfoil_peak_captured_strip_count": int(np.count_nonzero(xfoil_peak_captured)),
    "strip_csv_path": str(CSV_PATH),
}
REPORT_PATH.write_text(
    json.dumps(analysis_report, ensure_ascii=False, indent=2),
    encoding="utf-8",
)

print(f"Plot: {PNG_PATH}")
print(f"Strip data: {CSV_PATH}")
print(f"Natural-transition NeuralFoil, XTR = {NEURALFOIL_XTR_UPPER}/{NEURALFOIL_XTR_LOWER}")
if PLOT_XFOIL:
    print(f"Natural-transition XFOIL, XTR = {XFOIL_XTR_UPPER}/{XFOIL_XTR_LOWER}")
    print(f"XFOIL interior-peak strips: {int(np.count_nonzero(xfoil_peak_captured))}/{len(xfoil_peak_captured)}")
print(f"Alpha limit: {critical_alpha:.6f} deg")
print(
    "VLM Cl smoothing: Savitzky-Golay "
    f"window={limit_config['vlm_cl_smooth_window']}, "
    f"polyorder={limit_config['vlm_cl_smooth_polyorder']}, "
    "independent by VortexSheet"
)
print(f"Integrated VSPAERO wing CL: {wing_cl:.6f}")
print(f"Landing trim alpha: {balanced_alpha:.6f} deg")
print(f"Landing trim elevator: {balanced_elevator:.6f} deg")
print(f"Landing trim thrust: {balanced_thrust:.3f} N")
print(f"Landing trim propeller RPM: {balanced_rpm}")
print(f"Aircraft mass: {aircraft_mass_kg:.6f} kg")
print(f"Aircraft weight: {aircraft_weight_n:.3f} N")
print(f"Main-wing aerodynamic lift: {wing_lift_n:.3f} N")
print(
    f"Main-wing lift deficit (weight - wing lift; excludes tail/prop vertical force): "
    f"{lift_deficit_n:.3f} N"
)
print(f"Critical strip: sheet {critical_sheet}, trail vortex {rows[critical_index]['TrailVort']}")
