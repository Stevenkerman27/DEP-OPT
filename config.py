# Shared project configuration.

# Physical constants and units.
DENSITY = 1.225
MU = 1.85e-5
INCH_IN_M = 0.0254
G = 9.8
ULTIMATE_LOAD_FACTOR = 5.0
SAFETY_FACTOR = 1.5
MAX_SWEEP_LOCATION = 0.9954

# Structural and mass model.
MASS_PROP = {
    "S_density": 1.6,
    "CF_Strength": 450e6,
    "CF_rho": 2000,
    "fuse_mass": 1.0,
    "prop": [0.05, 0.05, 0.05, 0.12],
    "payload": 0.3,
}
D_MIN_MM = 5
D_MAX_MM = 20

# VSPAERO execution settings.
VSPAERO_CPU = 6
VSPAERO_WORKER_TIMEOUT_SECONDS = 900
VSPAERO_FARFIELD_FAST = 3
VSPAERO_FARFIELD_FINAL = 10
VSPAERO_WAKE_NODES_FAST = 16
VSPAERO_WAKE_NODES_FINAL = 32
VSPAERO_SOLVER_FAST = {
    "farfield": VSPAERO_FARFIELD_FAST,
    "wakenode": VSPAERO_WAKE_NODES_FAST,
}
VSPAERO_SOLVER_FINAL = {
    "farfield": VSPAERO_FARFIELD_FINAL,
    "wakenode": VSPAERO_WAKE_NODES_FINAL,
}
VSPAERO_CONTROL_SURFACE_SE_CONST_FLAG = 1
VSPAERO_MESH_QUALITY_MIN_EDGE = 1.0e-8
VSPAERO_MESH_QUALITY_MAX_ASPECT = 1.0e4

# Common aircraft geometry.
TAIL_POS = {"x": 0.615, "y": 0, "z": -0.03, "yr": 0}
AIRFOIL_CFG = {
    "filename": None,
    "Camber": 0.04,
    "CamberLoc": 0.4,
    "ThickChord": 0.12,
}
TAIL_AIRFOIL_CFG = {
    "filename": None,
    "Camber": 0,
    "CamberLoc": 0,
    "ThickChord": 0.12,
}
FLAP_CFG = {
    "name": "Flaperon",
    "c": 1,
    "Length_Start": 0.3,
    "Length_End": 0.25,
    "EtaStart": 0.8,
    "EtaEnd": 0.07,
}
ELEVATOR_CFG = {
    "name": "elevator",
    "c": 0,
    "Length_Start": 0.05,
    "Length_End": 0.05,
    "EtaStart": 0.727,
    "EtaEnd": 0,
}
TAIL_CFG = {"root": 0.168, "tip": 0.095, "span": 0.24}
PROP_DIAMETERS_INCH = {"lift": 8, "tip": 13}
FUSELAGE_HALF_WIDTH = 0.07

# Main optimization configuration.
AUTO_WING_POS = {"name": "Mainwing", "x": 0, "y": 0, "z": 0, "yr": 0}
AUTO_TAIL_POS = {"name": "HT", **TAIL_POS}
AUTO_DEF_CFG = {"elevator": -5, "Flaperon": -15}
INCLUDE_WEIGHT = 1
INCLUDE_TAKEOFF_LANDING = 1
CRUISE_SPEED = 15
LANDING_SPEED = 7
MAX_AOA = 10
CRUISE_CD0 = 0.08
LANDING_CD0 = 0.10
REFERENCE_CD0_S = 0.06
MESH_INTERVAL = 0.01

# Section-capacity and three-dimensional stall-limit analysis.
STALL_LIMIT_CONFIG = {
    "enabled": True,
    "short_chord_fraction": 0.5,
    "vlm_cl_smooth_window": 7,
    "vlm_cl_smooth_polyorder": 2,
    "flap_transition_width": 0.03,
    "flap_transition_segments": 16,
    "alpha_min": -8.0,
    "alpha_max": 18.0,
    "alpha_points": 61,
    "clmax_alpha_points": 19,
    "clmax_refine_half_width": 0.5,
    "clmax_refine_points": 9,
    "clmax_stability_cl_atol": 0.03,
    "clmax_stability_alpha_atol": 0.5,
    "clmax_fit_residual_atol": 0.005,
    "reynolds_cache_step": 5000.0,
    "root_iterations": 20,
    "cl_residual_atol": 1.0e-3,
    "n_crit": 9.0,
    "xtr_upper": 1.0,
    "xtr_lower": 1.0,
    "model_size": "xlarge",
}

XFOIL_CONFIG = {
    "executable": r"D:\3D\Projects\XFOIL6.99\xfoil.exe",
    "xtr_upper": 0.01,
    "xtr_lower": 0.55,
    "alpha_min_deg": 0.0,
    "alpha_max_deg": 18.0,
    "alpha_step": 0.25,
    "hinge_cache_step": 0.01,
    "max_iter": 200,
    "timeout_seconds": 60,
}

MAIN_WING_SHEET_COUNT = 2

SLSQP_ACC = 0.01
SLSQP_SENS_STEP = 0.1
SLSQP_MAXIT = 15
SLSQP_SENSITIVITY = {
    "span": 1.0,
    "Mean_chord": 1.0,
    "taper": 1.0,
    "wing_angle": 1.0,
    "thrust_ratio": 1.0,
    "thrust_ratio_landing": 1.0,
}

DESIGN_BOUNDS = {
    "span": (0.75, 1.0, 0.90),
    "Mean_chord": (0.12, 0.22, 0.20),
    "taper": (0.50, 1.0, 0.99),
    "wing_angle": (0.0, 2.0, 1.0),
    "thrust_ratio": (0.30, 0.82, 0.70),
    "thrust_ratio_landing": (0.30, 0.82, 0.40),
}

# Standalone evaluation configuration.
EVALUATION_WING_POS = {"name": "wing", "x": 0, "y": 0, "z": 0, "yr": 0}
EVALUATION_TAIL_POS = {"name": "HT", **TAIL_POS}
EVALUATION_DEF_CFG = {"elevator": -5, "Flaperon": 20}
EVALUATION_TYPICAL_SPEED = 12
EVALUATION_CD0 = 0.12
EVALUATION_MESH_INTERVAL = 0.0093
EVALUATION_KN = 0.1

EVALUATION_CONFIGURATIONS = {
    "opt": {
        "Mean_chord": 0.197,
        "taper": 0.59,
        "span": 1.0,
        "wing_angle": 2.0,
        "CG": 0.072748,
    },
    "base": {
        "Mean_chord": 0.24,
        "taper": 1.0,
        "span": 0.97,
        "wing_angle": 0.0,
        "CG": 0.072192137,
    },
}
