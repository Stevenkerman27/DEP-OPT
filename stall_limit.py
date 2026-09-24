"""NeuralFoil section capacity and VSPAERO LOD parsing primitives."""

from pathlib import Path
import math
import os
import sys

os.environ["PATH"] = (
    os.path.join(sys.prefix, "Library", "bin")
    + os.pathsep
    + os.environ["PATH"]
)
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import neuralfoil as nf
import numpy as np
from scipy.signal import savgol_filter

import config as project_config


def four_series_coordinates(camber, camber_location, thickness, points=201):
    """Return closed-trailing-edge NACA four-series coordinates."""
    if points < 20:
        raise ValueError("At least 20 airfoil points are required")
    if not 0.0 <= camber < 1.0 or not 0.0 < camber_location < 1.0:
        raise ValueError("Invalid Four-Series camber parameters")
    if not 0.0 < thickness < 1.0:
        raise ValueError("Invalid Four-Series thickness")

    x = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, points)))
    thickness_shape = 5.0 * thickness * (
        0.2969 * np.sqrt(x)
        - 0.1260 * x
        - 0.3516 * x ** 2
        + 0.2843 * x ** 3
        - 0.1036 * x ** 4
    )
    camber_line = np.where(
        x < camber_location,
        camber / camber_location ** 2 * (2.0 * camber_location * x - x ** 2),
        camber / (1.0 - camber_location) ** 2
        * ((1.0 - 2.0 * camber_location) + 2.0 * camber_location * x - x ** 2),
    )
    slope = np.where(
        x < camber_location,
        2.0 * camber / camber_location ** 2 * (camber_location - x),
        2.0 * camber / (1.0 - camber_location) ** 2 * (camber_location - x),
    )
    theta = np.arctan(slope)
    upper = np.column_stack((x - thickness_shape * np.sin(theta), camber_line + thickness_shape * np.cos(theta)))
    lower = np.column_stack((x + thickness_shape * np.sin(theta), camber_line - thickness_shape * np.cos(theta)))
    upper[0] = [0.0, 0.0]
    lower[-1] = [1.0, 0.0]
    upper[-1] = [1.0, 0.0]
    lower[0] = [0.0, 0.0]
    coordinates = np.vstack((upper[::-1], lower[1:]))
    coordinates[:, 0] = np.clip(coordinates[:, 0], 0.0, 1.0)
    return coordinates


def _bernstein_basis(x, count):
    basis = np.empty((len(x), count), dtype=float)
    order = count - 1
    for index in range(count):
        basis[:, index] = (
            math.comb(order, index)
            * x ** index
            * (1.0 - x) ** (order - index)
        )
    return basis


def _fit_kulfan_parameters(coordinates):
    """Fit the fixed NeuralFoil eight-weight Kulfan representation."""
    coordinates = np.asarray(coordinates, dtype=float).copy()
    coordinates[:, 0] = np.clip(coordinates[:, 0], 0.0, 1.0)
    leading_edge = int(np.argmin(coordinates[:, 0]))
    if leading_edge == 0 or leading_edge == len(coordinates) - 1:
        raise ValueError("Airfoil coordinates must contain upper and lower surfaces")
    upper = coordinates[:leading_edge + 1][::-1]
    lower = coordinates[leading_edge:]
    upper_x = upper[:, 0]
    lower_x = lower[:, 0]
    upper_y = upper[:, 1]
    lower_y = lower[:, 1]
    upper_basis = np.sqrt(upper_x)[:, None] * (1.0 - upper_x)[:, None] * _bernstein_basis(upper_x, 8)
    lower_basis = np.sqrt(lower_x)[:, None] * (1.0 - lower_x)[:, None] * _bernstein_basis(lower_x, 8)
    upper_le = upper_x * (1.0 - upper_x) ** 8.5
    lower_le = lower_x * (1.0 - lower_x) ** 8.5
    upper_te = upper_x / 2.0
    lower_te = -lower_x / 2.0
    upper_design = np.column_stack((upper_basis, np.zeros((len(upper), 8)), upper_le, upper_te))
    lower_design = np.column_stack((np.zeros((len(lower), 8)), lower_basis, lower_le, lower_te))
    design = np.vstack((upper_design, lower_design)).tolist()
    target = np.concatenate((upper_y, lower_y)).tolist()
    size = len(design[0])
    normal = [[0.0 for _ in range(size)] for _ in range(size)]
    right = [0.0 for _ in range(size)]
    for row, value in zip(design, target, strict=True):
        for i in range(size):
            right[i] += row[i] * value
            for j in range(size):
                normal[i][j] += row[i] * row[j]
    for pivot in range(size):
        pivot_row = max(range(pivot, size), key=lambda row: abs(normal[row][pivot]))
        if abs(normal[pivot_row][pivot]) < 1.0e-12:
            raise ValueError("Kulfan least-squares system is singular")
        normal[pivot], normal[pivot_row] = normal[pivot_row], normal[pivot]
        right[pivot], right[pivot_row] = right[pivot_row], right[pivot]
        pivot_value = normal[pivot][pivot]
        for column in range(pivot, size):
            normal[pivot][column] /= pivot_value
        right[pivot] /= pivot_value
        for row in range(size):
            if row == pivot:
                continue
            factor = normal[row][pivot]
            for column in range(pivot, size):
                normal[row][column] -= factor * normal[pivot][column]
            right[row] -= factor * right[pivot]
    solution = np.asarray(right, dtype=float)
    if solution[-1] < 0.0:
        reduced_design = [row[:-1] for row in design]
        reduced_size = size - 1
        normal = [[0.0 for _ in range(reduced_size)] for _ in range(reduced_size)]
        right = [0.0 for _ in range(reduced_size)]
        for row, value in zip(reduced_design, target, strict=True):
            for i in range(reduced_size):
                right[i] += row[i] * value
                for j in range(reduced_size):
                    normal[i][j] += row[i] * row[j]
        for pivot in range(reduced_size):
            pivot_row = max(range(pivot, reduced_size), key=lambda row: abs(normal[row][pivot]))
            if abs(normal[pivot_row][pivot]) < 1.0e-12:
                raise ValueError("Bounded Kulfan least-squares system is singular")
            normal[pivot], normal[pivot_row] = normal[pivot_row], normal[pivot]
            right[pivot], right[pivot_row] = right[pivot_row], right[pivot]
            pivot_value = normal[pivot][pivot]
            for column in range(pivot, reduced_size):
                normal[pivot][column] /= pivot_value
            right[pivot] /= pivot_value
            for row in range(reduced_size):
                if row == pivot:
                    continue
                factor = normal[row][pivot]
                for column in range(pivot, reduced_size):
                    normal[row][column] -= factor * normal[pivot][column]
                right[row] -= factor * right[pivot]
        solution = np.r_[np.asarray(right, dtype=float), 0.0]
    if not np.isfinite(solution).all():
        raise ValueError("Kulfan least-squares fit returned nonfinite parameters")
    return {
        "upper_weights": solution[:8],
        "lower_weights": solution[8:16],
        "leading_edge_weight": solution[16],
        "TE_thickness": solution[17],
    }


def _circular_fillet_surface(
    surface,
    hinge_point,
    hinge_y,
    angle,
    transition_width,
    transition_segments,
):
    sample_count = transition_segments * 16
    tangent_length = transition_width / 2.0
    fixed_x = np.linspace(
        hinge_point - transition_width,
        hinge_point,
        sample_count + 1,
    )
    moving_x = np.linspace(
        hinge_point,
        hinge_point + transition_width,
        sample_count + 1,
    )
    fixed = np.column_stack((
        fixed_x,
        np.interp(fixed_x, surface[:, 0], surface[:, 1]),
    ))
    moving_y = np.interp(moving_x, surface[:, 0], surface[:, 1])
    moving_dx = moving_x - hinge_point
    moving_dy = moving_y - hinge_y
    moving = np.column_stack((
        hinge_point + math.cos(angle) * moving_dx - math.sin(angle) * moving_dy,
        hinge_y + math.sin(angle) * moving_dx + math.cos(angle) * moving_dy,
    ))

    fixed_tangent = np.gradient(fixed, axis=0)
    moving_tangent = np.gradient(moving, axis=0)
    fixed_tangent /= np.linalg.norm(fixed_tangent, axis=1)[:, None]
    moving_tangent /= np.linalg.norm(moving_tangent, axis=1)[:, None]
    tangent_cross = (
        fixed_tangent[-1, 0] * moving_tangent[0, 1]
        - fixed_tangent[-1, 1] * moving_tangent[0, 0]
    )
    tangent_dot = float(np.dot(fixed_tangent[-1], moving_tangent[0]))
    turn_angle = math.atan2(tangent_cross, tangent_dot)
    if abs(turn_angle) < 1.0e-6 or abs(turn_angle) >= math.pi:
        raise ValueError("Cannot construct a circular flap fillet for parallel tangents")

    radius = tangent_length / math.tan(abs(turn_angle) / 2.0)
    normal_side = -math.copysign(1.0, tangent_cross)
    fixed_normal = normal_side * np.column_stack((
        fixed_tangent[:, 1],
        -fixed_tangent[:, 0],
    ))
    moving_normal = normal_side * np.column_stack((
        moving_tangent[:, 1],
        -moving_tangent[:, 0],
    ))
    fixed_offset = fixed + radius * fixed_normal
    moving_offset = moving + radius * moving_normal
    center_distance = np.linalg.norm(
        fixed_offset[:, None, :] - moving_offset[None, :, :],
        axis=2,
    )
    fixed_index, moving_index = np.unravel_index(
        int(np.argmin(center_distance)),
        center_distance.shape,
    )
    if center_distance[fixed_index, moving_index] > transition_width * 0.02:
        raise ValueError("Circular flap fillet does not meet both airfoil surfaces")

    center = 0.5 * (
        fixed_offset[fixed_index] + moving_offset[moving_index]
    )
    fixed_contact = fixed[fixed_index]
    moving_contact = moving[moving_index]
    start_angle = math.atan2(
        fixed_contact[1] - center[1],
        fixed_contact[0] - center[0],
    )
    end_angle = math.atan2(
        moving_contact[1] - center[1],
        moving_contact[0] - center[0],
    )
    if turn_angle < 0.0:
        arc_angle = -((start_angle - end_angle) % (2.0 * math.pi))
    else:
        arc_angle = (end_angle - start_angle) % (2.0 * math.pi)
    arc_angles = np.linspace(
        start_angle,
        start_angle + arc_angle,
        transition_segments + 1,
    )
    arc = center + radius * np.column_stack((
        np.cos(arc_angles),
        np.sin(arc_angles),
    ))
    arc[0] = fixed_contact
    arc[-1] = moving_contact

    fixed_branch_x = np.linspace(
        hinge_point - transition_width,
        fixed_x[fixed_index],
        transition_segments + 1,
    )
    fixed_branch = np.column_stack((
        fixed_branch_x,
        np.interp(fixed_branch_x, surface[:, 0], surface[:, 1]),
    ))
    moving_branch_x = np.linspace(
        moving_x[moving_index],
        hinge_point + transition_width,
        transition_segments + 1,
    )
    moving_branch_y = np.interp(
        moving_branch_x,
        surface[:, 0],
        surface[:, 1],
    )
    moving_branch_dx = moving_branch_x - hinge_point
    moving_branch_dy = moving_branch_y - hinge_y
    moving_branch = np.column_stack((
        hinge_point
        + math.cos(angle) * moving_branch_dx
        - math.sin(angle) * moving_branch_dy,
        hinge_y
        + math.sin(angle) * moving_branch_dx
        + math.cos(angle) * moving_branch_dy,
    ))
    fixed_leading = surface[surface[:, 0] < hinge_point - transition_width]
    moving_trailing = surface[surface[:, 0] > hinge_point + transition_width]
    moving_trailing_dx = moving_trailing[:, 0] - hinge_point
    moving_trailing_dy = moving_trailing[:, 1] - hinge_y
    moving_trailing = np.column_stack((
        hinge_point
        + math.cos(angle) * moving_trailing_dx
        - math.sin(angle) * moving_trailing_dy,
        hinge_y
        + math.sin(angle) * moving_trailing_dx
        + math.cos(angle) * moving_trailing_dy,
    ))
    return np.vstack((
        fixed_leading,
        fixed_branch,
        arc[1:-1],
        moving_branch,
        moving_trailing,
    ))


def _deflect_airfoil_coordinates(
    coordinates,
    hinge_point,
    deflection,
    transition_width,
    transition_segments,
):
    coordinates = np.asarray(coordinates, dtype=float).copy()
    if deflection == 0.0 or hinge_point >= 1.0:
        return coordinates, 1.0, 0.0
    if not 0.0 < hinge_point < 1.0:
        raise ValueError("Control-surface hinge point must lie inside the chord")
    if not 0.0 < transition_width <= min(hinge_point, 1.0 - hinge_point):
        raise ValueError("Flap transition width must fit on both sides of the hinge")
    if not isinstance(transition_segments, int) or transition_segments < 4:
        raise ValueError("Flap transition must contain at least four coordinate segments")

    leading_edge = int(np.argmin(coordinates[:, 0]))
    upper = coordinates[:leading_edge + 1][::-1]
    lower = coordinates[leading_edge:]
    upper_y = float(np.interp(hinge_point, upper[:, 0], upper[:, 1]))
    lower_y = float(np.interp(hinge_point, lower[:, 0], lower[:, 1]))
    hinge_y = 0.5 * (upper_y + lower_y)

    angle = -math.radians(deflection)
    upper = _circular_fillet_surface(
        upper,
        hinge_point,
        hinge_y,
        angle,
        transition_width,
        transition_segments,
    )
    lower = _circular_fillet_surface(
        lower,
        hinge_point,
        hinge_y,
        angle,
        transition_width,
        transition_segments,
    )
    coordinates = np.vstack((upper[::-1], lower[1:]))
    leading_edge_index = int(np.argmin(coordinates[:, 0]))
    leading_edge = coordinates[leading_edge_index].copy()
    trailing_edge = 0.5 * (coordinates[0] + coordinates[-1])
    chord_vector = trailing_edge - leading_edge
    chord_scale = float(np.linalg.norm(chord_vector))
    chord_angle = math.atan2(chord_vector[1], chord_vector[0])
    cos_angle = math.cos(chord_angle)
    sin_angle = math.sin(chord_angle)
    relative = coordinates - leading_edge
    coordinates[:, 0] = (cos_angle * relative[:, 0] + sin_angle * relative[:, 1]) / chord_scale
    coordinates[:, 1] = (-sin_angle * relative[:, 0] + cos_angle * relative[:, 1]) / chord_scale
    coordinates[leading_edge_index] = [0.0, 0.0]
    coordinates[0] = [1.0, 0.0]
    coordinates[-1] = [1.0, 0.0]
    return coordinates, chord_scale, -math.degrees(chord_angle)


def _airfoil_parameters(
    airfoil_cfg,
    flap_deflection=0.0,
    hinge_point=1.0,
    transition_width=project_config.STALL_LIMIT_CONFIG["flap_transition_width"],
    transition_segments=project_config.STALL_LIMIT_CONFIG["flap_transition_segments"],
    return_transform=False,
):
    if airfoil_cfg["filename"] is not None:
        raise NotImplementedError("Airfoil files are deferred until the Four-Series path is validated")
    coordinates = four_series_coordinates(
        airfoil_cfg["Camber"],
        airfoil_cfg["CamberLoc"],
        airfoil_cfg["ThickChord"],
    )
    coordinates, chord_scale, alpha_offset = _deflect_airfoil_coordinates(
        coordinates,
        hinge_point,
        flap_deflection,
        transition_width,
        transition_segments,
    )
    parameters = _fit_kulfan_parameters(coordinates)
    if return_transform:
        return parameters, chord_scale, alpha_offset
    return parameters


def _neuralfoil_polar_from_parameters(
    parameters,
    alpha,
    reynolds,
    limit_config,
    alpha_offset=0.0,
    reynolds_scale=1.0,
):
    alpha_values = np.atleast_1d(np.asarray(alpha, dtype=float))
    reynolds_values = np.broadcast_to(np.asarray(reynolds, dtype=float), alpha_values.shape)
    rows = []
    for alpha_value, reynolds_value in zip(alpha_values, reynolds_values, strict=True):
        result = nf.get_aero_from_kulfan_parameters(
            parameters,
            alpha=float(alpha_value + alpha_offset),
            Re=float(reynolds_value * reynolds_scale),
            n_crit=limit_config["n_crit"],
            xtr_upper=limit_config["xtr_upper"],
            xtr_lower=limit_config["xtr_lower"],
            model_size=limit_config["model_size"],
        )
        rows.append({key: float(np.asarray(value).reshape(-1)[0]) for key, value in result.items()})
    return {key: np.asarray([row[key] for row in rows], dtype=float) for key in rows[0]}


def neuralfoil_polar(airfoil_cfg, alpha, reynolds, limit_config=None):
    limit_config = project_config.STALL_LIMIT_CONFIG if limit_config is None else limit_config
    parameters = _airfoil_parameters(airfoil_cfg)
    return _neuralfoil_polar_from_parameters(parameters, alpha, reynolds, limit_config)


def _local_quadratic_peak(alpha, cl, index):
    if index == 0 or index == len(alpha) - 1:
        return None, None
    x0, x1, x2 = alpha[index - 1:index + 2]
    y0, y1, y2 = cl[index - 1:index + 2]
    slope_left = (y1 - y0) / (x1 - x0)
    slope_right = (y2 - y1) / (x2 - x1)
    quadratic = (slope_right - slope_left) / (x2 - x0)
    if quadratic >= 0.0:
        return None, None
    linear = slope_left - quadratic * (x0 + x1)
    vertex = float(-linear / (2.0 * quadratic))
    fitted = quadratic * alpha[index - 1:index + 2] ** 2 + linear * alpha[index - 1:index + 2]
    constant = y1 - quadratic * x1 ** 2 - linear * x1
    fitted += constant
    residual = float(np.max(np.abs(fitted - cl[index - 1:index + 2])))
    if not x0 <= vertex <= x2:
        return residual, None
    return residual, vertex


def predict_section_clmax(
    airfoil_cfg,
    reynolds,
    limit_config=None,
    flap_deflection=0.0,
    hinge_point=1.0,
):
    """Find the usable finite NeuralFoil lift-capacity peak."""
    limit_config = project_config.STALL_LIMIT_CONFIG if limit_config is None else limit_config
    alpha = np.linspace(
        limit_config["alpha_min"],
        limit_config["alpha_max"],
        limit_config["clmax_alpha_points"],
    )
    if reynolds <= 0.0 or len(alpha) < 3:
        raise ValueError("Invalid sectional capacity sampling or Reynolds number")

    parameters, chord_scale, alpha_offset = _airfoil_parameters(
        airfoil_cfg,
        flap_deflection=flap_deflection,
        hinge_point=hinge_point,
        transition_width=limit_config["flap_transition_width"],
        transition_segments=limit_config["flap_transition_segments"],
        return_transform=True,
    )
    polar = _neuralfoil_polar_from_parameters(
        parameters,
        alpha,
        reynolds,
        limit_config,
        alpha_offset=alpha_offset,
        reynolds_scale=chord_scale,
    )
    cl = np.asarray(polar["CL"], dtype=float).reshape(-1)
    confidence = np.asarray(polar["analysis_confidence"], dtype=float).reshape(-1)
    if cl.shape != alpha.shape or confidence.shape != alpha.shape:
        raise ValueError("NeuralFoil returned an invalid sectional polar shape")
    finite = np.isfinite(cl)
    indices = np.flatnonzero(finite)
    if len(indices) == 0:
        raise ValueError("NeuralFoil returned no finite sectional capacity")

    index = int(indices[np.argmax(cl[indices])])
    raw_cl = float(cl[index])
    raw_alpha = float(alpha[index])
    if raw_cl <= 0.0:
        raise ValueError("NeuralFoil returned a nonpositive sectional capacity")

    fit_residual, fit_alpha = _local_quadratic_peak(alpha, cl, index)
    status = "boundary_limited" if index in (0, len(alpha) - 1) else "grid_peak"
    if fit_alpha is not None:
        status = "stable_peak" if fit_residual <= limit_config["clmax_fit_residual_atol"] else "grid_sensitive_peak"
    usable = status in ("stable_peak", "boundary_limited")
    used_cl = raw_cl
    used_alpha = raw_alpha
    if fit_alpha is not None:
        refined_alpha = np.linspace(
            max(alpha[0], fit_alpha - limit_config["clmax_refine_half_width"]),
            min(alpha[-1], fit_alpha + limit_config["clmax_refine_half_width"]),
            limit_config["clmax_refine_points"],
        )
        refined = _neuralfoil_polar_from_parameters(
            parameters,
            refined_alpha,
            reynolds,
            limit_config,
            alpha_offset=alpha_offset,
            reynolds_scale=chord_scale,
        )
        refined_cl = np.asarray(refined["CL"], dtype=float).reshape(-1)
        refined_confidence = np.asarray(refined["analysis_confidence"], dtype=float).reshape(-1)
        refined_index = int(np.nanargmax(refined_cl))
        used_cl = float(refined_cl[refined_index])
        used_alpha = float(refined_alpha[refined_index])
        used_confidence = float(refined_confidence[refined_index])
    else:
        used_confidence = float(confidence[index])

    return {
        "valid": True,
        "usable": usable,
        "reason": "ok" if usable else status,
        "cl_max": used_cl,
        "cl_max_raw": raw_cl,
        "alpha_peak": used_alpha,
        "confidence": used_confidence,
        "peak_quality": status,
        "fit_residual": fit_residual,
        "reynolds": float(reynolds),
    }


def predict_strip_clmax(
    airfoil_cfg,
    reynolds,
    limit_config=None,
    flap_deflection=None,
    hinge_point=None,
):
    limit_config = project_config.STALL_LIMIT_CONFIG if limit_config is None else limit_config
    reynolds = np.asarray(reynolds, dtype=float)
    if reynolds.ndim != 1 or len(reynolds) == 0 or not np.isfinite(reynolds).all() or np.any(reynolds <= 0.0):
        raise ValueError("Expected one positive finite Reynolds number per strip")
    flap_deflection = (
        np.zeros_like(reynolds)
        if flap_deflection is None
        else np.broadcast_to(np.asarray(flap_deflection, dtype=float), reynolds.shape)
    )
    hinge_point = (
        np.ones_like(reynolds)
        if hinge_point is None
        else np.broadcast_to(np.asarray(hinge_point, dtype=float), reynolds.shape)
    )
    if not np.isfinite(flap_deflection).all() or not np.isfinite(hinge_point).all():
        raise ValueError("Flap geometry must be finite")
    step = limit_config["reynolds_cache_step"]
    if step <= 0.0:
        raise ValueError("Reynolds cache step must be positive")
    cache = {}
    states = []
    for value, deflection, hinge in zip(
        reynolds,
        flap_deflection,
        hinge_point,
        strict=True,
    ):
        key = (
            float(np.round(value / step) * step),
            float(deflection),
            float(hinge),
        )
        if key not in cache:
            cache[key] = predict_section_clmax(
                airfoil_cfg,
                key[0],
                limit_config,
                flap_deflection=key[1],
                hinge_point=key[2],
            )
        states.append(cache[key])
    usable = all(state["valid"] and state["usable"] for state in states)
    return {
        "valid": all(state["valid"] for state in states),
        "usable": usable,
        "reason": "ok" if usable else "neuralfoil_capacity_unusable",
        "unusable_strip_indices": [
            index for index, state in enumerate(states) if not state["usable"]
        ],
        "states": states,
        "cl_max": np.asarray([state["cl_max"] for state in states], dtype=float),
        "reynolds": reynolds,
    }


def read_lod_cases(path):
    """Read every OpenVSP LOD case, preserving condition and strip rows."""
    path = Path(path)
    cases = []
    condition = {}
    names = None
    rows = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "Sref_":
            if names is not None:
                cases.append({"condition": condition, "rows": rows})
            condition = {}
            names = None
            rows = []
        if len(fields) >= 2 and fields[0].endswith("_"):
            condition[fields[0]] = float(fields[1])
        elif fields[:3] == ["Iter", "VortexSheet", "TrailVort"]:
            names = fields
        elif names is not None and len(fields) == len(names) and fields[0].isdigit():
            rows.append(dict(zip(names, map(float, fields), strict=True)))
    if names is not None:
        cases.append({"condition": condition, "rows": rows})
    if not cases or any(not case["rows"] for case in cases):
        raise ValueError("LOD contains an empty strip solution")
    return cases


def lod_case_alpha(case):
    return float(case["condition"]["AoA_"])


def lod_case_wing_rows(case):
    return [row for row in case["rows"] if row["IsARotor"] == 0.0]


def lod_case_main_wing_rows(case, sheet_count=None):
    sheet_count = project_config.MAIN_WING_SHEET_COUNT if sheet_count is None else sheet_count
    rows = lod_case_wing_rows(case)
    sheets = sorted({row["VortexSheet"] for row in rows})
    if len(sheets) < sheet_count:
        raise ValueError("LOD does not contain the configured main-wing vortex sheets")
    main_sheets = set(sheets[:sheet_count])
    return [row for row in rows if row["VortexSheet"] in main_sheets]


def lod_case_reynolds(case, speed, density=None, viscosity=None):
    density = project_config.DENSITY if density is None else density
    viscosity = project_config.MU if viscosity is None else viscosity
    rows = lod_case_main_wing_rows(case)
    if speed <= 0.0 or density <= 0.0 or viscosity <= 0.0:
        raise ValueError("Speed, density and viscosity must be positive")
    return np.asarray([
        density * speed * row["V/Vref"] * row["Chord"] / viscosity
        for row in rows
    ], dtype=float)


def lod_case_capacity_excess(case, cl_max):
    rows = lod_case_main_wing_rows(case)
    cl_max = np.asarray(cl_max, dtype=float)
    cl = np.asarray([row["Cl"] for row in rows], dtype=float)
    if cl_max.shape != cl.shape or not np.isfinite(cl_max).all():
        raise ValueError("One finite Clmax value is required per main-wing strip")
    residual = cl - cl_max
    index = int(np.argmax(residual))
    return {
        "alpha": lod_case_alpha(case),
        "residual": float(residual[index]),
        "critical_strip": index,
        "cl": cl,
        "cl_max": cl_max,
    }


def smooth_lod_main_wing_cl(case, smooth_window, smooth_polyorder):
    rows = lod_case_main_wing_rows(case)
    cl = np.asarray([
        row["Cl"] / row["V/Vref"] ** 2
        for row in rows
    ], dtype=float)
    for sheet_id in sorted({int(row["VortexSheet"]) for row in rows}):
        sheet_indices = np.asarray([
            index for index, row in enumerate(rows)
            if int(row["VortexSheet"]) == sheet_id
        ], dtype=int)
        sheet_indices = sheet_indices[
            np.argsort([rows[index]["Yavg"] for index in sheet_indices])
        ]
        if len(sheet_indices) < smooth_window:
            raise ValueError(
                f"VortexSheet {sheet_id} has {len(sheet_indices)} strips; "
                f"at least {smooth_window} are required for smoothing"
            )
        sheet_span = np.asarray(
            [rows[index]["Yavg"] for index in sheet_indices],
            dtype=float,
        )
        uniform_span = np.linspace(sheet_span[0], sheet_span[-1], len(sheet_span))
        uniform_cl = np.interp(
            uniform_span,
            sheet_span,
            cl[sheet_indices],
        )
        uniform_cl = savgol_filter(
            uniform_cl,
            smooth_window,
            smooth_polyorder,
            mode="interp",
        )
        cl[sheet_indices] = np.interp(
            sheet_span,
            uniform_span,
            uniform_cl,
        )
    return rows, cl


def lod_main_wing_short_chord_mask(rows, short_chord_fraction):
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
        short_chord[sheet_indices] = (
            sheet_chord < short_chord_fraction * np.median(sheet_chord)
        )
    return short_chord


def solve_lod_stall_limit(cases, capacity_predictor, limit_config=None):
    """Find the first eligible strip crossing of NeuralFoil capacity."""
    limit_config = (
        project_config.STALL_LIMIT_CONFIG
        if limit_config is None
        else limit_config
    )
    smooth_window = limit_config["vlm_cl_smooth_window"]
    smooth_polyorder = limit_config["vlm_cl_smooth_polyorder"]
    short_chord_fraction = limit_config["short_chord_fraction"]
    ordered = sorted(cases, key=lod_case_alpha)
    scan = []
    case_rows = []
    case_cl = []
    case_capacity = []
    case_short_chord = []
    capacity_usable = True
    for case in ordered:
        rows, cl = smooth_lod_main_wing_cl(
            case,
            smooth_window,
            smooth_polyorder,
        )
        capacity = capacity_predictor(case)
        cl_max = np.asarray(capacity["cl_max"], dtype=float)
        if cl_max.shape != cl.shape or not np.isfinite(cl_max).all():
            raise ValueError("One finite Clmax value is required per main-wing strip")
        short_chord = lod_main_wing_short_chord_mask(
            rows,
            short_chord_fraction,
        )
        if len(rows) != len(cl_max):
            raise ValueError("NeuralFoil capacities do not match main-wing strips")
        if case_rows:
            row_keys = [
                (row["VortexSheet"], row["TrailVort"])
                for row in rows
            ]
            previous_row_keys = [
                (row["VortexSheet"], row["TrailVort"])
                for row in case_rows[-1]
            ]
            if row_keys != previous_row_keys:
                raise ValueError("Alpha-scan main-wing strip ordering changed between cases")
        states = capacity["states"]
        if len(states) != len(rows):
            raise ValueError("NeuralFoil states do not match main-wing strips")
        capacity_usable = capacity_usable and all(
            state["valid"] and state["usable"]
            for index, state in enumerate(states)
            if not short_chord[index]
        )
        case_rows.append(rows)
        case_cl.append(cl)
        case_capacity.append(cl_max)
        case_short_chord.append(short_chord)

    strip_crossings = []
    for case_index in range(len(ordered) - 1):
        left_alpha = lod_case_alpha(ordered[case_index])
        right_alpha = lod_case_alpha(ordered[case_index + 1])
        left_residual = case_cl[case_index] - case_capacity[case_index]
        right_residual = case_cl[case_index + 1] - case_capacity[case_index + 1]
        excluded = case_short_chord[case_index] | case_short_chord[case_index + 1]
        for strip_index in range(len(left_residual)):
            if (
                not excluded[strip_index]
                and left_residual[strip_index] < 0.0 <= right_residual[strip_index]
            ):
                alpha_limit = left_alpha + (right_alpha - left_alpha) * (
                    -left_residual[strip_index]
                    / (right_residual[strip_index] - left_residual[strip_index])
                )
                strip_crossings.append((
                    float(alpha_limit),
                    case_index,
                    strip_index,
                ))

    if not strip_crossings:
        start_residual = case_cl[0] - case_capacity[0]
        eligible = ~case_short_chord[0]
        reason = (
            "stall_limit_below_scan"
            if np.any(start_residual[eligible] >= 0.0)
            else "stall_limit_above_scan"
        )
        return {
            "alpha_limit": None,
            "bracket": None,
            "start": {"alpha": lod_case_alpha(ordered[0])},
            "end": {"alpha": lod_case_alpha(ordered[-1])},
            "valid": False,
            "usable": False,
            "reason": reason,
            "capacity_usable": capacity_usable,
        }
    strip_crossings.sort(key=lambda crossing: crossing[0])
    alpha_limit, crossing_case_index, critical_strip = strip_crossings[0]
    left_alpha = lod_case_alpha(ordered[crossing_case_index])
    right_alpha = lod_case_alpha(ordered[crossing_case_index + 1])
    usable = np.isfinite(alpha_limit) and capacity_usable
    reason = "ok" if capacity_usable else "neuralfoil_capacity_unusable"
    return {
        "alpha_limit": float(alpha_limit),
        "bracket": [float(left_alpha), float(right_alpha)],
        "critical_strip": critical_strip,
        "critical_strip_key": (
            case_rows[crossing_case_index][critical_strip]["VortexSheet"],
            case_rows[crossing_case_index][critical_strip]["TrailVort"],
        ),
        "below": {
            "alpha": left_alpha,
            "residual": float(
                case_cl[crossing_case_index][critical_strip]
                - case_capacity[crossing_case_index][critical_strip]
            ),
        },
        "above": {
            "alpha": right_alpha,
            "residual": float(
                case_cl[crossing_case_index + 1][critical_strip]
                - case_capacity[crossing_case_index + 1][critical_strip]
            ),
        },
        "valid": np.isfinite(alpha_limit),
        "usable": usable,
        "reason": reason,
        "capacity_usable": capacity_usable,
    }
