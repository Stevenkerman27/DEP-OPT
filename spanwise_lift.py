from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter

import stall_limit


def smooth_vlm_cl_by_sheet(rows, cl_values, smooth_window, smooth_polyorder):
    smoothed_cl = np.asarray(cl_values, dtype=float).copy()
    sheet_ids = sorted({int(row["VortexSheet"]) for row in rows})
    for sheet_id in sheet_ids:
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
            smoothed_cl[sheet_indices],
        )
        uniform_cl = savgol_filter(
            uniform_cl,
            smooth_window,
            smooth_polyorder,
            mode="interp",
        )
        smoothed_cl[sheet_indices] = np.interp(
            sheet_span,
            uniform_span,
            uniform_cl,
        )
    return smoothed_cl


def xfoil_clmax_by_strip(
    reynolds,
    flap_deflections,
    hinge_points,
    airfoil_config,
    limit_config,
    xfoil_config,
    output_dir,
):
    output_dir = Path(output_dir)
    airfoil_path = output_dir / "xfoil_spanwise_section.dat"
    polar_path = output_dir / "xfoil_spanwise_polar.pol"
    base_coordinates = stall_limit.four_series_coordinates(
        airfoil_config["Camber"],
        airfoil_config["CamberLoc"],
        airfoil_config["ThickChord"],
    )

    reynolds = np.asarray(reynolds, dtype=float)
    flap_deflections = np.broadcast_to(
        np.asarray(flap_deflections, dtype=float),
        reynolds.shape,
    )
    hinge_points = np.broadcast_to(
        np.asarray(hinge_points, dtype=float),
        reynolds.shape,
    )
    clmax = np.full(len(reynolds), np.nan, dtype=float)
    reynolds_xfoil = np.full(len(reynolds), np.nan, dtype=float)
    alpha_peak = np.full(len(reynolds), np.nan, dtype=float)
    peak_captured = np.zeros(len(reynolds), dtype=bool)
    alpha_limit_fallback = np.zeros(len(reynolds), dtype=bool)
    capacity_cache = {}

    for index, (reynolds_value, deflection, hinge) in enumerate(
        zip(reynolds, flap_deflections, hinge_points, strict=True)
    ):
        rounded_hinge = float(
            np.round(hinge / xfoil_config["hinge_cache_step"])
            * xfoil_config["hinge_cache_step"]
        )
        coordinates, chord_scale, alpha_offset = (
            stall_limit._deflect_airfoil_coordinates(
                base_coordinates,
                rounded_hinge,
                deflection,
                limit_config["flap_transition_width"],
                limit_config["flap_transition_segments"],
            )
        )
        current_reynolds = float(
            np.round(
                reynolds_value * chord_scale
                / limit_config["reynolds_cache_step"]
            )
            * limit_config["reynolds_cache_step"]
        )
        cache_key = (current_reynolds, float(deflection), rounded_hinge)
        requested_alpha_end = xfoil_config["alpha_max_deg"]

        if cache_key not in capacity_cache:
            np.savetxt(
                airfoil_path,
                coordinates,
                fmt="%.9f",
                header="NACA four-series section with flap",
                comments="",
            )
            commands = [
                f"LOAD {airfoil_path.name}",
                "PANE",
                "OPER",
                f"VISC {current_reynolds:.8g}",
                "VPAR",
                f"XTR {xfoil_config['xtr_upper']:.8g} "
                f"{xfoil_config['xtr_lower']:.8g}",
                f"N {limit_config['n_crit']:.8g}",
                "",
                f"ITER {xfoil_config['max_iter']}",
                "PACC",
                polar_path.name,
                "",
                f"ASEQ {xfoil_config['alpha_min_deg']:.8g} "
                f"{requested_alpha_end:.8g} "
                f"{xfoil_config['alpha_step']:.8g}",
                "PACC",
                "",
                "QUIT",
                "",
            ]
            polar_path.unlink(missing_ok=True)
            process = subprocess.run(
                [xfoil_config["executable"]],
                input="\n".join(commands),
                cwd=output_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=xfoil_config["timeout_seconds"],
                check=True,
            )
            polar_rows = []
            for line in polar_path.read_text(
                encoding="ascii",
                errors="replace",
            ).splitlines():
                fields = line.split()
                if len(fields) >= 7:
                    try:
                        polar_rows.append([float(value) for value in fields[:7]])
                    except ValueError:
                        continue
            polar = np.asarray(polar_rows, dtype=float).reshape((-1, 7))
            if len(polar) == 0 or not np.isfinite(polar[:, 1]).any():
                raise RuntimeError(
                    f"XFOIL returned no finite polar at Re={current_reynolds}\n"
                    f"{process.stdout}"
                )
            polar = polar[np.argsort(polar[:, 0])]
            peak_index = int(np.nanargmax(polar[:, 1]))
            peak_captured_value = (
                peak_index < len(polar) - 1
                and polar[-1, 0] - polar[peak_index, 0] >= 0.5
            )
            if peak_captured_value:
                capacity_cache[cache_key] = (
                    float(polar[peak_index, 1]),
                    float(polar[peak_index, 0]),
                    True,
                    False,
                )
            else:
                limit_index = int(
                    np.argmin(np.abs(polar[:, 0] - requested_alpha_end))
                )
                alpha_distance = abs(
                    polar[limit_index, 0] - requested_alpha_end
                )
                if alpha_distance > xfoil_config["alpha_step"] / 2.0 + 1.0e-8:
                    raise RuntimeError(
                        f"XFOIL polar did not reach the {xfoil_config['alpha_max_deg']:g} deg "
                        f"limit at Re={current_reynolds}; nearest alpha="
                        f"{polar[limit_index, 0]:.3f}"
                    )
                capacity_cache[cache_key] = (
                    float(polar[limit_index, 1]),
                    float(polar[limit_index, 0]),
                    False,
                    True,
                )

        (
            clmax[index],
            alpha_peak[index],
            peak_captured[index],
            alpha_limit_fallback[index],
        ) = capacity_cache[cache_key]
        reynolds_xfoil[index] = current_reynolds

    return {
        "clmax": clmax,
        "reynolds": reynolds_xfoil,
        "alpha_peak": alpha_peak,
        "peak_captured": peak_captured,
        "alpha_limit_fallback": alpha_limit_fallback,
    }


def plot_spanwise_lift_curves(
    span,
    vlm_cl,
    xfoil_clmax,
    output_path,
    title,
    semispan,
    reference_clmax=None,
    xfoil_transition=None,
    xfoil_alpha_limit=None,
):
    order = np.argsort(span)
    span = np.asarray(span, dtype=float)[order]
    vlm_cl = np.asarray(vlm_cl, dtype=float)[order]
    xfoil_clmax = np.asarray(xfoil_clmax, dtype=float)[order]
    full_span = np.r_[-span[::-1], span]
    full_vlm_cl = np.r_[vlm_cl[::-1], vlm_cl]
    full_xfoil_clmax = np.r_[xfoil_clmax[::-1], xfoil_clmax]

    figure, axis = plt.subplots(figsize=(10, 5.4), constrained_layout=True)
    axis.plot(
        full_span,
        full_vlm_cl,
        color="#1769aa",
        linewidth=2.0,
        label="VSPAERO local section Cl",
    )
    if reference_clmax is not None:
        reference_clmax = np.asarray(reference_clmax, dtype=float)[order]
        full_reference_clmax = np.r_[reference_clmax[::-1], reference_clmax]
        axis.plot(
            full_span,
            full_reference_clmax,
            color="#e87500",
            linewidth=2.0,
            label="Reference section Clmax",
        )
        axis.fill_between(
            full_span,
            full_reference_clmax,
            full_vlm_cl,
            where=(
                np.isfinite(full_vlm_cl)
                & np.isfinite(full_reference_clmax)
                & (full_vlm_cl >= full_reference_clmax)
            ),
            interpolate=True,
            color="#c43c39",
            alpha=0.22,
            label="Cl >= Clmax",
        )
    xfoil_label = "XFOIL section lift capacity"
    if xfoil_alpha_limit is not None:
        xfoil_label += f" (peak or Cl at {xfoil_alpha_limit:g} deg)"
    if xfoil_transition is not None:
        xfoil_label += (
            f" (XTR {xfoil_transition[0]:g}/{xfoil_transition[1]:g})"
        )
    axis.plot(
        full_span,
        full_xfoil_clmax,
        color="#238636",
        linewidth=2.0,
        label=xfoil_label,
    )
    axis.set_xlabel("Spanwise position y (m)")
    axis.set_ylabel("Section lift coefficient")
    axis.set_title(title)
    axis.set_xlim(-semispan, semispan)
    axis.set_ylim(bottom=0.0)
    axis.grid(color="#888888", linewidth=0.6, alpha=0.25)
    axis.legend(frameon=False, loc="best")
    figure.savefig(output_path, dpi=220)
    plt.close(figure)
