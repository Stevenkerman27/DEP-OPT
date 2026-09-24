import openvsp as vsp
import math
import numpy as np
import pandas as pd
import re, io
import os
import json
import subprocess
import sys
import tempfile
import shutil
import uuid
import prop
import config as project_config
import stall_limit
import warnings
# This is for vsp3.41!!!
density = project_config.DENSITY
mu = project_config.MU
inch_in_m = project_config.INCH_IN_M
max_sweeploc = project_config.MAX_SWEEP_LOCATION
case_name = "test"
file_name = "test.vsp3"
G = project_config.ULTIMATE_LOAD_FACTOR
SF = project_config.SAFETY_FACTOR
g = project_config.G

mass = 3
mass_prop = dict(project_config.MASS_PROP)
D_min = project_config.D_MIN_MM
D_max = project_config.D_MAX_MM

CPU = project_config.VSPAERO_CPU
VSPAERO_WORKER_TIMEOUT_SECONDS = project_config.VSPAERO_WORKER_TIMEOUT_SECONDS
far_1 = project_config.VSPAERO_FARFIELD_FAST
far_2 = project_config.VSPAERO_FARFIELD_FINAL
wakeN_1 = project_config.VSPAERO_WAKE_NODES_FAST
wakeN_2 = project_config.VSPAERO_WAKE_NODES_FINAL
solver_config0 = dict(project_config.VSPAERO_SOLVER_FAST)
solver_config1 = dict(project_config.VSPAERO_SOLVER_FINAL)
VSPAERO_WING_GEOM_SET = vsp.SET_FIRST_USER
VSPAERO_PROP_GEOM_SET = vsp.SET_FIRST_USER + 1
VSPAERO_FUSELAGE_GEOM_SET = vsp.SET_FIRST_USER + 2
VSPAERO_THICK_GEOM_SET = VSPAERO_FUSELAGE_GEOM_SET
# OpenVSP 3.51 keeps PROP actuator disks only when the propeller geometry is
# included in the thin geometry set.  A no-propeller solve must exclude them.
VSPAERO_NO_PROP_THIN_GEOM_SET = VSPAERO_WING_GEOM_SET
VSPAERO_PROP_THIN_GEOM_SET = vsp.SET_ALL
VSPAERO_NO_GEOM_SET = vsp.SET_NONE
vspaero_geometry_signature = None
def get_lod_df():
    path = None
    for candidate in (case_name + ".lod", case_name + "_DegenGeom.lod"):
        if os.path.exists(candidate):
            path = candidate
            break
    if path is None:
        return None

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.read().splitlines()
    if not lines:
        return None

    header_idx = None
    header_type = None
    for i, line in enumerate(lines):
        if line.strip().startswith("Iter"):
            header_idx = i
            header_type = "iter"
            break
        if re.match(r"\s*Wing\s+", line):
            header_idx = i
            header_type = "wing"
            break
    if header_idx is None:
        return None

    data_lines = []
    for line in lines[header_idx + 1:]:
        if not line.strip():
            continue
        tokens = line.strip().split()
        if not tokens:
            continue
        component = tokens[1] if header_type == "iter" and len(tokens) > 1 else tokens[0]
        if component == "2":
            break
        data_lines.append(line)
    if not data_lines:
        return None

    df = pd.read_csv(io.StringIO("\n".join([lines[header_idx]] + data_lines)),
                     sep=r"\s+", engine="python")
    required = ("Yavg", "Cl", "V/Vref", "Chord")
    if any(column not in df.columns for column in required):
        return None
    if "dArea" not in df.columns and "S" not in df.columns:
        return None

    for column in required:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    if "dArea" in df.columns:
        df["dArea"] = pd.to_numeric(df["dArea"], errors="coerce")
    else:
        df["S"] = pd.to_numeric(df["S"], errors="coerce")
        df["dArea"] = np.abs(np.diff(df["S"].to_numpy(float), prepend=0.0))
    df = df.dropna(subset=["Yavg", "Cl", "V/Vref", "Chord", "dArea"]).reset_index(drop=True)
    return df

def read_halfwing_cp():
    df = get_lod_df()
    if df is None:
        raise ValueError("无法获取 .lod 数据")

    y = df["Yavg"].to_numpy(float)
    Cl = df["Cl"].to_numpy(float)
    dA = df["dArea"].to_numpy(float)
    V = df["V/Vref"].to_numpy(float)
    w = Cl * dA * V ** 2
    y_cp = np.sum(w * y) / np.sum(w)

    return y_cp

def calculate_parasite_drag_from_lod(Vref, t_over_c):
    if t_over_c is None:
        return 0.0
    df = get_lod_df()
    if df is None:
        return 0.0

    V_local = df["V/Vref"].to_numpy(float) * Vref
    chord = df["Chord"].to_numpy(float)
    dA = df["dArea"].to_numpy(float)
    Re = np.maximum((density * V_local * chord) / mu, 1e5)
    cf = 0.455 / (np.log10(Re) ** 2.58)
    ff = 1.0 + 1.2 * t_over_c + 100.0 * (t_over_c ** 4)
    s_wet_local = 2.003 * dA * (1.0 + 0.25 * t_over_c)
    dp_local = 0.5 * density * V_local ** 2 * s_wet_local * cf * ff
    return np.sum(dp_local) * 2

def calculate_parasite_drag(V, chord, S_ref, t_over_c):
    if t_over_c is None:
        return 0.0, 0.0
    Re = max(density * V * chord / mu, 1e5)
    cf = 0.455 / (math.log10(Re) ** 2.58)
    ff = 1.0 + 1.2 * t_over_c + 100.0 * (t_over_c ** 4)
    s_wet = 2.003 * S_ref * (1.0 + 0.25 * t_over_c)
    cd_p = s_wet / S_ref * cf * ff
    return cd_p, 0.5 * density * V ** 2 * s_wet * cf * ff

def tb_wingbox_mass(
    n_ult: float,         # 极限载荷因子
    M_G: float,           # 全机设计质量 [kg]
    eta_cp: float,        # 压力中心比距 η_cp
    b_st: float,          # 半翼展 [m]
    rho: float,           # 材料密度 [kg/m^3]
    sigma_r: float,       # 材料许用应力 [Pa]
    eta_t: float = 0.8,   # 厚度效率系数 (0.6–0.9)
    Rin: float = 1.0,     # 弯矩分布修正系数
    k0: float = 0.36,     # Torenbeek 基准常数
) -> float:
    Rcant = 1.0  # 全悬臂翼
    term = (1.05 * Rcant / eta_t + 3.67)
    W_box = k0 * n_ult * Rin * M_G * g * eta_cp * b_st * (rho / sigma_r) * term
    return W_box

def tb_rib(t_r, t_t, rho, k, S):
    m = rho * k * S * (t_r + t_t)
    return m

def mass_simpwing(span, S, y_c, force, t = 0.001):
    sigma_allow = mass_prop["CF_Strength"]
    M = y_c * force * g  * SF * G
    solution = 0
    for D in range(D_min, D_max):
        D = D / 1000
        Di = (D - 2.0 * t)
        I = math.pi / 64.0 * (D**4 - Di**4)
        sigma_max = M * D / 2 / I
        if sigma_max <= sigma_allow:
            solution = D
            print("requires " + str(D) + "m of CF tube")
            break
    if solution == 0:
        raise ValueError("未找到满足要求的外径。")
    CF_mass = mass_prop["CF_rho"] * span * math.pi * (D **2 - (Di)**2) / 4 

    prop_mass_total = sum(mass_prop["prop"]) * 2
    wing_mass = CF_mass * 2 + prop_mass_total + S * mass_prop["S_density"]

    mass =  wing_mass + mass_prop["fuse_mass"] + mass_prop["payload"]
    return mass, wing_mass

def mass_sim_iter(span, S, y_cp, prop_pos, spanlist, chordlist,
                  tol=1e-3, max_iter=20):
    # ====== 第一次，用压力中心算 ======
    m_total, m_wing = mass_simpwing(span, S, y_cp, (mass_prop["fuse_mass"] + mass_prop["payload"]) / 2.0)

    for it in range(max_iter):
        # ====== 计算机翼重心 ======
        # 假设翼面质量在平面上均匀分布，按展向积分重心
        # spanlist 和 chordlist 为等长数组，分别为每段展向长度和平均弦长
        # 积分变量 y 取为每段中点
        y_seg = [sum(spanlist[:i]) + s/2 for i, s in enumerate(spanlist)]
        area_seg = [s * c for s, c in zip(spanlist, chordlist)]
        total_area = sum(area_seg)
        y_cg = sum(y*a for y, a in zip(y_seg, area_seg)) / total_area

        M_prop = 0.0
        for mi, yi in zip(mass_prop["prop"], prop_pos):
            F_i = mi * g * SF * G    # 乘上安全系数和冲击倍数
            M_prop += F_i * yi

        eq_force = m_wing + (M_prop / (y_cg * g * SF * G))
        # === 调用 mass_simpwing 重新计算翼盒 ===
        m_total_new, m_wing_new = mass_simpwing(span, S, y_cg, eq_force)
        # ====== 检查收敛 ======
        if abs(m_wing_new - m_wing) / m_wing < tol:
            print(f"收敛于第 {it+1} 次迭代: m_wing={m_wing_new:.4f} kg, y_cg={y_cg:.4f} m")
            return m_total_new, m_wing_new

        print(f"iter {it+1}: m_wing={m_wing_new:.4f} kg, Δ={(m_wing_new-m_wing):.4f}")
        m_total, m_wing = m_total_new, m_wing_new

    raise RuntimeError("迭代未收敛，请检查输入参数或强度裕度。")

def next_tess_value(x):
    n = max((x - 6) // 4 + 1, 0)
    return 5 + n * 4

def _rdp(points, eps):
    #Ramer–Douglas–Peucker poly线简化
    start, end = points[0], points[-1]
    def point_line_dist(pt):
        x0,y0 = pt
        x1,y1 = start; x2,y2 = end
        num = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)
        den = np.hypot(y2-y1, x2-x1)
        return num/den if den!=0 else 0
    dmax, idx = 0.0, 0
    for i in range(1, len(points)-1):
        d = point_line_dist(points[i])
        if d > dmax:
            dmax, idx = d, i
    if dmax > eps:
        left = _rdp(points[:idx+1], eps)
        right = _rdp(points[idx:], eps)
        return left[:-1] + right
    else:
        return [start, end]

def FFD_wing(chord, tip, span, twist, prop_pos, prop_D, ctrl, ctrl_dist, air_spd, eps=2e-3):
    from scipy.interpolate import CubicSpline

    # 1. 构造控制点
    ctrl_pts = [[0, chord]]
    for n in range(len(prop_pos)):
        dist = ctrl_dist * prop_D[n] / 2
        # 上洗侧
        loc = prop_pos[n] - dist
        chord_loc = np.interp(loc, [0, span], [chord, tip])
        lift_out = chord_loc * ctrl[n]
        ctrl_pts.append([loc, chord_loc + lift_out])
        # 下洗侧
        loc = prop_pos[n] + dist
        if loc < span:
            chord_loc = np.interp(loc, [0, span], [chord, tip])
            lift_in = chord_loc * ctrl[n]
            ctrl_pts.append([loc, chord_loc - lift_in])
    ctrl_pts.append([span, tip])
    ctrl_pts = np.array(ctrl_pts)

    # 2. 三次样条（起点斜率=0）
    spline = CubicSpline(
        ctrl_pts[:,0], ctrl_pts[:,1],
        bc_type=((1, 0.0), 'not-a-knot')
    )

    # 3. 密采样并简化
    t_dense = np.linspace(0, span, 500)
    curve = spline(t_dense)
    pts = list(zip(t_dense, curve))
    simp = _rdp(pts, eps)

    # 4. 拆成段长列表和弦长列表
    xs, ys = zip(*simp)
    # spanlist: 每段长度 = 相邻 x 差值
    spanlist  = [xs[i+1] - xs[i] for i in range(len(xs)-1)]
    # chordlist: 对应每个节点（段头和段尾）的弦长
    chordlist = list(ys)
    twistlist = np.zeros_like(spanlist)

    wing_S = 0
    for x in range(0, len(spanlist)):
        wing_S += spanlist[x]*(chordlist[x]+chordlist[x+1])
    Cl_target =  mass * g / (0.5 * density * air_spd**2 * wing_S)
    return spanlist, chordlist, twistlist, wing_S, Cl_target

def think_trapwing(root, tip, span,fuse_w, twist, air_spd):
    if (fuse_w == 0):
        spanlist = [span]
        chordlist = [root, tip]
        twistlist = [twist]
    else: 
        spanlist = [fuse_w, span - fuse_w]
        chordlist = [root, root, tip]
        twistlist = [0,twist]
    wing_S = 0
    for x in range(0, len(spanlist)):
        wing_S += spanlist[x]*(chordlist[x]+chordlist[x+1])
    Cl_target =  mass * g / (0.5 * density * air_spd**2 * wing_S)
    # Add curve
    return spanlist, chordlist, twistlist, wing_S, Cl_target

def think_2section_trapwing(chord, tip, span, Xsec, twist, air_spd):
    spanlist = [Xsec, span - Xsec]
    chordlist = [chord, chord, tip]
    twistlist = [twist, twist]
    wing_S = 0
    for x in range(0, len(spanlist)):
        wing_S += spanlist[x]*(chordlist[x]+chordlist[x+1])
    Cl_target =  mass * g / (0.5 * density * air_spd**2 * wing_S)
    # Add curve
    return spanlist, chordlist, twistlist, wing_S, Cl_target

def generate_elliptical_wing(chord_root: float, chord_tip: float, semispan: float, nSecs: int, air_spd):
    # 1) 计算椭圆半展长 a
    a = chord_root * semispan / math.sqrt(chord_root**2 - chord_tip**2)
    print("Length to be modelled: " +str(a))
    # 2) 在参数区间 [0, a] 上均匀取 nSecs+1 个“截面”，对应弦长线性插值
    idx = np.arange(nSecs + 1)
    chord_nodes = chord_root + (chord_tip - chord_root) * idx / nSecs
    # 3) 根据椭圆方程 y = a * sqrt(1 - (x/a)^2)，这里 x 对应 chord_nodes/chord_root
    pos_nodes = a * np.sqrt(np.clip(1 - (chord_nodes / chord_root)**2, 0.0, 1.0))
    # 4) 将节点按展长从根部（0）排序，计算各段展长
    order = np.argsort(pos_nodes)
    pos_nodes_sorted = pos_nodes[order]
    chord_nodes_sorted = chord_nodes[order]
    spanlist = np.diff(pos_nodes_sorted).tolist()
    chordlist = chord_nodes_sorted.tolist()
    twistlist = np.zeros_like(spanlist)
    wing_S = 0
    for x in range(0, len(spanlist)):
        wing_S += spanlist[x]*(chordlist[x]+chordlist[x+1])
    Cl_target =  mass * g / (0.5 * density * air_spd**2 * wing_S)
    return spanlist, chordlist, twistlist, wing_S, Cl_target

def refine_wing_mesh(chordlist: list, spanlist: list, prop_pos_list: list, prop_D_list: list, wave_amp: list, tess_int: float, air_spd: float):
    # 1) 原始节点
    pos_nodes = np.concatenate(([0.0], np.cumsum(spanlist)))
    chord_nodes = np.array(chordlist)
    wing_span = pos_nodes[-1]

    # 2) 构建 refine_intervals（带 margin）和 sine_intervals（未裁剪）
    margin = tess_int
    refine_intervals = []
    sine_intervals  = []
    for prop_pos, prop_D in zip(prop_pos_list, prop_D_list):
        r = prop_D / 2
        # 未裁剪的正弧波区间
        s0u, s1u = prop_pos - r, prop_pos + r
        sine_intervals.append((s0u, s1u))
        # 用于加密网格的区间（裁剪到翼展范围后，再加 margin）
        s0 = max(s0u - margin, 0.0)
        s1 = min(s1u + margin, wing_span)
        if s1 > s0:
            refine_intervals.append((s0, s1))

    if len(wave_amp) != len(sine_intervals):
        raise ValueError("wave_amp 列表长度必须与螺旋桨数量一致")

    # 3) 保留区间外的原始节点，但始终保留翼根(0.0)和翼尖(span)两个节点
    mask_keep = np.ones_like(pos_nodes, dtype=bool)
    for s0, s1 in refine_intervals:
        mask_keep &= (pos_nodes < s0) | (pos_nodes > s1)
    mask_keep[0] = True
    mask_keep[-1] = True
    kept_pos = pos_nodes[mask_keep]

    # 4) 在 refine_intervals 内生成新节点
    new_pts = []
    for s0, s1 in refine_intervals:
        n = int(np.floor((s1 - s0) / tess_int))
        if n > 0:
            new_pts.append(np.linspace(s0, s1, n + 2, endpoint=True)[1:-1])
    new_pos = np.concatenate(new_pts) if new_pts else np.array([], dtype=float)

    # 5) 合并 kept_pos, new_pos，并强制插入 sine_intervals 两端的裁剪后边界点
    boundary_pts = []
    all_pos = np.unique(
        np.sort(np.concatenate((kept_pos, new_pos, boundary_pts))))

    # ———————— 去重过密节点，只删内部点 ————————
    threshold = tess_int * 0.1
    diffs = np.diff(all_pos)

    keep_mask = np.ones_like(all_pos, dtype=bool)
    # 对每个差值，如果太小，就删掉后一个节点，但不要删第一个和最后一个
    for i, d in enumerate(diffs):
        if d < threshold and 0 < (i+1) < (len(all_pos)-1):
            keep_mask[i+1] = False

    # 强保留首尾节点
    keep_mask[0] = True
    keep_mask[-1] = True

    all_pos = all_pos[keep_mask]

    # 6) 插值原始弦长
    chord_all = np.interp(all_pos, pos_nodes, chord_nodes)
    # 7) 在 sine_intervals 未裁剪区间上叠加正弧波
    wave_amp = np.array(wave_amp) * chordlist[-1] #不超过翼尖长度
    for (s0u, s1u), amp in zip(sine_intervals, wave_amp):
        # 先裁剪出 all_pos 中真正在翼上的那段
        s0c, s1c = max(s0u, 0.0), min(s1u, wing_span)
        mask = (all_pos >= s0c) & (all_pos <= s1c)
        # 用未裁剪的长度 (s1u - s0u) 来算 theta
        theta = (all_pos[mask] - s0u) / (s1u - s0u) * 2 * np.pi
        chord_all[mask] += amp * np.sin(theta)
    # ————— 使用 RDP 简化点集 —————
    # 构造 (位置, 弦长) 对
    pts = list(zip(all_pos.tolist(), chord_all.tolist()))
    # eps 
    simplified = _rdp(pts, 0.001)
    # 拆回位置和弦长
    simp_pos, simp_chord = zip(*simplified)

    # 8) 零扭转，长度为简化后段数
    twist_all = np.zeros(len(simp_pos) - 1)

    # 9) 重新计算简化后的 spanlist 和 chordlist
    spanlist_refined  = np.diff(simp_pos).tolist()
    chordlist_refined = list(simp_chord)
    twistlist_refined = twist_all.tolist()

    # 10) 面积 & Cl_target（根据原始翼面积计算）
    wing_S    = sum(spanlist[i] * (chordlist[i] + chordlist[i+1]) for i in range(len(spanlist)))
    Cl_target = mass * g / (0.5 * density * air_spd**2 * wing_S)

    return spanlist_refined, chordlist_refined, twistlist_refined, wing_S, Cl_target

def ini_geom():
    global vspaero_geometry_signature
    vsp.ClearVSPModel()
    vspaero_geometry_signature = None
    vsp.SetSetName(VSPAERO_WING_GEOM_SET, "wing")
    vsp.SetSetName(VSPAERO_PROP_GEOM_SET, "prop")
    vsp.SetSetName(VSPAERO_FUSELAGE_GEOM_SET, "fuselage")

def create_wing(pos, spanlist, chordlist, twistlist, sweeploc, tess_int, airfoil_cfg, sub_cfg = []):
    #position
    xloc = pos["x"]
    yloc = pos["y"]
    zloc = pos["z"]
    yrot = pos["yr"]
    #airfoil
    cam = airfoil_cfg["Camber"]
    cam_loc = airfoil_cfg["CamberLoc"]
    thick = airfoil_cfg["ThickChord"]
    nSecs = len(spanlist) 
    # Add a wing
    wid = vsp.AddGeom( "WING", "" )
    vsp.SetSetFlag(wid, VSPAERO_WING_GEOM_SET, True)
    vsp.Update()
    wing_name = pos["name"]
    vsp.SetGeomName(wid, wing_name)
    #set posotion
    vsp.SetParmVal( wid, "X_Rel_Location", "XForm", xloc)
    vsp.SetParmVal( wid, "Y_Rel_Location", "XForm", yloc)
    vsp.SetParmValUpdate( wid, "Z_Rel_Location", "XForm", zloc)
    vsp.SetParmValUpdate( wid, "Y_Rel_Rotation", "XForm", yrot)

    # Set symmetry
    sym_parm_wing = vsp.FindParm(wid, "Sym_Planar_Flag", "Sym") 
    vsp.SetParmVal(sym_parm_wing, 0)     # no symmetry
    tessW = next_tess_value(int(chordlist[0] / tess_int * 2))
    vsp.SetParmVal( wid, "Tess_W", "Shape", tessW)
    print("Set tess_W to :" +str(tessW))
    vsp.SetParmValUpdate( wid, "TECluster", "WingGeom", 0.9)
    vsp.SetParmValUpdate( wid, "LECluster", "WingGeom", 0.8)
    #===== Insert A Couple More Sections =====//
    print("Creating sections")
    for i in range(1,nSecs+1):
        vsp.InsertXSec( wid, 1, vsp.XS_FOUR_SERIES)
    vsp.CutXSec( wid, 1 )
    vsp.SetParmVal( wid, "Root_Chord", "XSec_1", chordlist[0])
    #set airfoil
    if airfoil_cfg["filename"]:
        xsec_surf_id = vsp.GetXSecSurf(wid, 0)
        vsp.ChangeXSecShape(xsec_surf_id, 0, vsp.XS_FILE_AIRFOIL)
        src_xsec_id = vsp.GetXSec(xsec_surf_id, 0)
        vsp.ReadFileAirfoil(src_xsec_id, airfoil_cfg["filename"])
    else:
        vsp.SetParmVal( wid, "ThickChord", "XSecCurve_0", thick)
        vsp.SetParmVal( wid, "Camber", "XSecCurve_0", cam)
        vsp.SetParmVal( wid, "CamberLoc", "XSecCurve_0", cam_loc)
    vsp.Update()
    Npanel = 0
    for i in range (1, nSecs+1):
        sec_index = "XSec_" + str(i)
        curve_index = "XSecCurve_" + str(i)
        tessU = round(spanlist[i-1]/((chordlist[i] / chordlist[0]) * tess_int)) + 1 #经验
        #print(sec_index,curve_index)
        vsp.SetParmVal( wid, "SectTess_U", sec_index, tessU)
        vsp.SetParmVal( wid, "InCluster", sec_index, chordlist[i-1] / chordlist[i])
        vsp.SetDriverGroup(wid, i, vsp.SPAN_WSECT_DRIVER, vsp.ROOTC_WSECT_DRIVER, vsp.TIPC_WSECT_DRIVER)#控制分段长度，翼根弦长，翼尖弦长
        vsp.SetParmVal( wid, "Span", sec_index, spanlist[i-1])
        vsp.SetParmVal( wid, "Tip_Chord", sec_index, chordlist[i])
        vsp.SetParmVal( wid, "Twist", sec_index, twistlist[i-1])
        vsp.SetParmVal( wid, "Sweep_Location", sec_index, sweeploc)
        Npanel = Npanel + tessU * tessW
        #set airfoil
        if airfoil_cfg["filename"]:
            xsec_surf_id = vsp.GetXSecSurf(wid, i)
            vsp.ChangeXSecShape(xsec_surf_id, i, vsp.XS_FILE_AIRFOIL)
            src_xsec_id = vsp.GetXSec(xsec_surf_id, i)
            vsp.ReadFileAirfoil(src_xsec_id, airfoil_cfg["filename"])
        else:
            vsp.SetParmVal( wid, "ThickChord", curve_index, thick)
            vsp.SetParmVal( wid, "Camber", curve_index, cam)
            vsp.SetParmValUpdate( wid, "CamberLoc", curve_index, cam_loc)
        vsp.Update()
    print("Wing Added")

    #Add CS
    if sub_cfg:
        len_start = sub_cfg["Length_Start"] 
        len_end = sub_cfg["Length_End"]
        EtaEnd = sub_cfg["EtaEnd"]
        EtaStart = sub_cfg["EtaStart"]
        subsurf_id = vsp.AddSubSurf(wid, vsp.SS_CONTROL, 0)
        vsp.Update()
        # 列出这个 SubSurface 的所有 Parm ID，以及它们的 Name 和 Display Group
        parm_id_vec = vsp.GetSubSurfParmIDs(subsurf_id)
        print("Found %d parms on sub-surf %s" % (len(parm_id_vec), subsurf_id))  #Debug用
        for pid in parm_id_vec:
            name  = vsp.GetParmName(pid)
            group = vsp.GetParmDisplayGroupName(pid)
            
            print("  ParmID = %-20s   Name = %-10s   Group = %s" % (pid, name, group))
        #startu, endu = surface_U(spanlist, start_l, sum(spanlist) * len_sub)
        if sub_cfg["c"]:
            name_to_value = {"EtaFlag": 1, "EtaStart": EtaStart,"EtaEnd": EtaEnd, "Abs_Rel_Flag":1, 
                            "Length_C_Start": len_start, "Length_C_End": len_end,
                            "SE_Const_Flag": project_config.VSPAERO_CONTROL_SURFACE_SE_CONST_FLAG } #设置起始和末尾相对长度
        else:
            name_to_value = {"EtaFlag": 1, "EtaStart": EtaStart,"EtaEnd": EtaEnd, "Abs_Rel_Flag":0, 
                            "Length_Start": len_start, "Length_End": len_end,
                            "SE_Const_Flag": project_config.VSPAERO_CONTROL_SURFACE_SE_CONST_FLAG } #设置起始和末尾绝对长度
        for pid in parm_id_vec:
            name = vsp.GetParmName(pid)
            if name in name_to_value:
                vsp.SetParmVal(pid, name_to_value[name])
        vsp.Update()

        # 创建一个空的 VSPAERO 控制面组
        group_index = vsp.CreateVSPAEROControlSurfaceGroup()
        available_cs = vsp.GetAvailableCSNameVec( group_index )
        #print("Available control surfaces:", available_cs)
        # 构造前缀（在名称中通常为 "wing_name_"）
        prefix = wing_name + "_"
        # 遍历并收集所有以 prefix 开头的控制面索引
        indices = []
        # enumerate 从 0 开始计数，我们后面再 +1 转为 1-based
        for idx, cs_name in enumerate(available_cs):
            # 判断名称是否以指定前缀开头
            if cs_name.startswith(prefix):
                # 把对应的 1-based 索引加入结果列表
                indices.append(idx + 1)

        if len(indices) != 1:
            raise RuntimeError(
                f"Expected one VSPAERO control surface for {wing_name}, "
                f"found {indices}; available surfaces: {available_cs}"
            )
        # 同时收集这些控制面的名字
        added_cs_names = []
        for i in indices:
            # 记得 indices 是 1-based，所以要减 1 取列表元素
            added_cs_names.append(available_cs[i - 1])

        vsp.Update()
        # 给这个组命名
        cs_group_name = sub_cfg["name"]
        vsp.SetVSPAEROControlGroupName(cs_group_name, group_index)
        vsp.AddSelectedToCSGroup(indices, group_index)
        print(f"added {', '.join(added_cs_names)} to {cs_group_name}")
    vsp.Update()
    return Npanel

def place_prop(span, liftprop_Dia, tipprop_Dia, prop_ele, tess_int, fuse_w):
    clearance = 0.004
    # 计算螺旋桨数量
    Nprops = math.floor((span - fuse_w - 0.5 * tipprop_Dia * inch_in_m) /(liftprop_Dia * inch_in_m + clearance)) + 1
    # 计算布置间隔
    space = ((Nprops - 1) * liftprop_Dia + 0.5 * tipprop_Dia) * inch_in_m
    gap = (span - fuse_w - space) / Nprops
    # 坐标：均匀分布的升力桨 + 翼尖桨
    prop_pos = np.append(
        np.linspace(gap + 0.5 * liftprop_Dia * inch_in_m + fuse_w,
                    span - 0.5 * (tipprop_Dia + liftprop_Dia) * inch_in_m - gap,
                    Nprops - 1), span)
    
    # 直径：前面是升力桨，最后是翼尖桨
    prop_D = [liftprop_Dia] * (Nprops - 1) + [tipprop_Dia]

    prop_ele = [prop_ele] * (len(prop_pos)-1)
    prop_ele.append(0.002)

    #Add prop
    prop_id = []
    for i in range(0, len(prop_pos)):
        prop_id.append(vsp.AddGeom( "PROP", "" ))
        vsp.SetSetFlag(prop_id[-1], VSPAERO_PROP_GEOM_SET, True)
        sym_parm_prop = vsp.FindParm(prop_id[i], "Sym_Planar_Flag", "Sym")  
        vsp.SetParmValUpdate(sym_parm_prop, 0) #0 for none, 1 for XY, 2 for XZ
        vsp.SetParmVal( prop_id[i], "PropMode", "Design", vsp.PROP_DISK )
        vsp.SetParmVal( prop_id[i], "Diameter", "Design", prop_D[i] * inch_in_m )
        vsp.SetParmVal( prop_id[i], "X_Rel_Location", "XForm", -0.03)
        vsp.SetParmVal( prop_id[i], "Y_Rel_Location", "XForm", prop_pos[i] )
        vsp.SetParmVal( prop_id[i], "Z_Rel_Location", "XForm", prop_ele[i] )
        tessU = int(prop_D[i] * inch_in_m / tess_int * 0.9)
        vsp.SetParmVal( prop_id[i], "Tess_U", "Shape", tessU )
        vsp.SetParmVal( prop_id[i], "Tess_W", "Shape", tessU )
    vsp.Update()
    return prop_pos, prop_D, Nprops

def place_single_prop(prop_Dia, prop_pos, tess_int,  prop_ele = 0):
    prop_Dia = prop_Dia * inch_in_m
    prop_id = vsp.AddGeom( "PROP", "" )
    vsp.SetSetFlag(prop_id, VSPAERO_PROP_GEOM_SET, True)
    sym_parm_prop = vsp.FindParm(prop_id, "Sym_Planar_Flag", "Sym")  
    vsp.SetParmValUpdate(sym_parm_prop, 0) #0 for none, 1 for XY, 2 for XZ
    vsp.SetParmVal( prop_id, "PropMode", "Design", vsp.PROP_DISK )
    vsp.SetParmVal( prop_id, "Diameter", "Design", prop_Dia)
    vsp.SetParmVal( prop_id, "X_Rel_Location", "XForm", -0.03)
    vsp.SetParmVal( prop_id, "Y_Rel_Location", "XForm", prop_pos )
    vsp.SetParmVal( prop_id, "Z_Rel_Location", "XForm", prop_ele )
    tessU = int(prop_Dia / tess_int * 0.7)
    vsp.SetParmVal( prop_id, "Tess_U", "Shape", tessU )
    vsp.SetParmVal( prop_id, "Tess_W", "Shape", tessU )
    vsp.Update()

    return [prop_Dia]

def D0(cd0, spd, s):
    d0 = 0.5 * density * spd**2 * cd0 * s
    return d0

def solve_small_linear_system(matrix, rhs):
    matrix = np.asarray(matrix, dtype=float)
    rhs = np.asarray(rhs, dtype=float)
    n = len(rhs)
    if matrix.shape != (n, n):
        raise ValueError(f"Expected a square matrix with shape {(n, n)}, got {matrix.shape}")

    a = matrix.tolist()
    b = rhs.tolist()
    for column in range(n):
        pivot_row = max(range(column, n), key=lambda row: abs(a[row][column]))
        pivot = a[pivot_row][column]
        if abs(pivot) < 1.0e-12:
            raise ValueError(f"Singular Broyden Jacobian at column {column}")
        if pivot_row != column:
            a[column], a[pivot_row] = a[pivot_row], a[column]
            b[column], b[pivot_row] = b[pivot_row], b[column]

        for row in range(column + 1, n):
            factor = a[row][column] / a[column][column]
            a[row][column] = 0.0
            for j in range(column + 1, n):
                a[row][j] -= factor * a[column][j]
            b[row] -= factor * b[column]

    solution = [0.0] * n
    for row in range(n - 1, -1, -1):
        solution[row] = (
            b[row] - sum(a[row][j] * solution[j] for j in range(row + 1, n))
        ) / a[row][row]
    return np.asarray(solution, dtype=float)

def _set_vspaero_geometry_sets(analysis_name, thick_geom_set, thin_geom_set):
    vsp.SetIntAnalysisInput(analysis_name, "GeomSet", [thick_geom_set], 0)
    vsp.SetIntAnalysisInput(analysis_name, "ThinGeomSet", [thin_geom_set], 0)

def _read_vspgeom_mesh_quality(geom_file):
    lines = open(geom_file, "r", encoding="utf-8", errors="ignore").read().splitlines()
    face_error_count = sum(line.strip() == "# FACE ERROR" for line in lines)
    point_count, panel_count, symmetry = map(int, lines[2].split())
    points = np.asarray(
        [[float(value) for value in lines[3 + index].split()] for index in range(point_count)],
        dtype=float,
    )
    invalid_panel_count = 0
    zero_edge_panel_count = 0
    max_finite_aspect = 0.0
    min_edge = float("inf")
    for index in range(panel_count):
        fields = lines[4 + point_count + index].split()
        node_numbers = [int(value) for value in fields[1:5]]
        if any(number < 1 or number > point_count for number in node_numbers):
            invalid_panel_count += 1
            continue
        coordinates = points[np.asarray(node_numbers) - 1]
        edges = np.linalg.norm(
            np.roll(coordinates, -1, axis=0) - coordinates,
            axis=1,
        )
        current_min_edge = float(edges.min())
        min_edge = min(min_edge, current_min_edge)
        if current_min_edge <= project_config.VSPAERO_MESH_QUALITY_MIN_EDGE:
            zero_edge_panel_count += 1
            continue
        max_finite_aspect = max(max_finite_aspect, float(edges.max() / current_min_edge))
    return {
        "point_count": point_count,
        "panel_count": panel_count,
        "symmetry": symmetry,
        "face_error_count": face_error_count,
        "invalid_panel_count": invalid_panel_count,
        "zero_edge_panel_count": zero_edge_panel_count,
        "min_edge": min_edge,
        "max_finite_aspect": max_finite_aspect,
    }

def _validate_vspgeom_mesh(geom_file):
    quality = _read_vspgeom_mesh_quality(geom_file)
    print("VSPAERO mesh quality:", quality)
    if quality["face_error_count"] > 0:
        raise RuntimeError(f"VSPAERO geometry contains FACE ERROR records: {geom_file}")
    if quality["invalid_panel_count"] > 0:
        raise RuntimeError(f"VSPAERO geometry contains invalid panel nodes: {geom_file}")
    if quality["zero_edge_panel_count"] > 0:
        raise RuntimeError(f"VSPAERO geometry contains zero-length panel edges: {geom_file}")
    if quality["max_finite_aspect"] > project_config.VSPAERO_MESH_QUALITY_MAX_ASPECT:
        raise RuntimeError(
            f"VSPAERO geometry contains extreme sliver panels: "
            f"aspect={quality['max_finite_aspect']:.6g}, file={geom_file}"
        )
    return quality

def _run_vspaero_compute_geometry(thick_geom_set, thin_geom_set):
    compgeom_name = "VSPAEROComputeGeometry"
    print(compgeom_name)
    vsp.SetAnalysisInputDefaults(compgeom_name)
    vsp.SetIntAnalysisInput(compgeom_name, "Symmetry", [2], 0)
    _set_vspaero_geometry_sets(compgeom_name, thick_geom_set, thin_geom_set)
    vsp.PrintAnalysisInputs(compgeom_name)
    vsp.WriteVSPFile(file_name, vsp.SET_ALL)
    print("\tExecuting...")
    compgeom_resid = vsp.ExecAnalysis(compgeom_name)
    print("COMPLETE")
    vsp.PrintResults(compgeom_resid)
    try:
        geom_file = vsp.GetStringResults(compgeom_resid, "VSPGeomFileName")[0]
    except Exception:
        geom_file = vsp.GetStringResults(compgeom_resid, "DegenGeomFileName")[0]
    print("Generated VSPAERO geometry file:", geom_file)
    _validate_vspgeom_mesh(geom_file)
    return compgeom_resid

def _apply_actuator_disk_settings(prop_D, RPM, Ct, Cp):
    Nprops = len(prop_D)
    if Nprops == 0:
        return
    num_disks = vsp.GetNumActuatorDisks()
    if num_disks != Nprops:
        raise ValueError(f"Expected {Nprops} actuator disks, found {num_disks}")
    for i in range(Nprops):
        disk_id = vsp.FindActuatorDisk(i)
        disk_parm_ids = vsp.FindContainerParmIDs(disk_id)
        rpm_id = None
        ct_id = None
        cp_id = None
        for parm_id in disk_parm_ids:
            parm_name = vsp.GetParmName(parm_id)
            if parm_name == "RotorRPM":
                rpm_id = parm_id
            elif parm_name == "RotorCT":
                ct_id = parm_id
            elif parm_name == "RotorCP":
                cp_id = parm_id
        if rpm_id is None or ct_id is None or cp_id is None:
            raise ValueError(f"Actuator disk {i} is missing RotorRPM/RotorCT/RotorCP parameters")
        vsp.SetParmVal(rpm_id, RPM[i])
        vsp.SetParmVal(ct_id, Ct[i])
        vsp.SetParmVal(cp_id, Cp[i])

def _set_prop_shown(has_prop):
    for geom in vsp.FindGeoms():
        type_name = vsp.GetGeomTypeName(geom).upper()
        if type_name == "PROPELLER" or type_name == "PROP":
            vsp.SetSetFlag(geom, vsp.SET_SHOWN, has_prop)
    vsp.Update()

def _apply_control_surface_angles(angle):
    if not angle:
        return
    cs_group_container_id = vsp.FindContainer("VSPAEROSettings", 0)
    Num_cs = vsp.GetNumControlSurfaceGroups()
    for i in range(Num_cs):
        cs_name = vsp.GetVSPAEROControlGroupName(i)
        if cs_name in angle:
            grp_name = f"ControlSurfaceGroup_{i}"
            defl_parm = vsp.FindParm(cs_group_container_id, "DeflectionAngle", grp_name)
            vsp.SetParmValUpdate(defl_parm, angle[cs_name])

def _run_vspaero_sweep(CG, AlphaStart_input, AlphaEnd_input, AlphaNpts_input,
                       air_spd, wing_cfg, sol_config, angle,
                       thick_geom_set, thin_geom_set, has_propellers,
                       reynolds=None, wake_num_iter=None, prepare_model=True):
    wing_S = wing_cfg["wing_S"]
    bref = wing_cfg["bref"]
    cref = wing_cfg["cref"]
    if reynolds is None:
        Re = cref * air_spd * density / mu
    else:
        if not math.isfinite(reynolds) or reynolds <= 0.0:
            raise ValueError(f"reynolds must be positive and finite, got {reynolds}")
        Re = reynolds

    analysis_name = "VSPAEROSweep"
    print(analysis_name)
    vsp.SetAnalysisInputDefaults(analysis_name)
    _set_vspaero_geometry_sets(analysis_name, thick_geom_set, thin_geom_set)
    vsp.SetDoubleAnalysisInput(analysis_name, "AlphaEnd", [AlphaEnd_input], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "AlphaStart", [AlphaStart_input], 0)
    vsp.SetIntAnalysisInput(analysis_name, "AlphaNpts", [AlphaNpts_input], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "Xcg", [CG], 0)
    vsp.SetIntAnalysisInput(analysis_name, "NCPU", [CPU], 0)
    vsp.SetIntAnalysisInput(analysis_name, "Symmetry", [2], 0)
    if has_propellers:
        vsp.SetIntAnalysisInput(analysis_name, "PropBladesMode", [0], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "Sref", [wing_S], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "bref", [bref], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "cref", [cref], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "ReCref", [Re], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "Vinf", [air_spd], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "Vref", [air_spd], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "Rho", [density], 0)
    wake_num_iter = 5 if wake_num_iter is None else wake_num_iter
    if not isinstance(wake_num_iter, int) or wake_num_iter <= 0:
        raise ValueError(f"wake_num_iter must be a positive integer, got {wake_num_iter}")
    vsp.SetIntAnalysisInput(analysis_name, "WakeNumIter", [wake_num_iter], 0)
    vsp.SetIntAnalysisInput(analysis_name, "FarDistToggle", [1], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "FarDist", [sol_config["farfield"]], 0)
    vsp.SetIntAnalysisInput(analysis_name, "NumWakeNodes", [sol_config["wakenode"]], 0)
    if prepare_model:
        _apply_control_surface_angles(angle)
        vsp.Update()
    vsp.PrintAnalysisInputs(analysis_name)
    print("\n[INFO] 开始执行 VSPAEROSweep 分析...")
    res_id = vsp.ExecAnalysis(analysis_name)
    print("[INFO] 分析完成")
    return res_id

def _find_result_file(suffixes):
    for suffix in suffixes:
        path = case_name + suffix
        if os.path.exists(path):
            return path
    raise FileNotFoundError(f"Cannot find result file for case {case_name}: {suffixes}")

def _synchronize_control_surface_taglist():
    csf_path = case_name + ".csf"
    taglist_path = case_name + ".ControlSurfaces.taglist"
    if not os.path.exists(csf_path) or not os.path.exists(taglist_path):
        return

    with open(taglist_path, "r", encoding="utf-8", errors="ignore") as taglist_file:
        taglist_lines = taglist_file.read().splitlines()
    tag_count = int(taglist_lines[0].strip())
    tag_names = [line.strip() for line in taglist_lines[1:] if line.strip()]
    if len(tag_names) != tag_count:
        raise RuntimeError(f"Invalid control-surface taglist {taglist_path}")
    for tag_name in tag_names:
        if not os.path.exists(tag_name + ".tag"):
            raise RuntimeError(f"Missing control-surface tag file {tag_name}.tag")
    if not tag_names:
        raise RuntimeError(f"No control-surface tag files found for {case_name}")

    with open(csf_path, "r", encoding="utf-8", errors="ignore") as csf_file:
        csf_lines = csf_file.read().splitlines()
    vspaero_name_indices = [
        index for index, line in enumerate(csf_lines)
        if line.strip().startswith("VSPAERO Name:")
    ]
    tagfile_name_indices = [
        index for index, line in enumerate(csf_lines)
        if line.strip().startswith("Tagfile Name:")
    ]
    if len(vspaero_name_indices) != tag_count or len(tagfile_name_indices) != tag_count:
        raise RuntimeError(f"Control-surface count mismatch in {csf_path}")

    new_vspaero_names = [tag_name[len(case_name):] for tag_name in tag_names]

    for index, vspaero_name in zip(vspaero_name_indices, new_vspaero_names):
        csf_lines[index] = "VSPAERO Name: " + vspaero_name
    for index, tag_name in zip(tagfile_name_indices, tag_names):
        csf_lines[index] = "Tagfile Name: " + tag_name
    with open(csf_path, "w", encoding="utf-8", newline="\n") as csf_file:
        csf_file.write("\n".join(csf_lines) + "\n")

    vspaero_path = case_name + ".vspaero"
    if not os.path.exists(vspaero_path):
        return
    with open(vspaero_path, "r", encoding="utf-8", errors="ignore") as vspaero_file:
        vspaero_lines = vspaero_file.read().splitlines()
    vspaero_surface_indices = [
        index for index, line in enumerate(vspaero_lines)
        if "_Surf" in line and "_SS_CONT_" in line
    ]
    if len(vspaero_surface_indices) != tag_count:
        raise RuntimeError(f"Control-surface reference count mismatch in {vspaero_path}")
    for index, vspaero_name in zip(vspaero_surface_indices, new_vspaero_names):
        vspaero_lines[index] = vspaero_name
    with open(vspaero_path, "w", encoding="utf-8", newline="\n") as vspaero_file:
        vspaero_file.write("\n".join(vspaero_lines) + "\n")


def _read_vspaero_polar(drag_column="CDi"):
    polar_name = _find_result_file((".polar", "_DegenGeom.polar"))
    with open(polar_name, "r", encoding="utf-8", errors="ignore") as polar_file:
        lines = polar_file.read().splitlines()
    header_index = next(
        index for index, line in enumerate(lines)
        if line.strip().startswith("Beta ")
    )
    columns = lines[header_index].split()
    rows = [
        line.split()
        for line in lines[header_index + 1:]
        if len(line.split()) == len(columns)
    ]
    if not rows:
        raise ValueError(f"No polar rows are present in {polar_name}")
    df = pd.DataFrame(rows, columns=columns).apply(pd.to_numeric, errors="coerce")
    cl_column = "CLtot" if "CLtot" in df.columns else "CL"
    cm_column = "CMytot" if "CMytot" in df.columns else "CMy"
    if drag_column not in df.columns:
        if drag_column == "CDi" and "CDtot" in df.columns:
            drag_column = "CDtot"
        else:
            raise ValueError(f"Column {drag_column} is not present in {polar_name}")
    return df[cl_column].tolist(), df[drag_column].tolist(), df[cm_column].tolist()

def solve_stall_limit_from_lod(
    lod_path,
    air_spd,
    airfoil_cfg,
    limit_config=None,
    wing_cfg=None,
    flap_cfg=None,
    flap_deflection=0.0,
):
    """Solve the sectional stall limit from a completed multi-case LOD file."""
    limit_config = project_config.STALL_LIMIT_CONFIG if limit_config is None else limit_config
    cases = stall_limit.read_lod_cases(lod_path)
    if len(cases) < 2:
        raise RuntimeError(
            f"LOD contains {len(cases)} case; expected multiple alpha cases: {lod_path}"
        )
    if flap_deflection != 0.0 and (wing_cfg is None or flap_cfg is None):
        raise ValueError("Wing and flap geometry are required for a deflected-flap stall limit")

    def capacity_predictor(case):
        rows = stall_limit.lod_case_main_wing_rows(case)
        reynolds = stall_limit.lod_case_reynolds(
            case, air_spd, density=density, viscosity=mu
        )
        flap_angles = np.zeros(len(rows), dtype=float)
        hinge_points = np.ones(len(rows), dtype=float)
        if flap_deflection != 0.0:
            semispan = wing_cfg["bref"] / 2.0
            if semispan <= 0.0:
                raise ValueError("Wing semispan must be positive for flap capacity mapping")
            eta_start = flap_cfg["EtaStart"]
            eta_end = flap_cfg["EtaEnd"]
            length_start = flap_cfg["Length_Start"]
            length_end = flap_cfg["Length_End"]
            eta_min = min(eta_start, eta_end)
            eta_max = max(eta_start, eta_end)
            for index, row in enumerate(rows):
                eta = row["Yavg"] / semispan
                if eta_min <= eta <= eta_max:
                    span_fraction = (eta - eta_start) / (eta_end - eta_start)
                    flap_length = length_start + span_fraction * (length_end - length_start)
                    hinge_points[index] = 1.0 - flap_length
                    flap_angles[index] = -flap_deflection
        return stall_limit.predict_strip_clmax(
            airfoil_cfg,
            reynolds,
            limit_config=limit_config,
            flap_deflection=flap_angles,
            hinge_point=hinge_points,
        )

    result = stall_limit.solve_lod_stall_limit(
        cases,
        capacity_predictor,
        limit_config=limit_config,
    )
    result["lod_path"] = os.path.abspath(lod_path)
    result["case_count"] = len(cases)
    return result

def run_stall_limit_sweep(CG, air_spd, wing_cfg, airfoil_cfg,
                          angle=None, prop_D=None, RPM=None, Ct=None, Cp=None,
                          sol_config=None, limit_config=None,
                          isolated=True, flap_cfg=None, flap_deflection=0.0):
    """Run a real VSPAERO alpha sweep and solve the sectional stall limit."""
    limit_config = project_config.STALL_LIMIT_CONFIG if limit_config is None else limit_config
    sol_config = solver_config0 if sol_config is None else sol_config
    angle = {} if angle is None else angle
    prop_D = [] if prop_D is None else prop_D
    RPM = [] if RPM is None else RPM
    Ct = [] if Ct is None else Ct
    Cp = [] if Cp is None else Cp

    alpha_min = limit_config["alpha_min"]
    alpha_max = limit_config["alpha_max"]
    alpha_points = limit_config["alpha_points"]
    if alpha_points < 2:
        raise ValueError("STALL_LIMIT_CONFIG alpha_points must be at least 2")
    if alpha_max <= alpha_min:
        raise ValueError("STALL_LIMIT_CONFIG alpha_max must exceed alpha_min")

    aero_function = _runaero_isolated if isolated else runaero
    aero_function(
        CG, alpha_min, alpha_max, alpha_points, air_spd, wing_cfg, 0.0,
        sol_config, angle, prop_D, RPM, Ct, Cp,
    )

    lod_path = _find_result_file((".lod", "_DegenGeom.lod"))
    result = solve_stall_limit_from_lod(
        lod_path,
        air_spd,
        airfoil_cfg,
        limit_config=limit_config,
        wing_cfg=wing_cfg,
        flap_cfg=flap_cfg,
        flap_deflection=flap_deflection,
    )
    result["alpha_min"] = float(alpha_min)
    result["alpha_max"] = float(alpha_max)
    result["alpha_points"] = int(alpha_points)
    return result

def _runaero_isolated(CG, AlphaStart_input, AlphaEnd_input, AlphaNpts_input,
                      air_spd, wing_cfg, Cl_target, sol_config, angle,
                      prop_D, RPM, Ct, Cp, reynolds=None, wake_num_iter=None):
    source_model_path = os.path.abspath(file_name)
    worker_dir = os.path.join(os.getcwd(), f"vsp_worker_{uuid.uuid4().hex}")
    os.makedirs(worker_dir)
    worker_case_name = case_name
    worker_model_path = os.path.join(worker_dir, worker_case_name + ".vsp3")
    shutil.copy2(source_model_path, worker_model_path)
    request_path = os.path.join(worker_dir, "request.json")
    response_path = os.path.join(worker_dir, "response.json")
    log_path = os.path.join(worker_dir, "worker.log")

    request = {
        "case_name": worker_case_name,
        "file_name": worker_case_name + ".vsp3",
        "model_file": worker_model_path,
        "workdir": worker_dir,
        "CG": CG,
        "AlphaStart_input": AlphaStart_input,
        "AlphaEnd_input": AlphaEnd_input,
        "AlphaNpts_input": AlphaNpts_input,
        "air_spd": air_spd,
        "wing_cfg": wing_cfg,
        "Cl_target": Cl_target,
        "sol_config": sol_config,
        "angle": angle,
        "prop_D": np.asarray(prop_D, dtype=float).tolist(),
        "RPM": np.asarray(RPM, dtype=float).tolist(),
        "Ct": np.asarray(Ct, dtype=float).tolist(),
        "Cp": np.asarray(Cp, dtype=float).tolist(),
        "reynolds": reynolds,
        "wake_num_iter": wake_num_iter,
        "density": density,
    }
    with open(request_path, "w", encoding="utf-8") as request_file:
        json.dump(request, request_file, ensure_ascii=False)

    worker_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "infrastructure_worker.py")
    print(f"[VSPAERO worker] start case={case_name} dir={worker_dir}", flush=True)
    with open(log_path, "w", encoding="utf-8") as worker_log:
        startup_info = subprocess.STARTUPINFO()
        startup_info.dwFlags |= subprocess.STARTF_USESHOWWINDOW
        startup_info.wShowWindow = 0
        worker_process = subprocess.Popen(
            [sys.executable, worker_path, request_path, response_path],
            cwd=worker_dir,
            stdout=worker_log,
            stderr=subprocess.STDOUT,
            startupinfo=startup_info,
            creationflags=subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.CREATE_NO_WINDOW,
        )
        try:
            returncode = worker_process.wait(timeout=VSPAERO_WORKER_TIMEOUT_SECONDS)
        except subprocess.TimeoutExpired as error:
            subprocess.run(
                ["taskkill", "/PID", str(worker_process.pid), "/T", "/F"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=False,
            )
            raise TimeoutError(
                f"isolated VSPAERO worker exceeded {VSPAERO_WORKER_TIMEOUT_SECONDS} seconds; "
                f"worker directory: {worker_dir}; log: {log_path}"
            ) from error
    print(f"[VSPAERO worker] returncode={returncode} case={case_name}", flush=True)
    if returncode != 0:
        raise RuntimeError(
            f"isolated VSPAERO worker failed with exit code {returncode:#x}; "
            f"worker directory: {worker_dir}; log: {log_path}"
        )

    if not os.path.exists(response_path):
        raise RuntimeError(
            f"isolated VSPAERO worker exited successfully but produced no response; "
            f"worker directory: {worker_dir}; log: {log_path}"
        )
    with open(response_path, "r", encoding="utf-8") as response_file:
        response = json.load(response_file)
    for worker_filename in os.listdir(worker_dir):
        if worker_filename.startswith(worker_case_name):
            worker_result_path = os.path.join(worker_dir, worker_filename)
            if os.path.isfile(worker_result_path):
                parent_result_path = os.path.join(os.getcwd(), worker_filename)
                shutil.copy2(worker_result_path, parent_result_path)
    shutil.rmtree(worker_dir)
    return (
        response["drag"], response["alpha"], response["lift"], response["netdrag"],
        response["power"], response["Cl_list"], response["CMy"]
    )

def runaero(CG, AlphaStart_input, AlphaEnd_input, AlphaNpts_input, air_spd, wing_cfg, Cl_target, sol_config, angle=[], prop_D=[],RPM =[],Ct =[],Cp =[], reynolds=None, wake_num_iter=None):
    global vspaero_geometry_signature
    wing_S = wing_cfg["wing_S"]
    Nprops = len(prop_D)
    Re = wing_cfg["cref"] * air_spd * density / mu if reynolds is None else reynolds
    print("Re: " + str(Re))
    vsp.DeleteAllResults()
    _set_prop_shown(Nprops > 0)
    _apply_control_surface_angles(angle)
    vsp.Update()
    thin_geom_set = VSPAERO_PROP_THIN_GEOM_SET if Nprops > 0 else VSPAERO_NO_PROP_THIN_GEOM_SET
    geometry_signature = (Nprops > 0, Nprops, thin_geom_set)
    if geometry_signature != vspaero_geometry_signature:
        if Nprops > 0:
            _apply_actuator_disk_settings(prop_D, RPM, Ct, Cp)
        _run_vspaero_compute_geometry(VSPAERO_NO_GEOM_SET, thin_geom_set)
        _synchronize_control_surface_taglist()
        vspaero_geometry_signature = geometry_signature
    _run_vspaero_sweep(CG, AlphaStart_input, AlphaEnd_input, AlphaNpts_input,
                       air_spd, wing_cfg, sol_config, angle,
                       VSPAERO_NO_GEOM_SET, thin_geom_set,
                       Nprops > 0, reynolds, wake_num_iter, prepare_model=False)
    Cl_list, Cd_list, CMy_list = _read_vspaero_polar("CDi")

    thrust = 0
    power = 0
    for i in range(Nprops):
        thrust += prop_D[i] ** 4 * (RPM[i] / 60) ** 2 * Ct[i] * density * 2
        power += prop_D[i] ** 5 * (RPM[i] / 60) ** 3 * Cp[i] * density * 2

    if AlphaNpts_input >= 2:
        alpha_list = np.linspace(AlphaStart_input, AlphaEnd_input, AlphaNpts_input)
        alpha = np.interp(Cl_target, Cl_list, alpha_list)
        Cd_target = np.interp(Cl_target, Cl_list, Cd_list)
        CMy = np.interp(Cl_target, Cl_list, CMy_list)
        drag = 0.5 * Cd_target * density * wing_S * air_spd ** 2
        lift = 0.5 * Cl_target * density * wing_S * air_spd ** 2 + thrust * np.sin(alpha / 360 * 2 * np.pi)
    else:
        alpha = AlphaEnd_input
        drag = 0.5 * Cd_list[0] * density * wing_S * air_spd ** 2
        lift = 0.5 * Cl_list[0] * density * wing_S * air_spd ** 2 + thrust * np.sin(alpha / 360 * 2 * np.pi)
        CMy = CMy_list[0]
    net_drag = drag - thrust * np.cos(alpha / 360 * 2 * np.pi)
    return drag, alpha, lift, net_drag, power, Cl_list, CMy

def cal_cg(Kn, cfg, AOA, typ_speed):
    Mean_chord = cfg["cref"]
    CG1 = Mean_chord / 2 
    _, _, lift1, _, _, _, CMy0 = runaero(CG1, 0, 0, 1, typ_speed, cfg, solver_config0)
    _, _, lift2, _, _, _, CMy1 = runaero(CG1, AOA, AOA, 1, typ_speed, cfg, 0, solver_config0)
    dMda1 = (CMy1 - CMy0) / AOA

    CG2 = 0
    _, _, lift1, _, _, _, CMy0 = runaero(CG2, 0, 0, 1, typ_speed, cfg, 0, solver_config0)
    _, _, lift2, _, _, _, CMy1 = runaero(CG2, AOA, AOA, 1, typ_speed, cfg, 0, solver_config0)
    dMda2 = (CMy1 - CMy0) / AOA
    XNP = CG1 - (dMda1 * (CG1 - CG2))/(dMda1 - dMda2)

    CG = XNP - Kn * Mean_chord
    return CG

def make_broyden_config(prop_data, max_AOA, weight):
    return {
        "max_it": 4,
        "max_aero_calls": 6,
        "moment_scale": 0.1,
        "variable_scale": [1.0, weight, 5.0],
        "max_steps": [2.0, 0.8 * weight, 10.0],
        "alpha_bounds": [-1.0, max_AOA],
        "thrust_bounds": [0.0, 2.0],
        "elevator_bounds": [-25.0, 25.0],
        "lift_rel_tol": 0.01,
        "drag_rel_tol": 0.001,
        "moment_abs_tol": 1.0e-4,
        "final_tolerance_factor": 3.0,
        "broyden_damping": 0.85,
        "physics_diagonal": [0.10, -1.0, 0.8],
        "isolate_vspaero": True,
        "propdata": prop_data,
    }

def single_point(f_cond, geo_info, config):
    CG = geo_info["CG"]
    cfg = {"wing_S": geo_info["wing_S"], "bref": geo_info["bref"], "cref": geo_info["cref"]}
    spanlist = geo_info["spanlist"]
    span = geo_info["span"]
    chordlist = geo_info["chordlist"]
    wing_S = geo_info["wing_S"]
    def_cfg = geo_info["def_cfg"]
    prop_D = geo_info["prop_D"]
    prop_D_inch = geo_info["prop_D_inch"]
    prop_pos = geo_info["prop_pos"]
    global mass

    max_AOA = f_cond["max_AOA"]
    speed = f_cond["speed"]
    thrust_ratio = f_cond["TR"]
    d0_others = f_cond["d0_others"] if "d0_others" in f_cond else f_cond["d0"]
    t_over_c = geo_info["t_over_c"] if "t_over_c" in geo_info else None
    Cl_target = f_cond["Cl_target"]
    prop_data = config["propdata"]

    max_it = config["max_it"]
    max_aero_calls = config["max_aero_calls"]
    moment_scale = config["moment_scale"]
    alpha_scale = config["variable_scale"][0]
    elevator_scale = config["variable_scale"][2]
    max_steps = np.asarray(config["max_steps"], dtype=float)
    damping = config["broyden_damping"]
    physics_diagonal = np.asarray(config["physics_diagonal"], dtype=float)
    debug_log = []

    initial_elevator = def_cfg["elevator"]
    vsp.WriteVSPFile(file_name, vsp.SET_ALL)
    def_cfg["elevator"] = 0.0
    initial_aero_function = _runaero_isolated if config["isolate_vspaero"] else runaero
    if Cl_target == -1:
        drag, alpha, lift, netdrag, power, _, CMy = initial_aero_function(
            CG, max_AOA, max_AOA, 1, speed, cfg, 0, solver_config0, def_cfg,
            [], [], [], []
        )
    else:
        drag, alpha, lift, netdrag, power, _, CMy = initial_aero_function(
            CG, 0, max_AOA, 2, speed, cfg, Cl_target, solver_config0, def_cfg,
            [], [], [], []
        )

    d0_wing = calculate_parasite_drag_from_lod(speed, t_over_c)
    netdrag = netdrag + d0_others + d0_wing

    y_cp = read_halfwing_cp()
    mass, wing_mass = mass_sim_iter(span, wing_S, y_cp, prop_pos, spanlist, chordlist)
    weight = mass * g

    alpha_min = (
        min(config["alpha_bounds"][0], max_AOA)
        if Cl_target == -1
        else config["alpha_bounds"][0]
    )
    alpha_max = min(config["alpha_bounds"][1], max_AOA)
    thrust_min = config["thrust_bounds"][0]
    thrust_max = config["thrust_bounds"][1] * weight
    elevator_min = config["elevator_bounds"][0]
    elevator_max = config["elevator_bounds"][1]
    state_scale = np.array([alpha_scale, weight, elevator_scale], dtype=float)
    state_min = np.array([alpha_min, thrust_min, elevator_min], dtype=float)
    state_max = np.array([alpha_max, thrust_max, elevator_max], dtype=float)
    residual_scale = np.array([weight, weight, moment_scale], dtype=float)

    if Cl_target == -1:
        active_variables = np.array([1, 2], dtype=int)
        active_residuals = np.array([1, 2], dtype=int)
    else:
        active_variables = np.array([0, 1, 2], dtype=int)
        active_residuals = np.array([0, 1, 2], dtype=int)

    residual_tolerance = np.array([
        config["drag_rel_tol"],
        config["moment_abs_tol"] / moment_scale,
    ]) if Cl_target == -1 else np.array([
        config["lift_rel_tol"],
        config["drag_rel_tol"],
        config["moment_abs_tol"] / moment_scale,
    ])

    def project_state(z):
        state = np.asarray(z, dtype=float) * state_scale
        state = np.clip(state, state_min, state_max)
        return state / state_scale

    def evaluate_state(z, solver_config):
        z = project_state(z)
        alpha_state, thrust_state, elevator_state = z * state_scale
        def_cfg["elevator"] = elevator_state
        RPM, Ct, Cp = prop.equal_thrust(
            prop_data, thrust_state, speed, prop_D_inch, thrust_ratio
        )
        aero_function = _runaero_isolated if config["isolate_vspaero"] else runaero
        drag_state, _, lift_state, netdrag_state, power_state, _, CMy_state = aero_function(
            CG, alpha_state, alpha_state, 1, speed, cfg, 0, solver_config,
            def_cfg, prop_D, RPM, Ct, Cp
        )
        d0_wing_state = calculate_parasite_drag_from_lod(speed, t_over_c)
        netdrag_state = netdrag_state + d0_others + d0_wing_state
        residual_state = np.array([
            lift_state - weight,
            netdrag_state,
            CMy_state,
        ]) / residual_scale
        if Cl_target == -1:
            residual_state = residual_state[1:]
        return {
            "alpha": alpha_state,
            "thrust": thrust_state,
            "elevator": elevator_state,
            "RPM": RPM,
            "drag": drag_state,
            "lift": lift_state,
            "netdrag": netdrag_state,
            "power": power_state,
            "CMy": CMy_state,
            "residual": residual_state,
            "D0_wing": d0_wing_state,
        }

    initial_state = np.array([
        alpha,
        max(netdrag, thrust_min),
        initial_elevator,
    ], dtype=float)
    z = project_state(initial_state / state_scale)
    state = evaluate_state(z, solver_config0)
    aero_calls = 2
    max_steps_z = max_steps / state_scale
    if aero_calls > max_aero_calls:
        raise RuntimeError(f"Initial aerodynamic calls exceeded budget {max_aero_calls}")

    jacobian = np.zeros((len(active_variables), len(active_variables)), dtype=float)
    for row, residual_index in enumerate(active_residuals):
        variable_position = np.where(active_variables == residual_index)[0]
        if len(variable_position) != 1:
            raise ValueError(f"No physical variable is available for residual {residual_index}")
        jacobian[row, variable_position[0]] = physics_diagonal[residual_index]
    print(f"[Broyden] physical Jacobian initialized: {jacobian}", flush=True)

    converged = False
    final_solver_used = False
    for iteration in range(1, max_it + 1):
        residual = state["residual"]
        if np.all(np.abs(residual) <= residual_tolerance):
            converged = True
            break

        print(f"[Broyden] solve iteration={iteration} residual={residual}", flush=True)
        delta_active = solve_small_linear_system(jacobian, -residual)
        delta_z = np.zeros(3, dtype=float)
        delta_z[active_variables] = np.clip(
            damping * delta_active,
            -max_steps_z[active_variables],
            max_steps_z[active_variables],
        )

        trial_z = project_state(z + delta_z)
        final_solver_used = iteration == max_it
        trial_state = evaluate_state(trial_z, solver_config1 if final_solver_used else solver_config0)
        aero_calls += 1
        if aero_calls > max_aero_calls:
            raise RuntimeError(f"Aerodynamic call budget exceeded {max_aero_calls}")

        accepted_step = trial_z - z
        jacobian_step = trial_state["residual"] - residual - jacobian @ accepted_step[active_variables]
        denominator = np.dot(accepted_step[active_variables], accepted_step[active_variables])
        if denominator < 1.0e-14:
            raise RuntimeError("Broyden step collapsed to zero")
        jacobian = jacobian + np.outer(jacobian_step, accepted_step[active_variables]) / denominator

        z = trial_z
        state = trial_state
        debug_log.append({
            "Iteration": iteration,
            "AeroCalls": aero_calls,
            "Alpha": state["alpha"],
            "Lift": state["lift"],
            "Drag": state["drag"],
            "Thrust": state["thrust"],
            "Weight": weight,
            "CMy": state["CMy"],
            "Elevator": state["elevator"],
            "NetDrag": state["netdrag"],
            "LiftResidual": (state["lift"] - weight) / weight,
            "DragResidual": state["netdrag"] / weight,
            "MomentResidual": state["CMy"] / moment_scale,
            "ResidualNorm": np.linalg.norm(state["residual"], ord=2),
            "D0_wing": state["D0_wing"],
        })

    final_state = state
    if not final_solver_used:
        final_state = evaluate_state(z, solver_config1)
        aero_calls += 1
        if aero_calls > max_aero_calls:
            raise RuntimeError(f"Aerodynamic call budget exceeded {max_aero_calls}")
        if debug_log:
            debug_log[-1].update({
                "AeroCalls": aero_calls,
                "Alpha": final_state["alpha"],
                "Lift": final_state["lift"],
                "Drag": final_state["drag"],
                "Thrust": final_state["thrust"],
                "CMy": final_state["CMy"],
                "Elevator": final_state["elevator"],
                "NetDrag": final_state["netdrag"],
                "LiftResidual": (final_state["lift"] - weight) / weight,
                "DragResidual": final_state["netdrag"] / weight,
                "MomentResidual": final_state["CMy"] / moment_scale,
                "ResidualNorm": np.linalg.norm(final_state["residual"], ord=2),
                "D0_wing": final_state["D0_wing"],
            })
    final_residual_tolerance = residual_tolerance * config["final_tolerance_factor"]
    final_converged = np.all(np.abs(final_state["residual"]) <= final_residual_tolerance)
    if aero_calls > max_aero_calls:
        raise RuntimeError(f"Aerodynamic call budget exceeded {max_aero_calls}")
    state = final_state
    converged = final_converged

    if converged:
        print(f"Broyden converged in {len(debug_log)} iterations, {aero_calls} aero calls!")
    else:
        print(f"Warning: Broyden did not converge within {max_it} iterations ({aero_calls} aero calls).")

    mass_result = {"mass": mass, "wing_mass": wing_mass}
    filename = f"convergence_{case_name}_V{speed:.1f}.csv"
    pd.DataFrame(debug_log).to_csv(filename, index=False)
    print(f"Convergence log saved to {filename}")
    return (
        state["lift"], state["drag"], state["power"], state["alpha"],
        state["RPM"], state["thrust"], mass_result, state["elevator"]
    )

if __name__ == "__main__":
    print("testing")
    import os
    # 获取当前脚本所在目录
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # 创建 outputs 子文件夹（如果不存在）
    output_dir = os.path.join(script_dir, "outputs")
    os.makedirs(output_dir, exist_ok=True)
    # 切换当前工作目录
    os.chdir(output_dir)

    span = 0.7

    airfoil_cfg = {
    "filename": "EMX-07.dat",
    "Camber" : 0.04,
    "CamberLoc": 0.4,
    "ThickChord": 0.12
    }

    flap_cfg = {
        "name" : "flap",
        "c": 1,
        "Length_Start" : 0.3,
        "Length_End" : 0.25,
        "EtaStart" : 0.7,
        "EtaEnd": 0.05
    }

    ELE_cfg = {"name": "elevator", "c": 0,"Length_Start" : 0.05,"Length_End" : 0.05, "EtaStart" : 0.8, "EtaEnd": 0}

    wing_pos = {"name": "mainwing", "x":0, "y":0, "z":0, "yr": 0}
    wing_misc = {}
    tail_pos = {"name": "tail","x":0.5, "y":0, "z":0, "yr": 0}

    tess_int = 0.005
    cruise_spd = project_config.EVALUATION_TYPICAL_SPEED
    # 1) 生成初始网格
    #spans, chords, twists, wing_S, Cl_target = generate_elliptical_wing(chord_root=0.2, chord_tip=0.05, semispan=0.7, nSecs=7, air_spd = cruise_spd)
    spans, chords, twists, wing_S, Cl_target = think_trapwing(0.2, 0.05, span, 0.05, 0, cruise_spd)
    print("Cl is: " + str(Cl_target))
    print("S: " + str(wing_S))
    ini_geom()
    prop_pos, prop_D_inch, Nprops = place_prop(span, 8, 13,0, tess_int, 0)
    prop_D = np.array(prop_D_inch) * inch_in_m
    # 2) 在螺旋桨滑流区间细化
    spans, chords, twists, Wing_S, Cl_target = refine_wing_mesh(chords, spans, prop_pos, prop_D, [0.2, 0.2, 0.1], tess_int*3, cruise_spd)
    print("new Cl is: " + str(Cl_target))
    print("new S: " + str(Wing_S))

    create_wing(wing_pos, spans, chords, twists, 0, tess_int, airfoil_cfg, flap_cfg)
    vsp.Update()
    create_wing(tail_pos, [0.2], [0.2, 0.1], [0], max_sweeploc, tess_int, airfoil_cfg, ELE_cfg)
    vsp.WriteVSPFile(file_name, vsp.SET_ALL)
    print(f"模型已保存")
