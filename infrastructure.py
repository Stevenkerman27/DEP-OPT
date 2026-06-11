import openvsp as vsp
from scipy.interpolate import CubicSpline
import math
import numpy as np
import pandas as pd
import re, io
import prop
import warnings
# This is for vsp3.41!!!
density = 1.225
mu = 1.85e-5
inch_in_m = 0.0254
max_sweeploc = 0.9954
case_name = "test"
file_name = "test.vsp3"
G = 7
SF = 2
g = 9.8

mass = 3
mass_prop = {"S_density": 2, "CF_Strength": 400e6, "CF_rho": 2000, "fuse_mass": 1.2, "prop": [0.05,0.05,0.05,0.12], "payload": 0.3}
D_min = 5
D_max = 20

CPU = 6
far_1 = 3
far_2 = 10
wakeN_1 = 16
wakeN_2 = 32
solver_config0 = {"farfield":far_1, "wakenode":wakeN_1} # fast setup
solver_config1 = {"farfield":far_2, "wakenode":wakeN_2} # slow setup
def read_halfwing_cp():
    # 读取全文并定位表头
    path = case_name + ".lod"
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.read().splitlines()
    if not lines:
        raise ValueError(f"文件 {path} 为空。")

    header_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith("Iter"):
            header_idx = i
            break
    if header_idx is None:
        raise ValueError("未找到数据表头行（Iter...）。")

    # 从表头后读取主翼部分
    data_lines = []
    for line in lines[header_idx + 1:]:
        if not line.strip():
            continue

        tokens = line.strip().split()
        if len(tokens) < 2:
            continue
        
        # tokens[1] 是 VortexSheet (即 Component)
        if tokens[1] == '2':
            break

        data_lines.append(line)

    if not data_lines:
        raise ValueError("未找到主翼数据区段。")

    df = pd.read_csv(io.StringIO("\n".join([lines[header_idx]] + data_lines)),
                     sep=r"\s+", engine="python")

    for col in ("Yavg", "Cl", "dSpan", "V/Vref", "Chord"):
        if col not in df.columns:
            raise ValueError(f"缺少列 {col}")
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["Yavg", "Cl", "dSpan", "V/Vref", "Chord"]).reset_index(drop=True)

    # === 计算压力中心 ===
    y  = df["Yavg"].to_numpy(float)
    Cl = df["Cl"].to_numpy(float)
    dS = df["dSpan"].to_numpy(float)
    V = df["V/Vref"].to_numpy(float)
    chord = df["Chord"].to_numpy(float)

    w = Cl * dS * V ** 2
    y_cp = np.sum(w * y) / np.sum(w)

    return y_cp

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
    spanlist = [fuse_w, span - fuse_w]
    if (fuse_w == 0):
        chordlist = [root, tip]
        twistlist = [twist]
    else: 
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
    vsp.ClearVSPModel()
    vsp.SetSetName(vsp.SET_FIRST_USER, "wing")
    vsp.SetSetName(vsp.SET_FIRST_USER + 1, "prop")

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
    vsp.SetSetFlag(wid, vsp.SET_FIRST_USER, True)
    vsp.Update()
    vsp.SetSetFlag(wid, vsp.SET_FIRST_USER + 1, True)
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
        vsp.Update
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
                            "Length_C_Start": len_start, "Length_C_End": len_end, "SE_Const_Flag":0 } #设置起始和末尾相对长度
        else:
            name_to_value = {"EtaFlag": 1, "EtaStart": EtaStart,"EtaEnd": EtaEnd, "Abs_Rel_Flag":0, 
                            "Length_Start": len_start, "Length_End": len_end, "SE_Const_Flag":0 } #设置起始和末尾绝对长度
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
    # 坐标：翼尖桨 + 均匀分布的升力桨 (从翼尖向翼根排列)
    prop_pos = np.concatenate(([span],
        np.linspace(span - 0.5 * (tipprop_Dia + liftprop_Dia) * inch_in_m - gap,
                    gap + 0.5 * liftprop_Dia * inch_in_m + fuse_w,
                    Nprops - 1)))

    # 直径：第一个是翼尖桨，后面是升力桨
    prop_D = [tipprop_Dia] + [liftprop_Dia] * (Nprops - 1)

    prop_ele = [0.002] + [prop_ele] * (len(prop_pos)-1)

    #Add prop
    prop_id = []
    for i in range(0, len(prop_pos)):
        prop_id.append(vsp.AddGeom( "PROP", "" ))
        vsp.SetSetFlag(prop_id[-1], vsp.SET_FIRST_USER + 1, True)
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
    vsp.SetSetFlag(prop_id, vsp.SET_FIRST_USER + 1, True)
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

def runaero(CG, AlphaStart_input, AlphaEnd_input, AlphaNpts_input, air_spd, wing_cfg, Cl_target, sol_config, angle=[], prop_D=[],RPM =[],Ct =[],Cp =[]):
    wing_S = wing_cfg["wing_S"]
    bref   = wing_cfg["bref"]
    cref   = wing_cfg["cref"]
    Nprops = len(prop_D)
    Re = cref * air_spd * density / mu
    print("Re: " + str(Re))
    vsp.DeleteAllResults()
    
    # 控制 PROP 可见性
    geom_ids = vsp.FindGeoms()
    for geom in geom_ids:
        type_name = vsp.GetGeomTypeName(geom).upper()
        if type_name == "PROPELLER" or type_name == "PROP":
            vsp.SetSetFlag(geom, vsp.SET_SHOWN, Nprops > 0)
    vsp.Update()

    # Analysis: VSPAero Compute Geometry to Create Vortex Lattice DegenGeom File #
    compgeom_name = "VSPAEROComputeGeometry"
    print( compgeom_name )
    #Set defaults
    vsp.SetAnalysisInputDefaults( compgeom_name )
    # Analysis method
    vsp.SetIntAnalysisInput( compgeom_name, "Symmetry", [2], 0 )

    vsp.PrintAnalysisInputs( compgeom_name )
    vsp.WriteVSPFile(file_name, vsp.SET_ALL) # 必须在 ExecAnalysis 前保存
    print( "\tExecuting..." )
    compgeom_resid = vsp.ExecAnalysis( compgeom_name)
    print( "COMPLETE" )
    vsp.Update()

    # 提取作动盘并配置性能
    if Nprops > 0:
        num_disks = vsp.GetNumActuatorDisks()
        print(f"\n[INFO] Total Actuator Disks found in Analysis: {num_disks}")
        if num_disks == Nprops:
            for i in range(0, Nprops):
                disk_id = vsp.FindActuatorDisk(i)
                vsp.SetParmValUpdate( vsp.FindParm(disk_id, "RotorRPM", "Rotor"), RPM[i] )  
                vsp.SetParmValUpdate( vsp.FindParm(disk_id, "RotorCT",  "Rotor"),  Ct[i])
                vsp.SetParmValUpdate( vsp.FindParm(disk_id, "RotorCP",  "Rotor"),  Cp[i])
        else:
            print(f"Warning: Expected {Nprops} actuator disks, found {num_disks}")

    # Get & Display Results
    vsp.PrintResults( compgeom_resid )
    degen_name = vsp.GetStringResults(compgeom_resid, "DegenGeomFileName")[0]
    print("Generated DegenGeom file:", degen_name)

    # Analysis: VSPAEROSweep #
    analysis_name = "VSPAEROSweep"
    print( analysis_name )
    # Set Aero defaults
    vsp.SetAnalysisInputDefaults( analysis_name )
    vsp.SetDoubleAnalysisInput( analysis_name, "AlphaEnd", [AlphaEnd_input], 0 )
    vsp.SetDoubleAnalysisInput( analysis_name, "AlphaStart", [AlphaStart_input], 0 )
    vsp.SetIntAnalysisInput( analysis_name, "AlphaNpts", [AlphaNpts_input], 0 )
    vsp.SetDoubleAnalysisInput( analysis_name, "Xcg", [CG], 0 )
    vsp.SetIntAnalysisInput( analysis_name, "NCPU", [CPU], 0 )
    vsp.SetIntAnalysisInput( analysis_name, "Symmetry", [2], 0 ) # 2 for XZ symmetry
    vsp.SetIntAnalysisInput( analysis_name, "PropBladesMode", [0], 0 ) # Static mode for disks
    vsp.SetDoubleAnalysisInput( analysis_name, "Sref", [wing_S], 0 )
    vsp.SetDoubleAnalysisInput( analysis_name, "bref", [bref], 0 )
    vsp.SetDoubleAnalysisInput( analysis_name, "cref", [cref], 0 )
    vsp.SetDoubleAnalysisInput( analysis_name, "ReCref", [Re], 0 )
    vsp.SetDoubleAnalysisInput( analysis_name, "Vinf", [air_spd], 0 )
    vsp.SetDoubleAnalysisInput( analysis_name, "Vref", [air_spd], 0 )
    vsp.SetDoubleAnalysisInput( analysis_name, "Rho", [density], 0 )
    #convergence params
    vsp.SetIntAnalysisInput( analysis_name, "WakeNumIter", [5], 0 )
    vsp.SetIntAnalysisInput( analysis_name, "FarDistToggle", [1], 0 )
    vsp.SetDoubleAnalysisInput( analysis_name, "FarDist", [sol_config["farfield"]], 0 )
    vsp.SetIntAnalysisInput( analysis_name, "NumWakeNodes", [sol_config["wakenode"]], 0 )

    # CS Setting 
    if angle:
        cs_group_container_id = vsp.FindContainer("VSPAEROSettings", 0)
        Num_cs = vsp.GetNumControlSurfaceGroups()
        for i in range(0, Num_cs):
            cs_name = vsp.GetVSPAEROControlGroupName(i)
            if cs_name in angle:
                # 内部组名：ControlSurfaceGroup_0, ControlSurfaceGroup_1, …
                grp_name = f"ControlSurfaceGroup_{i}"
                defl_parm = vsp.FindParm(cs_group_container_id, "DeflectionAngle", grp_name)
                vsp.SetParmValUpdate(defl_parm, angle[cs_name])
    #vsp.WriteVSPFile("auto.vsp3", vsp.SET_ALL)

    vsp.Update()

    #设置完毕        
    print( "Edited input parameter: " )
    vsp.PrintAnalysisInputs( analysis_name )
    vsp.WriteVSPFile(file_name, vsp.SET_ALL)
    # 执行分析
    print("\n[INFO] 开始执行 VSPAEROSweep 分析...")
    res_id = vsp.ExecAnalysis(analysis_name)
    print("[INFO] 分析完成")

    # 可选：输出结果摘要 vsp.PrintResults(res_id)
    polar_name = case_name + ".polar"
    df = pd.read_fwf(polar_name, skiprows=2)
    # 提取需要的数据列
    Cl_list = df['CLtot'].tolist()
    Cd_list = df['CDi'].tolist() #only induced drag
    CMy_list = df['CMytot'].tolist()
    #计算功率，推力
    thrust = 0
    power = 0
    if (len(prop_D)):
        for i in range(Nprops):
            thrust = thrust + prop_D[i]**4 * (RPM[i]/60)**2 * Ct[i] * density * 2 
            power = power + prop_D[i]**5 * (RPM[i]/60)**3 * Cp[i] * density * 2 
    #计算气动力
    if (AlphaNpts_input >= 2):
        alpha_list = np.linspace(AlphaStart_input, AlphaEnd_input, AlphaNpts_input) 
        alpha = np.interp(Cl_target, Cl_list, alpha_list)
        Cd_target = np.interp(Cl_target, Cl_list, Cd_list)
        CMy = np.interp(Cl_target, Cl_list, CMy_list)
        drag = 0.5 * Cd_target * density * wing_S * air_spd**2
        lift = 0.5 * Cl_target * density * wing_S * air_spd**2 + thrust * np.sin(alpha/360 * 2 * np.pi)
    else:
        alpha = AlphaEnd_input
        drag = 0.5 * Cd_list[0] * density * wing_S * air_spd**2
        lift = 0.5 * Cl_list[0] * density * wing_S * air_spd**2 + thrust * np.sin(alpha/360 * 2 * np.pi)
        CMy = CMy_list[0]

    net_drag = drag - thrust* np.cos(alpha/360 * 2 * np.pi) #加上螺旋桨净阻力
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

def single_point(f_cond, geo_info, config):
    debug_log = []
    CG=geo_info["CG"]
    cfg = {"wing_S": geo_info["wing_S"], "bref": geo_info["bref"],"cref": geo_info["cref"]}
    spanlist = geo_info["spanlist"]
    span = geo_info["span"]
    chordlist = geo_info["chordlist"]
    wing_S = geo_info["wing_S"]
    def_cfg = geo_info["def_cfg"]
    prop_D = geo_info["prop_D"]
    prop_D_inch = geo_info["prop_D_inch"]
    prop_pos = geo_info["prop_pos"]
    global mass

    max_AOA =f_cond["max_AOA"]
    speed = f_cond["speed"]
    thrust_ratio = f_cond["TR"]
    d0 = f_cond["d0"]
    Cl_target = f_cond["Cl_target"]

    prop_data = config["propdata"]

    drag_maxit = config["max_it"]
    omega_drag = config["relax"][0]
    omega_lift = config["relax"][1]
    max_step = config["max_alpha_step"]

    ele_def0 = 0
    ele_def1 = def_cfg["elevator"]
    def_cfg["elevator"] = ele_def0
    if (Cl_target == -1):
        drag, alpha, lift, netdrag0, power, _, CMy0 = runaero(CG, max_AOA, max_AOA, 1, speed, cfg, 0, solver_config0, def_cfg)
    else:
        drag, alpha, lift, netdrag0, power, _, CMy0 = runaero(CG, 0, max_AOA, 2, speed, cfg, Cl_target, solver_config0, def_cfg)

    netdrag0 = netdrag0 + d0

    RPM, Ct, Cp = prop.equal_thrust(prop_data, netdrag0, speed, prop_D_inch, thrust_ratio) #大螺旋桨转速
    thrust1 = netdrag0
    #修正升力
    total_lift = lift + netdrag0 * np.sin(alpha/360 * 2 * np.pi)
    Cl_target = mass * g / total_lift * Cl_target
    #第二次气动
    
    def_cfg["elevator"] = ele_def1

    if (Cl_target == -1):
        drag, alpha, lift, netdrag1, power, _, CMy1 = runaero(CG, max_AOA, max_AOA, 1, speed, cfg, 0,
                                                              solver_config1, def_cfg, prop_D, RPM, Ct, Cp)
    else:
        drag, alpha, lift, netdrag1, power, _, CMy1 = runaero(CG, max(-1, alpha-1), min(max_AOA, alpha+1), 2, speed, cfg, Cl_target,
                                                              solver_config1, def_cfg, prop_D, RPM, Ct, Cp)

    netdrag1 = netdrag1 + d0
    #校准质量
    y_cp = read_halfwing_cp()
    mass, wing_mass = mass_sim_iter(span, wing_S, y_cp, prop_pos, spanlist, chordlist)
    
    thrust0 = 0 
    drag_tol = (drag + d0) * config["tol"][0]
    lift_tol = mass * g * config["tol"][1]
    mom_tol  =  CMy0 * config["tol"][2]
    # 初始化历史量
    alpha0, alpha1 = None, alpha
    L0,     L1     = None, lift
    lift_error0, lift_error1 = None, (lift - mass * g)   # 给第一次迭代提供可用的 lift 残差

    # --- 进入主循环 ---
    for i in range(drag_maxit):
        alpha_this_iter = alpha 

        den = (netdrag0 - netdrag1)
        if abs(den) < 1e-4:
            # 解析“牛顿步”：T_new = T_old + netdrag / cos(a)
            c = max(np.cos(np.deg2rad(alpha)), 0.05)  # 防奇异
            thrust_target = thrust1 + netdrag1 / c
        else:
            thrust_target = (netdrag0 * thrust1 - netdrag1 * thrust0) / den
        thrust2 = (1 - omega_drag) * thrust1 + omega_drag * thrust_target

        # 升降舵
        if (abs(CMy1) < mom_tol) or (abs(CMy0 - CMy1) < 1e-4):
            ele_def2 = ele_def1  # 保持不变
        else:
            ele_def2 = (CMy0 * ele_def1 - CMy1 * ele_def0) / (CMy0 - CMy1)

        # 应用输入值并进行当前迭代的气动计算 ---
        def_cfg["elevator"] = ele_def2
        RPM, Ct, Cp = prop.equal_thrust(prop_data, thrust2, speed, prop_D_inch, thrust_ratio)

        if alpha > max_AOA:
            warnings.warn(f"Stall may happen!! Alpha={alpha:.2f}", category=None, stacklevel=1)
        
        drag, _, lift2, netdrag2, power, _, CMy2 = runaero(CG, alpha, alpha, 1, speed, cfg, 0, solver_config1, def_cfg, prop_D, RPM, Ct, Cp)

        # 计算本次迭代产生的误差 ---
        netdrag2 = netdrag2 + d0
        lift_error2 = lift2 - (mass * g)

        debug_log.append({
            "Iteration": i + 1,
            "Alpha": alpha,
            "Lift": lift2,
            "Drag": drag,
            "Thrust": thrust2,
            "Weight": mass * g,
            "CMy": CMy2,
            "Elevator": ele_def2,
            "NetDrag": netdrag2
        })

        # 检查是否收敛 ---
        if (abs(netdrag2) < drag_tol) and (abs(lift_error2) < lift_tol):
            print(f"Converged in {i+1} iterations!")
            break

        # ====== 更新历史点 ======
        netdrag0, netdrag1 = netdrag1, netdrag2
        thrust0,  thrust1  = thrust1,  thrust2
        CMy0,     CMy1     = CMy1,     CMy2
        ele_def0, ele_def1 = ele_def1, ele_def2
        alpha0,   alpha1   = alpha1,   alpha
        lift_error0, lift_error1 = lift_error1, lift_error2
        L0, L1 = L1, lift2

        # 迎角
        if (Cl_target !=-1): #仅非升力约束下考虑
            if (lift_error1 is not None) and (abs(lift_error1) < lift_tol):
                # 保持 alpha 不变，继续让推力/力矩收敛
                alpha = alpha1
                continue

            # 用 lift 的割线斜率做一次牛顿步 + 动态松弛 + 步长限幅
            # 物理保护：设定物理斜率的上下限，防止推力耦合干扰导致步长过小或方向错误
            q_inf = 0.5 * density * speed**2
            dL_dalpha_min = 0.05 * wing_S * q_inf  # 每度约 5% CL 增量的保守底线
            dL_dalpha_max = 0.15 * wing_S * q_inf  # 每度约 15% CL 增量的物理顶线
            
            if i < 1 or (alpha0 is None) or (L0 is None) or abs(alpha1 - alpha0) < 1e-4:
                # 历史不足：用比例法的小步修正，但限制斜率上限以防步长过小
                S_alpha_guess = abs(lift2) / max(abs(alpha1), 1e-3)
                S_alpha_guess = np.clip(S_alpha_guess, dL_dalpha_min, dL_dalpha_max)
                delta_alpha = (mass * g - lift2) / max(S_alpha_guess, 1e-4)
            else:
                dL_dalpha = (L1 - L0) / (alpha1 - alpha0)
                # 如果估算的斜率不物理（太小、为负、或因数据点过近导致极大），强制限制在物理区间
                dL_dalpha = np.clip(dL_dalpha, dL_dalpha_min, dL_dalpha_max)
                
                delta_alpha = (mass * g - lift2) / dL_dalpha

            # 步长限幅 + 抑制震荡
            delta_alpha = np.clip(delta_alpha, -max_step, max_step)
            alpha_target = alpha1 + delta_alpha

            alpha = (1.0 - omega_lift) * alpha1 + omega_lift * alpha_target
        else:
            alpha = max_AOA

    else:
        print("Warning: Did not converge within the maximum number of iterations.")

    mass_result = {"mass":mass, "wing_mass":wing_mass}
    if debug_log:
        speed = f_cond.get("speed", 0)
        filename = f"convergence_{case_name}_V{speed:.1f}.csv"
        pd.DataFrame(debug_log).to_csv(filename, index=False)
        print(f"Convergence log saved to {filename}")
    return lift2, drag, power, alpha_this_iter, RPM, thrust2, mass_result, ele_def2

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
    cruise_spd = 12
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