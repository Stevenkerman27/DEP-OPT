import os
import re
import numpy as np
import matplotlib.pyplot as plt

prop_choice = {8:"8x4", 9:"9x6", 10:"10x7", 11:"11x7", 12:"12x8", 13:"13x8", 86.6:"86.6"}

def read_apce_grouped(data_dir):
    # 匹配 apce_<size>_…_<4 位 RPM>.txt，并捕获 size（如 10x5）和 rpm
    pattern = re.compile(
    r'^apce_'
    r'(\d+(?:\.\d+)?x\d+(?:\.\d+)?).*_'  # 这里让“数字.数字”变成可选
    r'(\d{4})\.txt$'
)
    grouped = {}

    for fname in os.listdir(data_dir): #对每一个文件
        m = pattern.match(fname)
        if not m:
            continue
        size, rpm = m.group(1), int(m.group(2)) #捕获正则表达式中的尺寸和转速
        m = re.match(r'(\d+(?:\.\d+)?)x(\d+(?:\.\d+)?)', size)
        dia = int(m.group(1))
        path = os.path.join(data_dir, fname)
        # loadtxt：假设文件中前 3 列依次为 J, CT, CP
        arr = np.genfromtxt(path,names=True, dtype=None, encoding=None)
        # 然后按名字取列
        J  = arr['J']
        U = J * dia * 0.0254 * rpm/60 #来流速度
        CT = arr['CT']
        CP = arr['CP']
        eta = arr['eta']

        grouped.setdefault(size, []).append((rpm, U, CT, CP, eta))#多个表示相同螺旋桨尺寸的数据合并
    return grouped

def pick_propdata(prop_data, req_dia, req_speed, req_RPM):
    # 1) 选尺寸
    try:
        size = prop_choice[req_dia]
    except KeyError:
        raise ValueError(f"未知直径 {req_dia}，prop_choice 中无对应条目")

    grouped = prop_data
    if size not in grouped:
        raise ValueError(f"尺寸 {size} 在数据中不存在")

    curves = grouped[size]

    curves_sorted = sorted(curves, key=lambda x: x[0])
    # rpm_list = [rpm for rpm, *_ in curves_sorted] # 这行现在不是必需的

    CT_at_speed = []
    CP_at_speed = []
    eta_at_speed = []
    rpm_list_valid = []
    
    for rpm, U_arr, CT_arr, CP_arr, eta_arr in curves_sorted:
        # 1) 排序
        sort_idx = np.argsort(U_arr)
        Us   = U_arr[sort_idx]
        CTs  = CT_arr[sort_idx]
        CPs  = CP_arr[sort_idx]
        etas = eta_arr[sort_idx]

        if not (Us[0] <= req_speed <= Us[-1]):
            continue

        # 2) 在速度维度插值
        CT_s = np.interp(req_speed, Us, CTs)
        CP_s = np.interp(req_speed, Us, CPs)
        eta_s = np.interp(req_speed, Us, etas)

        CT_at_speed.append(CT_s)
        CP_at_speed.append(CP_s)
        eta_at_speed.append(eta_s)
        rpm_list_valid.append(rpm)

    # 5) 在 rpm 维度插值，得最终 CT, CP, eta
    CT_out = np.interp(req_RPM,  rpm_list_valid, CT_at_speed)
    CP_out = np.interp(req_RPM,  rpm_list_valid, CP_at_speed)
    eta_out = np.interp(req_RPM, rpm_list_valid, eta_at_speed)

    # 根据公式 K_Q = C_P / (2 * pi)
    KQ_out = CP_out / (2 * np.pi)

    return CT_out, CP_out, eta_out, KQ_out

def read_xrotor_grouped(data_dir):
    """
    读取 X-Rotor 数据（如 xr86.6_70.txt），并返回与 read_apce_grouped 相同格式的字典：
      { size_key: [ (rpm, U_arr, CT_arr, CP_arr, eta_arr), ... ] } <--- 返回值现在是5元组
    prop_choice 用来指定哪些直径要读取（key: float(dia) -> 任意value）。
    """
    xr_pattern = re.compile(r'^xr(\d+(?:\.\d+)?)[_.]\d+(?:\.txt)?$')
    grouped = {}

    for fname in os.listdir(data_dir):
        m = xr_pattern.match(fname)
        if not m:
            continue

        dia = float(m.group(1))
        size_key = prop_choice.get(dia)
        if size_key is None:
            continue

        path = os.path.join(data_dir, fname)
        with open(path, 'r') as f:
            lines = f.readlines()

        # 找到包含 rpm 和 V 的列头行
        header_idx = next((i for i, L in enumerate(lines) if 'rpm' in L and 'V' in L), None)
        if header_idx is None:
            raise RuntimeError(f"No header with 'rpm' and 'V' in {fname}")

        cols = re.split(r'\s+', lines[header_idx].strip())
        idx_rpm = cols.index('rpm')
        idx_V   = cols.index('V')
        idx_CT  = cols.index('CT')
        idx_CP  = cols.index('CP')
        idx_eta = cols.index('eff') # <--- 1. 新增：找到 'eff' 列的索引

        # 读取所有行到列表
        rpm_list, V_list, CT_list, CP_list = [], [], [], []
        eta_list = [] # <--- 2. 新增：为 eta 创建一个列表

        for row in lines[header_idx+1:]:
            parts = re.split(r'\s+', row.strip())
            # <--- 3. 修改：更新最大索引检查，确保包含 idx_eta
            if len(parts) <= max(idx_rpm, idx_V, idx_CT, idx_CP, idx_eta):
                continue
            try:
                rpm_list.append(float(parts[idx_rpm]))
                V_list.append(   float(parts[idx_V]))
                CT_list.append(  float(parts[idx_CT]))
                CP_list.append(  float(parts[idx_CP]))
                eta_list.append( float(parts[idx_eta])) # <--- 4. 新增：读取 eta 数据
            except ValueError:
                continue

        # 转为 NumPy 数组
        rpm_arr = np.array(rpm_list)
        V_arr   = np.array(V_list)
        CT_arr  = np.array(CT_list)
        CP_arr  = np.array(CP_list)
        eta_arr = np.array(eta_list) # <--- 5. 新增：转换 eta 列表

        # 对每个唯一 rpm 值分组，将对应的子数组作为一条曲线
        for rpm_val in np.unique(rpm_arr):
            mask = (rpm_arr == rpm_val)
            U_sub  = V_arr[mask]
            CT_sub = CT_arr[mask]
            CP_sub = CP_arr[mask]
            eta_sub = eta_arr[mask] # <--- 5. 新增：提取 eta 子数组
            
            # 按 U 排序，确保后续插值正常
            sort_idx = np.argsort(U_sub)
            U_sorted  = U_sub[sort_idx]
            CT_sorted = CT_sub[sort_idx]
            CP_sorted = CP_sub[sort_idx]
            eta_sorted = eta_sub[sort_idx] # <--- 5. 新增：排序 eta 子数组

            # <--- 5. 修改：附加 5 元素元组
            grouped.setdefault(size_key, []).append(
                (float(rpm_val), U_sorted, CT_sorted, CP_sorted, eta_sorted)
            )

    return grouped

def equal_thrust(data, thrust, req_speed, prop_property, thrust_ratio = 1):
    rho = 1.225              # 空气密度 [kg/m³]
    thrust_side = thrust/2   # 单侧机翼总推力
    if (thrust_ratio != 1):
        large_dia  = prop_property[0]
        small_dias = prop_property[1:]
    else:
        small_dias = []
        large_dia = prop_property[0]
    
    # 目标推力分配
    T_large_target  = thrust_side * thrust_ratio
    T_smalls_target = thrust_side - T_large_target

    # 计算给定 rpm 时的小螺旋桨总推力
    def smalls_residual(rpm_small):
        n = rpm_small/60.0
        Tsum = 0.0
        for d_inch in small_dias:
            CT, CP, eta, KQ = pick_propdata(data, d_inch, req_speed, rpm_small)
            D = d_inch * 0.0254
            Tsum += CT * rho * n**2 * D**4
        return Tsum - T_smalls_target

    if (thrust_ratio != 1):
        # 二分法求小螺旋桨 rpm
        lo, hi = 0.0, 7000.0
        for _ in range(20):
            mid = 0.5*(lo+hi)
            if smalls_residual(mid) >= 0:
                hi = mid
            else:
                lo = mid
        rpm_small = 0.5*(lo+hi)
    else: 
        rpm_small = 0

    # 计算给定 rpm 时的大螺旋桨推力残差
    def large_residual(rpm_large):
        nL = rpm_large/60.0
        CT_L, CP, eta, KQ = pick_propdata(data, large_dia, req_speed, rpm_large)
        D_L = large_dia * 0.0254
        T_L = CT_L * rho * nL**2 * D_L**4
        return T_L - T_large_target

    # 二分法求大螺旋桨 rpm
    lo, hi = 0.0, 7000.0
    for _ in range(20):
        mid = 0.5*(lo+hi)
        if large_residual(mid) >= 0:
            hi = mid
        else:
            lo = mid
    rpm_large = 0.5*(lo+hi)

    rpms = [rpm_large] + [rpm_small] * len(small_dias)
    # 一次性拿到所有 (CT, CP)
    results = []
    for d, rpm in zip(prop_property, rpms):
        CT_val, CP_val, _, _ = pick_propdata(data, d, req_speed, rpm)
        
        results.append((CT_val, CP_val))
    # 解包成两个列表
    CT, CP = map(list, zip(*results))
    return rpms, CT, CP


if __name__ == '__main__':
    rho = 1.225
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    data = read_apce_grouped('data')
    print(pick_propdata(data,13,7.5,3000)) # CT CP eta KQ
    print("13inch prop data printed")
    print(pick_propdata(data,8,7.5,5000)) # CT CP eta KQ
    print("8inch prop data printed")
    #print(equal_thrust(data, 2, 12, [8,8,8,12], 0.5))
    print('reading xr data')
    xr_data = read_xrotor_grouped('data')
    data.update(xr_data)
    #print(pick_propdata(xr_data, 86.6, 70, 5000))
    #print('xr data read')
    sizes_to_plot = [ '8x4', '13x8']
    mymarkers = ['o', 's', 'v']
    linestyles  = ['-', '--', '-.']
    styleindex = [0,1,2]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8,10), sharex=True)
    for size, mystyle in zip(sizes_to_plot, styleindex):
        if size not in data:
            print(f"Warning: no data for size {size!r}")
            continue

        # data[size] 是一个列表，元素为 (rpm, U_array, CT_array, CP_array)
        for rpm, U, CT, CP, eta in data[size]:
            m = re.match(r'(\d+)x\d+', size)
            dia = int(m.group(1))
            T = rho * CT * (dia**4) * (rpm/60)**2
            P = rho * CP * (dia**5) * (rpm/60)**3
            ax1.plot(U, CT,marker=mymarkers[mystyle],linestyle=linestyles[mystyle],label=f'{size}, {rpm} RPM')
            ax2.plot(U, eta,marker=mymarkers[mystyle],linestyle=linestyles[mystyle],label=f'{size}, {rpm} RPM')
            
    ax1.set_ylabel('CT')
    ax1.set_title('CT vs Airspeed')
    ax1.set_ylim(bottom=-0.01)
    ax1.grid(True)
    ax2.set_xlabel('Airspeed (U)')
    ax2.set_ylabel('eta')
    ax2.set_title('eta vs Airspeed')
    ax2.set_ylim(bottom=0)
    ax2.grid(True)
    handles, labels = ax2.get_legend_handles_labels()
    fig.legend(handles, labels)
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.show()