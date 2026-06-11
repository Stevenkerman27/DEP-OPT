# OpenVSP 3.50+ Python API 调用注意事项 (特别是 VSPAERO 与作动盘)

在从较早版本的 OpenVSP（如 3.2x, 3.4x）升级到 3.50+ 时，其底层数据结构和 API（特别是针对 VSPAERO 和 Actuator Disk 作动盘）发生了许多关键性的变化。以下是编写自动化脚本时必须遵守的规范和踩坑记录。

## 1. 几何体变更与作动盘激活
*   **组件名称变更**：旧版中的 `PROPELLER` 几何体类型已被废弃，现在统一使用 `PROP`。
*   **激活作动盘**：必须显式地将 `PROP` 组件的设计模式更改为作动盘模式：
    ```python
    prop_id = vsp.AddGeom( "PROP", "" )
    vsp.SetParmVal( prop_id, "PropMode", "Design", vsp.PROP_DISK ) # 2.0 即 PROP_DISK
    ```

## 2. 集合过滤 (GeomSet 与 ThinGeomSet)
*   **避免手动覆盖集合配置**：在 OpenVSP 3.50+ 中，VSPAERO 提取作动盘等特殊容器的逻辑与 `SET_SHOWN` (可见集) 强绑定。
*   如果在使用 `vsp.SetIntAnalysisInput` 时强制覆盖了 `GeomSet` 或 `ThinGeomSet`（例如设为 `SET_NONE` 和 `SET_FIRST_USER`），会导致分析引擎在提取结构时“漏掉”作动盘。
*   **最佳实践**：
    1. 调用 `vsp.SetAnalysisInputDefaults(analysis_name)` 初始化默认值（默认 `GeomSet` 为 `SET_NONE`，`ThinGeomSet` 为 `SET_SHOWN`）。
    2. **不修改**集合分析输入参数。
    3. 如果需要过滤掉某些几何体（如只分析机翼，不分析机身），请直接通过修改几何体的**可见性**来控制：
       ```python
       vsp.SetSetFlag(geom_id, vsp.SET_SHOWN, False) # 从分析中剔除
       vsp.SetSetFlag(geom_id, vsp.SET_SHOWN, True)  # 纳入分析
       ```

## 3. 作动盘性能参数注入时机
*   **分离的容器**：在 3.50+ 中，作动盘的性能参数（如 `RotorRPM`, `RotorCT`, `RotorCP`）不再挂载于 `PROP` 几何体下，而是被独立放置在一个名为 Actuator Disk 的特殊容器内。
*   **注入时机**：这个容器只有在运行了 `VSPAEROComputeGeometry` 之后才会被正确初始化和暴露。如果在几何计算前注入，参数将无法挂载并会在计算时被默认值覆盖（RPM = -2000, CT = 0.4, CP = 0.6）。
*   **正确流程**：
    ```python
    vsp.ExecAnalysis("VSPAEROComputeGeometry")
    vsp.Update()
    
    num_disks = vsp.GetNumActuatorDisks()
    if num_disks > 0:
        disk_id = vsp.FindActuatorDisk(0) # 获取第一个作动盘
        vsp.SetParmValUpdate( vsp.FindParm(disk_id, "RotorRPM", "Rotor"), 5000 )
        vsp.SetParmValUpdate( vsp.FindParm(disk_id, "RotorCT",  "Rotor"), 0.5 )
        vsp.SetParmValUpdate( vsp.FindParm(disk_id, "RotorCP",  "Rotor"), 0.5 )
    ```

## 4. 滑流力矩 (Swirl) 与 PropBladesMode、对称性
*   **关于 PropBladesMode**：对于纯作动盘 (`PROP_DISK` 模式)，VSPAERO 的 `PropBladesMode` (0=Static, 1=Unsteady, 2=Pseudo-Steady) 其实是**不起作用/冗余的**。作动盘作为一个纯数学的边界条件，其后洗（Downwash）由 `RotorCT` 决定，而其滑流旋转（Swirl）完全且直接由 `RotorCP` (功率/扭矩系数) 决定。
    *   *如果你想在作动盘上完全关闭滑流扭转效应（只留轴向加速），你应该将 `RotorCP` 设为 `0`。*
*   **不要开启对称性**：如果作动盘位于 $Y=0$ 对称面上，且你开启了气动对称平面（`Symmetry = 1` 或 `2`），那么作动盘产生的单边旋转滑流在计算全机气动力矩（`CMxtot`, `CMztot`）时会被镜像的相反力矩在数学上精确抵消为 `0.0`。
*   **解决方案**：要分析单桨滑流对滚转和偏航的影响，必须关闭对称性：
    ```python
    vsp.SetIntAnalysisInput("VSPAEROSweep", "Symmetry", [0], 0) 
    ```

## 5. 分析文件的命名与 VSP3 保存
*   在执行 `VSPAEROComputeGeometry` 或 `VSPAEROSweep` 等分析时，VSPAERO 需要基于当前的工作文件名来输出气动结果文件（如 `.polar`, `.history`, `.lod`）。
*   **坑点**：如果你的模型只在内存中创建/修改，而在执行分析前**没有**调用 `WriteVSPFile` 将其保存为 `.vsp3` 文件，VSPAERO 无法获取正确的文件前缀名，它会将所有的输出文件全部堆到目录下的 `Unnamed.polar`、`Unnamed.history` 中。
*   这不仅会让输出结果杂乱，还会导致脚本去读取指定的 `.polar` 文件时，读到的是旧的、残留的错误数据。
*   **正确流程**：
    ```python
    # 配置好分析的所有输入参数
    vsp.Update()
    
    # 必须在 ExecAnalysis 前保存，以便 VSPAERO 识别输出文件名
    vsp.WriteVSPFile("my_project_name.vsp3", vsp.SET_ALL) 
    
    res_id = vsp.ExecAnalysis("VSPAEROSweep")
    # 此时结果会被正确写入 my_project_name.polar 中
    ```
