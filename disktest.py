import openvsp as vsp
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(script_dir, "outputs")
os.makedirs(output_dir, exist_ok=True)
os.chdir(output_dir)

vsp.ClearVSPModel()

# 1. Create Rectangular Wing
wid = vsp.AddGeom("WING", "") 
vsp.SetGeomName(wid, "TestWing")
vsp.SetParmVal(wid, "Root_Chord", "XSec_1", 0.5)
vsp.SetParmVal(wid, "Tip_Chord", "XSec_1", 0.5)
vsp.SetParmVal(wid, "Span", "XSec_1", 2.0)
vsp.SetParmVal(wid, "Sweep", "XSec_1", 0.0)

# Set Mesh Density
vsp.SetParmVal(wid, "Tess_W", "Shape", 20)
vsp.SetParmVal(wid, "SectTess_U", "XSec_1", 30)

# Turn off wing symmetry
sym_parm_wing = vsp.FindParm(wid, "Sym_Planar_Flag", "Sym") 
vsp.SetParmVal(sym_parm_wing, 0) # 0 for none

# Add Control Surface (Flap)
# Using AddSubSurf
subsurf_id = vsp.AddSubSurf(wid, vsp.SS_CONTROL, 0)
vsp.Update()
parm_id_vec = vsp.GetSubSurfParmIDs(subsurf_id)
# EtaFlag=1, Abs_Rel_Flag=1 for relative length, Length_C_Start, Length_C_End, EtaStart, EtaEnd
name_to_value = {
    "EtaFlag": 1, "EtaStart": 0.8, "EtaEnd": 0.2, 
    "Abs_Rel_Flag": 1, "Length_C_Start": 0.25, "Length_C_End": 0.25
}
for pid in parm_id_vec:
    name = vsp.GetParmName(pid)
    if name in name_to_value:
        vsp.SetParmVal(pid, name_to_value[name])
vsp.Update()

# Create Control Surface Group
group_index = vsp.CreateVSPAEROControlSurfaceGroup()
available_cs = vsp.GetAvailableCSNameVec(group_index)
prefix = "TestWing_"
indices = []
for idx, cs_name in enumerate(available_cs):
    if cs_name.startswith(prefix):
        indices.append(idx + 1)
vsp.Update()
cs_group_name = "test_flap"
vsp.SetVSPAEROControlGroupName(cs_group_name, group_index)
vsp.AddSelectedToCSGroup(indices, group_index)
vsp.Update()

# 2. Add Propellers
# Prop 1
prop1_id = vsp.AddGeom("PROP", "")
vsp.SetGeomName(prop1_id, "Prop1")
vsp.SetParmValUpdate(vsp.FindParm(prop1_id, "Sym_Planar_Flag", "Sym"), 0)
vsp.SetParmVal(prop1_id, "PropMode", "Design", vsp.PROP_DISK)
vsp.SetParmVal(prop1_id, "Diameter", "Design", 0.3)
vsp.SetParmVal(prop1_id, "X_Rel_Location", "XForm", -0.1)
vsp.SetParmVal(prop1_id, "Y_Rel_Location", "XForm", 0.5)

# Prop 2
prop2_id = vsp.AddGeom("PROP", "")
vsp.SetGeomName(prop2_id, "Prop2")
vsp.SetParmValUpdate(vsp.FindParm(prop2_id, "Sym_Planar_Flag", "Sym"), 0)
vsp.SetParmVal(prop2_id, "PropMode", "Design", vsp.PROP_DISK)
vsp.SetParmVal(prop2_id, "Diameter", "Design", 0.5)
vsp.SetParmVal(prop2_id, "X_Rel_Location", "XForm", -0.1)
vsp.SetParmVal(prop2_id, "Y_Rel_Location", "XForm", 1.5)

vsp.Update()

# Write VSP3 file before ComputeGeometry
vsp.WriteVSPFile("disktest.vsp3", vsp.SET_ALL)

# 3. Compute Geometry
compgeom_name = "VSPAEROComputeGeometry"
vsp.SetAnalysisInputDefaults(compgeom_name)
vsp.SetIntAnalysisInput(compgeom_name, "Symmetry", [2], 0) # 2 for XZ
print("\n\tExecuting ComputeGeometry...")
vsp.ExecAnalysis(compgeom_name)
print("COMPLETE")

# 4. Configure Actuator Disks
vsp.Update()
num_disks = vsp.GetNumActuatorDisks()
print(f"\n[INFO] Total Actuator Disks found in Analysis: {num_disks}")

if num_disks == 2:
    disk1_id = vsp.FindActuatorDisk(0)
    vsp.SetParmValUpdate(vsp.FindParm(disk1_id, "RotorRPM", "Rotor"), 3000)
    vsp.SetParmValUpdate(vsp.FindParm(disk1_id, "RotorCT",  "Rotor"), 0.5)
    vsp.SetParmValUpdate(vsp.FindParm(disk1_id, "RotorCP",  "Rotor"), 0.5)

    disk2_id = vsp.FindActuatorDisk(1)
    vsp.SetParmValUpdate(vsp.FindParm(disk2_id, "RotorRPM", "Rotor"), 5000)
    vsp.SetParmValUpdate(vsp.FindParm(disk2_id, "RotorCT",  "Rotor"), 0.8)
    vsp.SetParmValUpdate(vsp.FindParm(disk2_id, "RotorCP",  "Rotor"), 0.8)
else:
    print("Warning: Expected 2 actuator disks, found ", num_disks)

# 5. Sweep Analysis
analysis_name = "VSPAEROSweep"
vsp.SetAnalysisInputDefaults(analysis_name)
vsp.SetDoubleAnalysisInput(analysis_name, "AlphaEnd", [5], 0)
vsp.SetDoubleAnalysisInput(analysis_name, "AlphaStart", [0], 0)
vsp.SetDoubleAnalysisInput(analysis_name, "Vinf", [12], 0)
vsp.SetDoubleAnalysisInput(analysis_name, "Vref", [12], 0)
vsp.SetDoubleAnalysisInput(analysis_name, "Sref", [2.0], 0) # Sref 
vsp.SetIntAnalysisInput(analysis_name, "AlphaNpts", [2], 0)
vsp.SetIntAnalysisInput(analysis_name, "NCPU", [8], 0)
vsp.SetIntAnalysisInput(analysis_name, "Symmetry", [2], 0) # 2 for XZ symmetry
vsp.SetIntAnalysisInput(analysis_name, "PropBladesMode", [0], 0)

# Set Deflection
angle = {"test_flap": 20}
cs_group_container_id = vsp.FindContainer("VSPAEROSettings", 0)
Num_cs = vsp.GetNumControlSurfaceGroups()
for i in range(0, Num_cs):
    cs_name = vsp.GetVSPAEROControlGroupName(i)
    if cs_name in angle:
        grp_name = f"ControlSurfaceGroup_{i}"
        defl_parm = vsp.FindParm(cs_group_container_id, "DeflectionAngle", grp_name)
        vsp.SetParmValUpdate(defl_parm, angle[cs_name])

vsp.Update()
vsp.WriteVSPFile("disktest.vsp3", vsp.SET_ALL)
print("\n[INFO] 开始执行 VSPAEROSweep 分析...")
vsp.ExecAnalysis(analysis_name)
print("[INFO] 分析完成")
