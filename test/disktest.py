import openvsp as vsp
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(script_dir, "outputs")
os.makedirs(output_dir, exist_ok=True)
os.chdir(output_dir)

vsp.ClearVSPModel()

# 1. Primary Wing (Shown)
wid1 = vsp.AddGeom( "WING", "" ) 
vsp.SetGeomName(wid1, "MainWing")
vsp.SetParmVal( wid1, "Root_Chord", "XSec_1", 0.2)
vsp.SetParmVal( wid1, "Tip_Chord", "XSec_1", 0.1)
vsp.SetParmVal( wid1, "Span", "XSec_1", 0.2)
# SET_SHOWN is True by default for new geoms, but let's be explicit
vsp.SetSetFlag(wid1, vsp.SET_SHOWN, True)

# 2. Primary Propeller as Actuator Disk (Shown)
prop_id1 = vsp.AddGeom( "PROP", "" )
vsp.SetGeomName(prop_id1, "MainProp")
vsp.SetParmVal( prop_id1, "PropMode", "Design", vsp.PROP_DISK )
vsp.SetParmVal( prop_id1, "Diameter", "Design", 0.2 )
vsp.SetParmVal( prop_id1, "X_Rel_Location", "XForm", -0.05)
vsp.SetSetFlag(prop_id1, vsp.SET_SHOWN, True)

# 3. Secondary Wing (Hidden)
wid2 = vsp.AddGeom( "WING", "" ) 
vsp.SetGeomName(wid2, "HiddenWing")
vsp.SetParmVal( wid2, "Root_Chord", "XSec_1", 0.5)
vsp.SetParmVal( wid2, "Span", "XSec_1", 1.0)
vsp.SetSetFlag(wid2, vsp.SET_SHOWN, False) # Hide it

# 4. Secondary Propeller as Actuator Disk (Hidden)
prop_id2 = vsp.AddGeom( "PROP", "" )
vsp.SetGeomName(prop_id2, "HiddenProp")
vsp.SetParmVal( prop_id2, "PropMode", "Design", vsp.PROP_DISK )
vsp.SetParmVal( prop_id2, "Diameter", "Design", 1.0 )
vsp.SetParmVal( prop_id2, "X_Rel_Location", "XForm", 2.0)
vsp.SetSetFlag(prop_id2, vsp.SET_SHOWN, False) # Hide it

vsp.Update()

# Calculate actual wing area for Sref (Only MainWing matters if filtering works)
# Area = Span * (Root + Tip) / 2 = 0.2 * (0.2 + 0.1) / 2 = 0.03 (per side)
# Total Area = 0.06
wing_S = 0.06

# Run ComputeGeometry
compgeom_name = "VSPAEROComputeGeometry"
vsp.SetAnalysisInputDefaults( compgeom_name )

print( "\n\tExecuting ComputeGeometry..." )
vsp.ExecAnalysis( compgeom_name)
print( "COMPLETE" )

# Configure Actuator Disk parameters AFTER ComputeGeometry
vsp.Update()
num_disks = vsp.GetNumActuatorDisks()
print(f"\n[INFO] Total Actuator Disks found in Analysis: {num_disks}")
print(f"       (Expected: 1, because the second prop is hidden)")

if num_disks > 0:
    for i in range(num_disks):
        disk_id = vsp.FindActuatorDisk(i)
        print(f"       Configuring Actuator Disk: {disk_id}")
        vsp.SetParmValUpdate( vsp.FindParm(disk_id, "RotorRPM", "Rotor"), 5000 )  
        vsp.SetParmValUpdate( vsp.FindParm(disk_id, "RotorCT",  "Rotor"),  0.5)
        vsp.SetParmValUpdate( vsp.FindParm(disk_id, "RotorCP",  "Rotor"),  0.5)
else:
    print("Warning: No actuator disks found to configure!")

# Run VSPAEROSweep
analysis_name = "VSPAEROSweep"
vsp.SetAnalysisInputDefaults( analysis_name )
vsp.SetDoubleAnalysisInput( analysis_name, "AlphaEnd", [2], 0 )
vsp.SetDoubleAnalysisInput( analysis_name, "AlphaStart", [0], 0 )
vsp.SetDoubleAnalysisInput( analysis_name, "Vinf", [12], 0 )
vsp.SetDoubleAnalysisInput( analysis_name, "Vref", [12], 0 )
vsp.SetDoubleAnalysisInput( analysis_name, "Sref", [wing_S], 0 )
vsp.SetIntAnalysisInput( analysis_name, "AlphaNpts", [2], 0 )
vsp.SetIntAnalysisInput( analysis_name, "NCPU", [8], 0 )
vsp.SetIntAnalysisInput( analysis_name, "Symmetry", [0], 0 )
vsp.SetIntAnalysisInput( analysis_name, "PropBladesMode", [0], 0 )

vsp.Update()
vsp.WriteVSPFile("disktest.vsp3", vsp.SET_ALL)
print("\n[INFO] 开始执行 VSPAEROSweep 分析...")
res_id = vsp.ExecAnalysis(analysis_name)
print("[INFO] 分析完成")