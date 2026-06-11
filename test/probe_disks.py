import openvsp as vsp
import os

vsp.ClearVSPModel()
prop_id = vsp.AddGeom( "PROP", "" )
vsp.SetParmVal( prop_id, "PropMode", "Design", vsp.PROP_DISK )
vsp.SetParmVal( prop_id, "Diameter", "Design", 0.2 )
vsp.SetSetFlag(prop_id, vsp.SET_FIRST_USER + 0, True)
vsp.Update()

compgeom_name = "VSPAEROComputeGeometry"
vsp.SetAnalysisInputDefaults( compgeom_name )
vsp.SetIntAnalysisInput(compgeom_name, 'ThinGeomSet', [vsp.SET_FIRST_USER + 0], 0) 

print("Before DeleteAllResults: NumActuatorDisks =", vsp.GetNumActuatorDisks())
print("Before SetAnalysisInputDefaults: NumActuatorDisks =", vsp.GetNumActuatorDisks())
vsp.SetAnalysisInputDefaults( compgeom_name )
print("After SetAnalysisInputDefaults: NumActuatorDisks =", vsp.GetNumActuatorDisks())

analysis_name = "VSPAEROSweep"
vsp.SetAnalysisInputDefaults( analysis_name )
vsp.SetIntAnalysisInput(analysis_name, "ActuatorDiskFlag", [1], 0)
vsp.Update()

print("After setting ActuatorDiskFlag: NumActuatorDisks =", vsp.GetNumActuatorDisks())

for i in range(vsp.GetNumActuatorDisks()):
    print(f"Disk {i} ID: {vsp.FindActuatorDisk(i)}")
