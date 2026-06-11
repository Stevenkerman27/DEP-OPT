import openvsp as vsp
import os

vsp.ClearVSPModel()
prop_id = vsp.AddGeom( "PROP", "" )
vsp.SetParmVal( prop_id, "PropMode", "Design", 2.0 ) # PROP_DISK
vsp.Update()

print(f"\n[INFO] PROP Geom ID: {prop_id}")
for p_id in vsp.FindContainerParmIDs(prop_id):
    print(f"    Geom Parm -> Group: {vsp.GetParmGroupName(p_id):15} | Name: {vsp.GetParmName(p_id):25} | Value: {vsp.GetParmVal(p_id)}")

num_disks = vsp.GetNumActuatorDisks()
if num_disks > 0:
    disk_id = vsp.FindActuatorDisk(0)
    print(f"\n[INFO] Actuator Disk ID: {disk_id}")
    groups = vsp.FindContainerGroupNames(disk_id)
    for g in groups:
        print(f"  Group: {g}")
    for p_id in vsp.FindContainerParmIDs(disk_id):
        print(f"    Rotor Parm: {vsp.GetParmName(p_id)} (Group: {vsp.GetParmGroupName(p_id)}) = {vsp.GetParmVal(p_id)}")

# Check for hidden "Driver" parameter in Rotor group
print("\n[INFO] Searching for 'Driver' parameters in Actuator Disk...")
for p_id in vsp.FindContainerParmIDs(disk_id):
    if "Driver" in vsp.GetParmName(p_id):
        print(f"  FOUND: {vsp.GetParmName(p_id)} in Group {vsp.GetParmGroupName(p_id)}")

# Try to find a parameter named "PropellerDriver" even if not listed
try:
    p = vsp.FindParm(disk_id, "PropellerDriver", "Rotor")
    if p:
        print(f"  FOUND hidden parm: PropellerDriver = {vsp.GetParmVal(p)}")
except:
    pass
