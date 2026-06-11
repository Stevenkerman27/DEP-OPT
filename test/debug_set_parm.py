import openvsp as vsp
import os

vsp.ClearVSPModel()
prop_id = vsp.AddGeom( "PROP", "" )
vsp.SetParmVal( prop_id, "PropMode", "Design", 2.0 ) # PROP_DISK
vsp.Update()

num_disks = vsp.GetNumActuatorDisks()
if num_disks > 0:
    disk_id = vsp.FindActuatorDisk(0)
    print(f"Actuator Disk ID: {disk_id}")
    
    rpm_parm = vsp.FindParm(disk_id, "RotorRPM", "Rotor")
    print(f"RotorRPM Parm ID: {rpm_parm}")
    
    if rpm_parm:
        print(f"Current Value: {vsp.GetParmVal(rpm_parm)}")
        vsp.SetParmValUpdate(rpm_parm, 5555.0)
        print(f"New Value after SetParmValUpdate: {vsp.GetParmVal(rpm_parm)}")
        
        # Check if it persists after a global Update
        vsp.Update()
        print(f"Value after global Update: {vsp.GetParmVal(rpm_parm)}")
    else:
        print("FAILED to find RotorRPM parameter")
else:
    print("FAILED to find Actuator Disk")
