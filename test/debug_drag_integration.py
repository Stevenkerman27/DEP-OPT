import numpy as np
import sys
import os

# Add parent directory to sys.path to ensure infrastructure can be imported
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.append(os.getcwd())

import infrastructure as opb

def main():
    # Parameters from evaluate.py (base config)
    Mean_chord = 0.197
    taper = 0.59
    span = 1
    ThickChord = 0.12
    V = 15.0 # typ_speed
    
    # Wing Area S = Mean_chord * span (for rectangular wing as in evaluate.py base config)
    wing_S = Mean_chord * span 
    
    print(f"Integration Check for evaluate.py parameters:")
    print(f"  Velocity: {V} m/s")
    print(f"  Mean Chord: {Mean_chord} m")
    print(f"  Span: {span} m")
    print(f"  S_ref: {wing_S:.4f} m^2")
    print(f"  t/c: {ThickChord}")
    
    cdp, dp = opb.calculate_parasite_drag(V, Mean_chord, wing_S, ThickChord)
    print(f"Calculated CDp: {cdp:.6f}, Dp: {dp:.6f} N")

if __name__ == "__main__":
    main()
