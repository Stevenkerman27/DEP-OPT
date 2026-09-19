# Parasitic Drag Estimation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement a physics-based parasitic drag estimation function in `infrastructure.py` and verify it with tests.

**Architecture:** Add a standalone calculation function to `infrastructure.py` that utilizes global constants (`density`, `mu`) and provides a Raymer/Prandtl-Schlichting based drag estimate.

**Tech Stack:** Python, NumPy.

---

### Task 1: Add Drag Calculation Function to Infrastructure

**Files:**
- Modify: `infrastructure.py`
- Create: `test/test_drag_calculation.py`

- [ ] **Step 1: Write the failing test**
Create `test/test_drag_calculation.py` with a test case for a standard wing.

```python
import infrastructure as opb
import numpy as np

def test_calculate_parasite_drag():
    # Standard case: V=50m/s, chord=0.2m, S_ref=0.5m^2, t/c=0.12
    V = 50.0
    chord = 0.2
    S_ref = 0.5
    t_over_c = 0.12
    
    cdp = opb.calculate_parasite_drag(V, chord, S_ref, t_over_c)
    
    # Re = 1.225 * 50 * 0.2 / 1.85e-5 = 662162
    # log10(Re) = 5.8209
    # Cf = 0.455 / (5.8209^2.58) = 0.455 / 93.35 = 0.004874
    # FF = 1 + 1.2*0.12 + 100*(0.12^4) = 1.144 + 0.020736 = 1.1647
    # Swet = 2.003 * 0.5 * (1 + 0.25*0.12) = 1.0015 * 1.03 = 1.0315
    # CDp = 0.004874 * 1.1647 * 1.0315 / 0.5 = 0.0117
    
    assert np.isclose(cdp, 0.0117, rtol=1e-2)
```

- [ ] **Step 2: Run test to verify it fails**
Run: `python test/test_drag_calculation.py`
Expected: `AttributeError: module 'infrastructure' has no attribute 'calculate_parasite_drag'`

- [ ] **Step 3: Implement the function in `infrastructure.py`**
Add the function to the end of `infrastructure.py`.

```python
def calculate_parasite_drag(V, chord, S_ref, t_over_c):
    """
    Calculates the parasitic drag coefficient (CDp) for a wing.
    
    Args:
        V (float): Velocity [m/s]
        chord (float): Reference chord (MAC) [m]
        S_ref (float): Reference area [m^2]
        t_over_c (float): Thickness-to-chord ratio (e.g., 0.12)
        
    Returns:
        float: Parasitic drag coefficient CDp
    """
    if V <= 0:
        return 0.0
        
    # Re = rho * V * c / mu
    Re = (density * V * chord) / mu
    
    # Cf = 0.455 / (log10(Re)^2.58) (Prandtl-Schlichting)
    cf = 0.455 / (np.log10(Re)**2.58)
    
    # FF = 1 + 1.2(t/c) + 100(t/c)^4 (Raymer Wing)
    ff = 1.0 + 1.2 * t_over_c + 100.0 * (t_over_c**4)
    
    # Swet = 2.003 * Sref * (1 + 0.25 * t/c) (Raymer)
    s_wet = 2.003 * S_ref * (1.0 + 0.25 * t_over_c)
    
    # CDp = (Cf * FF * Q * Swet) / Sref, assuming Q=1.0
    cdp = (cf * ff * s_wet) / S_ref
    
    return cdp
```

- [ ] **Step 4: Run test to verify it passes**
Run: `python test/test_drag_calculation.py`
Expected: `PASS` (or no output if using plain asserts and it succeeds)

- [ ] **Step 5: Commit changes**
```bash
git add infrastructure.py test/test_drag_calculation.py
git commit -m "feat: add parasitic drag estimation function to infrastructure"
```

---

### Task 2: Integration Example in Evaluate (Documentation/Verification)

**Files:**
- Create: `test/debug_drag_integration.py`

- [ ] **Step 1: Create a script to demonstrate integration with evaluate.py parameters**
This script simulates how `evaluate.py` would use the new function.

```python
import infrastructure as opb

# Parameters from evaluate.py (base config)
Mean_chord = 0.24
taper = 1
span = 0.97
ThickChord = 0.12
V = 12.0 # typ_speed
wing_S = Mean_chord * span # Approximation for rectangular wing

cdp = opb.calculate_parasite_drag(V, Mean_chord, wing_S, ThickChord)
print(f"Calculated CDp for base config at {V} m/s: {cdp:.6f}")

# Check if it's within expected range for this aircraft scale
assert 0.005 < cdp < 0.030
```

- [ ] **Step 2: Run the script**
Run: `python test/debug_drag_integration.py`
Expected: Output showing the calculated CDp and successful assertion.

- [ ] **Step 3: Commit**
```bash
git add test/debug_drag_integration.py
git commit -m "test: add integration debug script for parasitic drag"
```
