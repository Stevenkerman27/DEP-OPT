# Design Doc: Parasitic Drag Estimation for Wings

## 1. Overview
Implement a physics-based parasitic drag estimation module for wings using skin friction ($C_f$) and form factor ($FF$) methods. This replaces hardcoded drag constants with calculations based on wing geometry and flight conditions.

## 2. Requirements
- **Formula for $C_f$:** Prandtl-Schlichting (Turbulent Flow).
- **Formula for $FF$:** Raymer's Wing Form Factor.
- **Form Factor Inputs:** Thickness-to-chord ratio ($t/c$), ignoring sweep and twist.
- **Integration:** Must be accessible from `evaluate.py` via `infrastructure.py`.
- **Constants:** Use global `density` and `mu` from `infrastructure.py`.

## 3. Architecture
The implementation will be housed in `infrastructure.py` to maintain a single source of truth for aerodynamic calculations.

### 3.1 Function Specification: `calculate_parasite_drag`
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
```

### 3.2 Mathematical Formulas
1. **Reynolds Number ($Re$):**
   $$Re = \frac{\rho \cdot V \cdot \text{chord}}{\mu}$$
   *Variables $\rho$ and $\mu$ are pulled from `infrastructure.py` globals.*

2. **Skin Friction Coefficient ($C_f$):**
   $$C_f = \frac{0.455}{(\log_{10} Re)^{2.58}}$$

3. **Form Factor ($FF$):**
   $$FF = 1 + 1.2(t/c) + 100(t/c)^4$$

4. **Wetted Area ($S_{wet}$):**
   $$S_{wet} = 2.003 \cdot S_{ref} \cdot (1 + 0.25 \cdot t/c)$$

5. **Parasitic Drag Coefficient ($C_{D,p}$):**
   $$C_{D,p} = \frac{C_f \cdot FF \cdot S_{wet}}{S_{ref}}$$

## 4. Implementation Plan
1. **Infrastructure Update:** Add `calculate_parasite_drag` to `infrastructure.py`.
2. **Global Constants:** Ensure `density` and `mu` are correctly referenced.
3. **Verification:** Create `test/test_drag_calculation.py` to validate results.

## 5. Testing & Validation
- **Case 1:** $Re \approx 10^6, t/c = 0.12, S_{ref} = 1.0$. Expected $C_{D,p} \approx 0.008 - 0.012$.
- **Integration Test:** Verify `evaluate.py` can call the function and receive a valid float.
