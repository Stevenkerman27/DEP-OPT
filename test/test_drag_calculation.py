import numpy as np
import sys
import os

# Add parent directory to sys.path to ensure infrastructure can be imported
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.append(os.getcwd())

import infrastructure as opb

def test_calculate_parasite_drag():
    # Standard case: V=50m/s, chord=0.2m, S_ref=0.5m^2, t/c=0.12
    V = 50.0
    chord = 0.2
    S_ref = 0.5
    t_over_c = 0.12
    
    try:
        cdp, dp = opb.calculate_parasite_drag(V, chord, S_ref, t_over_c)
    except AttributeError:
        print("Failing as expected: calculate_parasite_drag not found")
        return
    except TypeError:
        print("Return type mismatch: expected tuple (cdp, dp)")
        return

    # Theoretical values for check:
    # Re = 1.225 * 50 * 0.2 / 1.85e-5 = 662162
    # log10(Re) = 5.8209
    # Cf = 0.455 / (5.8209^2.58) = 0.455 / 93.35 = 0.004874
    # FF = 1 + 1.2*0.12 + 100*(0.12^4) = 1.144 + 0.020736 = 1.1647
    # Swet = 2.003 * 0.5 * (1 + 0.25*0.12) = 1.0015 * 1.03 = 1.0315
    # CDp = 0.004874 * 1.1647 * 1.0315 / 0.5 = 0.0117
    # Dp = 0.5 * 1.225 * 50^2 * 0.5 * 0.0117 = 8.953125
    
    print(f"Calculated CDp: {cdp}, Dp: {dp} N")
    assert np.isclose(cdp, 0.0117, rtol=2e-2)
    assert np.isclose(dp, 8.95, rtol=2e-2)
    print("Test Passed!")

if __name__ == "__main__":
    test_calculate_parasite_drag()
