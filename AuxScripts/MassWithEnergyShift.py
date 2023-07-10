import numpy as np
import sys
sys.path.append('C:\\work\\Science\\BINP\\Kaon Mass Measure')
from AuxScripts.EnergyShift import Shift

    
def GetMass(delta_E: float | None = None):
    M_RCNC = np.array([497.577, 497.568, 497.533, 497.554, 497.544, 497.57, 497.589, 497.568, 497.604])
    err_M = np.array([0.036, 0.014, 0.011, 0.008, 0.007, 0.01, 0.013, 0.014, 0.029])
    # err_M = np.array([0.028, 0.01, 0.004, 0.004, 0.002, 0.002, 0.004, 0.007, 0.027])
    RC_data = np.array([0.1, 0.099, 0.089, 0.08, 0.071, 0.116, 0.191, 0.334, 1.453])
    E = np.array([504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 513.864])
    M_vis = np.round(M_RCNC + RC_data, 3)

    shift = Shift(M_vis, err_M, RC_data, E, delta_E)

    delta_E, deltaE_err = shift.Get_Shift()

    M_fit, err = shift.Get_Mass_fit()
    print("Delta E =", delta_E, "+/-", deltaE_err, "MeV")
    print("M_fit =", M_fit, "+/-", err, "MeV")
    print("M_new =", list(shift.Get_corrected_mass()))
    
    print("RC_shifted =", list(np.round(shift.Get_corrected_RC(), 4)))