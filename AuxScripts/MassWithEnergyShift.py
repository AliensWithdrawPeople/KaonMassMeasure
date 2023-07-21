import numpy as np
import sys
sys.path.append('C:\\work\\Science\\BINP\\Kaon Mass Measure')
from AuxScripts.EnergyShift import Shift

    
def GetMass(M_RCNC: np.ndarray, err_M: np.ndarray, RC_data: np.ndarray, E: np.ndarray, delta_E: float | None = None):
    M_vis = np.round(M_RCNC + RC_data, 3)

    shift = Shift(M_vis, err_M, RC_data, E, delta_E)

    delta_E, deltaE_err = shift.Get_Shift()

    M_fit, err = shift.Get_Mass_fit()
    print("Delta E =", delta_E, "+/-", deltaE_err, "MeV")
    print("M_fit =", M_fit, "+/-", err, "MeV")
    print("M_new =", list(shift.Get_corrected_mass()))
    
    print("RC_shifted =", list(np.round(shift.Get_corrected_RC(), 4)))