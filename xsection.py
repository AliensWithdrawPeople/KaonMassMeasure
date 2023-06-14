import numpy as np
from scipy import integrate
from scipy.interpolate import make_interp_spline, BSpline 

class Xsection:
    def __init__(self, energies: np.ndarray, xsection: np.ndarray) -> None:
        self.energies = energies
        self.xsection = xsection
        self._xsection_spline = make_interp_spline(self.energies, self.xsection)
        self._energy_range = np.linspace(500, 520, 1000)
        self._vis_xsection = self._xsection_spline(self._energy_range)
        print(self._energy_range)
        pass
    
    def Get_xsection(self, energy: float)->float:
        return float(np.interp(energy, self._energy_range, self._vis_xsection))