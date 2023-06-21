import numpy as np
from scipy.interpolate import make_interp_spline, BSpline 

class Xsection:
    def __init__(self, energies: np.ndarray, xsection: np.ndarray, interpolation_range: tuple[float, float]) -> None:
        self.energies = energies
        self.xsection = xsection
        
        data = np.array([energies, xsection]).T 
        data2 = data[data[:,0].argsort()]
        self._xsection_spline = make_interp_spline(data2[:, 0], data2[:, 1], k=2)
        
        self._energy_range = np.linspace(interpolation_range[0], interpolation_range[1], 10000)
        self._vis_xsection = self._xsection_spline(self._energy_range)
        pass
    
    def Get_xsection(self, energy: float)->float:
        return float(np.interp(energy, self._energy_range, self._vis_xsection))
    
    def Get_Data(self):
        return self._energy_range, self._vis_xsection