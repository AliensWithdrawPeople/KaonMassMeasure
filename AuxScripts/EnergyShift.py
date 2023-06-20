import numpy as np
from scipy.interpolate import make_interp_spline
from scipy import optimize as opt

def interpolate(x: np.ndarray, y: np.ndarray, energy_range: tuple) -> tuple[np.ndarray, np.ndarray]:
    data = np.array([x, y]).T
    data2 = data[data[:, 0].argsort()]
    spline = make_interp_spline(data2[:, 0], data2[:, 1], k=2)

    xs = np.linspace(energy_range[0], energy_range[1], 10000)
    return xs, spline(xs)

class Shift:
    def __init__(self, M_vis: np.ndarray, err_M: np.ndarray, RC_data: np.ndarray, E: np.ndarray) -> None:
        self._M_vis = M_vis
        self._err_M = err_M
        self._RC_data = RC_data
        self._E = E
        
        _, self._RC_interp = interpolate(E, RC_data, (504, 515))
        self._energy_range, self._M_vis_interp = interpolate(E, M_vis, (504, 515))
        
        self._fit_res, self._fit_err = opt.curve_fit(lambda x, mean, delta: mean +  self.Get_RC(x + delta), E, M_vis, p0=[497.611, 0], sigma=err_M, absolute_sigma=True)
        delta_E = np.round(self._fit_res[1], 3)
        self._M_new = self.Get_M_vis(E) - self.Get_RC(E + delta_E)
        self._M_fit, self._err = opt.curve_fit(lambda x, a: a, E, self._M_new, sigma=err_M)
        pass
    
    def Get_RC(self, energy: float | np.ndarray) -> float | np.ndarray:
        return np.interp(energy, self._energy_range, self._RC_interp)

    def Get_M_vis(self, energy: float | np.ndarray) -> float | np.ndarray:
        return np.interp(energy, self._energy_range, self._M_vis_interp)
    
    def Get_Shift(self)->tuple[float, float]:
        return (np.round(self._fit_res[1], 3), np.round(np.sqrt(np.diag(self._fit_err))[1], 3))
    
    def Get_Mass_fit(self)->tuple[float, float]:
        return (np.round(self._M_fit[0], 3), np.round(np.sqrt(np.diag(self._err))[0], 3))
    
    def Get_corrected_mass(self)->np.ndarray:
        return np.round(self._M_new, 3)