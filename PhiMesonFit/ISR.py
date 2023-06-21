import numpy as np
from scipy.interpolate import make_interp_spline 
from scipy import integrate

from xsection import Xsection


pi: float = np.pi
alpha: float = 1 / 137.035999
mK: float = 497.614  
"""Kaon mass (MeV)"""
me: float = 0.511 
"""Electron mass (MeV)"""
mPi: float = 139.570
"""Charged pion mass (Mev)"""

   
energies: list = [1004.066, 1010.466, 1012.955, 1015.068, 1016.105, 1017.155, 1017.156, 1018.046, 1019.118,  1019.214, 1019.421, 1019.902, 
                1021.222, 1021.309, 1022.078,  1022.744, 1023.264,  1025.320, 1027.956, 1029.090, 1033.907, 1040.028,  1049.864,  1050.862,  1059.947]
cross_sections: list = [6.87, 42.16, 96.74, 219.53, 366.33, 628.15, 624.76, 996.62, 1413.65,  1433.05, 1434.84, 1341.91,  833.20, 
                        807.54, 582.93,  443.71, 377.77,  199.26, 115.93, 96.96, 50.12, 31.27,  16.93,  17.47,  12.09]
xsec = Xsection(np.array(energies), np.array(cross_sections), (2 * mK, max(energies)))
    
# Taken from Kozyrev's article (https://doi.org/10.48550/arXiv.1604.02981)
SigmaBornExternal = xsec.Get_xsection

class Solver:
    def __init__(self, external_xsec, energies: np.ndarray, vis_xsec: np.ndarray, spread: np.ndarray) -> None:
        self.ens = energies
        self.spread = spread
        self.ext = external_xsec
        self.vis = vis_xsec
        self.born = external_xsec
                
        self.RC = np.array(len(energies) * [1])
        """ eps^i(s) * (1 + delta^i(s))"""
        
        self.__calculate__()
        pass
    
    def __calculate__(self) -> None:
        prev_RC = 0
        self.RC = self.step()
        counter = 0
        while(max(np.abs(np.subtract(self.RC, prev_RC) )) > 1e-3):
            if(counter > 100):
                raise OverflowError("There was a problem with ISR correction calculation: Too many iterations!")
            prev_RC = self.RC
            self.RC = self.step()
            counter = counter + 1
        
    def Get_ISR_correction(self)->np.ndarray:
        return self.RC
    
    def step(self) -> np.ndarray:
        RC = np.array([self.SigmaVis(energy=E, spread=sigma_E) / self.born(E) for (E, sigma_E) in list(zip(self.ens, self.spread))])
        self.born = Xsection(self.ens, self.vis / self.RC, (2 * mK, max(self.ens))).Get_xsection
        return RC
        
    
    def SigmaVis(self, energy: float, spread: float):
        return 1 / np.sqrt(2 * pi * spread**2) * integrate.quad(lambda E: np.exp(-(E - energy)**2 / (2 * spread**2)) * \
                              integrate.quad(lambda x: Solver.F(x, E*E) *  self.born(E * np.sqrt((1 - x)) ), 0, 1 - 2 * mK * mK / E / E), -np.inf, np.inf)

    @staticmethod
    def F(x: float, s: float)->float:
        """ F is radiator (see Fadin-Kuraev and Achasov's presentation)

        Args:
            x (float): variable of integration in integral for corrected cross-section; For x << 1 x is energy carried away by photons and pairs.
            s (float): Mandelstam variable, s = E_cm ^ 2

        Returns:
            float: value of F at x for certain energy. 
        """
        # L = ln(s / m_e^2)
        L: float = np.log(s / me / me)
        # beta = 2*alpha/pi * (L-1)
        beta: float = 2 * alpha / pi *(L - 1)

        p1: float = beta * x**(beta - 1.) * (1 + alpha / pi * (pi**2 / 3 - 1./2.) + 3./4. * beta - 1./24. * beta**2 * (L / 3. + 2 * pi ** 2 - 37./4.) )
        p2: float = beta * (1. - x /2.)

        p3: float = 4. * (2. - x) * np.log(1. / x)
        p4: float = 1. / x * (1. + 3. * (1.-x)**2) * np.log(1./(1.-x))

        # Photon part
        phPart: float = p1  - p2 + beta * beta / 8. * (p3 + p4 - 6. + x)

        E: float = np.sqrt(s / 4.)

        # Real e+e- part.
        if(x - 2 * me / E > 0):
            p5: float = 1 / 6. / x * (x - 2 * me / E)**beta
            p6: float = (np.log(s * x * x / me / me) - 5./3.)**2.
            p7: float = beta/3.*(np.log(s * x**2 / (me**2)) - 5./3)
            p8: float = L*L/2. * (2./3.*(1-(1-x)**3) / (1.-x) - (2.-x)*np.log(1./(1.-x)) + x/2.)

            elPart: float = (alpha / pi)**2 * (p5 * p6 * (2 - 2*x + x*x + p7) + p8)
        else:
            elPart: float = 0

        return phPart + elPart
    
def SigmaBorn(s: float, m_phi: float, gamma_phi: float, B_phi_gamma_phi_ee: float)->float:
    m_rho = 775.26
    m_omega = 782.66
    gamma_rho = 1 # todo
    gamma_omega = 1 # todo
    # B_VKK
    B_rho = 1 # todo
    B_omega = 1 # todo
    # K^0 momentum
    def p_K(s: float)->float: 
        return np.emath.sqrt(s - 4 * 497.614**2)
    
    def D(s: float, m_V: float, Gamma_V)->complex:
        return m_V**2 - s - 1j * Gamma_V(s)
    
    def Gamma_rho(s: float):
        return gamma_rho * m
    
    
    g_phi_gamma = np.emath.sqrt(3 * m_phi**3 * B_phi_gamma_phi_ee/ 4 / pi / alpha) 
    g_phi_KK = np.emath.sqrt(6 * pi * m_phi**2 * gamma_phi / p_K(m_omega**2)**3)
    g_phi = g_phi_gamma * g_phi_KK
    
    g_rho = np.emath.sqrt(3 * m_rho**3 * 7.04 / 4 / pi / alpha) * (g_phi_KK / np.sqrt(2))
    g_omega = np.emath.sqrt(3 * m_omega**3 * 0.60 / 4 / pi / alpha) * (-g_phi_KK / np.sqrt(2))
    
    rho_part = g_rho / D(s, m_rho, Gamma_rho)