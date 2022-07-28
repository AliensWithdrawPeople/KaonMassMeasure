from typing import List
import numpy as np
from scipy import integrate
from scipy.interpolate import make_interp_spline
from scipy.misc import derivative
from scipy.optimize import root_scalar, RootResults

pi: float = np.pi
alpha: float = 1 / 137.035999
mK: float = 497.614  
"""Kaon mass (MeV)"""
me: float = 0.511 
"""Electron mass (MeV)"""
mPi: float = 139.570
"""Charged pion mass (Mev)"""

def massFunc(s: float, psi: float)->float:
    eta: float = (1 - 0.9999*0.9999) / (1 + 0.9999*0.9999)
    # pions' beta ^2
    b: float = 1 - 4 * mPi**2 / (s / 4)
    return (s / 4 * (1 - (1 + (1 - eta**2)**0.5 * np.cos(psi)) * (1 - (1 - eta**2 * b)**0.5 )/ eta / eta ))**0.5

def massNC(s: float, psi: float, sigmaPsi: float)->float:
    return massFunc(s, psi) - sigmaPsi**2 / 2 * derivative(lambda x: massFunc(s, x), psi, 1e-5, 2, order=9)

def FormFactor(s: float)->complex:
    MPhi = 1.01919e+03
    WPhi = 4.19280e+00
    Smax = 2.83674e-02
    par3 = -2.10120e-02
    par4 = -1.95319e-02
    par5 = 4.70921e-03
    par6 = 1.49144e+02
    par7 = 1.48529e-02
    par8 = 1.25388e-02
    par9 = -3.89597e-03
    #double Mks = 497.614;

    Maphi1680: float = 1680
    Mrho2150: float = 2150
    Momega1420: float = 1425
    Momega1650: float = 1670
    WPhi1680: float = 170.
    Womega1420: float = 400.
    Womega1650: float = par6
    Warho2150: float = 500.
    aphi: complex = MPhi*MPhi/complex(s-MPhi*MPhi,np.sqrt(s) * WPhi)
    aphi1680: complex = Maphi1680*Maphi1680/complex(s-Maphi1680*Maphi1680,np.sqrt(s)*WPhi1680)
    aomega1420: complex = Momega1420*Momega1420/complex(s-Momega1420*Momega1420, -np.sqrt(s)*Womega1420)
    aomega1650: complex = Momega1650*Momega1650/complex(Momega1650*Momega1650-s, -np.sqrt(s)*Womega1650)
    arho2150: complex = Mrho2150*Mrho2150/complex(Mrho2150*Mrho2150-s, -np.sqrt(s)*Warho2150)
    ATot: complex = Smax*aphi + complex(par3,par4) + par5*aomega1420 + par7*aomega1650 + par8*aphi1680 - par9*arho2150
    return np.sqrt(1 / alpha)*ATot

def SigmaBorn(s:float)->float:
    # Old version (incorrect)
    #return 204.5 * alpha * alpha * pow(1 - 4 * mK * mK / s, 3./2.) * 2 / 3 * pi / s * pow(abs(FormFactor(s)), 2.)
    energies: List = [1004.066, 1010.466, 1012.955, 1015.068, 1016.105, 1017.155, 1017.156, 1018.046, 1019.118,  1019.214, 1019.421, 1019.902, 
                    1021.222, 1021.309, 1022.078,  1022.744, 1023.264,  1025.320, 1027.956, 1029.090, 1033.907, 1040.028,  1049.864,  1050.862,  1059.947]
    cross_sections: List = [6.87, 42.16, 96.74, 219.53, 366.33, 628.15, 624.76, 996.62, 1413.65,  1433.05, 1434.84, 1341.91,  833.20, 
                            807.54, 582.93,  443.71, 377.77,  199.26, 115.93, 96.96, 50.12, 31.27,  16.93,  17.47,  12.09]
    data = np.array([energies, cross_sections]).T 
    data2 = data[data[:,0].argsort()]
    spline = make_interp_spline(data2[:, 0], data2[:, 1], k=2)
    x = np.linspace(mK-10, max(energies), 10000)
    y = spline(x)
    y = np.where(x<=2*mK, 0, y)
    return np.interp(s, (x**2), y)

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

def epsCalc(s: float, massUpperLimit: float)->float:
    res: RootResults = root_scalar(lambda psi: massFunc(s, psi) - massUpperLimit, xtol=1e-5, x0=2, x1=3)
    print("Root finding is", str(res.converged) + "; number of iteration =", res.function_calls)
    psiCr: float = res.root
    enLim: float = (mK**2 - 4 * (mPi * np.cos(psiCr / 2))**2)**0.5 / np.sin(psiCr / 2)
    if res.converged:
        print("Energy limit =", enLim)
    # (2*enLim - 2 * mK) / mK
    return 1 - 4 * enLim**2 / s if res.converged else 1

def SigmaCorrected(s:float, eps=1)->float:
    return integrate.quad(lambda x: SigmaBorn(s * (1-x)) * F(x, s), 0, eps, epsabs = 1e-6, epsrel=1e-4, limit=500)[0]

def GetMassCorrected(s: float, psi: float, sigmaPsi: float, eps: float)->float:
    """ Compute mass of Ks with radiative and nonlinearity correction.

    Args:
        s (float): Mandelstam variable, s = E_cm ^ 2
        psi (float): critical angle
        sigmaPsi (float): standard deviation of critical angle

    Returns:
        float: corrected mass of Ks
    """        
    conv: float = integrate.quad(lambda x: massNC(s * (1-x), psi, sigmaPsi) * SigmaBorn(s * (1-x)) * F(x, s), 0, eps, epsabs = 1e-6, epsrel=1e-4, limit=500)[0]
    return conv / SigmaCorrected(s, eps)

energy: float =  2 * 509
s: float = energy**2
psi: float =   2.63568
sigmaPsi: float = 0.0164407
massUpperLimit: float = 505
maxPhotonEnergy: float = 10 

# eps: float = epsCalc(s, massUpperLimit)
# eps: float = 2 * maxPhotonEnergy / energy
# print("Corrected Mass = ", GetMassCorrected(s, psi, sigmaPsi, eps), "eps =", eps)
# print("Mass = ", massFunc(s, psi) )

# {2.73353, 2.65747, 2.63507, 2.61432, 2.5988, 2.56547}
# {505, 508, 509, 510, 511, 514}
# {497.652, 497.626, 497.640, 497.628, 497.668, 498.036}
# {0.004, 0.005, 0.005, 0.004, 0.006, 0.013}