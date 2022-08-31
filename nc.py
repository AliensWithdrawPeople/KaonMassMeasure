import numpy as np
from scipy import integrate
from scipy.misc import derivative


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

energy: float =  2 * 504.895
s: float = energy**2
psi: float = 2.73226
sigmaPsi: float = 0.0164407
print("deltaM = ", massNC(s, psi, sigmaPsi))