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

def responseFunction(psi: float, avgPsi: float, sigma: float)->float:
    return np.exp(-(psi - avgPsi)**2 / (2 * sigma**2)) / (2 * pi)**0.5 / sigma

def massNC(s: float, psi: float, sigmaPsi: list)->float:
    return integrate.quad(lambda x: massFunc(s, x) * responseFunction(x, psi, sigmaPsi) , 0, pi, epsabs = 1e-6, epsrel=1e-4, limit=500)[0]

energy: float =  2 * 510
s: float = energy**2
psi: float = 2.61185  
sigmaPsi: float = 0.0164407
print("deltaM = ", massFunc(s, psi) - massNC(s, psi, sigmaPsi))