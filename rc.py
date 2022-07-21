import numpy as np
from scipy import integrate

pi:float = np.pi
alpha:float = 1 / 137.035999
kaonMass:float = 497.614  
electronMass:float = 0.511
a:float = pi * pi / 6 - 1 / 4

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

    Maphi1680:float = 1680
    Mrho2150:float = 2150
    Momega1420:float = 1425
    Momega1650:float = 1670
    WPhi1680:float = 170.
    Womega1420:float = 400.
    Womega1650:float = par6
    Warho2150:float = 500.
    aphi:complex = MPhi*MPhi/complex(s-MPhi*MPhi,np.sqrt(s) * WPhi)
    aphi1680:complex = Maphi1680*Maphi1680/complex(s-Maphi1680*Maphi1680,np.sqrt(s)*WPhi1680)
    aomega1420:complex = Momega1420*Momega1420/complex(s-Momega1420*Momega1420, -np.sqrt(s)*Womega1420)
    aomega1650:complex = Momega1650*Momega1650/complex(Momega1650*Momega1650-s, -np.sqrt(s)*Womega1650)
    arho2150:complex = Mrho2150*Mrho2150/complex(Mrho2150*Mrho2150-s, -np.sqrt(s)*Warho2150)
    ATot:complex = Smax*aphi + complex(par3,par4) + par5*aomega1420 + par7*aomega1650 + par8*aphi1680 - par9*arho2150
    return np.sqrt(1 / alpha)*ATot

def SigmaBorn(s:float)->float:
    return 204.5 * alpha * alpha * pow(1 - 4 * kaonMass * kaonMass / s, 3./2.) * 2 / 3 * pi / s * pow(abs(FormFactor(s)), 2.)

def F(x: float, s: float)->float:
    # L = ln(s / m_e^2)
    L:float = np.log(s / electronMass / electronMass)
    # beta = 2*alpha/pi * (L-1)
    beta:float = 2 * alpha / pi *(L - 1)
    # Brackets before x^(beta - 1) in F.
    par1:float = 1 + alpha / pi * (pi * pi / 3 - 1./2.) + 3./4. * beta - 1./24. * beta * beta *(L / 3. + 2 * pi * pi - 37./4.)

    f:float = beta * pow(x, beta - 1.) * par1  - beta * (1. - x /2.) + beta * beta / 8. * (4. * (2. - x) * np.log(1. / x) + 1. / x * (1. + 3. * (1.-x)*(1.-x)) * np.log(1./(1.-x)) - 6. + x)

    E:float = np.sqrt(s / 4.)
    # In case of real e+e- pairs are not banned.
    epemPart:float = 0
    if(x - 2 * electronMass / E > 0): 
        epemPart =  alpha * alpha / pi / pi * \
        ( 1 / 6. / x * pow(x - 2 * electronMass / E, beta) * pow(np.log(s * x * x / electronMass / electronMass) -5./3. , 2.) * \
        (2 - 2*x + x*x + beta/3.*(np.log(s * x * x / electronMass /electronMass) -5./3)) + \
        L*L/2. * (2./3.*(1-pow(1-x, 3)) / (1.-x) - (2.-x)*np.log(1./(1.-x)) +x/2.) )
    return f + epemPart

def SigmaCorrected(s:float)->float:
    sigma:float = integrate.quad(lambda x: pow(x, -0.85), 0, 0.04)[0]
    return sigma

print(SigmaCorrected(1010*1010))