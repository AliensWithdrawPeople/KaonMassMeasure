#include "TF1.h"
#include "TF2.h"
#include "TTree.h"
#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <ctime> 
#include "TRandom.h"

class RadCor
{
private:
    TF1 *massFunc;
    TF1 *D;

    Float_t E;
    Float_t pRatio;
    Float_t psi;
    double s;
    double sigmaPsi = 1.64407e-02;

    double L;
    // Maybe: 1 / alpha = 137.071999.
    constexpr static double alpha = 1 / 137.035999;  
    constexpr static double kaonMass = 497.614;   
    constexpr static double electronMass = 0.511;
    constexpr static double pi = 3.14159265358979323846;
    // Relevant to D function.
    constexpr static double a = pi * pi / 6 - 1/4;

    /*
    Born cross-section of e+e- -> KsKl.
    For definition go to https://arxiv.org/abs/1604.02981v3 or https://doi.org/10.1142/S0217751X92001423.
    Look here: https://cmd.inp.nsk.su/websvn/filedetails.php?repname=Cmd3Sim&path=%2Ftrunk%2Fgenerator%2Fradcor%2Fsrc%2FTKnFormFactor.C
    */
    Float16_t Sigma0(const double &s);
    // Return FormFactor of e+e- -> KsKl.
    std::complex<double> FormFactor(const double &s);
    Float16_t DEval(const double &z);
public:
    RadCor(Float_t energy, Float_t ksdpsi, Float_t momentumRatio, bool isCrAngleMethod = true);
    double RadCorEval();
    ~RadCor();
};

RadCor::RadCor(Float_t energy, Float_t ksdpsi, Float_t momentumRatio, bool isCrAngleMethod = true)
{
    E = energy; pRatio = momentumRatio;
    psi = ksdpsi;

    if(isCrAngleMethod)
    { 
        massFunc = new TF1("mass1", "[0] * TMath::Sqrt(1 - (1 - 4 * 139.57 * 139.57 / [0] / [0]) * cos(x / 2) * cos(x / 2))");
    }
    else
    {
        massFunc = new TF1("MassLnY", "sqrt([0] * [0] * (1 - (1 + sqrt(1 - [1] *[1]) * cos(x))*(1 - sqrt(1 - [1] * [1] * (1 - 4 * 139.57 * 139.57 / [0] / [0])))/ [1] / [1] ))");
        massFunc->SetParameter(1, (1 - pRatio*pRatio) / (1 + pRatio*pRatio));
    }
    
    s = 4 * E * E;
    // Auxilary function 

    L = 2 * TMath::Log(2*E / electronMass);
    // par1 = b=2*alpha/pi * (L-1)
    double par1 = 2 * alpha/pi*(L - 1.);
    // par2 = 2*m_e/E
    double par2 = 2. * electronMass / E;

    /*
    Structure D-function 
    D(x,s) = D_gamma + D_e+e-. 
    For definition go to https://arxiv.org/abs/hep-ph/9703456v1.
    */
    /*
    D = new TF1("D func", [&](double *x, double *p){
        double Dgamma = par1/2. * pow((1-x[0]), par1/2.-1.) * (1. + 3.*par1/8. + par1*par1/16.*(9./8.- pi * pi/3.)) 
        - par1/4.*(1+x[0]) + par1*par1/32.*(-4.*(1+x[0])*TMath::Log(1.-x[0]) - (1+3.*x[0]*x[0])/(1.-x[0])*TMath::Log(x[0])-5.-x[0]);

        double Dee = 0;
        if(1-x[0]-par2 > 0)
        {
            Dee = par1/2*pow(1-x[0], par1/2.-1.)*(-par1*par1/288.*(2.*L - 15.)) + 
            pow(alpha / pi, 2.) * ( 1./12./(1.-x[0]) * pow(1.-x[0]-par2, par1/2.)*pow(2. * TMath::Log((1.-x[0]) / par2) - 5./3., 2.) *
            (1.+x[0] * x[0] + par1/6.*(2. * TMath::Log((1.- x[0]) / par2) - 5./3.)) + L*L/4.*(2./3.*(1-pow(x[0], 3.)) / x[0] + (1.-x[0]) / 2. +
            (1. + x[0])*TMath::Log(x[0])) );
        }
        return Dgamma + Dee;}, 1e-6, 1-1e-6, 1);
    */
}

RadCor::~RadCor()
{
    delete massFunc;
    delete D;
}

Float16_t RadCor::DEval(const double &z)
{
    Float16_t x = 1 - z;
    Float16_t b = 2 * alpha/pi*(L - 1.);
    Float16_t b2 = 0.5 * b;
    Float16_t D0   = 1. + 3./8.* b + b * b /16. * (9./8. - pi * pi / 3.);

    Float16_t D = b2*pow(x,b2-1) * D0 - 0.5*b2*(1+z) - b2 * b2 / 8 * (4*( 1 + z)*log(x)+(1+3*z*z)/x * log1p(-x) + 5 + z);
    return D;
}

std::complex<double> RadCor::FormFactor(const double &s)
{
    double MPhi = 1.01919e+03;
    double WPhi = 4.19280e+00;
    double Smax = 2.83674e-02;
    double par3 = -2.10120e-02;
    double par4 = -1.95319e-02;
    double par5 = 4.70921e-03;
    double par6 = 1.49144e+02;
    double par7 = 1.48529e-02;
    double par8 = 1.25388e-02;
    double par9 = -3.89597e-03;
    //double Mks = 497.614;

    double Maphi1680 = 1680;
    //double Mrho =775.49;
    //double Mrho1450 = 1465;
    //double Mrho1450b = 1450;
    //double Mrho1570 = 1570;
    //double Mrho1700 = 1720;
    //double Mrho1900 = 1900;
    double Mrho2150 = 2150;
    //double Momega = 782.65;
    double Momega1420 = 1425;
    double Momega1650 = 1670;
    double WPhi1680 = 170.;
    double Womega1420 = 400.;
    double Womega1650 = par6;
    double Warho2150 = 500.;
    std::complex<double> aphi = MPhi*MPhi/std::complex<double>(s-MPhi*MPhi,sqrt(s)*WPhi);
    std::complex<double> aphi1680 = Maphi1680*Maphi1680/std::complex<double>(s-Maphi1680*Maphi1680,sqrt(s)*WPhi1680);
    std::complex<double> aomega1420 = Momega1420*Momega1420/std::complex<double>(s-Momega1420*Momega1420, -sqrt(s)*Womega1420);
    std::complex<double> aomega1650 = Momega1650*Momega1650/std::complex<double>(Momega1650*Momega1650-s, -sqrt(s)*Womega1650);
    std::complex<double> arho2150 = Mrho2150*Mrho2150/std::complex<double>(Mrho2150*Mrho2150-s, -sqrt(s)*Warho2150);
    std::complex<double> ATot = Smax*aphi + std::complex<double> (par3,par4) + par5*aomega1420 + par7*aomega1650 + par8*aphi1680 - par9*arho2150 ;

    return sqrt(1 / alpha)*ATot;// maybe alpha^-1 = 137.071999
}

Float16_t RadCor::Sigma0(const double &s)
{
    auto ff = FormFactor(s);
    return 204.5 * alpha * alpha * pow(1 - 4 * kaonMass * kaonMass / s, 3./2.) * 2 / 3 * TMath::Pi() / s * abs(ff * std::conj(ff)) * 1e6;
}

double RadCor::RadCorEval()
{
    double beta = sqrt(1 - kaonMass * kaonMass / E / E);
    double b = -1 + 0.5*(1-beta)/beta*L + 1/beta*TMath::Log((1+beta)/2) + 
                0.5*(1+beta*beta)/beta * (-TMath::DiLog(-(1-beta)/(1+beta)) + TMath::DiLog((1-beta)/(1+beta)) - pi*pi/12 + 
                L*TMath::Log((1+beta)/2)-2*L*TMath::Log(beta) + 1.5*pow(TMath::Log((1+beta)/2), 2) -0.5*pow(TMath::Log(beta), 2) - 
                3*TMath::Log(beta)*TMath::Log((1+beta)/2) + L + 2*TMath::Log((1+beta)/2) );

    /* 
    Krc(s, x1, x2);
    p[0] -  auxilary parameter: p[0] == 0 stands for normalization factor, 
                                p[0] == 1 stands for convolution of mass func and Krc.  
    For definition go to https://arxiv.org/abs/hep-ph/9703456v1.
    */
    TF2 Krc("K RadCor", [&](double* x, double* p) { 
        massFunc->SetParameter(0, E * sqrt((1.-x[0]) * (1.-x[1])) );
        if(1 - 4 * kaonMass * kaonMass / ( s*(1.-x[0]) * (1.-x[1]) ) >= 0)
        { 
            return (fabs(p[0] - 1) < 0.1 ? (massFunc->Eval(psi)- sigmaPsi * sigmaPsi / 2 * massFunc->Derivative2(psi) ) * 1.004836  : 1) * // (massFunc->Eval(psi)- sigmaPsi * sigmaPsi / 2 * massFunc->Derivative2(psi) )
            (1 + 2 * alpha / pi * (1 + a + b)) * DEval(x[0]) * DEval(x[1]) * Sigma0(s*(1-x[0]) * (1-x[1]));
        }
        else
        { return 0.0; }
        }, 0, 1, 0, 1, 1);

    Krc.SetParameter(0, 0);

    double xmin = 1e-10;//0.017;
    double xmax = 1 - 4*kaonMass*kaonMass/s; //0.05;
    Double_t N = Krc.Integral(xmin, xmax, xmin, xmax, 1e-10);
    Krc.SetParameter(0, 1);
    double rc = Krc.Integral(xmin, xmax, xmin, xmax, 1e-10)/N;
    massFunc->SetParameter(0, E);
    return rc ; 
}


int radCor()
{
    auto rc = new RadCor(510, 2.61516, 1., true);
    // FullRec 2body 497.602 +- 0.003 MeV
    // CrAngle 2body 497.623 +- 0.007 MeV
    // FullRec MCGPJ 497.724 +- 0.003 MeV 
    // CrAngle MCGPJ 497.726 +- 0.006 MeV psi = 2.61516
    // Exp 497.605 +- 0.004 psi = 2.62255
    
    // CrAngle 497.601 +- 0.004 MeV
    // FullRec 497.603 +- 0.003 MeV

    // Toy MC
    // E =      514 (0.018)         511 (0.0155)        510 ()              509                     508                 505
    // psi =    2.56835 (2.5659)    2.60012 (2.59945)   2.61516 (2.61438)   2.63627 (2.63566)       2.65877 (2.65833)   2.73436 (2.7339)
    // Mass0 =  499.183 +- 0.017    497.982 +- 0.008    497.726 +- 0.006    497.728 +- 0.007        497.744 +- 0.009    497.762 +- 0.007
    // dMass =                      -0.434205           -0.0908659          -0.0995038              -0.187103           -0.198127
 

    auto revMassFunc = new TF1("revMassFunc", "2*TMath::ACos(TMath::Sqrt((1 - x * x / [0] / [0])/(1-4*139.57 * 139.57 / [0] / [0])))");
    revMassFunc->SetParameter(0, 505.);
    std::cout << "!!!!!!!!!!! " << revMassFunc->Eval(497.762) << std::endl;
    std::cout << rc->RadCorEval() << std::endl;
    return 0;
}