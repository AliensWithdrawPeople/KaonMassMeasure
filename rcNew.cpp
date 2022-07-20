#include "TF1.h"
#include <complex>
#include <cmath>

class RadCorNew
{
private:
    constexpr static double alpha = 1 / 137.035999;  
    constexpr static double kaonMass = 497.614;   
    constexpr static double electronMass = 0.511;
    constexpr static double pi = 3.14159265358979323846;
    constexpr static double a = pi * pi / 6 - 1 / 4;

    // s = 4 * e_beam^2
    double s;
    // Critical angle
    double psi;
    // L = ln(s / m_e^2)
    double L;
    // beta = 2*alpha/pi * (L-1)
    double beta;
    // Brackets before x^(beta - 1) in F.
    double par1;

    // Return FormFactor of e+e- -> KsKl.
    std::complex<double> FormFactor(const double &s);
    double SigmaBorn(const double &s);
    // Fadin-Kuraev RadCor. Page 12.
    double F(const double &x, const double &s);
    // Fadin-Kuraev RadCor. Formula 27.
    double SigmaCorrected(const double &s);
public:
    RadCorNew(double energy, double criticalAngle);
    ~RadCorNew();
    double GetSigmaBorn(const double &s)
    { return SigmaBorn(s); }
    double GetMassCorrected();
};

RadCorNew::RadCorNew(double energy, double criticalAngle)
{
    s = 4 * energy * energy;
    psi = criticalAngle;
    L = std::log(s / electronMass / electronMass);
    beta = 2 * alpha / pi *(L - 1);
    par1 = 1 + alpha / pi * (pi * pi / 3 - 1./2.) + 3./4. * beta - 1./24. * beta * beta *(L / 3. + 2 * pi * pi -37./4.); 
}

RadCorNew::~RadCorNew()
{
}

double RadCorNew::GetMassCorrected()
{
    return SigmaCorrected(s);
}

std::complex<double> RadCorNew::FormFactor(const double &s)
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

double RadCorNew::SigmaBorn(const double &s)
{
    auto ff = FormFactor(s);
    return 209.011 *alpha * alpha * std::pow(1 - 4 * kaonMass * kaonMass / s, 3./2.) * 2 / 3 * TMath::Pi() / s * std::abs(ff * std::conj(ff));
}

double RadCorNew::F(const double &x, const double &s)
{
    double f = beta * std::pow(x, beta - 1) * par1 - beta * (1 - x /2.) +
    beta * beta / 8. * (4 * (2 - x) * std::log(1 / x) + 1 / x * (1 + 3 * (1-x)*(1-x)) * std::log(1/(1-x)) - 6 + x);

    double E = std::sqrt(s) / 2;
    // In case of real e+e- pairs are not banned.
    double epemPart = 0;
    if(x - 2 * electronMass / E > 0)
    { 
        epemPart =  alpha * alpha / pi / pi * 
        (1 / 6. / x * std::pow(x - 2 * electronMass / E, beta) * std::pow(std::log(s * x * x / electronMass / electronMass) -5./3., 2) * 
        (2 - 2*x + x*x + beta/3.*(std::log(s * x * x / electronMass /electronMass) -5./3)) + 
        L*L/2. * (2./3.*(1-std::pow(1-x, 3)) / (1-x) - (2-x)*std::log(1/(1-x)) +x/2.) );
    }
    return f + epemPart;
}

double RadCorNew::SigmaCorrected(const double &s)
{
    // Photon energy upper bound 
    double eps = (std::sqrt(s) - 2 * kaonMass) / kaonMass;
    auto massFunc = new TF1("MassLnY", "sqrt([0] * [0] * (1 - (1 + sqrt(1 - [1] *[1]) * cos(x))*(1 - sqrt(1 - [1] * [1] * (1 - 4 * 139.57 * 139.57 / [0] / [0])))/ [1] / [1] ))");
    massFunc->SetParameter(1, (1 - 0.9999*0.9999) / (1 + 0.9999*0.9999));
    auto sigma = new TF1("sigma with radcor", 
    [&](double*x, double *p)
    { 
        // Not sure about sqrt(s * x[0]). It must be convolution of massFunc and cross section. 
        massFunc->SetParameter(0, sqrt(s * x[0]));
        return SigmaBorn(s * (1 - x[0])) * F(x[0], s); //fabs(p[0] - 1) < 0.1 ? massFunc->Eval(psi) : 1) * 
    }, 0, eps, 1);
    sigma->SetParameter(0, 1);
    double tmp = sigma->Integral(0, eps);
    sigma->SetParameter(0, 0);
    return sigma->Integral(0., eps);
}

void rcNew()
{
    auto rc = new RadCorNew(505, 2.61468);
    //std::cout << "Sigma = " << rc->GetMassCorrected() << std::endl;
    double energy = 505;
    double s = 4 * energy * energy;
    std::cout << "SigmaBorn = " << rc->GetSigmaBorn(s) << std::endl;
}