#include "TF1.h"
#include "TF2.h"
#include "TTree.h"
#include <vector>
#include <complex.h>
#include <algorithm>
#include <chrono>
#include <ctime> 

class RadCor
{
private:
    TF1 *massFunc;
    TF1 *Li2;
    TF1 *D;

    Float_t E;
    Float_t pRatio;
    Float_t psi;
    Float_t pK0;
    // Maybe: 1 / alpha = 137.071999.
    static const double alpha = 1 / 137.035999;  
    static const double kaonMass = 497.611;   
    static const double electronMass = 0.511;
 

    double kaonMass; 
    double electronMass;
    double s;

    // Born cross-section.
    double Sigma0(const double &s);
    // Return FormFactor of e+e- -> KsKl.
    std::complex<double> FormFactor(const double &s);
public:
    RadCor(Float_t K0momentum, Float_t energy, Float_t ksdpsi, Float_t momentumRatio);
    double RadCorEval();
    ~RadCor();
};

RadCor::RadCor(Float_t K0momentum, Float_t energy, Float_t ksdpsi, Float_t momentumRatio)
{
    E = energy; pRatio = momentumRatio;
    psi = ksdpsi; pK0 = K0momentum;

    s = 4 * E * E;
    // Auxilary function 
    Li2 = new TF1("Li2", "-TMath::Log(1-t)/t");
    /* 
    D(x,s) = D_gamma + D_e+e-. 
    [0] - b=2*alpha/pi * (L-1); [1] - L; [2] - 2*m_e/E.
    For definition go to https://arxiv.org/abs/hep-ph/9703456v1.
    */
    D = new TF1("D func", "pow([0]/2*(1-x), [0]/2-1) * ( 1+3*[0]/8+[0]*[0]/16*(9/8-TMath::Pi()*TMath::Pi()/3)) - [0]/4*(1+x)+ \
        [0]*[0]/32*(4*(1+x)*TMath::Log(1/(1-x)) + (1+3*x*x)/(1-x)*TMath::Log(1/x)-5-x) + \
        pow([0]/2*(1-x), [0]/2-1)*(-[0]*[0]/288*(2*[1] - 15)) + pow((1/137/TMath::Pi()), 2) * ( 1/12/(1-x) * \
        pow((1-x-[2]), [0]/2)*(pow(pow(TMath::Log((1-x), 2) /[2]/[2])-5/3), 2) *(1+pow(x, 2) +[0]/6*(pow(TMath::Log((1-x), 2) /[2]/[2])-5/3)) +\
        [1]*[1]/4*(2/3*(1-x^3)/x+1/2*(1-x) + (1+x)*TMath::Log(x)) ) * ((1-x-[2] > 0) ? 1 : 0)");

    /*
    Cross-section of e+e- -> KsKl.
    For definition go to https://arxiv.org/abs/1604.02981v3 or https://doi.org/10.1142/S0217751X92001423.
    Look here: https://cmd.inp.nsk.su/websvn/filedetails.php?repname=Cmd3Sim&path=%2Ftrunk%2Fgenerator%2Fradcor%2Fsrc%2FTKnFormFactor.C
    */
    Sigma0();
}

RadCor::~RadCor()
{
    delete massFunc;
    delete Li2;
    delete D;
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

double RadCor::Sigma0(const double &s)
{
    auto ff = FormFactor(s);
    return alpha * alpha * (1 - 4 * kaonMass * kaonMass / s) * TMath::Pi() / s * fabs(ff * std::conj(ff));
}

double RadCor::RadCorEval()
{
    double beta_ = sqrt(1 - kaonMass * kaonMass / E / E);
    double L = 2 * TMath::Log(2*E / electronMass);
    double a = TMath::Pi()*TMath::Pi() / 6 - 1/4;
    double b = -1 + (1/beta_-1)*L/2 + 1/beta_*TMath::Log((1+beta_)/2) + 
                (1/beta_ + beta_)/2 *(-Li2.Integral(0,(beta_-1)/(1+beta_)) + Li2.Integral(0, (1-beta_)/(1+beta_)) - TMath::Pi()*TMath::Pi()/12 + 
                L*TMath::Log((1+beta_)/2)-2*L*TMath::Log(beta_) + 1.5*pow(TMath::Log((1+beta_)/2), 2) -0.5*pow(TMath::Log(beta_), 2) - 
                3*TMath::Log(beta_)*TMath::Log((1+beta_)/2) + L + 2*TMath::Log((1+beta_)/2) );


    /* 
    Krc(s, x1, x2);
    p[0] -  auxilary parameter: p[0] == 0 stands for normalization factor, 
                                p[0] == 1 stands for convolution of mass func and Krc.  
    For definition go to https://arxiv.org/abs/hep-ph/9703456v1.
    */
    TF2 Krc("K RadCor", [&](double* x, double*p) { 
        massFunc->SetParameters(E * (1-x[0]) * (1-x[1]), pRatio);
        double bForDFormula = 2/137/TMath::Pi()*(L - 1);
        D.SetParameters(bForDFormula, L, 2 * electronMass / E);
        return (fabs(p[0] - 1) < 0.1 ? massFunc->Eval(psi) : 1) * (1 + 2/137/TMath::Pi() * (1+a+b)) * D.Eval(x[0]) * D.Eval(x[1]) * Sigma0(s*(1-x[0]) * (1-x[1])); 
        }, 0, 1, 0, 1, 1);
    Krc.SetParameter(0, 0);

    Double_t N = Krc.Integral(0., 1., 0., 1.);
    
    Krc.SetParameter(0, 1);
    return Krc.Integral(0., 1., 0., 1.) / N;
}























double EnergyHandler::RadCor(TF1 *massFunc, Float_t E, Float_t pRatio, Float_t psi)
{
    double kaonMass = 497.611; double electronMass = 0.511;
    double s = 4 * E * E;
    double beta_ = sqrt(1 - kaonMass * kaonMass / E / E);
    TF1 Li2("Li2", "-TMath::Log(1-t)/t");
    double L = 2 * TMath::Log(2*E / electronMass);
    double a = TMath::Pi()*TMath::Pi() / 6 - 1/4;
    double b = -1 + (1/beta_-1)*L/2 + 1/beta_*TMath::Log((1+beta_)/2) + 
                (1/beta_ + beta_)/2 *(-Li2.Integral(0,(beta_-1)/(1+beta_)) + Li2.Integral(0, (1-beta_)/(1+beta_)) - TMath::Pi()*TMath::Pi()/12 + 
                L*TMath::Log((1+beta_)/2)-2*L*TMath::Log(beta_) + 1.5*pow(TMath::Log((1+beta_)/2), 2) -0.5*pow(TMath::Log(beta_), 2) - 
                3*TMath::Log(beta_)*TMath::Log((1+beta_)/2) + L + 2*TMath::Log((1+beta_)/2) );

    /* 
    D(x,s) = D_gamma + D_e+e-. 
    [0] - b=2*alpha/pi * (L-1); [1] - L; [2] - 2*m_e/E.
    For definition go to https://arxiv.org/abs/hep-ph/9703456v1.
    */
    TF1 D("D func", "pow([0]/2*(1-x), [0]/2-1) * ( 1+3*[0]/8+[0]*[0]/16*(9/8-TMath::Pi()*TMath::Pi()/3)) - [0]/4*(1+x)+ \
        [0]*[0]/32*(4*(1+x)*TMath::Log(1/(1-x)) + (1+3*x*x)/(1-x)*TMath::Log(1/x)-5-x) + \
        pow([0]/2*(1-x), [0]/2-1)*(-[0]*[0]/288*(2*[1] - 15)) + pow((1/137/TMath::Pi()), 2) * ( 1/12/(1-x) * \
        pow((1-x-[2]), [0]/2)*(pow(pow(TMath::Log((1-x), 2) /[2]/[2])-5/3), 2) *(1+x^2+[0]/6*(pow(TMath::Log((1-x), 2) /[2]/[2])-5/3)) +\
        [1]*[1]/4*(2/3*(1-x^3)/x+1/2*(1-x) + (1+x)*TMath::Log(x)) ) * ((1-x-[2] > 0) ? 1 : 0) ");

    /*
    Cross-section of e+e- -> KsKl.
    For definition go to https://arxiv.org/abs/1604.02981v3 or https://doi.org/10.1142/S0217751X92001423.
    TO-DO!!!
    */
    TF1 sigma0("sigma0", "1");

    // Krc(s, x1, x2); 
    // For definition go to https://arxiv.org/abs/hep-ph/9703456v1.
    TF2 Krc("K RadCor", [&](double* x, double*p) { 
        massFunc->SetParameters(E * (1-x[0]) * (1-x[1]), pRatio);
        double bForDFormula = 2/137/TMath::Pi()*(L - 1);
        D.SetParameters(bForDFormula, L, 2 * electronMass / E);
        //sigma0.SetParameters(); // TO-DO!!!
        return (fabs(p[0] - 1) < 0.1 ? massFunc->Eval(psi) : 1) * (1 + 2/137/TMath::Pi() * (1+a+b)) * D.Eval(x[0]) * D.Eval(x[1]) * sigma0(s*(1-x[0]) * (1-x[1])); 
        }, 0, 1, 0, 1, 1);
    Krc.SetParameter(0, 0);
    Double_t N = Krc.Integral(0., 1., 0., 1.);
    
    Krc.SetParameter(0, 1);
    return Krc.Integral(0., 1., 0., 1.) / N;
}