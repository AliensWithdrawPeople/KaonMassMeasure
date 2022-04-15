#include "TF1.h"
#include "TF2.h"
#include "TTree.h"
#include <vector>
#include <algorithm>
#include <chrono>
#include <ctime> 

class RadCor
{
private:
    TF1 *massFunc;
    TF1 *Li2;
    TF1 *D;
    TF1 *sigma0;

    Float_t E;
    Float_t pRatio;
    Float_t psi;
    Float_t pK0;

    double kaonMass; 
    double electronMass;
    double s;

    void sigma0Def();
public:
    RadCor(Float_t K0momentum, Float_t energy, Float_t ksdpsi, Float_t momentumRatio);
    double RadCorEval();
    ~RadCor();
};

RadCor::RadCor(Float_t K0momentum, Float_t energy, Float_t ksdpsi, Float_t momentumRatio)
{
    E = energy; pRatio = momentumRatio;
    psi = ksdpsi; pK0 = K0momentum;

    kaonMass = 497.611; 
    electronMass = 0.511;
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
        pow((1-x-[2]), [0]/2)*(pow(pow(TMath::Log((1-x), 2) /[2]/[2])-5/3), 2) *(1+x^2+[0]/6*(pow(TMath::Log((1-x), 2) /[2]/[2])-5/3)) +\
        [1]*[1]/4*(2/3*(1-x^3)/x+1/2*(1-x) + (1+x)*TMath::Log(x)) ) * ((1-x-[2] > 0) ? 1 : 0) ");

    /*
    Cross-section of e+e- -> KsKl.
    For definition go to https://arxiv.org/abs/1604.02981v3 or https://doi.org/10.1142/S0217751X92001423.
    TO-DO!!!
    */
    sigma0Def();
}

RadCor::~RadCor()
{
    delete massFunc;
    delete Li2;
    delete D;
    delete sigma0;
}

void RadCor::sigma0Def()
{
    sigma0 = new TF1("sigma0", "1");

    // {rho, omega, phi}
    std::vector<double> mass = {775.26, 782.63, 1019.461};
    // Ã, Mev. {{rho->ee, omega->ee, phi->ee}, {rho, omega, phi}}
    std::vector<std::vector<double>> width = {{7.04e-3, 0.60e-3, 1.4447}, {147.4, 8.68, 4.249}};
    std::vector<double> branch;
    std::vector<double> Dv = {};
    std::vector<double> coupConstGamma;
    std::vector<double> coupConstKK;

    double alpha = 1/137;

    for(int i = 1; i < 3; i++)
    {
        coupConstGamma.push_back(sqrt(3 * mass[i] * mass[i] * mass[i] * width[0][i] / (4 * TMath::Pi() * alpha))); 
        coupConstKK.push_back( sqrt(6*TMath::Pi() * mass[i] * mass[i] * width[1][i] * branch[i] / () ) );
    }

    
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