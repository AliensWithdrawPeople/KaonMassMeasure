#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"

int auxFunc()
{
    double pRatio = 0.99999;
    double psi = 2.56547;

    double energy = 514;
    double mass =  499.118;
    double massK = 497.614;


    auto massFunc = new TF1("MassLnY", "sqrt([0] * [0] * (1 - (1 + sqrt(1 - [1] *[1]) * cos(x))*(1 - sqrt(1 - [1] * [1] * (1 - 4 * 139.57 * 139.57 / [0] / [0])))/ [1] / [1] ))");
    massFunc->SetParameter(1, (1 - pRatio*pRatio) / (1 + pRatio*pRatio));

    //auto tmpFunc = new TF1("MassLnY", 
    //"sqrt(x * x * (1 - (1 + sqrt(1 - [1] *[1]) * cos([0]))*(1 - sqrt(1 - [1] * [1] * (1 - 4 * 139.57 * 139.57 / x / x)))/ [0] / [0] ))", 480, 520);
    auto tmpFunc = new TF1("MassLnY", "sqrt(x * x * (1 - (1 + sqrt(1 - [1] *[1]) * cos([0]))*(1 - sqrt(1 - [1] * [1] * (1 - 4 * 139.57 * 139.57 / x / x)))/ [1] / [1] ))");
    tmpFunc->SetParameters(psi, (1 - pRatio*pRatio) / (1 + pRatio*pRatio));
    massFunc->SetParameter(0, energy);
    double Emin = tmpFunc->GetX(490, 450, 550, 1e-5, 10000);
    double Emax = tmpFunc->GetX(505, 450, 550, 1e-5, 10000);
    std::cout<< Emin << " : " << Emax << std::endl;


    auto revMassFunc = new TF1("revMassFunc", "2*TMath::ACos(TMath::Sqrt((1 - x * x / [0] / [0])/(1-4*139.57 * 139.57 / [0] / [0])))");
    revMassFunc->SetParameter(0, energy);
    std::cout << "At E = " << energy << " and M = " << mass << " Psi = " << revMassFunc->Eval(mass) << std::endl;


    //std::vector<Float_t> vM = {497.645, 497.64, 497.648, 497.633, 497.676, 497.742};
    //std::vector<Float_t> vMerr = {0.0049, 0.0053, 0.0055, 0.0047, 0.0072, 0.0148};

    //std::vector<Float_t> vM = {497.628, 497.639, 497.637, 497.637, 497.6684, 497.726};
    //std::vector<Float_t> vMerr = {0.004, 0.005,  0.005, 0.005, 0.007, 0.015};

    // MCGPJ without radcor
    std::vector<Float_t> vM = {497.771, 497.728, 497.714, 497.735,  497.953, 499.159};
    std::vector<Float_t> vMerr = {5.16393e-03, 5.26006e-03, 5.44217e-03, 4.68768e-03, 7.08306e-03, 1.47878e-02};

    std::vector<Float_t> vE = {505, 508, 509, 510, 511, 514};
    std::vector<Float_t> zeroes(vM.size(), 0.0);
    TGraphErrors gr(vM.size(), vE.data(), vM.data(), zeroes.data(), vMerr.data());

    // 2body Gen
    std::vector<Float_t> vM1 = {497.604, 497.615, 497.608, 0, 497.617,  497.612};
    std::vector<Float_t> vMerr1 = {4.44086e-03, 5.90347e-03, 6.28059e-03, 0, 5.80474e-03, 6.01228e-03};
    TGraphErrors gr1(vM1.size(), vE.data(), vM1.data(), zeroes.data(), vMerr1.data());
    gr1.SetMarkerColor(kRed);

    auto massKline = new TLine(505, massK, 514, massK);
    massKline->SetLineColor(kBlue);
    massKline->SetLineWidth(2);
    
    gr.DrawClone("AP");
    gr1.DrawClone("Same P");
    massKline->DrawClone("Same");
    
    

    std::vector<Float_t> vSigma = {1.37127e-02, 1.44228e-02, 1.45462e-02, 1.46740e-02, 1.72791e-02, 1.84825e-02,  2.11118e-02, 2.26650e-02};
    std::vector<Float_t> vSigmaErr = {2.41164e-04, 2.59689e-04, 2.51916e-04, 3.13132e-04, 3.32211e-04, 4.86097e-04, 4.64792e-04, 6.67834e-04};
    std::vector<Float_t> vWidth = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4};
    std::vector<Float_t> zeroes2(vSigma.size(), 0.0);
    TGraphErrors gr2(vSigma.size(), vWidth.data(), vSigma.data(), zeroes2.data(), vSigmaErr.data());

    // gr2.DrawClone("AP");

    return 0;
}