#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"

int auxFunc()
{
    double pRatio = 0.99999;
    double psi = 2.56547;

    double energy = 509;
    double mass =  497.729;
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

    std::vector<Float_t> vM = {497.652, 497.626, 497.640, 497.628, 497.668, 498.036};
    std::vector<Float_t> vMerr = {0.004, 0.005, 0.005, 0.004, 0.006, 0.013};

    std::vector<Float_t> vE = {505, 508, 509, 510, 511, 514};
    std::vector<Float_t> zeroes(vM.size(), 0.0);
    TGraphErrors gr(vM.size(), vE.data(), vM.data(), zeroes.data(), vMerr.data());
    auto gr1 = new TGraphErrors();
    gr1->AddPoint(510, 497.619);
    gr1->SetPointError(0, 0, 0.005);
    gr1->SetMarkerColor(kGreen);

    auto massKline = new TLine(505, massK, 514, massK);
    massKline->SetLineColor(kBlue);
    massKline->SetLineWidth(2);

    auto massKline2 = new TLine(505, 497.619, 514, 497.619);
    massKline2->SetLineColor(kGreen);
    massKline2->SetLineWidth(2);
    gr.DrawClone();
    massKline->DrawClone("Same");
    massKline2->DrawClone("Same");
    gr1->DrawClone("Same");
    
    

    std::vector<Float_t> vPsi = {2.73532, 2.73491, 2.73476, 2.73462, 2.73463};
    std::vector<Float_t> vPsiErr = {4.07562e-04, 3.35361e-04, 2.88769e-04, 2.57586e-04, 2.35667e-04};
    std::vector<Float_t> vWidth = {0.1, 0.15, 0.2, 0.25, 0.3};
    std::vector<Float_t> zeroes2(vPsi.size(), 0.0);
    TGraphErrors gr2(vPsi.size(), vWidth.data(), vPsi.data(), zeroes2.data(), vPsiErr.data());

    //gr2.DrawClone();

    return 0;
}