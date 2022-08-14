#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"

int auxFunc()
{
    double pRatio = 0.99999;
    double psi = 2.56547;

    double energy = 514;
    double mass =  499.069;
    double massK = 497.611;


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
    std::vector<Float_t> vM = {497.704, 497.687, 497.669, 497.735,  497.927, 499.078};
    std::vector<Float_t> vMerr = {6.81036e-03, 7.11921e-03, 7.26057e-03, 7.75657e-03, 9.52740e-03, 2.02138e-02};

    std::vector<Float_t> vE = {505, 508, 509, 510, 511, 514};
    std::vector<Float_t> zeroes(vM.size(), 0.0);
    TGraphErrors gr(vM.size(), vE.data(), vM.data(), zeroes.data(), vMerr.data());

    // MCGPJ with radcor and resolution correction
    // std::vector<Float_t> vMRad = { 497.613, 497.604, 497.610, 497.615, 497.608, 497.609};
    // std::vector<Float_t> vMerrRad = {3.32994e-03, 4.84249e-03, 4.76159e-03, 1.01080e-02, 5.09326e-03, 5.65773e-03};
    std::vector<Float_t> vMRad = { 497.596, 497.594, 497.595, 497.601, 497.633, 497.631};
    std::vector<Float_t> vMerrRad = {6.81036e-03, 7.11921e-03, 7.26057e-03, 7.75657e-03, 9.52740e-03, 2.02138e-02};
    TGraphErrors grRCNC(vMRad.size(), vE.data(), vMRad.data(), zeroes.data(), vMerrRad.data());
    grRCNC.SetMarkerColor(kGreen);


    // 2body Gen without resolution correction
    std::vector<Float_t> vM2b = {497.585, 497.594, 497.599, 497.594, 497.615, 497.603};
    std::vector<Float_t> vMerr2b = {5.01664e-03, 6.86874e-03, 7.28854e-03, 1.02356e-02, 7.93818e-03, 6.48984e-03};
    TGraphErrors gr2b(vM2b.size(), vE.data(), vM2b.data(), zeroes.data(), vMerr2b.data());
    gr2b.SetMarkerColor(kGreen);

    // 2body Gen with resolution correction
    std::vector<Float_t> vM2bCorr = {497.604, 497.615, 497.609, 497.603, 497.622,  497.613};
    std::vector<Float_t> vMerr2bCorr = {4.44086e-03, 6.99183e-03, 6.28059e-03, 1.02157e-02, 7.56118e-03, 6.01228e-03};
    TGraphErrors gr2bCorr(vM2bCorr.size(), vE.data(), vM2bCorr.data(), zeroes.data(), vMerr2bCorr.data());
    gr2bCorr.SetMarkerColor(kRed);

    auto massKline = new TLine(505, massK, 514, massK);
    massKline->SetLineColor(kBlue);
    massKline->SetLineWidth(2);
    
    // gr.DrawClone("AP");
    // gr2bCorr.DrawClone("Same P");
    // gr2b.DrawClone("Same AP");
    // gr2bCorr.DrawClone("Same P");
    grRCNC.DrawClone("same AP");
    massKline->DrawClone("Same");
    
    

    std::vector<Float_t> vSigma = {1.37127e-02, 1.44228e-02, 1.45462e-02, 1.46740e-02, 1.72791e-02, 1.84825e-02,  2.11118e-02, 2.26650e-02};
    std::vector<Float_t> vSigmaErr = {2.41164e-04, 2.59689e-04, 2.51916e-04, 3.13132e-04, 3.32211e-04, 4.86097e-04, 4.64792e-04, 6.67834e-04};
    std::vector<Float_t> vWidth = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4};
    std::vector<Float_t> zeroes2(vSigma.size(), 0.0);
    TGraphErrors gr2(vSigma.size(), vWidth.data(), vSigma.data(), zeroes2.data(), vSigmaErr.data());

    // gr2.DrawClone("AP");

    return 0;
}