#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"

int auxFunc()
{
    double pRatio = 0.99999;
    double psi = 2.56547;

    double energy = 514;
    double mass =  499.109;
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


    std::vector<Float_t> vE = {505, 508, 509, 510, 511, 514};
    std::vector<Float_t> zeroes(6, 0.0);

    // Full reconstruction
    // Avg energy radcor
    std::vector<Float_t> vMRad2 = { 497.597, 497.593, 497.597, 497.583, 497.617, 497.642};
    std::vector<Float_t> vMRad2NC = {497.609, 497.617, 497.608, 497.611, 497.617, 497.612};
    std::vector<Float_t> vMerrRad2 = {0.007, 0.009, 0.008, 0.010, 0.010, 0.019};
    TGraphErrors grRC(vMRad2.size(), vE.data(), vMRad2.data(), zeroes.data(), vMerrRad2.data());
    TGraphErrors grRCNC(vMRad2NC.size(), vE.data(), vMRad2NC.data(), zeroes.data(), vMerrRad2.data());
    // grRC2.SetMarkerColor(kRed);
    
    // RMS
    // std::vector<Float_t> vM2 = { 497.597, 497.600, 497.597, 497.587, 497.606, 497.581};
    std::vector<Float_t> vM2 = {497.620, 497.618, 497.615, 497.610, 497.633, 497.704};
    std::vector<Float_t> vMerr2 = {0.007, 0.009, 0.008, 0.010, 0.010, 0.019};
    TGraphErrors grRCNCnew(vM2.size(), vE.data(), vM2.data(), zeroes.data(), vMerr2.data());
    grRCNCnew.SetMarkerColor(kBlack);

    std::vector<Float_t> vDeltaMRad = {94, 61, 53, 90, 320, 1504};
    std::vector<Float_t> vDeltaMRadErr = { 0.008, 0.007, 0.004, 0.010, 0.013, 0.017};
    TGraphErrors grRCdeltaM(vDeltaMRad.size(), vE.data(), vDeltaMRad.data(), zeroes.data(), vDeltaMRadErr.data());
    // grRCdeltaM.SetMarkerColor(kGreen);

    // deltaM_NC
    std::vector<Float_t> vDeltaM_NC_RMS = {0.023, 0.018, 0.018, 0.023, 0.027, 0.123}; //RMS sigma
    std::vector<Float_t> vDeltaM_NC_Fit = { 0.012, 0.017, 0.011, 0.024, 0.011, 0.031}; //Fit sigma
    std::vector<Float_t> vDeltaM_NC_FitErr = {0.0006, 0.0004, 0.0006, 0.0010, 0.0007, 0.0009};
    std::vector<Float_t> vDeltaM_NC_RMSerr = {0.0006, 0.0004, 0.0006, 0.0010, 0.0021, 0.0050};
    TGraphErrors grDeltaM_NC_RMS(vDeltaM_NC_RMS.size(), vE.data(), vDeltaM_NC_RMS.data(), zeroes.data(), vDeltaM_NC_RMSerr.data());
    TGraphErrors grDeltaM_NC_Fit(vDeltaM_NC_Fit.size(), vE.data(), vDeltaM_NC_Fit.data(), zeroes.data(), vDeltaM_NC_FitErr.data());

    // Critical angle
    std::vector<Float_t> vMCrAngle = {497.685, 497.695, 497.67, 497.728, 497.914, 499.079};
    std::vector<Float_t> vMCrAngleRCNC_Fit = {497.596, 497.615, 497.612, 497.621, 497.598, 497.582}; //Fit sigma
    std::vector<Float_t> vMCrAngleRCNC_RMS = {497.608, 497.643, 497.608, 497.634, 497.646, 497.8}; //RMS sigma

    std::vector<Float_t> vMCrAngleErr = {0.010, 0.013, 0.013, 0.014, 0.015, 0.034};
    std::vector<Float_t> vMCrAngleRCNC_FitErr = {0.014, 0.016, 0.016, 0.018, 0.021, 0.057};
    std::vector<Float_t> vMCrAngleRCNC_RMSErr = {0.010, 0.013, 0.013, 0.013, 0.015, 0.034};

    TGraphErrors grMCrAngle(vMCrAngle.size(), vE.data(), vMCrAngle.data(), zeroes.data(), vMCrAngleErr.data());
    TGraphErrors grMCrAngleRCNC_Fit(vMCrAngleRCNC_Fit.size(), vE.data(), vMCrAngleRCNC_Fit.data(), zeroes.data(), vMCrAngleRCNC_FitErr.data());
    TGraphErrors grMCrAngleRCNC_RMS(vMCrAngleRCNC_RMS.size(), vE.data(), vMCrAngleRCNC_RMS.data(), zeroes.data(), vMCrAngleRCNC_RMSErr.data());

    auto massKline = new TLine(505, massK, 514, massK);
    massKline->SetLineColor(kBlue);
    massKline->SetLineWidth(2);

    grMCrAngle.SetMarkerColor(kRed);
    grMCrAngleRCNC_RMS.SetMarkerColor(kRed);
    grMCrAngleRCNC_Fit.SetMarkerColor(kBlack);
    grRCNCnew.SetMarkerColor(kRed);

    // grRCNC.DrawClone("AP");
    // grRCNCnew.DrawClone("AP same");
    grDeltaM_NC_Fit.DrawClone("AP");
    grDeltaM_NC_RMS.DrawClone("P same");
    // massKline->DrawClone("Same");

    std::vector<Float_t> vSigma = {1.37127e-02, 1.44228e-02, 1.45462e-02, 1.46740e-02, 1.72791e-02, 1.84825e-02,  2.11118e-02, 2.26650e-02};
    std::vector<Float_t> vSigmaErr = {2.41164e-04, 2.59689e-04, 2.51916e-04, 3.13132e-04, 3.32211e-04, 4.86097e-04, 4.64792e-04, 6.67834e-04};
    std::vector<Float_t> vWidth = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4};
    std::vector<Float_t> zeroes2(vSigma.size(), 0.0);
    TGraphErrors grSigma(vSigma.size(), vWidth.data(), vSigma.data(), zeroes2.data(), vSigmaErr.data());

    std::vector<Float_t> vAngle = { 2.73092, 2.65629, 2.63373, 2.62245, 2.61362, 2.59637, 2.55372};
    std::vector<Float_t> vEtmp = {505, 508, 509, 509.527, 510, 511, 514};
    TGraph grAngleVsE(vEtmp.size(), vEtmp.data(), vAngle.data());
    std::cout << "Psi angle vs Energy correlation factor = " << grAngleVsE.GetCorrelationFactor() << std::endl;

    std::vector<Float_t> vMPDG = {497.742, 497.661, 497.625, 497.634, 497.583, 497.607};
    std::vector<Float_t> vMPDGerr = {0.085, 0.033, 0.031, 0.024, 0.021, 0.017};
    std::vector<Float_t> vExpNum = {505, 506, 507, 508, 509, 510};
    std::vector<Float_t> zeroes8(vMPDG.size(), 0.0);
    TGraphErrors grMPDG(vMPDG.size(), vExpNum.data(), vMPDG.data(), zeroes8.data(), vMPDGerr.data());

    // grMPDG.DrawClone("same AP");
    // massKline->DrawClone("Same");

    // gr2.DrawClone("AP");

    return 0;
}