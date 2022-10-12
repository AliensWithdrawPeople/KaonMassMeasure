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


    //std::vector<Float_t> vM = {497.645, 497.64, 497.648, 497.633, 497.676, 497.742};
    //std::vector<Float_t> vMerr = {0.0049, 0.0053, 0.0055, 0.0047, 0.0072, 0.0148};

    //std::vector<Float_t> vM = {497.628, 497.639, 497.637, 497.637, 497.6684, 497.726};
    //std::vector<Float_t> vMerr = {0.004, 0.005,  0.005, 0.005, 0.007, 0.015};

    // MCGPJ without radcor
    std::vector<Float_t> vM = {497.701, 497.680, 497.660, 497.701,  497.937, 499.109};
    std::vector<Float_t> vMerr = {7e-03, 10e-03, 8e-03, 8e-03, 1.02194e-02, 1.8e-02};

    std::vector<Float_t> vE = {505, 508, 509, 510, 511, 514};
    std::vector<Float_t> zeroes(vM.size(), 0.0);
    TGraphErrors gr(vM.size(), vE.data(), vM.data(), zeroes.data(), vMerr.data());

    // MCGPJ with radcor and resolution correction
    // kaon energy = sum of pions energies
    // std::vector<Float_t> vMRad = { 497.613, 497.604, 497.610, 497.615, 497.608, 497.609};
    // std::vector<Float_t> vMerrRad = {3.32994e-03, 4.84249e-03, 4.76159e-03, 1.01080e-02, 5.09326e-03, 5.65773e-03};
    
    // Fadin-Kuraev radcor
    std::vector<Float_t> vMRad = { 497.596, 497.601, 497.595, 497.603, 497.633, 497.631};
    std::vector<Float_t> vMerrRad = {6.81036e-03, 7.11921e-03, 7.26057e-03, 8.3e-03, 9.52740e-03, 2.02138e-02};
    TGraphErrors grRCNC(vMRad.size(), vE.data(), vMRad.data(), zeroes.data(), vMerrRad.data());
    grRCNC.SetMarkerColor(kGreen);

    // Avg energy radcor

    // Old Kl cut: (dPhi < 0.5 || dPhi > 2 * TMath::Pi() - 0.5) && phen0[j] > 40
    // std::vector<Float_t> vMRad2 = { 497.601, 497.604, 497.599, 497.611, 497.616, 497.650};
    // std::vector<Float_t> vMerrRad2 = {0.008, 0.005, 0.011, 0.011, 0.017, 0.028};

    // New Kl cut: (dPhi < 1 || dPhi > 2 * TMath::Pi() - 1) && fabs(dTheta) < 1 && phen0[j] > 40
    // std::vector<Float_t> vMRad2 = { 497.603, 497.613, 497.606, 497.616, 497.626, 497.616};
    // std::vector<Float_t> vMerrRad2 = {0.009, 0.010, 0.011, 0.011, 0.012, 0.023};

    std::vector<Float_t> vMRad2 = { 497.597, 497.593, 497.597, 497.583, 497.617, 497.642};
    std::vector<Float_t> vMRad2NC = { 497.609, 497.610, 497.608, 497.607, 497.626, 497.673};
    // RMS sigma { 497.619, 497.619, 497.614, 497.616, 497.638, 497.759}
    // Fit sigma { 497.609, 497.610, 497.608, 497.607, 497.26, 497.673}
    std::vector<Float_t> vMerrRad2 = {0.007, 0.008, 0.008, 0.008, 0.010, 0.024};
    // std::vector<Float_t> vMerrRad2 = {0.009, 0.010, 0.011, 0.011, 0.012, 0.023};
    // std::vector<Float_t> vMerrRad2 = {0, 0, 0, 0, 0, 0};
    TGraphErrors grRC2(vMRad2.size(), vE.data(), vMRad2.data(), zeroes.data(), vMerrRad2.data());
    TGraphErrors grRCNC2(vMRad2NC.size(), vE.data(), vMRad2NC.data(), zeroes.data(), vMerrRad2.data());
    // grRC2.SetMarkerColor(kRed);
                                
    std::vector<Float_t> vM2 = { 497.602, 497.613, 497.606, 497.616, 497.626, 497.616};
    std::vector<Float_t> vMerr2 = {0, 0, 0, 0, 0, 0};
    TGraphErrors grRCNC2new(vM2.size(), vE.data(), vM2.data(), zeroes.data(), vMerr2.data());
    grRCNC2new.SetMarkerColor(kBlack);

    std::vector<Float_t> vDeltaMRad = { 0.104, 0.087, 0.063, 0.118, 0.320, 1.467};
    std::vector<Float_t> vDeltaMRadErr = { 0.008, 0.007, 0.004, 0.010, 0.013, 0.017};
    TGraphErrors grRCdeltaM(vDeltaMRad.size(), vE.data(), vDeltaMRad.data(), zeroes.data(), vDeltaMRadErr.data());
    // grRCdeltaM.SetMarkerColor(kGreen);

    // deltaM_NC
    std::vector<Float_t> vDeltaM_NC_RMS = { 0.022, 0.022, 0.017, 0.033, 0.021, 0.117}; //RMS sigma
    std::vector<Float_t> vDeltaM_NC_Fit = { 0.012, 0.017, 0.011, 0.024, 0.011, 0.031}; //Fit sigma
    std::vector<Float_t> vDeltaM_NC_FitErr = {0.0006, 0.0004, 0.0006, 0.0010, 0.0007, 0.0009};
    std::vector<Float_t> vDeltaM_NC_RMSerr = {0.0006, 0.0004, 0.0006, 0.0010, 0.0021, 0.0050};
    TGraphErrors grDeltaM_NC_RMS(vDeltaM_NC_RMS.size(), vE.data(), vDeltaM_NC_RMS.data(), zeroes.data(), vDeltaM_NC_RMSerr.data());
    TGraphErrors grDeltaM_NC_Fit(vDeltaM_NC_Fit.size(), vE.data(), vDeltaM_NC_Fit.data(), zeroes.data(), vDeltaM_NC_FitErr.data());

    // Avg energy radcor without Kl cut
    std::vector<Float_t> vMRad3 = { 497.609, 497.609, 497.610, 497.600, 497.633, 497.678};
    std::vector<Float_t> vMerrRad3 = {0, 0, 0, 0, 0, 0};
    TGraphErrors grRCNC3(vMRad3.size(), vE.data(), vMRad3.data(), zeroes.data(), vMerrRad3.data());
    grRCNC3.SetMarkerColor(kRed);

    {
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
    }

    auto massKline = new TLine(505, massK, 514, massK);
    massKline->SetLineColor(kBlue);
    massKline->SetLineWidth(2);
    {
    // Without cuts
    std::vector<Float_t> vERadMean = { 0.104652, 0.0806505, 0.069, 0.0754306, 0.110864, 0.30721, 1.5044};
    std::vector<Float_t> vERadMeanErr = {0.00126345, 0.00108585, 0.000956731, 0.000943914, 0.00117893, 0.00218002, 0.00639725};
    
    // With cuts
    std::vector<Float_t> vERadMean1 = { 0.088976, 0.0805728, 0.0665978, 0.070706, 0.106681, 0.313561, 1.4743};
    std::vector<Float_t> vERadMeanErr1 = {0.0025272, 0.00253048, 0.00200995, 0.00195715, 0.00254587, 0.0048777, 0.0141468};

    std::vector<Float_t> vE2 = {505, 508, 509, 509.527, 510, 511, 514};
    std::vector<Float_t> zeroes3(vERadMean.size(), 0.0);
    TGraphErrors grDeltaERad(vERadMean.size() - 1, vE2.data(), vERadMean.data(), zeroes3.data(), vERadMeanErr.data());
    TGraphErrors grDeltaERad1(vERadMean1.size() - 1, vE2.data(), vERadMean1.data(), zeroes3.data(), vERadMeanErr1.data());
    // grDeltaERad1.SetMarkerColor(kRed);
    // grDeltaERad.DrawClone("AP");
    // grDeltaERad1.DrawClone("same P");
    }
    
    // gr.DrawClone("AP");
    // gr2b.DrawClone("Same AP");
    // gr2bCorr.DrawClone("Same P");
    // grRCNC3.DrawClone("same AP");
    // grRCNC2new.DrawClone("AP");
    // grRC2.DrawClone("AP");
    grRCNC2.DrawClone("AP");
    massKline->DrawClone("Same");

    // grDeltaM_NC_RMS.SetMarkerColor(kRed);
    // grDeltaM_NC_RMS.DrawClone("AP");
    // grDeltaM_NC_Fit.DrawClone("P same");

    // grRCdeltaM.DrawClone("AP");
    
    

    std::vector<Float_t> vSigma = {1.37127e-02, 1.44228e-02, 1.45462e-02, 1.46740e-02, 1.72791e-02, 1.84825e-02,  2.11118e-02, 2.26650e-02};
    std::vector<Float_t> vSigmaErr = {2.41164e-04, 2.59689e-04, 2.51916e-04, 3.13132e-04, 3.32211e-04, 4.86097e-04, 4.64792e-04, 6.67834e-04};
    std::vector<Float_t> vWidth = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4};
    std::vector<Float_t> zeroes2(vSigma.size(), 0.0);
    TGraphErrors gr2(vSigma.size(), vWidth.data(), vSigma.data(), zeroes2.data(), vSigmaErr.data());

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