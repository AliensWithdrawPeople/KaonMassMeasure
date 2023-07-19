#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TCanvas.h"

int MassVsE()
{
    TCanvas c("canv", "canv");
    std::vector<Float_t> zeroes(100, 0.0);
/*
******************************
* Prev results:
******************************
*/ 
    std::vector<Float_t> vMPDG = {497.742, 497.661, 497.625, 497.634, 497.583, 497.607};
    std::vector<Float_t> vMPDGerr = {0.085, 0.033, 0.031, 0.024, 0.021, 0.017};
    std::vector<Float_t> vExpNum = {505, 506, 507, 508, 509, 510};
    TGraphErrors grMPDG(vMPDG.size(), vExpNum.data(), vMPDG.data(), zeroes.data(), vMPDGerr.data());
    // grMPDG.DrawClone("same AP");

/*
******************************
* Mass vs E_beam in Exp:
******************************
*/ 

    // std::vector<Float_t> vM_Exp = {497.545, 497.566, 497.539, 497.547, 497.544, 497.571, 497.579, 497.566, 497.6};
    // std::vector<Float_t> vM_Exp = {497.632, 497.608, 497.585, 497.596, 497.604, 497.619, 497.62, 497.597, 497.635};
    std::vector<Float_t> vM_Exp = {497.571, 497.549, 497.527, 497.539, 497.538, 497.55, 497.545, 497.517, 497.542};

    std::vector<Float_t> vM_Exp2 = {497.577, 497.568, 497.533, 497.554, 497.544, 497.57, 497.589, 497.568, 497.604};
    std::vector<Float_t> vM_Exp_vis = {497.677, 497.667, 497.622, 497.634, 497.615, 497.686, 497.78, 497.902, 499.057};

    std::vector<Float_t> vMerrExp = {0.036, 0.014, 0.011, 0.008, 0.007, 0.01, 0.013, 0.014, 0.029};

    std::vector<Float_t> vRC = {497.652, 497.647, 497.638, 497.627, 497.627, 497.679, 497.761, 497.912, 499.053};
    // std::vector<Float_t> vE = {505, 508, 508.5, 509, 509.5, 510, 510.5, 511, 514};
    std::vector<Float_t> vE = {504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 513.864};
    

    TGraphErrors grRCNC_exp(vM_Exp.size(), vE.data(), vM_Exp.data(), zeroes.data(), vMerrExp.data());
    TGraphErrors grRCNC_exp2(vM_Exp.size(), vE.data(), vM_Exp2.data(), zeroes.data(), vMerrExp.data());
    TGraphErrors grNC_exp(vM_Exp.size(), vE.data(), vM_Exp_vis.data(), zeroes.data(), vMerrExp.data());

    TGraphErrors grRC(vM_Exp.size(), vE.data(), vRC.data(), zeroes.data(), zeroes.data());

    grRCNC_exp.SetTitle("Black -- with Energy controll, Blue -- without");
    grRCNC_exp.GetXaxis()->SetTitle("E_{beam}, MeV");
    grRCNC_exp.GetYaxis()->SetTitle("M^{(FullRec)}_{RC NC RecCor}, #frac{MeV}{c^{2}}");
    grRCNC_exp.GetYaxis()->SetTitleOffset(1.2);
    grRCNC_exp.SetName("EnControl");
    grRCNC_exp.SetMarkerSize(2);
    grRCNC_exp.Fit("pol0", "ME");
    
    grRCNC_exp2.SetName("NoEnControl");
    grRCNC_exp2.SetMarkerColor(kBlue);
    grRCNC_exp2.SetLineColor(kBlue);
    grRCNC_exp2.SetMarkerStyle(22);
    grRCNC_exp2.SetMarkerSize(2);
    grRCNC_exp2.Fit("pol0", "ME+");

    grRC.SetName("RC");
    grRC.SetMarkerColor(kBlue);
    grRC.SetLineColor(kBlue);
    grRC.SetMarkerStyle(22);
    grRC.SetMarkerSize(2);

    grRCNC_exp.DrawClone("AP");

    // grNC_exp.DrawClone("AP");
    // grRC.DrawClone("P same");
    // grRCNC_exp2.DrawClone("AP");

/*
******************************
* Mass vs E_beam in MC:
******************************
*/ 
    // std::vector<Float_t> vM_MC_NoSmearing = {497.590, 497.608, 497.614, 497.615, 497.605, 497.618, 497.601, 497.635, 497.600};
    // std::vector<Float_t> vM_MC_NoSmearingErr = {0.011, 0.013, 0.010, 0.009, 0.009, 0.009, 0.010, 0.010, 0.034};

    // std::vector<Float_t> vM_MC_WithSmearing = {497.627, 497.604, 497.612, 497.603, 497.606, 497.611, 497.607, 497.613,  497.609};
    // std::vector<Float_t> vM_MC_WithSmearingErr = {0.008, 0.009 , 0.003, 0.006, 0.006, 0.003, 0.006, 0.003, 0.023};
    
    // std::vector<Float_t> vM_MC_WrongEnergy = {497.575, 497.582, 497.584, 497.586, 497.588, 497.593};
    // std::vector<Float_t> vM_MC_WrongEnergyErr = {0.002, 0.002, 0.003, 0.003, 0.003, 0.003};
    // std::vector<Float_t> vE_short = {508.5, 509, 509.5, 510, 510.5, 511};

    // TGraphErrors grMC_NoSmearing(vE.size(), vE.data(), vM_MC_NoSmearing.data(), zeroes.data(), vM_MC_NoSmearingErr.data());
    // TGraphErrors grMC_WithSmearing(vE.size(), vE.data(), vM_MC_WithSmearing.data(), zeroes.data(), vM_MC_WithSmearingErr.data());
    // TGraphErrors grMC_WrongEnergy(vE_short.size(), vE_short.data(), vM_MC_WrongEnergy.data(), zeroes.data(), vM_MC_WrongEnergyErr.data());

    // TLine massKline(504.5, 497.614, 514.5, 497.614);
    // massKline.SetLineStyle(2);
    // massKline.SetLineWidth(4);
    // massKline.SetLineColor(kRed);

    // grMC_NoSmearing.SetTitle("Black -- without Energy Smearing, Blue -- with Energy Smearing, Red line --- mass in generator");
    // grMC_NoSmearing.GetXaxis()->SetTitle("E_{beam}, MeV");
    // grMC_NoSmearing.GetYaxis()->SetTitle("M^{(FullRec)}_{RC NC}, #frac{MeV}{c^{2}}");
    // grMC_NoSmearing.GetYaxis()->SetTitleOffset(1.2);
    // grMC_NoSmearing.SetName("NoSmear");
    // grMC_NoSmearing.SetMarkerSize(2);
    // grMC_NoSmearing.Fit("pol0", "ME");
    
    // grMC_WithSmearing.SetName("Smear");
    // grMC_WithSmearing.SetMarkerColor(kBlue);
    // grMC_WithSmearing.SetLineColor(kBlue);
    // grMC_WithSmearing.SetMarkerStyle(22);
    // grMC_WithSmearing.SetMarkerSize(2);
    // grMC_WithSmearing.Fit("pol0", "ME+");

    // grMC_WrongEnergy.SetName("Smear");
    // grMC_WrongEnergy.SetMarkerColor(kBlue);
    // grMC_WrongEnergy.SetLineColor(kBlue);
    // grMC_WrongEnergy.SetMarkerStyle(21);
    // grMC_WrongEnergy.SetMarkerSize(2);
    // grMC_WrongEnergy.Fit("pol0", "ME+");

    // grMC_WrongEnergy.DrawClone("AP");


    // grMC_NoSmearing.DrawClone("AP");
    // grMC_WithSmearing.DrawClone("AP");
    // massKline.DrawClone("same");

/*
******************************
* Mass vs E_beam in EXP for Cowboys and Sailors:
******************************
*/

    // std::vector<Float_t> vM_Exp_cowboy = {497.747, 497.659, 497.615, 497.605, 497.621, 497.691, 497.765, 497.899, 499.158};
    // std::vector<Float_t> vMerr_Exp_cowboy = {0.0549, 0.0214, 0.0117, 0.0117, 0.0096, 0.0099, 0.0112, 0.0145, 0.0464};
    // TGraphErrors grNC_exp_cowboy(vM_Exp_cowboy.size(), vE.data(), vM_Exp_cowboy.data(), zeroes.data(), vMerr_Exp_cowboy.data());

    // grNC_exp_cowboy.SetTitle("Black -- Cowboy, Blue -- Sailor");
    // grNC_exp_cowboy.GetXaxis()->SetTitle("E_{beam}, MeV");
    // grNC_exp_cowboy.GetYaxis()->SetTitle("M^{(FullRec)}_{NC RecCor}, #frac{MeV}{c^{2}}");
    // grNC_exp_cowboy.GetYaxis()->SetTitleOffset(1.2);
    // grNC_exp_cowboy.SetName("Cowboy");
    // grNC_exp_cowboy.SetMarkerSize(1);

    // std::vector<Float_t> vM_Exp_sailor = {497.777, 497.644, 497.625, 497.63, 497.624, 497.682, 497.764, 497.885};
    // std::vector<Float_t> vMerr_Exp_sailor = {0.0547, 0.0211, 0.0118, 0.0118, 0.0097, 0.0098, 0.0112, 0.0144};
    // TGraphErrors grNC_exp_sailor(vM_Exp_sailor.size(), vE.data(), vM_Exp_sailor.data(), zeroes.data(), vMerr_Exp_sailor.data());

    // grNC_exp_sailor.SetTitle("Black -- Cowboy, Blue -- Sailor");
    // grNC_exp_sailor.GetXaxis()->SetTitle("E_{beam}, MeV");
    // grNC_exp_sailor.GetYaxis()->SetTitle("M^{(FullRec)}_{NC RecCor}, #frac{MeV}{c^{2}}");
    // grNC_exp_sailor.GetYaxis()->SetTitleOffset(1.2);
    // grNC_exp_sailor.SetName("Sailor");
    // grNC_exp_sailor.SetMarkerColor(kBlue);
    // grNC_exp_sailor.SetMarkerSize(1);

    // grNC_exp_cowboy.DrawClone("AP");
    // grNC_exp_sailor.DrawClone("P Same");

    
/*
******************************
* Mass_Cowboys - Mass_ Sailorsvs E_beam in EXP:
******************************
*/
    // std::vector<Float_t> vDiffM_Exp_cowboy_sailor = {-0.03, 0.015, -0.01, -0.025, -0.003, 0.009, 0.001, 0.014};
    // std::vector<Float_t> vDiffMerr_Exp_cowboy_sailor = {0.077499, 0.030053, 0.016617, 0.016617, 0.013647, 0.01393, 0.015839, 0.020436};
    // TGraphErrors grNC_Exp_cowboy_sailor_diff(vDiffM_Exp_cowboy_sailor.size(), vE.data(), vDiffM_Exp_cowboy_sailor.data(), zeroes.data(), vDiffMerr_Exp_cowboy_sailor.data());
    // grNC_Exp_cowboy_sailor_diff.SetTitle("Cowboy Sailor Difference");
    // grNC_Exp_cowboy_sailor_diff.GetXaxis()->SetTitle("E_{beam}, MeV");
    // grNC_Exp_cowboy_sailor_diff.GetYaxis()->SetTitle("#Delta M, #frac{MeV}{c^{2}}");
    // grNC_Exp_cowboy_sailor_diff.GetYaxis()->SetTitleOffset(1.2);
    // grNC_Exp_cowboy_sailor_diff.SetName("CowboySailorDiff");
    // grNC_Exp_cowboy_sailor_diff.SetMarkerSize(1);

    // grNC_Exp_cowboy_sailor_diff.DrawClone("AP");
    
    c.DrawClone();
    return 0;
}