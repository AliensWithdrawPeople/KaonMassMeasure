#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLatex.h"

int MassVsE()
{
    TCanvas c("canv", "canv");
    std::vector<Float_t> zeroes(100, 0.0);
/*
******************************
* Prev results:
******************************
*/ 
    std::vector<Float_t> vMPDG = {0,    497.742, 497.661, 497.625, 497.583, 497.607, 497.580, 0};
    std::vector<Float_t> vMPDGerr = {0, 0.085,   0.033,   0.031,   0.021,   0.017,   0.034, 0};
    std::vector<Float_t> vExpNum = {500, 505, 506, 507, 508, 509, 510, 520};
    TGraphErrors grMPDG(vMPDG.size(), vExpNum.data(), vMPDG.data(), zeroes.data(), vMPDGerr.data());
    grMPDG.GetYaxis()->SetTitle("M_{K^{0}}, MeV/c^{2}");
    grMPDG.GetYaxis()->CenterTitle(true);
    grMPDG.GetXaxis()->SetTickLength(0.);
    grMPDG.GetXaxis()->SetNdivisions(0);
    TLatex *latex1 = new TLatex(grMPDG.GetX()[1], grMPDG.GetY()[1],"       CMD '85"); 
    TLatex *latex2 = new TLatex(grMPDG.GetX()[2], grMPDG.GetY()[2],"       CMD '87"); 
    TLatex *latex3 = new TLatex(grMPDG.GetX()[3], grMPDG.GetY()[3],"       NA48 '02"); 
    TLatex *latex4 = new TLatex(grMPDG.GetX()[4], grMPDG.GetY()[4],"       KLOE '07"); 
    TLatex *latex5 = new TLatex(grMPDG.GetX()[5], grMPDG.GetY()[5],"       CLEO-c data '14"); 
    TLatex *latex6 = new TLatex(grMPDG.GetX()[6], grMPDG.GetY()[6],"       This work"); 
    latex1->SetTextAlign(10);
    latex2->SetTextAlign(10);
    latex3->SetTextAlign(10);
    latex4->SetTextAlign(10);
    latex5->SetTextAlign(10);
    latex6->SetTextAlign(10);

    // latex1->SetTextAngle(90);
    // latex2->SetTextAngle(90);
    // latex3->SetTextAngle(90);
    // latex4->SetTextAngle(90);
    // latex5->SetTextAngle(90);
    // latex6->SetTextAngle(90);

    grMPDG.GetListOfFunctions()->Add(latex1); 
    grMPDG.GetListOfFunctions()->Add(latex2); 
    grMPDG.GetListOfFunctions()->Add(latex3); 
    grMPDG.GetListOfFunctions()->Add(latex4); 
    grMPDG.GetListOfFunctions()->Add(latex5); 
    grMPDG.GetListOfFunctions()->Add(latex6); 

    grMPDG.GetXaxis()->SetRangeUser(505, 510);
    grMPDG.GetYaxis()->SetRangeUser(497.54, 497.8);
    grMPDG.DrawClone("AP");

/*
******************************
* Mass vs E_beam in Exp:
******************************
*/ 
    std::vector<Float_t> vE = {504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 511.444, 513.864};

    // 0.3 > |\theta_{K_S} - \pi/2|; M_{avg} = 497.580 \pm 0.004 MeV; Events = 319198
    std::vector<Float_t> vM_Exp = {497.529, 497.591, 497.576, 497.572, 497.562, 497.586, 497.599, 497.583, 497.588, 497.603};
    // With phi width = 4.5 MeV
    // std::vector<Float_t> vM_Exp = {497.529, 497.585, 497.565, 497.563, 497.565, 497.587, 497.587, 497.592, 497.605, 497.5987};
    std::vector<Float_t> vMerrExp = {0.046, 0.018, 0.012, 0.011, 0.01, 0.01, 0.011, 0.015, 0.019, 0.053};

    // 0.4 > |\theta_{K_S} - \pi/2| > 0.3; M_{avg} = 497.564 \pm 0.006 MeV; Events = 93320
    std::vector<Float_t> vM_Exp_relaxed = {497.624, 497.544, 497.561, 497.550, 497.565, 497.563, 497.582, 497.552, 497.602, 497.611};
    std::vector<Float_t> vMerrExp_relaxed = {0.066, 0.031, 0.018, 0.017, 0.012, 0.013, 0.017, 0.025, 0.036, 0.107};
    
    // 0.3 > |\theta_{K_S} - \pi/2| > 0.2; M_{avg} = 497.571 \pm 0.006 MeV; Events = 101771
    std::vector<Float_t> vM_Exp_tightened = {497.466, 497.574, 497.553, 497.563, 497.563, 497.587, 497.593, 497.583, 497.52, 497.522};
    std::vector<Float_t> vMerrExp_tightened = {0.07, 0.03, 0.017, 0.016, 0.012, 0.012, 0.016, 0.023, 0.032, 0.096};

    TGraphErrors grRCNC_exp(vM_Exp.size(), vE.data(), vM_Exp.data(), zeroes.data(), vMerrExp.data());
    TGraphErrors grRCNC_exp_relaxed(vM_Exp_relaxed.size(), vE.data(), vM_Exp_relaxed.data(), zeroes.data(), vMerrExp_relaxed.data());
    TGraphErrors grRCNC_exp_tightened(vM_Exp_tightened.size(), vE.data(), vM_Exp_tightened.data(), zeroes.data(), vMerrExp_tightened.data());

    std::vector<Float_t> vM_Exp_NoEnergyCorr = {497.529, 497.591, 497.576, 497.572, 497.562, 497.586, 497.599, 497.583, 497.588, 497.607};
    std::vector<Float_t> vMerrExp_NoEnergyCorr = {0.046, 0.018, 0.012, 0.011, 0.01, 0.01, 0.011, 0.015, 0.019, 0.054};
    TGraphErrors grRCNC_exp_NoEnergyCorr(vM_Exp_NoEnergyCorr.size(), vE.data(), vM_Exp_NoEnergyCorr.data(), 
                                        zeroes.data(), vMerrExp_NoEnergyCorr.data());

    grRCNC_exp.SetTitle("Black -- with Energy correction, Blue -- without");
    grRCNC_exp.GetXaxis()->SetTitle("E_{beam}, MeV");
    grRCNC_exp.GetYaxis()->SetTitle("M^{(FullRec)}_{RC NC RecCor}, #frac{MeV}{c^{2}}");
    grRCNC_exp.GetYaxis()->SetTitleOffset(1.2);
    grRCNC_exp.SetName("EnControl");
    grRCNC_exp.SetMarkerSize(1.2);
    grRCNC_exp.Fit("pol0", "ME");

    grRCNC_exp_relaxed.SetName("KsTheta_relaxed");
    grRCNC_exp_relaxed.SetMarkerSize(1.2);
    grRCNC_exp_relaxed.SetMarkerStyle(22);
    grRCNC_exp_relaxed.SetMarkerColor(kBlue);

    grRCNC_exp_tightened.SetName("KsTheta_tightened");
    grRCNC_exp_tightened.SetMarkerSize(1.2);
    grRCNC_exp_tightened.SetMarkerStyle(23);
    grRCNC_exp_tightened.SetMarkerColor(kMagenta);

    grRCNC_exp_NoEnergyCorr.SetName("NoEnControl");
    grRCNC_exp_NoEnergyCorr.SetMarkerSize(1.5);
    grRCNC_exp_NoEnergyCorr.SetMarkerStyle(22);
    grRCNC_exp_NoEnergyCorr.SetMarkerColor(kBlue);
    grRCNC_exp_NoEnergyCorr.Fit("pol0", "ME");
    grRCNC_exp_NoEnergyCorr.SetLineColor(kBlue);
    grRCNC_exp_NoEnergyCorr.GetFunction("pol0")->SetLineColor(kBlue);

    grRCNC_exp.DrawClone("AP");
    grRCNC_exp_relaxed.DrawClone("same P");
    // grRCNC_exp_tightened.DrawClone("same P");
    // grRCNC_exp_NoEnergyCorr.DrawClone("same P");

/*
******************************
* Mass vs E_beam in MC:
******************************
*/ 
    std::vector<Float_t> vM_MC_NoSmearing = {497.64, 497.657, 497.647, 497.642, 497.642, 497.642, 497.65, 497.642, 497.647, 497.633,};
    std::vector<Float_t> vM_MC_NoSmearingErr = {0.00331564, 0.00494999, 0.00522022, 0.00300056, 0.00274137, 0.0028886, 0.00562059, 0.00583294, 0.00588136, 0.00558898};

    // std::vector<Float_t> vM_MC_WithSmearing = {497.611, 497.619, 497.605, 497.604, 497.606, 497.604, 497.611, 497.608, 497.612, 497.604};
    std::vector<Float_t> vM_MC_WithSmearing = {497.623, 497.631, 497.618, 497.616, 497.619, 497.616, 497.623, 497.619, 497.625, 497.619};
    std::vector<Float_t> vM_MC_WithSmearingErr = {0.00343265, 0.00515374, 0.00547423, 0.00315212, 0.00283681, 0.00300401, 0.00584891, 0.00605414, 0.00610191, 0.00575733};
    
    TGraphErrors grMC_NoSmearing(vE.size(), vE.data(), vM_MC_NoSmearing.data(), zeroes.data(), vM_MC_NoSmearingErr.data());
    TGraphErrors grMC_WithSmearing(vE.size(), vE.data(), vM_MC_WithSmearing.data(), zeroes.data(), vM_MC_WithSmearingErr.data());

    TLine massKline(504.5, 497.614, 514.5, 497.614);
    massKline.SetLineStyle(2);
    massKline.SetLineWidth(4);
    massKline.SetLineColor(kRed);

    grMC_NoSmearing.SetTitle("Black -- without Corrections, Blue -- with Corrections, Red line --- mass in generator");
    grMC_NoSmearing.GetXaxis()->SetTitle("E_{beam}, MeV");
    grMC_NoSmearing.GetYaxis()->SetTitle("M^{(FullRec)}_{RC NC}, #frac{MeV}{c^{2}}");
    grMC_NoSmearing.GetYaxis()->SetTitleOffset(1.2);
    grMC_NoSmearing.SetName("NoSmear");
    grMC_NoSmearing.SetMarkerSize(2);
    grMC_NoSmearing.Fit("pol0", "ME");
    
    grMC_WithSmearing.SetTitle("Black -- without Corrections, Blue -- with Corrections, Red line --- mass in generator");
    grMC_WithSmearing.GetXaxis()->SetTitle("E_{beam}, MeV");
    grMC_WithSmearing.GetYaxis()->SetTitle("M^{(FullRec)}_{RC NC}, #frac{MeV}{c^{2}}");
    grMC_WithSmearing.GetYaxis()->SetTitleOffset(1.2);
    grMC_WithSmearing.SetName("Smear");
    grMC_WithSmearing.SetMarkerColor(kBlue);
    grMC_WithSmearing.SetLineColor(kBlue);
    grMC_WithSmearing.SetMarkerStyle(22);
    grMC_WithSmearing.SetMarkerSize(2);
    grMC_WithSmearing.Fit("pol0", "ME+");
    grMC_WithSmearing.GetFunction("pol0")->SetLineColor(kBlue);


    // grMC_NoSmearing.DrawClone("AP");
    // grMC_WithSmearing.DrawClone("AP");
    // massKline.DrawClone("same");

/*
******************************
* Mass vs E_beam in EXP for Cowboys and Sailors:
******************************
*/

    std::vector<Float_t> vM_Exp_cowboy = {497.466, 497.601, 497.564, 497.567, 497.575, 497.606, 497.593, 497.602, 497.577, 497.538};
    std::vector<Float_t> vMerr_Exp_cowboy = {0.0537, 0.0239, 0.0144, 0.0136, 0.0111, 0.0116, 0.0139, 0.0199, 0.0258, 0.0724};
    TGraphErrors grNC_exp_cowboy(vM_Exp_cowboy.size(), vE.data(), vM_Exp_cowboy.data(), zeroes.data(), vMerr_Exp_cowboy.data());

    grNC_exp_cowboy.SetTitle("Black -- Cowboy, Blue -- Sailor");
    grNC_exp_cowboy.GetXaxis()->SetTitle("E_{beam}, MeV");
    grNC_exp_cowboy.GetYaxis()->SetTitle("M_{K_{S}}, #frac{MeV}{c^{2}}");
    grNC_exp_cowboy.GetYaxis()->SetTitleOffset(1.2);
    grNC_exp_cowboy.SetName("Cowboy");
    grNC_exp_cowboy.SetMarkerSize(1);

    std::vector<Float_t> vM_Exp_sailor = {497.552, 497.575, 497.568, 497.562, 497.559, 497.581, 497.58, 497.56, 497.596, 497.673};
    std::vector<Float_t> vMerr_Exp_sailor = {0.0616, 0.0235, 0.0144, 0.0137, 0.0111, 0.0116, 0.0141, 0.0187, 0.0262, 0.0761};
    TGraphErrors grNC_exp_sailor(vM_Exp_sailor.size(), vE.data(), vM_Exp_sailor.data(), zeroes.data(), vMerr_Exp_sailor.data());

    grNC_exp_sailor.SetTitle("Black -- Cowboy, Blue -- Sailor");
    grNC_exp_sailor.GetXaxis()->SetTitle("E_{beam}, MeV");
    grNC_exp_sailor.GetYaxis()->SetTitle("MM_{K_{S}}, #frac{MeV}{c^{2}}");
    grNC_exp_sailor.GetYaxis()->SetTitleOffset(1.2);
    grNC_exp_sailor.SetName("Sailor");
    grNC_exp_sailor.SetMarkerColor(kBlue);
    grNC_exp_sailor.SetMarkerSize(1);

    // grNC_exp_cowboy.DrawClone("AP");
    // grNC_exp_sailor.DrawClone("P Same");

    
/*
******************************
* Mass_Cowboys - Mass_Sailors vs E_beam in EXP:
******************************
*/
    std::vector<Float_t> vDiffM_Exp_cowboy_sailor = {-0.086, 0.026, -0.004, 0.005, 0.016, 0.025, 0.013, 0.042, -0.019, -0.135};
    std::vector<Float_t> vDiffMerr_Exp_cowboy_sailor = {0.0811, 0.032, 0.0177, 0.0165, 0.0121, 0.013, 0.0171, 0.0254, 0.0354, 0.1046};
    TGraphErrors grNC_Exp_cowboy_sailor_diff(vDiffM_Exp_cowboy_sailor.size(), vE.data(), vDiffM_Exp_cowboy_sailor.data(), zeroes.data(), vDiffMerr_Exp_cowboy_sailor.data());
    grNC_Exp_cowboy_sailor_diff.SetTitle("Cowboy Sailor Difference");
    grNC_Exp_cowboy_sailor_diff.GetXaxis()->SetTitle("E_{beam}, MeV");
    grNC_Exp_cowboy_sailor_diff.GetYaxis()->SetTitle("#Delta M, #frac{MeV}{c^{2}}");
    grNC_Exp_cowboy_sailor_diff.GetYaxis()->SetTitleOffset(1.2);
    grNC_Exp_cowboy_sailor_diff.SetName("CowboySailorDiff");
    grNC_Exp_cowboy_sailor_diff.SetMarkerSize(1);

    // grNC_Exp_cowboy_sailor_diff.DrawClone("AP");
    
    c.DrawClone();
    return 0;
}