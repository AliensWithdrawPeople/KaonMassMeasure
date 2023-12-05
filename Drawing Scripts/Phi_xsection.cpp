#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TFile.h"
#include "TCanvas.h"
#include <TLegend.h>

int Phi_xsection()
{
    TCanvas canv("canv");
/************ phi -> KsKl cross section vs energy ************/
    std::vector<double> zeroes(100, 0.0);

    std::vector<double> energy = {1.001859, 1.00596, 1.0096, 1.015724, 1.016808, 1.017914, 1.019056, 1.019912, 1.020916, 1.02207, 1.022896, 1.027728, 1.0338, 1.03979, 1.0498, 1.06};
    std::vector<double> energySpread = {0.0003492, 0.0003625, 0.0003626, 0.0003937, 0.0003473, 0.0003556, 0.0003819, 0.0004158, 0.0004071, 0.0003791, 0.0003914, 0.0003672, 0.00036, 0.00041, 0.000417, 0.000378};

    // Old FF
    // std::vector<double> xsec_vis = {2.52409, 8.54869, 22.92775, 219.98564, 378.05288, 657.75841, 986.31979, 973.88074, 778.15612, 518.83165, 398.13079, 145.33342, 75.84839, 49.53809, 32.36666, 23.08197};
    // std::vector<double> xsec_vis_err = {0.143951, 0.166218, 0.54809, 1.360034, 1.128901, 2.044079, 2.033836, 2.35264, 2.098364, 1.980744, 1.932598, 1.130711, 0.799666, 0.682457, 0.506897, 0.455273};

    // New FF: Const bkg fit
    std::vector<double> xsec_vis2 = {2.55493, 8.57642, 22.91284, 220.21836, 378.98208, 657.72612, 984.25605, 974.97838, 778.57844, 518.9548, 397.77018, 145.34911, 75.96123, 49.65799, 32.40645, 23.18225};
    std::vector<double> xsec_vis_err2 = {0.144931, 0.167673, 0.552109, 1.487758, 1.418073, 2.096495, 1.879297, 2.051185, 2.637251, 2.267091, 2.161583, 1.197029, 0.838923, 0.703349, 0.517876, 0.462314};
    auto grVisXsec2 = new TGraphErrors(energy.size(), energy.data(), xsec_vis2.data(), energySpread.data(), xsec_vis_err2.data());

    // New FF: MC + pol2 fit
    std::vector<double> xsec_vis = {2.51468, 8.67344, 23.48487, 224.97181, 389.29359, 674.34957, 1016.52104, 999.00898, 795.56775, 532.48244, 409.77972, 149.8607, 78.52898, 51.27216, 33.31665, 23.73516};
    std::vector<double> xsec_vis_err = {0.143783, 0.168631, 0.559077, 1.506623, 1.44527, 2.131825, 1.926259, 2.088932, 2.681391, 2.307578, 2.203054, 1.217629, 0.854318, 0.715368, 0.525421, 0.467961};
    
    TFile *top = new TFile("C:/work/Science/BINP/Kaon Mass Measure/PhiMesonFit/vcs_sig_mc_bkg_pol2.root", "recreate");

    auto grVisXsec = new TGraphErrors(energy.size(), energy.data(), xsec_vis.data(), energySpread.data(), xsec_vis_err.data());
    grVisXsec->GetYaxis()->SetTitle("#sigma_{vis}, nb");
    grVisXsec->GetXaxis()->SetTitle("E_{avg}, MeV");
    grVisXsec->SetName("vcs");
    grVisXsec->SetMarkerColor(kRed);
    grVisXsec2->SetMarkerColor(kBlue);
    grVisXsec->SetMarkerSize(0.9);
    grVisXsec->DrawClone("AP");
    grVisXsec2->DrawClone("P same");
    
    grVisXsec->Write();
    top->Write();
    top->Save();

/************ Efficiency vs energy ************/
    // Outdated (old Kaon FF in MCGPJ) efficiency.
    std::vector<double> efficiency_old = {0.25922, 0.25138, 0.24396, 0.23382, 0.23255, 0.23099, 0.23027, 0.22983, 0.22779, 0.22737, 0.22637, 0.22331, 0.21986, 0.21707, 0.21214, 0.20391};
    std::vector<double> efficiency_err_old = {0.00044, 0.000436, 0.000431, 0.00045, 0.000424, 0.000423, 0.000423, 0.000451, 0.000298, 0.000421, 0.00042, 0.000437, 0.000416, 0.000414, 0.000531, 0.000524};
    
    auto grEff_old = new TGraphErrors(energy.size(), energy.data(), efficiency_old.data(), zeroes.data(), efficiency_err_old.data());
    grEff_old->SetMarkerColor(kRed);
    grEff_old->SetName("grEff_old");

    // Up-to-date efficiency.
    std::vector<double> efficiency = {0.25945, 0.25114, 0.2444, 0.2339, 0.23235, 0.23133, 0.23111, 0.22996, 0.22806, 0.22772, 0.22709, 0.22395, 0.22032, 0.21754, 0.21277, 0.2047};
    std::vector<double> efficiency_err = {0.00048, 0.00047, 0.00047, 0.00049, 0.00049, 0.00028, 0.00025, 0.00025, 0.00047, 0.00047, 0.00047, 0.00044, 0.00061, 0.0006, 0.00059, 0.00058};
    
    auto grEff = new TGraphErrors(energy.size(), energy.data(), efficiency.data(), zeroes.data(), efficiency_err.data());
    grEff->SetName("grEff");

    grEff->GetYaxis()->SetTitleOffset(1.2);
    grEff->GetYaxis()->SetTitle("#varepsilon_{MC}");
    grEff->GetXaxis()->SetTitle("E_{cms}, GeV");
    grEff->SetName("eff");
    grEff->SetMarkerColor(kBlue);
    grEff->SetTitle("Selection efficiency of e^{+}e^{-}#rightarrowK_{S}K_{L} events");

    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->AddEntry(grEff_old,"With old Kaon FF in MCGPJ","ep");
    legend->AddEntry(grEff,"With new Kaon FF in MCGPJ","ep");

    // grEff->DrawClone("AP");
    // grEff_old->DrawClone("P same");
    // legend->DrawClone();

/************ Background vs energy ************/
    // Const bkg; chi2 fit
    // std::vector<double> bckg = {0.5448, 0.32112, 0.61351, 0.71244, 1.23553, 1.67868, 2.71171, 2.57525, 2.08286, 1.31758, 1.41573, 0.88321, 0.8125,  0.76288, 0.71542, 0.84304};
    // std::vector<double> bckg_err = {0.03429, 0.01627, 0.04449, 0.0373, 0.03013, 0.04185, 0.03296, 0.03633, 0.04488, 0.04353, 0.0528,  0.04186, 0.03907, 0.03983, 0.03503, 0.03965};
    
    // Const bkg, only right sideband (range = (540, 576)); NLL fit
    // std::vector<double> bckg = {0.69951, 0.47594, 0.42957, 0.55441, 0.86295, 0.77104, 1.04199, 1.10664, 1.00505, 0.88549, 0.81114, 1.00279, 1.04139, 1.02398, 1.2634, 1.2652};
    // std::vector<double> bckg_err = {0.03887, 0.01981, 0.03722, 0.0329, 0.02517, 0.02833, 0.02039, 0.02377, 0.03114, 0.03567, 0.03994, 0.04461, 0.04424, 0.04616, 0.04659, 0.0486};

    // MC + bkg fit
    std::vector<double> bckg = {0.60188, 0.40148, 0.39382, 0.23414, 0.23994, 0.21656, 0.29378, 0.20309, 0.1284, 0.35398, 0.41952, 0.56358, 0.7418, 0.77178, 0.96882, 1.05519};
    std::vector<double> bckg_err = {0.03891, 0.02062, 0.04516, 0.05381, 0.02114, 0.09759, 0.00026, 0.00369, 0.04822, 0.04436, 0.00083, 0.04753, 0.05373, 0.05456, 0.05113, 0.05295};

    TGraphErrors grBckg(energy.size(), energy.data(), bckg.data(), zeroes.data(), bckg_err.data());
    grBckg.GetYaxis()->SetTimeOffset(1.2);
    grBckg.GetYaxis()->SetTitle("#sigma_{bkg}, nb");
    grBckg.GetXaxis()->SetTitle("E_{cms}, GeV");
    grBckg.SetName("bckg");
    grBckg.SetMarkerColor(kBlue);
    grBckg.DrawClone("AP");
    
    canv.DrawClone();
    return 0;
}