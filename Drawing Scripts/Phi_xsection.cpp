#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TFile.h"
#include "TCanvas.h"

int Phi_xsection()
{
/************ phi -> KsKl cross section vs energy ************/
    std::vector<double> zeroes(100, 0.0);

    std::vector<double> energy = {1.001859, 1.00596, 1.0096, 1.015724, 1.016808, 1.017914, 1.019056, 1.019912, 1.020916, 1.02207, 1.022896, 1.027728, 1.0338, 1.03979, 1.0498, 1.06};
    std::vector<double> energySpread = {0.0003492, 0.0003625, 0.0003626, 0.0003937, 0.0003473, 0.0003556, 0.0003819, 0.0004158, 0.0004071, 0.0003791, 0.0003914, 0.0003672, 0.00036, 0.00041, 0.000417, 0.000378};

    // std::vector<double> xsec_vis = {2.53683, 8.59876, 22.92589, 220.26869, 379.10793, 660.65136, 991.66518, 981.29475, 784.45118, 521.50047, 400.65176, 146.01895, 76.14629, 49.99045, 32.65257, 23.26833};
    std::vector<double> xsec_vis = {2.54701, 8.63326, 23.01787, 221.15247, 380.62902, 663.30208, 995.64402, 985.23198, 787.59862, 523.59288, 402.25929, 146.60482, 76.45181, 50.19102, 32.78358, 23.36169};
    std::vector<double> xsec_vis_err = {0.14526, 0.16786, 0.55025, 1.37041, 1.14811, 2.14582, 2.2168, 2.58212, 2.20359, 2.04632, 1.9765, 1.1414, 0.8061, 0.69147, 0.51343, 0.46079};

    TFile *top = new TFile("vcs_16pts_NewEff_ver4.root", "recreate");

    auto grVisXsec = new TGraphErrors(energy.size(), energy.data(), xsec_vis.data(), energySpread.data(), xsec_vis_err.data());
    grVisXsec->GetYaxis()->SetTitle("#sigma_{vis}, nb");
    grVisXsec->GetXaxis()->SetTitle("E_{avg}, MeV");
    grVisXsec->SetName("vcs");
    grVisXsec->SetMarkerColor(kRed);
    grVisXsec->SetMarkerSize(0.9);
    grVisXsec->DrawClone("AP");
    
    grVisXsec->Write();
    top->Write();
    top->Save();


/************ Efficiency vs energy ************/
    TCanvas canv("canv");

    // Up-to-date efficiency.
    // std::vector<double> efficiency = {0.25922, 0.25138, 0.24396, 0.23382, 0.23255, 0.23099, 0.23027, 0.2285, 0.22779, 0.22737, 0.22637, 0.23285, 0.21986, 0.21707, 0.21214, 0.20391};
    std::vector<double> efficiency = {0.25922, 0.25138, 0.24396, 0.23382, 0.23255, 0.23099, 0.23027, 0.22983, 0.22779, 0.22737, 0.22637, 0.22331, 0.21986, 0.21707, 0.21214, 0.20391};
    std::vector<double> efficiency_err = {0.00044, 0.000436, 0.000431, 0.00045, 0.000424, 0.000423, 0.000423, 0.000451, 0.000298, 0.000421, 0.00042, 0.000437, 0.000416, 0.000414, 0.000531, 0.000524};
    TGraphErrors grEff(energy.size(), energy.data(), efficiency.data(), zeroes.data(), efficiency_err.data());

    // Old data for E = 1.010, 1.016, 1.024 GeV (smaller data size). Now points have correct errors.
    std::vector<double> efficiency_old = {0.26591, 0.2577, 0.2517, 0.24325, 0.24028, 0.23805, 0.23718, 0.23539, 0.23547, 0.23485, 0.23438, 0.23102, 0.22891, 0.22683, 0.22242, 0.21371};
    std::vector<double> efficiency_old_err = {0.000444, 0.000439, 0.001371, 0.001355, 0.000429, 0.000428, 0.000427, 0.000426, 0.000426, 0.000426, 0.000426, 0.001331, 0.000422, 0.000421, 0.000541, 0.000533};
    TGraphErrors grEff_old(energy.size(), energy.data(), efficiency_old.data(), zeroes.data(), efficiency_old_err.data());

    // With additional points (E = 1.025, 1.026 GeV).
    std::vector<double> energy_bigger = {1.001859, 1.00596, 1.0096, 1.015724, 1.016808, 1.017914, 1.019056, 1.019912, 1.020916, 1.02207, 1.022896, 1.025, 1.026, 1.027728, 1.0338, 1.03979, 1.0498, 1.06};
    std::vector<double> efficiency_bigger = {0.26591, 0.2577, 0.25025, 0.2409, 0.24028, 0.23805, 0.23718, 0.23539, 0.23547, 0.23485, 0.23438, 0.23306, 0.23334, 0.23285, 0.22891, 0.22683, 0.22242, 0.21371};
    std::vector<double> efficiency_bigger_err = {0.000444, 0.000439, 0.000435, 0.000455, 0.000429, 0.000428, 0.000427, 0.000426, 0.000426, 0.000426, 0.000426, 0.000423, 0.000423, 0.000669, 0.000422, 0.000421, 0.000541, 0.000533};
    TGraphErrors grEff_bigger(energy_bigger.size(), energy_bigger.data(), efficiency_bigger.data(), zeroes.data(), efficiency_bigger_err.data());

    grEff.GetYaxis()->SetTitle("#varepsilon_{MC}");
    grEff.GetXaxis()->SetTitle("E_{cms}, GeV");
    grEff.SetName("eff");
    grEff.SetMarkerColor(kBlue);
    grEff.SetMarkerSize(0.9);
    // grEff.SetTitle(" Red -- old, Blue -- new (difference in points E = 1.010, 1.016, 1.024 GeV)");
    grEff.SetTitle("Blue -- new, Magenta -- with additional points (E = 1.025, 1.026 GeV)");

    grEff_old.GetYaxis()->SetTitle("#varepsilon_{MC}");
    grEff_old.GetXaxis()->SetTitle("E_{cms}, GeV");
    grEff_old.SetName("eff_old");
    grEff_old.SetMarkerColor(kRed);
    grEff_old.SetMarkerSize(0.9);

    grEff_bigger.GetYaxis()->SetTitle("#varepsilon_{MC}");
    grEff_bigger.GetXaxis()->SetTitle("E_{cms}, GeV");
    grEff_bigger.SetName("eff_bigger");
    grEff_bigger.SetMarkerColor(kMagenta);
    grEff_bigger.SetMarkerSize(0.9);
    
    grEff.DrawClone("AP");
    // grEff_old.DrawClone("P same");
    // grEff_bigger.DrawClone("P same");

    canv.DrawClone();


/************ Background vs energy ************/
    // std::vector<double> bckg = {0.55096, 0.33513, 0.63147, 0.88467, 0.8514, 0.83545, 0.95334};
    // std::vector<double> bckg_err = {0.09305, 0.04184, 0.14596, 0.10435, 0.1095, 0.0959, 0.10704};
    // std::vector<double> energyShort = {1.001859, 1.00596, 1.0096, 1.0338, 1.03979, 1.0498, 1.06};
    // TGraphErrors grBckg(energyShort.size(), energyShort.data(), bckg.data(), zeroes.data(), bckg_err.data());
    // grBckg.GetYaxis()->SetTitle("N_{bckg} / L, nb");
    // grBckg.GetXaxis()->SetTitle("E_{cms}, GeV");
    // grBckg.SetName("bckg");
    // grBckg.SetMarkerColor(kBlue);
    // grBckg.SetMarkerSize(1.);
    // grBckg.DrawClone("AP");

    return 0;
}