#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TFile.h"

int Phi_xsection()
{
    std::vector<double> zeroes(100, 0.0);

    std::vector<double> energy = {1.001859, 1.00596, 1.0096, 1.015724, 1.016808, 1.017914, 1.019056, 1.019912, 1.020916, 1.02207, 1.022896, 1.027728};
    std::vector<double> energySpread = {0.0003492, 0.0003625, 0.0003626, 0.0003937, 0.0003473, 0.0003556, 0.0003819, 0.0004158, 0.0004071, 0.0003791, 0.0003914, 0.0003672};

    // std::vector<double> energy = {1.015724, 1.016808, 1.017914, 1.019056, 1.019912, 1.020916};
    // std::vector<double> energySpread = {0.0003937, 0.0003473, 0.0003556, 0.0003819, 0.0004158, 0.0004071};


    // std::vector<double> xsec_vis = {219.29452, 379.46311, 661.0057, 991.93678, 987.60974, 783.37357};
    // std::vector<double> xsec_vis_err = {1.33853, 1.12915, 2.10807, 2.17229, 2.4543, 2.47224};


    std::vector<double> xsec_vis = {2.4649, 8.38149, 22.0801, 219.29452, 379.46311, 661.0057, 991.93678, 987.60974, 783.37357, 522.4324, 402.48297, 149.54976};
    std::vector<double> xsec_vis_err = {0.14137, 0.16368, 0.53163, 1.33853, 1.12915, 2.10807, 2.17229, 2.4543, 2.47224, 2.01235, 1.94623, 1.14083};


    std::vector<double> efficiency = {0.26591, 0.2577, 0.2517, 0.24325, 0.24028, 0.23805, 0.23718, 0.23539, 0.23547, 0.23485, 0.23438, 0.22889};
    // std::vector<double> efficiency = {0.25902, 0.25123, 0.24398, 0.24325, 0.24028, 0.23805, 0.23718, 0.23539, 0.23547, 0.23485, 0.23438, 0.22889};
    std::vector<double> efficiency_err = {0.000438, 0.000434, 0.000429, 0.000429, 0.000427, 0.000426, 0.000425, 0.000424, 0.000424, 0.000424, 0.000424, 0.00042};
    TFile *top = new TFile("vcs.root", "recreate");

    auto grVisXsec = new TGraphErrors(energy.size(), energy.data(), xsec_vis.data(), energySpread.data(), xsec_vis_err.data());
    grVisXsec->GetYaxis()->SetTitle("#sigma_{vis}, nb");
    grVisXsec->GetXaxis()->SetTitle("E_{avg}, MeV");
    grVisXsec->SetName("vcs");
    grVisXsec->SetMarkerColor(kRed);
    grVisXsec->SetMarkerSize(1.5);
    grVisXsec->DrawClone("AP");
    
    grVisXsec->Write();
    top->Write();
    top->Save();

    TGraphErrors grEff(energy.size(), energy.data(), efficiency.data(), zeroes.data(), efficiency_err.data());
    grEff.GetYaxis()->SetTitle("#varepsilon_{MC}");
    grEff.GetXaxis()->SetTitle("E_{cms}, GeV");
    grEff.SetName("eff");
    grEff.SetMarkerColor(kBlue);
    grEff.SetMarkerSize(1.5);
    grEff.DrawClone("AP");

    return 0;

    // [p0]/([p4]*(x-[p5])*(x-[p5]))/(1-[p1]*exp([p2]*(x-[p3])))
}