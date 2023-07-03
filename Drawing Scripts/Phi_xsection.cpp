#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TFile.h"

int Phi_xsection()
{
    std::vector<double> zeroes(100, 0.0);

    // std::vector<double> energy = {500.9293, 502.9799, 504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 511.4479, 513.864};
    // std::vector<double> energy = {1001.859, 1005.96, 1009.6, 1015.724, 1016.808, 1017.914, 1019.056, 1019.912, 1020.916, 1022.07, 1022.896, 1027.728};
    // std::vector<double> energySpread = {0.3492, 0.3625, 0.3626, 0.3937, 0.3473, 0.3556, 0.3819, 0.4158, 0.4071, 0.3791, 0.3914, 0.3672};

    std::vector<double> energy = {1.001859, 1.00596, 1.0096, 1.015724, 1.016808, 1.017914, 1.019056, 1.019912, 1.020916, 1.02207, 1.022896, 1.027728};
    std::vector<double> energySpread = {0.0003492, 0.0003625, 0.0003626, 0.0003937, 0.0003473, 0.0003556, 0.0003819, 0.0004158, 0.0004071, 0.0003791, 0.0003914, 0.0003672};

    // std::vector<double> xsec_vis = {0.568, 1.6317, 3.0119, 37.9855, 68.9547, 121.7992, 180.229, 178.9924, 146.8741, 90.1184, 70.3458, 27.91};
    // std::vector<double> xsec_vis_err = {0.033, 0.036, 0.082, 0.273, 0.237, 0.393, 0.305, 0.345, 0.424, 0.369, 0.378, 0.238};

    std::vector<double> xsec_vis = {5.6382, 11.0821, 17.8353, 200.678, 360.1538, 625.5156, 908.8247, 900.9489, 734.6828, 444.511, 347.8403, 135.666};
    std::vector<double> xsec_vis_err = {0.328, 0.248, 0.483, 1.441, 1.237, 2.017, 1.54, 1.735, 2.123, 1.821, 1.87, 1.156};

    std::vector<double> efficiency = {0.100742, 0.14724, 0.168875, 0.189286, 0.191459, 0.194718, 0.19831, 0.198671, 0.199915, 0.202736, 0.202236, 0.205726};
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

    TGraphErrors grEff(energy.size(), energy.data(), efficiency.data(), zeroes.data(), zeroes.data());
    grEff.GetYaxis()->SetTitle("#eps");
    grEff.GetXaxis()->SetTitle("E_{avg}, MeV");
    grEff.SetName("eff");
    grEff.SetMarkerColor(kBlue);
    grEff.SetMarkerSize(1.5);
    // grEff.DrawClone("AP");

    return 0;
}