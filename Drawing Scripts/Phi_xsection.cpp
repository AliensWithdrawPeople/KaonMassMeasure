#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TAxis.h"

int Phi_xsection()
{
    std::vector<double> zeroes(100, 0.0);

    std::vector<double> energy = {500.9293, 502.9799, 504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 511.4479, 513.864};
    std::vector<double> energyErr = {0.004, 0.002, 0.005, 0.004, 0.004, 0.004};

    std::vector<double> xsec_vis = {0.568, 1.6317, 3.0119, 37.9855, 68.9547, 121.7992, 180.229, 178.9924, 146.8741, 90.1184, 70.3458, 27.91};
    std::vector<double> xsec_vis_err = {0.033, 0.036, 0.082, 0.273, 0.237, 0.393, 0.305, 0.345, 0.424, 0.369, 0.378, 0.238};

    std::vector<double> efficiency = {0.100742, 0.14724, 0.168875, 0.189286, 0.191459, 0.194718, 0.19831, 0.198671, 0.199915, 0.202736, 0.202236, 0.205726};

    TGraphErrors grVisXsec(energy.size(), energy.data(), xsec_vis.data(), zeroes.data(), xsec_vis_err.data());
    grVisXsec.GetYaxis()->SetTitle("#sigma_{vis}, nb");
    grVisXsec.GetXaxis()->SetTitle("E_{avg}, MeV");
    grVisXsec.SetName("VisXsec");
    grVisXsec.SetMarkerColor(kRed);
    grVisXsec.SetMarkerSize(1.5);
    grVisXsec.DrawClone("AP");


    TGraphErrors grEff(energy.size(), energy.data(), efficiency.data(), zeroes.data(), zeroes.data());
    grEff.GetYaxis()->SetTitle("#eps");
    grEff.GetXaxis()->SetTitle("E_{avg}, MeV");
    grEff.SetName("eff");
    grEff.SetMarkerColor(kBlue);
    grEff.SetMarkerSize(1.5);
    // grEff.DrawClone("AP");

    return 0;
}