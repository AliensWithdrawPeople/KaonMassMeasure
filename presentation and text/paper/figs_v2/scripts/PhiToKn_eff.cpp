#include <TFile.h>
#include <TAxis.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLine.h>



void PhiToKn_eff()
{
    std::vector<double> zeroes(100, 0.0);
    std::vector<double> energy = {1.001859, 1.00596, 1.0096, 1.015724, 1.016808, 1.017914, 1.019056, 1.019912, 1.020916, 1.02207, 1.022896, 1.027728, 1.0338, 1.03979, 1.0498, 1.06};
    std::vector<double> efficiency = {0.25945, 0.25114, 0.2444, 0.2339, 0.23235, 0.23133, 0.23111, 0.22996, 0.22806, 0.22772, 0.22709, 0.22395, 0.22032, 0.21754, 0.21277, 0.2047};
    std::vector<double> efficiency_err = {0.00048, 0.00047, 0.00047, 0.00049, 0.00049, 0.00028, 0.00025, 0.00025, 0.00047, 0.00047, 0.00047, 0.00044, 0.00061, 0.0006, 0.00059, 0.00058};
    
    auto grEff = new TGraphErrors(energy.size(), energy.data(), efficiency.data(), zeroes.data(), efficiency_err.data());


    grEff->GetYaxis()->SetRangeUser(0.198, 0.275);
    grEff->GetYaxis()->SetNdivisions(505);
    grEff->GetYaxis()->SetTitle("#varepsilon_{MC}");
    grEff->GetXaxis()->SetTitle("E_{cms}, GeV");

    TCanvas canvas("canv", "", 800, 600);
    TLatex lat; 
    double x = 0.2;
    double y = 0.85;
    lat.SetNDC();
    lat.SetTextFont(72);
    double delx = 0.130*696*gPad->GetWh()/(472*gPad->GetWw());

    grEff->DrawClone("AP");

    lat.DrawLatex(x, y, "CMD-3");
    canvas.DrawClone();

    canvas.SaveAs("C:/work/Science/KaonMassMeasure/presentation and text/paper/figs_v2/PhiToKn_efficiency.pdf");
}