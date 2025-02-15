#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLatex.h"

void Mass_Exp()
{    
    TCanvas canvas("canv", "", 800, 600);

    std::vector<Float_t> zeroes(10, 0.0);
    std::vector<Float_t> energies = {504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 511.444, 513.864};
    std::vector<Float_t> mass = {497.535, 497.588, 497.577, 497.575, 497.563, 497.588, 497.602, 497.586, 497.59, 497.609};
    std::vector<Float_t> mass_err = {0.046, 0.018, 0.012, 0.011, 0.01, 0.01, 0.011, 0.015, 0.019, 0.053};

    TGraphErrors gr_mass(energies.size(), energies.data(), mass.data(), zeroes.data(), mass_err.data());

    gr_mass.GetXaxis()->SetTitle("E_{beam}, MeV");
    gr_mass.GetYaxis()->SetTitle("M, MeV");
    gr_mass.GetYaxis()->SetTitleOffset(1.7);
    gr_mass.GetYaxis()->SetNdivisions(505);
    gr_mass.Fit("pol0", "ME+");

    TLatex lat; 
    double x = 0.18;
    double y = 0.85;
    lat.SetNDC();
    lat.SetTextFont(72);
    double delx = 0.130*696*gPad->GetWh()/(472*gPad->GetWw());

    
    gr_mass.DrawClone("AP");
    lat.DrawLatex(x, y, "#splitline{#font[72]{CMD-3} #chi^{2}/ndf = 9.5/9}{M_{K^{0}_{S}} = 497.582 #pm 0.004 MeV}");
    canvas.DrawClone();

    canvas.SaveAs("C:/work/Science/KaonMassMeasure/presentation and text/paper/figs_v2/Mass_Exp.pdf");
}