#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLatex.h"

void workflow_check_mc()
{    
    TCanvas canvas("canv", "", 800, 600);

    std::vector<Float_t> zeroes(10, 0.0);
    std::vector<Float_t> energies = {504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 511.444, 513.864};
    std::vector<Float_t> mass = {497.612, 497.616, 497.614, 497.611, 497.614, 497.613, 497.622, 497.616, 497.611, 497.62};
    std::vector<Float_t> mass_err = {0.00286283, 0.00309153, 0.00299324, 0.00149237, 0.00145347, 0.00173928, 0.00433594, 0.00606485, 0.00743935, 0.0160282};

    TGraphErrors gr_mass(energies.size(), energies.data(), mass.data(), zeroes.data(), mass_err.data());

    TLine gen_mass(504.5, 497.614, 514.5, 497.614);
    gen_mass.SetLineStyle(2);
    // gen_mass.SetLineWidth(4);

    gr_mass.GetXaxis()->SetTitle("E_{beam}, MeV");
    gr_mass.GetYaxis()->SetTitle("M, MeV");
    gr_mass.GetYaxis()->SetTitleOffset(1.7);
    gr_mass.GetYaxis()->SetNdivisions(505);
    gr_mass.Fit("pol0", "ME+");
    // gr_mass.GetFunction("pol0")->SetLineStyle(kBlue);

    TLatex lat; 
    double x = 0.18;
    double y = 0.85;
    lat.SetNDC();
    lat.SetTextFont(72);
    double delx = 0.130*696*gPad->GetWh()/(472*gPad->GetWw());

    
    gr_mass.DrawClone("AP");
    gen_mass.DrawClone("same");
    lat.DrawLatex(x, y, "#splitline{#font[72]{CMD-3} #chi^{2}/ndf = 8.2/9}{M_{K^{0}_{S}} = 497.613 #pm 0.001 MeV}");
    canvas.DrawClone();

    canvas.SaveAs("C:/work/Science/KaonMassMeasure/presentation and text/paper/figs_v2/workflow_check_mc.pdf");
}