#include <TFile.h>
#include <TAxis.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLine.h>
#include <RooPlot.h>


void inv_mass_fit()
{
    auto inv_mass = TFile::Open("./m_inv_505.root")->Get<RooPlot>("M_inv_frame");

    TCanvas canvas("canv", "", 800, 600);
    TLatex lat; 
    double x = 0.2;
    double y = 0.85;
    lat.SetNDC();
    lat.SetTextFont(72);
    double delx = 0.130*696*gPad->GetWh()/(472*gPad->GetWw());

    inv_mass->GetYaxis()->SetTitle("Events/1.9 MeV");
    inv_mass->DrawClone("col");

    lat.DrawLatex(x, y, "CMD-3");
    canvas.SetLogy();
    canvas.DrawClone();

    canvas.SaveAs("C:/work/Science/KaonMassMeasure/presentation and text/paper/figs_v2/M_inv_fit.pdf");
}