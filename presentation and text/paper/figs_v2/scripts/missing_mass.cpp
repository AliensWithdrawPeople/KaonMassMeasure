#include <TFile.h>
#include <TAxis.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLine.h>


void missing_mass()
{
    auto missing_mass = TFile::Open("./exp509.5.root")->Get<TH1D>("hMissingMass");

    TCanvas canvas("canv", "", 1000, 1000);
    TLatex lat; 
    double x = 0.2;
    double y = 0.85;
    lat.SetNDC();
    lat.SetTextFont(72);
    double delx = 0.130*696*gPad->GetWh()/(472*gPad->GetWw());

    // hMlnY->GetYaxis()->SetRangeUser(485, 508);
    // hMlnY->GetXaxis()->SetRangeUser(-0.7, 0.7);

    missing_mass->Rebin(4);

    // missing_mass->GetYaxis()->SetNdivisions(505);
    missing_mass->GetYaxis()->SetRangeUser(0.1, 3e5);
    missing_mass->GetXaxis()->SetRangeUser(60, 580);

    missing_mass->GetXaxis()->SetTitle("M_{miss}, MeV");
    missing_mass->GetYaxis()->SetTitle("Events/4 MeV");
    missing_mass->DrawClone("");

    TLine cut(350, 0.1, 350, 1e4);

    cut.SetLineWidth(4);
    cut.SetLineStyle(2);
    cut.SetLineColor(kRed);

    cut.DrawClone("same");

    lat.DrawLatex(x, y, "CMD-3");

    canvas.SetLogy();
    canvas.DrawClone();

    canvas.SaveAs("C:/work/Science/KaonMassMeasure/presentation and text/paper/figs_v2/Missing_mass.pdf");
}