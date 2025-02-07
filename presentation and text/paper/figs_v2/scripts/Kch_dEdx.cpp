#include <TFile.h>
#include <TAxis.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TTree.h>
#include <TF1.h>


void Kch_dEdx()
{
    TCanvas canvas("canv", "", 800, 600);

    auto tree = TFile::Open("./Kch_dEdx.root")->Get<TTree>("kChargedTree");
    auto hist = new TH2D("hist", "", 100, 40, 400, 100, 0, 3e4);
    tree->Draw("tdedx:tptot >> hist");

    TLine cut1(60, 7e3, 150, 7e3);
    TLine cut2(150, 7e3, 150, 3e4);
    TLine cut3(60, 7e3, 60, 3e4);
    cut1.SetLineStyle(2);
    cut2.SetLineStyle(2);
    cut3.SetLineStyle(2);

    cut1.SetLineWidth(3);
    cut2.SetLineWidth(3);
    cut3.SetLineWidth(3);

    TLatex lat; 
    double x = 0.25;
    double y = 0.85;
    lat.SetNDC();
    lat.SetTextFont(72);
    double delx = 0.130*696*gPad->GetWh()/(472*gPad->GetWw());

    hist->GetYaxis()->SetTitle("dE/dx, a.u");
    hist->GetXaxis()->SetTitle("P, MeV");
    hist->GetYaxis()->SetMaxDigits(3);


    hist->DrawClone("P");
    cut1.DrawClone("same");
    cut2.DrawClone("same");
    cut3.DrawClone("same");

    lat.DrawLatex(x, y, "CMD-3");
    canvas.DrawClone();

    canvas.SaveAs("C:/work/Science/KaonMassMeasure/presentation and text/paper/figs_v2/Kch_dEdx.pdf");
}