#include <TFile.h>
#include <TAxis.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLine.h>


void M_vs_lnY_pfx_MC509()
{
    auto hMlnY = TFile::Open("./Hists_MC509_noCorr.root")->Get<TH2D>("hMlnY");
    hMlnY->GetYaxis()->SetRangeUser(485, 508);
    hMlnY->GetXaxis()->SetRangeUser(-0.7, 0.7);
    auto hMlnY_pfx = hMlnY->ProfileX();

    TCanvas canvas("canv", "", 1000, 1000);
    TLatex lat; 
    double x = 0.2;
    double y = 0.85;
    lat.SetNDC();
    lat.SetTextFont(72);
    double delx = 0.130*696*gPad->GetWh()/(472*gPad->GetWw());


    hMlnY_pfx->Rebin(15);
    hMlnY_pfx->GetYaxis()->SetRangeUser(496.5, 498);
    hMlnY_pfx->GetXaxis()->SetRangeUser(-0.6, 0.6);
    hMlnY_pfx->Fit("pol0", "SQLME", "", -0.3, 0.3);

    hMlnY_pfx->GetYaxis()->SetNdivisions(505);
    hMlnY_pfx->GetXaxis()->SetTitle("ln(Y)");
    hMlnY_pfx->GetYaxis()->SetTitle("M, MeV");
    hMlnY_pfx->DrawClone();

    TLine left_line(-0.3, 496.5, -0.3, 498);
    TLine right_line(0.3, 496.5, 0.3, 498);

    left_line.SetLineWidth(4);
    right_line.SetLineWidth(4);
    left_line.SetLineStyle(2);
    right_line.SetLineStyle(2);
    left_line.SetLineColor(kRed);
    right_line.SetLineColor(kRed);

    left_line.DrawClone("same");
    right_line.DrawClone("same");

    lat.DrawLatex(x, y, "CMD-3");
    canvas.DrawClone();

    canvas.SaveAs("C:/work/Science/KaonMassMeasure/presentation and text/paper/figs_v2/M_vs_lnY_pfx_MC509.pdf");
}