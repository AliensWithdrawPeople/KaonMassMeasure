#include <TFile.h>
#include <TAxis.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLine.h>


void M_vs_lnY_MC509()
{
    auto hMlnY = TFile::Open("./Hists_MC509_noCorr.root")->Get<TH2D>("hMlnY");

    TCanvas canvas("canv", "", 1000, 1000);
    TLatex lat; 
    double x = 0.2;
    double y = 0.85;
    lat.SetNDC();
    lat.SetTextFont(72);
    double delx = 0.130*696*gPad->GetWh()/(472*gPad->GetWw());

    hMlnY->Rebin2D(4, 4);
    hMlnY->GetYaxis()->SetRangeUser(485, 508);
    hMlnY->GetXaxis()->SetRangeUser(-0.7, 0.7);

    hMlnY->GetYaxis()->SetNdivisions(505);
    hMlnY->GetXaxis()->SetTitle("ln(Y)");
    hMlnY->GetYaxis()->SetTitle("M, MeV");
    hMlnY->DrawClone("col");

    TLine upper_line(-0.3, 505, 0.3, 505);
    TLine lower_line(-0.3, 490, 0.3, 490);
    TLine left_line(-0.3, 490, -0.3, 505);
    TLine right_line(0.3, 490, 0.3, 505);

    upper_line.SetLineWidth(4);
    lower_line.SetLineWidth(4);
    left_line.SetLineWidth(4);
    right_line.SetLineWidth(4);

    upper_line.SetLineStyle(2);
    lower_line.SetLineStyle(2);
    left_line.SetLineStyle(2);
    right_line.SetLineStyle(2);

    upper_line.SetLineColor(kRed);
    lower_line.SetLineColor(kRed);
    left_line.SetLineColor(kRed);
    right_line.SetLineColor(kRed);

    upper_line.DrawClone("same");
    lower_line.DrawClone("same");
    left_line.DrawClone("same");
    right_line.DrawClone("same");

    lat.DrawLatex(x, y, "CMD-3");
    canvas.DrawClone();

    canvas.SaveAs("C:/work/Science/KaonMassMeasure/presentation and text/paper/figs_v2/M_vs_lnY_MC509.pdf");
}