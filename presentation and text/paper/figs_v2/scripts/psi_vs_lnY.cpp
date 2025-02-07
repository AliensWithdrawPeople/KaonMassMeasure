#include <TFile.h>
#include <TAxis.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLine.h>
#include <RooPlot.h>


void psi_vs_lnY()
{
    auto hist = TFile::Open("./Hists_Exp510.5_psi_vs_lnY.root")->Get<TH2D>("hPsilnY");
    hist->Rebin2D(1, 8);

    TCanvas canvas("canv", "", 800, 600);
    TLatex lat; 
    double x = 0.2;
    double y = 0.85;
    lat.SetNDC();
    lat.SetTextFont(72);
    double delx = 0.130*696*gPad->GetWh()/(472*gPad->GetWw());

    hist->GetYaxis()->SetRangeUser(2.4, 3.1);

    hist->GetYaxis()->SetTitle("#psi, rad");
    hist->GetXaxis()->SetTitle("ln(Y)");
    hist->Draw("scat");

    lat.DrawLatex(x, y, "CMD-3");
    canvas.DrawClone();

    canvas.SaveAs("C:/work/Science/KaonMassMeasure/presentation and text/paper/figs_v2/psi_vs_lnY.pdf");
}