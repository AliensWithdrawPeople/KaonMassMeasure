#include <TFile.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLine.h>
#include <RooPlot.h>


void kaon_energy_spectrum()
{
    auto hist = TFile::Open("./Hists_MC514_std_noEsmear.root")->Get<TH1D>("hEnergySpectrum");
    hist->Rebin(8);

    TCanvas canvas("canv", "", 800, 600);
    TLatex lat; 
    double x = 0.2;
    double y = 0.85;
    lat.SetNDC();
    lat.SetTextFont(72);
    double delx = 0.130*696*gPad->GetWh()/(472*gPad->GetWw());

    hist->GetXaxis()->SetRangeUser(505, 515);
    hist->GetYaxis()->SetTitle("Events/0.08 MeV");
    hist->GetXaxis()->SetTitle("E_{K}, MeV");
    hist->Draw();

    lat.DrawLatex(x, y, "CMD-3");
    canvas.SetLogy();
    canvas.DrawClone();

    canvas.SaveAs("C:/work/Science/KaonMassMeasure/presentation and text/paper/figs_v2/kaon_energy_spectrum.pdf");
}