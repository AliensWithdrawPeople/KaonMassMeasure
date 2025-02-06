#include <TFile.h>
#include <TAxis.h>
#include <TF1.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TGraphErrors.h>

void energy_drift()
{
    // Setup: start
    auto min_y = 509.65;
    auto max_y = 510.18;

    auto min_x = 61360;
    auto max_x = 61880;
    // Setup: end
    auto raw_canvas = TFile::Open("./energy_drift.root")->Get<TCanvas>("c1");

    auto emeas = (TGraphErrors* )raw_canvas->GetListOfPrimitives()->FindObject("grEmeas");
    auto comptonMeanEnergy = (TGraphErrors* )raw_canvas->GetListOfPrimitives()->FindObject("grComptonMeanEnergy");
    auto KchEnergy = (TGraphErrors* )raw_canvas->GetListOfPrimitives()->FindObject("grKchEnergy");
    
    emeas->GetXaxis()->SetTitle("N_{run}");
    emeas->GetYaxis()->SetTitle("E, MeV");
    emeas->GetXaxis()->SetRangeUser(min_x, max_x);
    emeas->GetYaxis()->SetRangeUser(min_y, max_y);

    emeas->SetMarkerStyle(20);
    KchEnergy->SetMarkerStyle(4);
    comptonMeanEnergy->SetLineWidth(4);

    TF1 emeas_approx1("energy_approx1", "364.312 / (x[0] - 61332.6) / (x[0] - 61332.6) + 509.925", 61380, 61560);
    TF1 energy_approx1("energy_approx1", "364.312 / (x[0] - 61332.6) / (x[0] - 61332.6) + 509.925 - 0.2", 61380, 61560);

    TF1 emeas_approx2("energy_approx2", "598.053 / (x[0] - 61479.6) / (x[0] - 61479.6) + 509.933", 61560, 61689);
    TF1 energy_approx2("energy_approx2", "598.053 / (x[0] - 61479.6) / (x[0] - 61479.6) + 509.933 - 0.200", 61560, 61689);

    TF1 emeas_approx3("energy_approx3", "0.0101 * sin(0.0397 * (x[0]-61704.6)) + 509.960", 61689, 61856);
    TF1 energy_approx3("energy_approx3", "0.0101 * sin(0.0397 * (x[0]-61704.6)) + 509.960 - 0.2", 61689, 61856);

    emeas_approx1.SetLineStyle(2);
    emeas_approx2.SetLineStyle(2);
    emeas_approx3.SetLineStyle(2);
    
    energy_approx1.SetLineStyle(9);
    energy_approx2.SetLineStyle(9);
    energy_approx3.SetLineStyle(9);

    emeas_approx1.SetLineWidth(3);
    energy_approx1.SetLineWidth(3);
    emeas_approx2.SetLineWidth(3);
    energy_approx2.SetLineWidth(3);
    emeas_approx3.SetLineWidth(3);
    energy_approx3.SetLineWidth(3);

    TCanvas canvas("canv", "", 800, 600);

    emeas->DrawClone("AP");
    KchEnergy->DrawClone("P same");
    comptonMeanEnergy->DrawClone("same");

    emeas_approx1.DrawClone("same");
    energy_approx1.DrawClone("same");

    emeas_approx2.DrawClone("same");
    energy_approx2.DrawClone("same");

    emeas_approx3.DrawClone("same");
    energy_approx3.DrawClone("same");

    TLatex lat; 
    double x = 0.2;
    double y = 0.85;
    lat.SetNDC();
    lat.SetTextFont(72);
    double delx = 0.130*696*gPad->GetWh()/(472*gPad->GetWw());

    lat.DrawLatex(x, y, "CMD-3");
    // canvas.SetLogy();
    canvas.DrawClone();

    canvas.SaveAs("C:/work/Science/KaonMassMeasure/presentation and text/paper/figs_v2/Energy_drift_510MeV.pdf");
}