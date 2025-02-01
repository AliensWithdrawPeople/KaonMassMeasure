#include <TFile.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLine.h>

void Ks_mom()
{
    // E = 509 MeV
    // auto y_label = "Events/0.5 MeV";
    // auto min_y = 5e1;
    // auto max_y = 5e5;
    // bool is_log_scale = true;
    // int rebin_factor = 1;
    // auto min_mom = 1.06958e+02 - 5 * 4.28612e+00;
    // auto max_mom = 1.06958e+02 + 5 * 4.28612e+00;

    // E = 517 MeV
    auto y_label = "Events/2 MeV";
    auto min_y = 0;
    auto max_y = 2.4e3;
    bool is_log_scale = false;
    int rebin_factor = 4;
    auto min_mom = 85;
    auto max_mom = 1.39491e+02 + 5 * 4.76301e+00;
    auto hKsMom = TFile::Open("./Kn517.root")->Get<TH1D>("hKsMom");

    TCanvas canvas("canv", "", 1000, 1000);
    TLatex lat; 
    double x = 0.2;
    double y = 0.85;
    lat.SetNDC();
    lat.SetTextFont(72);
    double delx = 0.130*696*gPad->GetWh()/(472*gPad->GetWw());

    hKsMom->Rebin(rebin_factor);

    hKsMom->GetXaxis()->SetTitle("P_{K^{0}_{S}}, MeV");
    hKsMom->GetYaxis()->SetTitle(y_label);
    hKsMom->GetListOfFunctions()->Remove(hKsMom->GetFunction("gaus"));

    hKsMom->GetYaxis()->SetMaxDigits(2);
    hKsMom->GetYaxis()->SetNdivisions(505);
    hKsMom->GetYaxis()->SetRangeUser(min_y, max_y);
    hKsMom->GetXaxis()->SetRangeUser(47., 255.);

    hKsMom->DrawClone();

    TLine left_line(min_mom, min_y, min_mom, 0.8 * max_y);
    TLine right_line(max_mom, min_y, max_mom, 0.8 * max_y);

    left_line.SetLineWidth(4);
    right_line.SetLineWidth(4);
    left_line.SetLineStyle(2);
    right_line.SetLineStyle(2);
    left_line.SetLineColor(kRed);
    right_line.SetLineColor(kRed);

    left_line.DrawClone("same");
    right_line.DrawClone("same");

    lat.DrawLatex(x, y, "CMD-3");

    if(is_log_scale)
    { canvas.SetLogy(); }
    canvas.DrawClone();

    canvas.SaveAs("C:/work/Science/KaonMassMeasure/presentation and text/paper/figs_v2/Ks_Mom_Exp517.pdf");
}