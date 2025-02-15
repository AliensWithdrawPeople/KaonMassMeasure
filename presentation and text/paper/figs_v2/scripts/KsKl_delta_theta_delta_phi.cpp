#include <TFile.h>
#include <TAxis.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLine.h>


void KsKl_delta_theta_delta_phi()
{
    // Configuration: start
    auto is_MC = false;
    auto save = true;
    
    auto hist = TFile::Open("./exp509.5.root")->Get<TH2D>("hdThetadPhi");
    if(is_MC)
    { hist = TFile::Open("./MC509.5_KSKL_dTheta_dPhi.root")->Get<TH2D>("hdThetadPhi"); }
    // Configuration: stop
    
    TCanvas canvas("canv", "", 1000, 1000);
    TLatex lat; 
    double x = 0.2;
    double y = 0.85;
    lat.SetNDC();
    lat.SetTextFont(72);
    lat.SetTextColor(kWhite);
    double delx = 0.130*696*gPad->GetWh()/(472*gPad->GetWw());

    hist->Rebin2D(6, 6);
    hist->GetYaxis()->SetRangeUser(-1.5, 1.5);
    hist->GetXaxis()->SetRangeUser(1.5, 4.5);

    hist->GetYaxis()->SetNdivisions(505);
    hist->GetXaxis()->SetTitle("#Delta #varphi, rad");
    hist->GetYaxis()->SetTitle("#Delta #theta, rad");
    // hist->GetZaxis()->SetTitle("Entries/0.03#times0.03 rad^{2}");
    hist->GetZaxis()->SetMaxDigits(2);
    hist->DrawClone("colz");

    TLine upper_line(3.14 - 1, 1, 3.14 + 1, 1);
    TLine lower_line(3.14 - 1, -1, 3.14 + 1, -1);
    TLine left_line(3.14 - 1, -1, 3.14 - 1, 1);
    TLine right_line(3.14 + 1, -1, 3.14 + 1, 1);

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

    if(is_MC)
    { lat.DrawLatex(x, y, "CMD-3 Simulation"); }
    else
    { lat.DrawLatex(x, y, "CMD-3"); }

    canvas.SetRightMargin(0.15);
    canvas.SetTopMargin(0.07);
    canvas.DrawClone();

    if(save)
    {
        if(is_MC)
        { canvas.SaveAs("C:/work/Science/KaonMassMeasure/presentation and text/paper/figs_v2/KsKl_delta_theta_delta_phi_MC509.5.pdf"); }
        else
        { canvas.SaveAs("C:/work/Science/KaonMassMeasure/presentation and text/paper/figs_v2/KsKl_delta_theta_delta_phi_EXP509.5.pdf"); }
    }
}