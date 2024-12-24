#include <TCanvas.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TLine.h>
#include <TFile.h>
#include <TProfile.h>

#include <vector>
#include <memory>
#include <iostream>
#include <string>


int Mass_vs_lnY()
{
    auto filename = [](std::string point) { return "C:/work/Science/BINP/Kaon Mass Measure/hists/MC/Hists_MC" + point + ".root"; };
    auto filename_exp = [](std::string point) { return "C:/work/Science/BINP/Kaon Mass Measure/hists/Exp/Hists_Exp" + point + ".root"; };
    std::vector<std::string> points = {"505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514"};

    auto canv = new TCanvas("canv","c",800,600);
    canv->Divide(2, 2);

    // for(int i = 0; i < points.size(); i++)
    for(int i = 4; i < 8; i++)
    {
        auto pad = canv->cd(i - 4 + 1);
        std::cout << i << " : " << points[i] << std::endl;
        auto file = TFile::Open(filename(points[i]).c_str());
        auto hist = file->Get<TProfile>("hMlnYpfx");
        // hist->SetTitle(("E_{beam} = " + points[i] + " MeV, Black -- MC, Blue -- Exp").c_str());
        hist->SetTitle(("E_{beam} = " + points[i] + " MeV").c_str());
        hist->Fit("pol0", "SMQE", "", -0.27, 0.27);
        // hist->SetMarkerStyle(1);
        hist->GetYaxis()->SetRangeUser(497.5, 497.8);
        hist->GetXaxis()->SetRangeUser(-0.35, 0.35);

        hist->GetYaxis()->SetTitle("M_{K_{S}}, MeV/c^{2}   ");
        hist->GetYaxis()->SetTitleOffset(1);
        hist->GetXaxis()->SetTitle("ln(Y)");
        hist->DrawClone();

        auto file_exp = TFile::Open(filename_exp(points[i]).c_str());
        auto hist_exp = file_exp->Get<TProfile>("hMlnYpfx");
        hist_exp->SetTitle(("E_{beam} = " + points[i] + " MeV, EXP").c_str());
        hist_exp->SetMarkerStyle(23);
        hist_exp->SetMarkerColor(kBlue);
        hist_exp->GetYaxis()->SetRangeUser(497.3, 497.8);
        hist_exp->GetXaxis()->SetRangeUser(-0.6, 0.6);
        // hist_exp->DrawClone("same");
    }
    return 0;
}