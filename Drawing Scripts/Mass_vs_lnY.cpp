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
    std::vector<std::string> points = {"505", "508", "508.5", "509", "509.5", "510", "511", "511.5", "514"};

    TLine lower_bound(-0.8, 490, 0.8, 490);
    TLine higher_bound(-0.8, 505, 0.8, 505);
    lower_bound.SetLineColor(kBlue);
    higher_bound.SetLineColor(kBlue);

    auto canv = new TCanvas("canv","c",800,600);
    canv->Divide(3, 3);

    for(int i = 0; i < points.size(); i++)
    {
        auto pad = canv->cd(i + 1);
        std::cout << i << " : " << points[i] << std::endl;
        auto file = TFile::Open(filename(points[i]).c_str());
        auto hist = file->Get<TH2D>("hMassVsKsTheta")->ProfileX(("pfx" + points[i]).c_str());
        // auto hist = file->Get<TH2D>("hMlnY");
        hist->SetTitle(("E_{beam} = " + points[i] + " MeV, Black -- MC, Blue -- Exp").c_str());
        hist->Rebin(15);
        // hist->Fit("pol0", "SQMEL", "", -0.3, 0.3);
        // hist->SetMarkerStyle(1);
        hist->GetYaxis()->SetRangeUser(497.3, 497.8);
        hist->GetXaxis()->SetRangeUser(-0.6, 0.6);

        hist->GetYaxis()->SetTitle("M_{K_{S}}, MeV/c^{2}   ");
        hist->GetXaxis()->SetTitle("#theta_{K_{S}} - #pi/2, rad");
        
        // hist->GetYaxis()->SetTitleOffset(0.9);
        hist->DrawClone();
        // lower_bound.Draw("same");
        // higher_bound.Draw("same");


        auto file_exp = TFile::Open(filename_exp(points[i]).c_str());
        auto hist_exp = file_exp->Get<TH2D>("hMassVsKsTheta")->ProfileX(("pfx_exp" + points[i]).c_str());
        hist_exp->SetTitle(("E_{beam} = " + points[i] + " MeV, EXP").c_str());
        hist_exp->Rebin(15);
        hist_exp->SetMarkerStyle(23);
        hist_exp->SetMarkerColor(kBlue);
        hist_exp->GetYaxis()->SetRangeUser(497.3, 497.8);
        hist_exp->GetXaxis()->SetRangeUser(-0.6, 0.6);
        hist_exp->DrawClone("same");

    }
    canv->DrawClone();
    return 0;
}