#include <TCanvas.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TLine.h>
#include <TFile.h>
#include <TProfile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TF1.h>

#include <vector>
#include <memory>
#include <iostream>
#include <string>

int Mass_MC_Exp_DrawAll()
{
    auto filename_MC = [](std::string point) { return "C:/work/Science/BINP/Kaon Mass Measure/hists/MC/Hists_MC" + point + ".root"; };
    auto filename_Exp = [](std::string point) { return "C:/work/Science/BINP/Kaon Mass Measure/hists/Exp/Hists_Exp" + point + ".root"; };
    std::vector<std::string> points = {"505", "508", "508.5", "509", "509.5", "510", "511", "511.5", "514"};

    auto canv = new TCanvas("canv","c",800,600);
    canv->Divide(3, 3);
    for(int i = 0; i < points.size(); i++)
    {
        auto pad = canv->cd(i + 1);

        auto file = TFile::Open(filename_MC(points[i]).c_str());
        auto hist = file->Get<TH2D>("hMassVsKsTheta");
        auto pfx_MC = hist->ProfileX(("MC_pfx" + points[i]).c_str());
        pfx_MC->Rebin(5);
        pfx_MC->GetXaxis()->SetRangeUser(-0.6, 0.6);
        pfx_MC->GetYaxis()->SetRangeUser(497.5, 497.7);
        pfx_MC->SetTitle(("E_{beam} = " + points[i] + " MeV, Black -- MC, Blue -- Exp").c_str());

        pfx_MC->GetYaxis()->SetTitle("M_{K^{0}_{S}}, MeV/c^{2}");
        pfx_MC->GetXaxis()->SetTitle("#theta_{K^{0}_{S}} - #pi/2, rad");
        auto func = new TF1("func_MC", "pol0", -0.3, 0.3);


        auto file_Exp = TFile::Open(filename_Exp(points[i]).c_str());
        auto hist_exp = file_Exp->Get<TH2D>("hMassVsKsTheta");
        auto pfx_Exp = hist_exp->ProfileX(("Exp_pfx" + points[i]).c_str());
        pfx_Exp->Rebin(5);
        pfx_Exp->GetXaxis()->SetRangeUser(-0.6, 0.6);
        pfx_Exp->GetYaxis()->SetRangeUser(497.5, 497.7);

        pfx_Exp->SetTitle(("E_{beam} = " + points[i] + " MeV, Black -- MC, Blue -- Exp").c_str());
        pfx_Exp->SetMarkerColor(kBlue);
        pfx_Exp->GetYaxis()->SetTitle("M_{K^{0}_{S}}, MeV/c^{2}");
        pfx_Exp->GetXaxis()->SetTitle("#theta_{K^{0}_{S}} - #pi/2, rad");
        auto func_exp = new TF1("func_Exp", "pol0", -0.3, 0.3);
        func_exp->SetLineColor(kBlue);
        pfx_MC->Fit(func, "", "goff", -0.3, 0.3);
        pfx_Exp->Fit(func_exp, "", "goff", -0.3, 0.3);

        pfx_Exp->DrawClone("");
        pfx_MC->DrawClone("same");
        func->DrawClone("same");
        func_exp->DrawClone("same");
    }

    canv->DrawClone();
    return 0;
}