#include <TCanvas.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TLine.h>
#include <TFile.h>
#include <TProfile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

#include <vector>
#include <memory>
#include <iostream>
#include <string>


int PhiXsec_DrawAll()
{
    auto filename = [](std::string point) { return "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/PhiXSection/Kn" + point + ".root"; };
    std::vector<std::string> points = {"501", "503", "505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514", "517", "520", "525", "530"};

    auto canv = new TCanvas("canv","c",800,600);
    canv->Divide(4, 4);
    for(int i = 0; i < points.size(); i++)
    // for(int i = 0; i < 16; i++)
    {
        // auto pad = canv->cd(i + 1 - 8);
        auto pad = canv->cd(i + 1);
        pad->SetLogy();
        std::cout << i << " : " << points[i] << std::endl;
        auto file = TFile::Open(filename(points[i]).c_str());
        auto hist = file->Get<TH1D>("hMass");
        hist->Rebin(3);
        hist->SetTitle(("E_{beam} = " + points[i] + " MeV").c_str());

        hist->GetYaxis()->SetTitleOffset(1.2);
        hist->GetYaxis()->SetTitle("Events/3 MeV");
        hist->GetXaxis()->SetTitle("M^{(inv)}_{K_{S}}, MeV/c^{2}");

        auto res2 = hist->Fit("pol0", "SQMEL0", "", 540, 578);
        TLine right_bkg(540, res2->Parameter(0), 580, res2->Parameter(0));
        right_bkg.SetLineColor(kBlue);
        right_bkg.SetLineStyle(2);

        // auto mean_bkg_lvl = (res1->Parameter(0) + res2->Parameter(0)) / 2.;
        auto mean_bkg_lvl = res2->Parameter(0);
        TLine mean_bkg(420, mean_bkg_lvl, 580, mean_bkg_lvl);
        mean_bkg.SetLineColor(kRed);
        hist->SetTitle(("Energy = " + points[i] + "MeV").c_str());
        hist->DrawClone("PE0");
        right_bkg.DrawClone("same");
        mean_bkg.DrawClone("same");
    }
    canv->DrawClone();
    return 0;
}