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

int DeltaThetaMC_DrawAll()
{
    auto filename = [](std::string point) { return "C:/work/Science/BINP/Kaon Mass Measure/hists/MC/Hists_MC" + point + ".root"; };
    std::vector<std::string> points = {"505", "508", "508.5", "509", "509.5", "505", "510", "510.5", "511", "511.5", "514", "510"};

    auto canv = new TCanvas("canv","c",800,600);
    canv->Divide(3, 3);
    // for(int i = 0; i < points.size(); i++)
    for(int i = 0; i < points.size(); i++)
    {
        auto pad = canv->cd(i + 1);
        std::cout << i << " : " << i+1 << " : " << points[i] << " : " << filename(points[i]) << std::endl;
        auto file = TFile::Open(filename(points[i]).c_str());
        auto hist = file->Get<TH2D>("hThetaDiffRecGenVsTheta_PionPos");
        auto pfx = hist->ProfileX(("pfx" + points[i]).c_str());
        pfx->Rebin(16);
        pfx->GetXaxis()->SetRangeUser(0.5, 2.6);
        pfx->SetTitle(("For #pi^{-}, E_{beam} = " + points[i] + " MeV").c_str());

        pfx->GetYaxis()->SetTitle("#theta^{(rec)}_{#pi^{-}} - #theta^{(gen)}_{#pi^{-}}, rad");
        pfx->GetXaxis()->SetTitle("#theta_{#pi^{-}}, rad");
        TF1 func("func", "pol1", 0.6, 2.5);
        func.SetParameter(0, -0.007);
        func.SetParameter(1, 0.005);
        // auto res1 = pfx->Fit("func", "SQMEL", "", 0.7, 2.4);
        // auto res2 = pfx->Fit("func", "SQMEL", "", 0.7, 2.4);
        pfx->DrawClone("");
    }
    canv->DrawClone();
    return 0;
}