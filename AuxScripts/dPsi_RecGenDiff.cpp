#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TF1.h"
#include "TProfile.h"

#include <vector>
#include <memory>
#include <iostream>


int dPsi_RecGenDiff()
{
    std::vector<std::string> points = {"501", "503", "505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514", "517", "520", "525", "530"};
    std::string path = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/MC/KsKl_Smeared/";
    auto filename = [&path](std::string point) { return path + "MC" + point + ".root"; };

    std::vector<double> constPar = {};
    std::vector<double> constParErr = {};
    std::vector<double> constY = {};
    std::vector<double> constYerr = {};

    std::vector<double> diff = {};
    std::vector<double> diffErr = {};
    std::vector<double> energy = {505, 508, 508.5, 509, 509.5, 510, 510.5, 511, 511.5, 514, 517, 520, 525, 530};
    std::vector<double> energyErr(100, 0.0);

    TF1 formula("formula", "pol2", -2, 2, "");
    formula.SetParameters(0., 0., 0.0001);
    formula.FixParameter(1, 0.);

    auto file1  = new TFile("dPsi_RecDiff_out.root", "recreate");
    int counter = 0;
    for (const auto &point : points)
    {
        auto file = std::unique_ptr<TFile>(TFile::Open(filename(point).c_str()));
        auto ksTr = std::unique_ptr<TTree>(file->Get<TTree>("ksTree"));
        TH2D hPsiRecGenDiff("hPsiRecGenDiff", "", 250, -1.6, 1.6, 250, -1, 1);
        TH2D hY_RecGenDiff("hY_RecGenDiff", "", 250, -1.6, 1.6, 250, -1, 1);
        ksTr->Draw("ksdpsi_gen - ksdpsi : kstheta - TMath::Pi()/2 >> hPsiRecGenDiff", 
                    "nhitPos > 10 && nhitNeg > 10 && 1.1 < piThetaPos && piThetaPos < TMath::Pi() - 1.1 && 1.1 < piThetaNeg && piThetaNeg < TMath::Pi() - 1.1", "goff");
        ksTr->Draw("Y - Y_gen : kstheta - TMath::Pi()/2 >> hY_RecGenDiff", 
                    "nhitPos > 10 && nhitNeg > 10 && 1.1 < piThetaPos && piThetaPos < TMath::Pi() - 1.1 && 1.1 < piThetaNeg && piThetaNeg < TMath::Pi() - 1.1", "goff");
        
        auto psiPfx = hPsiRecGenDiff.ProfileX();
        psiPfx->SetName(("psiPfx" + std::to_string(counter)).c_str());
        psiPfx->SetTitle(("PsiRecGenDiff pfx" + point).c_str());
        TFitResultPtr resPsi = psiPfx->Fit("formula", "SQME", "goff", -0.5, 0.5);
        std::cout << point << " psi chi2 / ndf = " << resPsi->Chi2() << "/" << resPsi->Ndf() << std::endl;

        psiPfx->SetName(("psiPfx" + std::to_string(counter)).c_str());
        psiPfx->SetTitle(("PsiRecGenDiff pfx" + point).c_str());


        auto YPfx = hY_RecGenDiff.ProfileX();
        YPfx->Rebin(5);
        YPfx->SetName(("Ypfx" + std::to_string(counter)).c_str());
        YPfx->SetTitle(("YRecGenDiff pfx" + point).c_str());
        TFitResultPtr resY = YPfx->Fit("pol0", "SQME", "goff", -0.5, 0.5);
        std::cout << point << " Y chi2 / ndf = " << resY->Chi2() << "/" << resY->Ndf() << std::endl;

        YPfx->SetName(("Ypfx" + std::to_string(counter)).c_str());
        YPfx->SetTitle(("YRecGenDiff pfx" + point).c_str());



        constPar.push_back(resPsi->Parameter(0));
        diff.push_back(resPsi->Parameter(2));

        constParErr.push_back(resPsi->ParError(0));
        diffErr.push_back(resPsi->ParError(2));


        constY.push_back(resY->Parameter(0));
        constYerr.push_back(resY->ParError(0));

        std::cout << point << " MeV delta ksdpsi = " << resPsi->Parameter(0) << " + " << resPsi->Parameter(2) << " * x * x" << std::endl; 
        std::cout << point << " MeV delta Y = " << resY->Parameter(0) << "\n\n" << std::endl; 

        file1->cd();
        psiPfx->Write();
        YPfx->Write();

        counter++;
    }
    
    TGraphErrors grDiff(energy.size(), energy.data(), diff.data(), energyErr.data(), diffErr.data());
    TGraphErrors grConstPar(energy.size(), energy.data(), constPar.data(), energyErr.data(), constParErr.data());

    TGraphErrors grConstY(energy.size(), energy.data(), constY.data(), energyErr.data(), constYerr.data());

    grDiff.SetName("grDiff");
    grConstPar.SetName("grConstPar");

    grConstY.SetName("grConstY");

    grDiff.GetXaxis()->SetTitle("E_{beam}, MeV");
    grDiff.GetYaxis()->SetTitle("ksdpsi_gen - ksdpsi par[2], rad");
    grDiff.GetYaxis()->SetTitleOffset(1.);

    grConstPar.GetXaxis()->SetTitle("E_{beam}, MeV");
    grConstPar.GetYaxis()->SetTitle("ksdpsi_gen - ksdpsi par[0], rad");
    grConstPar.GetYaxis()->SetTitleOffset(1.);

    grConstY.GetXaxis()->SetTitle("E_{beam}, MeV");
    grConstY.GetYaxis()->SetTitle("Y - Y_gen");
    grConstY.GetYaxis()->SetTitleOffset(1.);

    grDiff.Write();
    grConstPar.Write();
    grConstY.Write();

    grDiff.DrawClone();

    file1->Close();

    return 0;
}