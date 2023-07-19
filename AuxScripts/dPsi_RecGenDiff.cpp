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
    auto filename = [&path](std::string point) { return (path + "MC" + point + ".root").c_str(); };

    std::vector<double> constPar = {};
    std::vector<double> constParErr = {};

    std::vector<double> diff = {};
    std::vector<double> diffErr = {};
    std::vector<double> energy = {501, 503, 505, 508, 508.5, 509, 509.5, 510, 510.5, 511, 511.5, 514, 517, 520, 525, 530};
    std::vector<double> energyErr(100, 0.0);

    TF1 formula("formula", "[0] + [1] * x * x", -2, 2, "");
    formula.SetParameters(0., 0.0001);
    auto file1  = new TFile("dPsi_RecDiff_out.root", "recreate");
    int counter = 0;
    for (const auto &point : points)
    {
        auto file = std::unique_ptr<TFile>(TFile::Open(filename(point)));
        auto ksTr = std::unique_ptr<TTree>(file->Get<TTree>("ksTree"));
        TH2D hist("hist", "", 250, -1.6, 1.6, 250, -1, 1);
        ksTr->Draw("ksdpsi_gen - ksdpsi : kstheta - TMath::Pi()/2 >> hist", "fabs(log(Y)) < 0.3 && nhitPos > 10 && nhitNeg > 10", "goff");
        auto pfx = hist.ProfileX();
    
        pfx->SetName(("pfx" + std::to_string(counter)).c_str());
        pfx->SetTitle(("pfx" + point).c_str());
        TFitResultPtr res = pfx->Fit("formula", "SQME", "goff", -1, 1);

        std::cout << point << " chi2 / ndf = " << res->Chi2() << "/" << res->Ndf() << std::endl;

        constPar.push_back(res->Parameter(0));
        diff.push_back(res->Parameter(1));

        constParErr.push_back(res->ParError(0));
        diffErr.push_back(res->ParError(1));

        file1->cd();
        pfx->Write();

        counter++;
    }
    
    TGraphErrors grDiff(energy.size(), energy.data(), diff.data(), energyErr.data(), diffErr.data());
    TGraphErrors grConstPar(energy.size(), energy.data(), constPar.data(), energyErr.data(), constParErr.data());

    grDiff.SetName("grDiff");
    grConstPar.SetName("grConstPar");

    grDiff.GetXaxis()->SetTitle("E_{beam}, MeV");
    grDiff.GetYaxis()->SetTitle("ksdpsi_gen - ksdpsi par[2], rad");
    grDiff.GetYaxis()->SetTitleOffset(1.);

    grConstPar.GetXaxis()->SetTitle("E_{beam}, MeV");
    grConstPar.GetYaxis()->SetTitle("ksdpsi_gen - ksdpsi par[0], rad");
    grConstPar.GetYaxis()->SetTitleOffset(1.);

    grDiff.Write();
    grConstPar.Write();

    grDiff.DrawClone();

    file1->Close();

    return 0;
}