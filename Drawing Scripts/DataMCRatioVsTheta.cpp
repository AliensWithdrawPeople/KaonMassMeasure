#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <iostream>

int DataMCRatioVsTheta(std::string point)
{
    auto file1 = TFile::Open(("C:/work/Science/BINP/Kaon Mass Measure/tr_ph/PhiXSection/Kn" + point + ".root").c_str());
    auto file2 = TFile::Open(("C:/work/Science/BINP/Kaon Mass Measure/tr_ph/PhiXSection/MC/MC_Kn" + point + ".root").c_str());
    
    auto tr_EXP = file1->Get<TTree>("Kn");
    auto tr_MC = file2->Get<TTree>("Kn_MC");


    auto hExp = new TH1D("hExp", "", 80, 0, 3.2);
    auto hMC = new TH1D("hMC", "", 80, 0, 3.2);

    tr_EXP->Draw("piThetaPos >> hExp", "", "goff");
    tr_MC->Draw("piThetaPos >> hMC", "", "goff");

    hMC->Scale(hExp->Integral(int(1. / hExp->GetBinWidth(0)), int(2.14 / hExp->GetBinWidth(0))) / 
                hMC->Integral(int(1. / hMC->GetBinWidth(0)), int(2.14 / hMC->GetBinWidth(0))));
    
    hExp->Divide(hMC);

    auto res = hExp->Fit("pol0", "SQEM", "", 1.1, 2.0);
    auto histMean = hExp->Integral(int(1. / hExp->GetBinWidth(0)), int(2.14 / hExp->GetBinWidth(0))) / (int(2.14 / hExp->GetBinWidth(0)) - int(1. / hExp->GetBinWidth(0)));

    auto missingIntegral = hExp->Integral(int(1. / hExp->GetBinWidth(0)), int(1.1 / hExp->GetBinWidth(0))) ;
                            // + hExp->Integral(int(2.01 / hExp->GetBinWidth(0)), int(2.138 / hExp->GetBinWidth(0)));
    
    auto width =  int(1.1 / hExp->GetBinWidth(0)) - int(1. / hExp->GetBinWidth(0));// + int(2.138 / hExp->GetBinWidth(0)) - int(2.01 / hExp->GetBinWidth(0));
    auto diff = (res->Parameter(0) * width - missingIntegral);
    std::cout << diff << std::endl; 
    std::cout << hExp->GetBinContent(int(1.1 / hExp->GetBinWidth(0))) << std::endl; 
    std::cout << hExp->GetBinCenter(int(1.1 / hExp->GetBinWidth(0))) << std::endl; 

    hExp->Draw();
    // hDataMC->Draw();

    return 0;
}