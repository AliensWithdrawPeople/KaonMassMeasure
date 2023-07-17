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

int PiTheta()
{
    std::string path = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/MC/KsKl_Smeared/";
    const auto filename = [&path](std::string point) {return path + "MC" + point + ".root"; };
    
    auto file = TFile::Open(filename("510").c_str());
    auto ksTr = file->Get<TTree>("ksTree");

    auto hPiTheta1 = std::make_unique<TH2D>(TH2D("hPiTheta1", "#theta_{#pi^{+}} vs #theta_{#pi^{-}} @ #theta_{K_{S}} #approx #pi/2", 200, 0.5, 2.5, 200, 0.5, 2.5));
    auto hPiTheta2 = std::make_unique<TH2D>(TH2D("hPiTheta2", "#theta_{#pi^{+}} vs #theta_{#pi^{-}} @ #theta_{K_{S}} #approx #pi/4", 200, 0.5, 2.5, 200, 0.5, 2.5));
    auto hPiTheta3 = std::make_unique<TH2D>(TH2D("hPiTheta3", "#theta_{#pi^{+}} vs #theta_{#pi^{-}} @ #theta_{K_{S}} #approx 3#pi/4", 200, 0.5, 2.5, 200, 0.5, 2.5));

    ksTr->Draw("piThetaPos : piThetaNeg >> hPiTheta1", "fabs(kstheta - TMath::Pi()/2) < 0.1", "goff");
    ksTr->Draw("piThetaPos : piThetaNeg >> hPiTheta2", "fabs(kstheta - TMath::Pi()/4) < 0.1", "goff");
    ksTr->Draw("piThetaPos : piThetaNeg >> hPiTheta3", "fabs(kstheta - 3 * TMath::Pi()/4) < 0.1", "goff");

    hPiTheta1->DrawClone("col");
    hPiTheta2->DrawClone("col");
    hPiTheta3->DrawClone("col");
    return 0;
}