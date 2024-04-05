#include "TF1.h"
#include "TH1.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLatex.h"

int RC_EnergyCorrection()
{
    TCanvas *c = new TCanvas("c1","hists with different scales",600,400);
    gStyle->SetOptStat(kFALSE);
    std::vector<Float_t> zeroes(100, 0.0);
/*
******************************
* deltaE^{(RC)} vs E_beam:
******************************
*/
    std::vector<Float_t> vE = {504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 513.864};
    std::vector<Float_t> vE0 = {504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 511.444, 513.864};
    
    std::vector<Float_t> vEcorr_smeared = {0.105504, 0.0760349, 0.069423, 0.0598971, 0.0769491, 0.118673, 0.197142, 0.336072, 0.467132, 1.54063};
    std::vector<Float_t> vEcorrErr_smeared = {0.00258291, 0.00285727, 0.00271313, 0.00137076, 0.00132589, 0.00158087, 0.00393903, 0.0054766, 0.00680275, 0.0145707};

    TGraphErrors grEcorr_smeared(vEcorr_smeared.size(), vE0.data(), vEcorr_smeared.data(), zeroes.data(), vEcorrErr_smeared.data());

    std::vector<Float_t> vEcorr = {0.0984153, 0.0744702, 0.0639822, 0.0584108, 0.0760487, 0.119492, 0.201323, 0.338991, 0.463166, 1.54571};
    std::vector<Float_t> vEcorrErr = {0.00194074, 0.0022625, 0.0020729, 0.00108828, 0.00104003, 0.00122872, 0.00306742, 0.00423206, 0.00521031, 0.0112342};
    TGraphErrors grEcorr(vEcorr.size(), vE0.data(), vEcorr.data(), 
                                        zeroes.data(), vEcorrErr.data());

    grEcorr_smeared.SetTitle("Black -- with Energy smearing, Blue -- without");
    grEcorr_smeared.GetXaxis()->SetTitle("E_{beam}, MeV");
    grEcorr_smeared.GetYaxis()->SetTitle("#DeltaE^{(RC)}, MeV");
    grEcorr_smeared.GetYaxis()->SetTitleOffset(1.2);
    grEcorr_smeared.SetName("Ecorr_smeared");
    grEcorr_smeared.SetMarkerSize(1.5);

    grEcorr.SetName("Ecorr");
    grEcorr.SetMarkerSize(1.5);
    grEcorr.SetMarkerStyle(22);
    grEcorr.SetMarkerColor(kBlue);

    auto h2 = new TH2F("h", "Axes", 100, 504 , 512, 100, 0.03, 0.6);
    h2->SetTitle("Black -- with Energy smearing, Blue -- without");
    h2->GetXaxis()->SetTitle("E_{beam}, MeV");
    h2->GetYaxis()->SetTitle("#DeltaE^{(RC)}, MeV");

    h2->Draw();
    c->Update();

    TGaxis *axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
         gPad->GetUxmax(), gPad->GetUymax(), 0.03, 0.6 * 1.0748, 510, "+L");
    axis->SetTitle("#DeltaM^{(RC)}, MeV ");
    axis->SetTitleFont(42);
    axis->SetTitleSize(0.05);
    axis->SetTitleOffset(0.7);
    axis->Draw();

    grEcorr_smeared.DrawClone("P same");
    grEcorr.DrawClone("P same");
    std::cout << gPad->GetUxmax() << " : " << gPad->GetUymin() << std::endl;

    return 0;
}
