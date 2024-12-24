#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLatex.h"

int sigma_theta_data_mc_diff()
{
    TCanvas c("canv", "canv");
    std::vector<Float_t> zeroes(100, 0.0);

    std::vector<Float_t> vE = {507.862, 508.404, 508.957, 509.528, 509.956};

    std::vector<Float_t> sigma_diff = {30, 15, 24, 3, 13};
    std::vector<Float_t> sigma_diff_err = {13, 8, 8, 7, 8};

    TGraphErrors gr_DataMC_diff(sigma_diff.size(), vE.data(), sigma_diff.data(), zeroes.data(), sigma_diff_err.data());
    gr_DataMC_diff.SetTitle("EXP-MC diff");
    gr_DataMC_diff.GetXaxis()->SetTitle("E_{beam}, MeV");
    gr_DataMC_diff.GetYaxis()->SetTitle("#delta#sigma^{(#theta)}_{M}/#sigma^{(#theta)}_{M}, %");
    gr_DataMC_diff.SetName("DataMC_diff");
    gr_DataMC_diff.SetMarkerSize(1);

    gr_DataMC_diff.DrawClone("AP");
    
    c.DrawClone();
    return 0;
}