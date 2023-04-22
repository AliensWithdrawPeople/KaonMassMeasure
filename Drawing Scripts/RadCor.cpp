#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TAxis.h"

int RadCor()
{
    std::vector<Float_t> zeroes(100, 0.0);
    std::vector<Float_t> RC = {0.073, 0.075, 0.073, 0.069, 0.064, 0.112, 0.186, 0.328, 1.51};
    std::vector<Float_t> RC_smeared = {0.1, 0.099, 0.089, 0.08, 0.086, 0.116, 0.191, 0.334, 1.453};
    std::vector<Float_t> vE = {505, 508, 508.5, 509, 509.5, 510, 510.5, 511, 514};

    TGraph grRC(RC.size(), vE.data(), RC.data());
    TGraph grRC_smearing(RC.size(), vE.data(), RC_smeared.data());

    grRC_smearing.SetMarkerColor(kBlue);
    grRC_smearing.SetMarkerSize(1.5);
    grRC_smearing.SetMarkerStyle(22);

    grRC.SetMarkerSize(1.5);
    grRC.SetMarkerStyle(23);

    grRC.SetTitle("ISR correction (black -- without energy smearing, blue  -- with smearing)");
    grRC.GetXaxis()->SetTitle("E, MeV");
    grRC.GetYaxis()->SetTitle("#DeltaM^{(RC)}, #frac{MeV}{c^{2}}");

    grRC.DrawClone("AP");
    grRC_smearing.DrawClone("P same");
    return 0;
}