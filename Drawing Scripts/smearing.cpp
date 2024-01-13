#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TAxis.h"

int smearing()
{
    std::vector<double> zeroes(100, 0.0);

    std::vector<double> energy = {510.316, 510.274, 510.266, 510.272, 510.270, 510.270};
    std::vector<double> energyErr = {0.004, 0.002, 0.005, 0.004, 0.004, 0.004};

    std::vector<double> energyUnCutted = {510.317,  510.276, 510.270, 510.270, 510.269,  510.269};
    std::vector<double> energyUnCuttedErr = {0.003, 0.001, 0.003, 0.002, 0.002, 0.002};

    std::vector<double> points = {1, 5, 7, 9, 11, 13};

    TGraphErrors grSmearCut(points.size(), points.data(), energy.data(), zeroes.data(), energyErr.data());
    TGraphErrors grSmearNoCut(points.size(), points.data(), energyUnCutted.data(), zeroes.data(), energyUnCuttedErr.data());

    grSmearCut.SetTitle("E_{mean} = 510.458 MeV, #sigma_{beam} = 0.288 MeV (Black -- Without Cuts, Red -- With Cuts)");
    grSmearCut.GetXaxis()->SetTitle("number of points");
    grSmearCut.GetYaxis()->SetTitle("E_{avg}, MeV");
    grSmearCut.SetName("WithCuts");
    grSmearCut.SetMarkerColor(kRed);
    grSmearCut.SetMarkerSize(1.5);
    grSmearCut.DrawClone("AP");
    
    grSmearNoCut.SetName("WithoutCuts");
    grSmearNoCut.SetMarkerStyle(22);
    grSmearNoCut.SetMarkerSize(1.5);
    grSmearNoCut.DrawClone("P Same");

/*
******************************
* M(with smearing) - M(without smearing) vs E_beam in MC:
******************************
*/ 
    
    std::vector<double> vE = {504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 511.444, 513.864};
    std::vector<double> vDeltaM_smearing = {-0.007, -0.029, -0.018, -0.018, -0.001, 0.001, 0.001, 0.033, -0.011, 0.007};
    std::vector<double> vDeltaM_smearing_err = {0.008, 0.011, 0.012, 0.006, 0.006, 0.006, 0.013, 0.015, 0.017, 0.03};
    TGraphErrors grDeltaM_Smearing(vE.size(), vE.data(), vDeltaM_smearing.data(), zeroes.data(), vDeltaM_smearing_err.data());

    grDeltaM_Smearing.SetTitle("M^{(MC)}_{with energy smearing} - M^{(MC)}_{without energy smearing} vs E_{beam}");
    grDeltaM_Smearing.GetYaxis()->SetTitle("#DeltaM^{(MC)}_{energy smearing}, MeV/c^{2}");
    grDeltaM_Smearing.GetXaxis()->SetTitle("E_{avg}, MeV");
    grDeltaM_Smearing.SetMarkerColor(kRed);
    grDeltaM_Smearing.SetMarkerSize(1.5);
    grDeltaM_Smearing.DrawClone("AP");

    return 0;
}