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

    return 0;
}