#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TF1.h" 

using std::vector;

int KchEnergyShift()
{
    vector<double> zeroes(100, 0.0);

    vector<double> E = {505, 508, 508.5, 508.5, 508.5, 509, 509.5, 509.5, 509.5, 510, 510, 510, 510.5, 510.5, 511, 511, 514};
    vector<double> Pavg = {106.339, 119.779, 121.882, 121.882, 121.882, 123.952, 125.989, 125.989, 125.989, 127.996, 127.996, 127.996, 129.974, 129.974, 131.924, 131.924, 143.10};
    vector<double> energyShift = {5.567, 4.228, 4.075, 4.051, 4.078, 3.907, 3.753, 3.746, 3.778, 3.722, 3.711, 3.729, 3.684, 3.681, 3.559, 3.561, 3.106};
    vector<double> energyShiftErr = {0.031, 0.005, 0.023, 0.024, 0.017,  0.010, 0.009, 0.026, 0.009, 0.009, 0.013, 0.009, 0.013, 0.008, 0.022, 0.009, 0.011};
    vector<double> energyShiftMC_NoFSR = {5.443, 4.144,  3.922, 3.922, 3.922, 3.770, 3.673, 3.673, 3.673, 3.597, 3.597, 3.597, 3.487, 3.487, 3.300, 3.300, 2.860};
    vector<double> energyShiftMC_FSR = {5.462, 4.122, 3.913, 3.913, 3.913, 3.757, 3.670, 3.670, 3.670, 3.605, 3.605, 3.605, 3.484, 3.484, 3.314, 3.314, 2.879};

    vector<double> dataMCratio = {1.02, 1.03, 1.04, 1.04, 1.04, 1.04, 1.02, 1.02, 1.03, 1.03, 1.03, 1.03, 1.06, 1.06, 1.07, 1.07, 1.08};
    vector<double> dataMCratioErr = {0.008, 0.002, 0.009, 0.008, 0.017, 0.004, 0.003, 0.01, 0.003, 0.004, 0.005, 0.004, 0.005, 0.003, 0.009, 0.004, 0.005};

    TGraphErrors grShiftExp(Pavg.size(), Pavg.data(), energyShift.data(), zeroes.data(), energyShiftErr.data());
    TGraphErrors grShiftMC_NoFSR(Pavg.size(), Pavg.data(), energyShiftMC_NoFSR.data(), zeroes.data(), zeroes.data());
    TGraphErrors grShiftMC_FSR(Pavg.size(), Pavg.data(), energyShiftMC_FSR.data(), zeroes.data(), zeroes.data());

    TGraphErrors grDataMC_ratio(Pavg.size(), Pavg.data(), dataMCratio.data(), zeroes.data(), dataMCratioErr.data());
    grDataMC_ratio.GetXaxis()->SetTitle("P_{avg}, #frac{MeV}{c}");
    grDataMC_ratio.GetYaxis()->SetTitle("data/MC");

    std::string fromStr = 
        "0.3890836 * pow(1+(493.677/x)*(493.677/x), 1.343142/2) + 0.5472272 + 0.3864033 * 0.3864033 *" 
        "(sqrt(x*x + 493.677*493.677) - 506.3635)/((sqrt(x*x + 493.677*493.677) - 506.3635)*(sqrt(x*x + 493.677*493.677) - 506.3635) + 4.410*4.410/16) +"
        "0.04314388 * exp(-(x-110.9493)*(x-110.9493)/2/2.074735/2.074735)";
    
    TF1 func("cmd2Cal", fromStr.c_str(), 100, 150, "");
    func.SetLineColor(kRed);

    grShiftExp.SetTitle("Kch Energy shift (Black -- exp, Blue -- MC, Red -- CMD2 preprint)");
    grShiftExp.GetXaxis()->SetTitle("P_{avg}, #frac{MeV}{c}");
    // grShiftExp.GetXaxis()->SetTitle("E_{mean}, MeV");
    grShiftExp.GetYaxis()->SetTitle("#DeltaE, MeV");
    grShiftMC_FSR.SetLineColor(kBlue);
    grShiftMC_FSR.SetMarkerColor(kBlue);

    grShiftExp.DrawClone("AP");
    grShiftMC_FSR.DrawClone("P same");
    func.DrawClone("same");

    grDataMC_ratio.DrawClone("AP");

    return 0;
}