#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1D.h"
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
    vector<double> Pavg2 = {106.339, 119.779, 121.882, 123.952, 125.989, 127.996, 129.974, 131.924, 143.105};

/*
******************************
* Energy shift vs Pavg:
******************************
*/ 
    vector<double> energyShift = {5.567, 4.228, 4.075, 4.051, 4.078, 3.907, 3.753, 3.746, 3.778, 3.722, 3.711, 3.729, 3.684, 3.681, 3.559, 3.561, 3.106};
    vector<double> energyShiftErr = {0.031, 0.005, 0.023, 0.024, 0.017,  0.010, 0.009, 0.026, 0.009, 0.009, 0.013, 0.009, 0.013, 0.008, 0.022, 0.009, 0.011};
    vector<double> energyShiftMC_NoFSR = {5.443, 4.144,  3.922, 3.922, 3.922, 3.770, 3.673, 3.673, 3.673, 3.597, 3.597, 3.597, 3.487, 3.487, 3.300, 3.300, 2.860};
    vector<double> energyShiftMC_FSR = {5.462, 4.122, 3.913, 3.913, 3.913, 3.757, 3.670, 3.670, 3.670, 3.605, 3.605, 3.605, 3.484, 3.484, 3.314, 3.314, 2.879};
    vector<double> meanDelta = {5.477, 4.224, 4.067, 3.914, 3.757, 3.717, 3.684, 3.562, 3.104};
    vector<double> meanDeltaErr = {0.03, 0.009, 0.008, 0.009, 0.004, 0.005, 0.007, 0.009, 0.012};

    TGraphErrors grShiftExp(Pavg.size(), Pavg.data(), energyShift.data(), zeroes.data(), energyShiftErr.data());
    TGraphErrors grShiftMC_NoFSR(Pavg.size(), Pavg.data(), energyShiftMC_NoFSR.data(), zeroes.data(), zeroes.data());
    TGraphErrors grShiftMC_FSR(Pavg.size(), Pavg.data(), energyShiftMC_FSR.data(), zeroes.data(), zeroes.data());
    TGraphErrors grMeanDelta(Pavg2.size(), Pavg2.data(), meanDelta.data(), zeroes.data(), meanDeltaErr.data());

    std::string fromStr = 
        "0.3890836 * pow(1+(493.677/x)*(493.677/x), 1.343142/2) + 0.5472272 + 0.3864033 * 0.3864033 *" 
        "(sqrt(x*x + 493.677*493.677) - 506.3635)/((sqrt(x*x + 493.677*493.677) - 506.3635)*(sqrt(x*x + 493.677*493.677) - 506.3635) + 4.410*4.410/16) +"
        "0.04314388 * exp(-(x-110.9493)*(x-110.9493)/2/2.074735/2.074735)";
    
    TF1 CMD2Preprint("cmd2Cal", fromStr.c_str(), 100, 150, "");
    CMD2Preprint.SetLineColor(kRed);

    grShiftExp.SetTitle("Kch Energy shift (Black -- exp, Blue -- MC, Magenta -- E^{Kch}_{avg} - compton_mean, Red -- CMD2 preprint)");
    grShiftExp.GetXaxis()->SetTitle("P_{avg}, #frac{MeV}{c}");
    // grShiftExp.GetXaxis()->SetTitle("E_{mean}, MeV");
    grShiftExp.GetYaxis()->SetTitle("#DeltaE, MeV");

    grShiftMC_FSR.SetLineColor(kBlue);
    grShiftMC_FSR.SetMarkerColor(kBlue);

    grMeanDelta.SetMarkerColor(kMagenta);

    grShiftExp.DrawClone("AP");
    grShiftMC_FSR.DrawClone("P same");
    grMeanDelta.DrawClone("P same");
    CMD2Preprint.DrawClone("same");

/*
******************************
* EnergyShift Data - MC vs Pavg:
******************************
*/ 
    vector<double> dataMCdiff = {0.11, 0.11, 0.17, 0.15, 0.17, 0.15, 0.08, 0.08, 0.11, 0.12, 0.11, 0.12, 0.2, 0.2, 0.25, 0.25, 0.23};
    vector<double> dataMCdiffErr = {0.044, 0.007, 0.035, 0.03, 0.068, 0.014, 0.013, 0.037, 0.013, 0.013, 0.018, 0.013, 0.018, 0.011, 0.031, 0.013, 0.016};

    vector<double> meanDelta_MCdiff = {0.02, 0.1, 0.15, 0.16, 0.09, 0.11, 0.2, 0.25, 0.23};
    vector<double> meanDelta_MCdiffErr = {0.042, 0.013, 0.011, 0.013, 0.006, 0.007, 0.01, 0.013, 0.017};

    TGraphErrors grDataMC_diff(Pavg.size(), Pavg.data(), dataMCdiff.data(), zeroes.data(), dataMCdiffErr.data());
    TGraphErrors grMeanDelta_MC_diff(Pavg2.size(), Pavg2.data(), meanDelta_MCdiff.data(), zeroes.data(), meanDelta_MCdiffErr.data());
    grMeanDelta_MC_diff.SetMarkerColor(kBlue);
    grMeanDelta_MC_diff.SetName("meanDelta_MCratio");

    grDataMC_diff.SetName("DataMCdiff");
    grDataMC_diff.SetTitle("Black - data-MC, Blue - (E^{(K^{#pm})}_{avg} - compton_mean)-MC");
    grDataMC_diff.GetXaxis()->SetTitle("P_{avg}, #frac{MeV}{c}");
    grDataMC_diff.GetYaxis()->SetTitle("data-MC, MeV");
    grDataMC_diff.DrawClone("AP");
    grMeanDelta_MC_diff.DrawClone("P same");

/*
******************************
* (Data - MC) - fit_curve_vals(pol2: 8.19e-5 * mom * mom - 0.01488 * mom + 0.714251) 
******************************
*/ 
    vector<double> dataMC_fitCurve_diff = {0.05, -0.0, 0.05, 0.03, 0.06, 0.02, -0.06, -0.06, -0.03, -0.03, -0.05, -0.03, 0.04, 0.03, 0.07, 0.07, -0.04};
    TGraphErrors grDataMC_fitCurve_diff(Pavg.size(), Pavg.data(), dataMC_fitCurve_diff.data(), zeroes.data(), dataMCdiffErr.data());
    TH1D hDataMC_fitCurve_diff("hist", "hist", 5, -0.1, 0.1);
    for(auto diff : dataMC_fitCurve_diff)
    { hDataMC_fitCurve_diff.Fill(diff); }

    grDataMC_fitCurve_diff.DrawClone("AP");
    hDataMC_fitCurve_diff.DrawClone();

/*
******************************
* EnergyShift vs M correlation
******************************
*/ 
    vector<double> M = {497.525, 497.561, 497.54, 497.54, 497.54, 497.553, 497.548, 497.548, 497.548, 497.573, 497.573, 497.573, 497.586, 497.586, 497.565, 497.565, 497.609};
    vector<double> M_err = {0.042, 0.013, 0.011, 0.011, 0.011, 0.013, 0.006, 0.006, 0.006, 0.007, 0.007, 0.007, 0.01, 0.01, 0.013,  0.013, 0.017};
    // M - M_mean (pol0); M_mean (pol0) = 497.561 +/- 0.004 MeV
    vector<double> M_diff = {-0.037, -0.001, -0.022, -0.022, -0.022, -0.009, -0.014, -0.014, -0.014, 0.011, 0.011, 0.011, 0.024, 0.024, 0.003, 0.003, 0.047};
    vector<double> M_diff_err = {0.036, 0.014, 0.011, 0.011, 0.011, 0.009, 0.008, 0.008, 0.008, 0.01, 0.01, 0.01, 0.013, 0.013, 0.014, 0.014, 0.029};


    TGraphErrors grDataMC_diff_vs_M(M_diff.size(), M_diff.data(), dataMCdiff.data(), M_diff_err.data(), dataMCdiffErr.data());

    grDataMC_diff_vs_M.SetTitle(("M_{mean} = 497.561 #pm 0.004 MeV, CorrelationFactor = " + std::to_string(grDataMC_diff_vs_M.GetCorrelationFactor())).c_str());
    grDataMC_diff_vs_M.GetXaxis()->SetTitle("M - M_{mean}, MeV");
    grDataMC_diff_vs_M.GetYaxis()->SetTitle("#DeltaE Data - MC, MeV");

    // grDataMC_diff_vs_M.DrawClone("AP");

    return 0;
}