#include <TFile.h>
#include <TGraphErrors.h>
#include <TRatioPlot.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TRoot.h>

#include <vector>


int DrawBcs()
{
    std::string filename = "bcs_RunLumi_v9.root";
    std::string fullFilename = "C:/work/Science/BINP/Kaon Mass Measure/PhiMesonFit/" + filename; 
    auto file = TFile::Open(fullFilename.c_str());
    auto bcs = file->Get<TGraphErrors>("bcs");
    auto f_bcs = file->Get<TF1>("f_bcs");

    int N_points = bcs->GetN();
    double x = 0;
    double y = 0;
    double f_y = 0;
    double yErr = 0;

    std::vector<double> xs = {};

    std::vector<double> ratios = {};
    std::vector<double> ratioErrs = {};

    std::vector<double> residuals = {};
    std::vector<double> residualErrs = {};
    for(int i = 0; i < N_points; i++)
    {
        x = bcs->GetPointX(i);
        y = bcs->GetPointY(i);
        yErr = bcs->GetErrorY(i);
        f_y = f_bcs->Eval(x);

        xs.push_back(x);
        ratios.push_back((y - f_y) / f_y);
        residuals.push_back(y - f_y);

        ratioErrs.push_back(yErr / f_y);
        residualErrs.push_back(yErr);
    }

    std::vector<double> zeroes(100, 0.0);
    TGraphErrors grRatio(xs.size(), xs.data(), ratios.data(), zeroes.data(), ratioErrs.data());
    TGraphErrors grResidual(xs.size(), xs.data(), residuals.data(), zeroes.data(), residualErrs.data());

    grRatio.SetName("Ratio");
    grRatio.SetTitle("Ratio");
    grResidual.SetName("Residual");
    grResidual.SetTitle("Residual");

    grResidual.GetYaxis()->SetTitle("#sigma^{(data)}_{born} - #sigma^{(fit)}_{born}, nb");
    grResidual.GetXaxis()->SetTitle("E_{cms}, GeV");

    grRatio.GetYaxis()->SetTitle("#sigma^{(data)}_{born} / #sigma^{(fit)}_{born}");
    grRatio.GetXaxis()->SetTitle("E_{cms}, GeV");

    grRatio.DrawClone("ap");
    grResidual.DrawClone("ap");

    return 0;
}