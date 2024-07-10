#include <TFile.h>
#include <TGraphErrors.h>
#include <TRatioPlot.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TRoot.h>
#include <TRandom.h>

#include <vector>


int DrawBcs()
{
    // std::string filename = "bcs_bkg_pol0.root";
    std::string filename = "bcs_kskl_track_rec_eff.root";
    std::string fullFilename = "C:/work/Science/BINP/Kaon Mass Measure/PhiMesonFit/bcs/" + filename; 
    auto file = TFile::Open(fullFilename.c_str());
    auto vcs = file->Get<TGraphErrors>("vcs");
    auto f_vcs = file->Get<TF1>("f_vcs");

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

    grRatio.GetYaxis()->SetTitle("(data - fit) / fit");
    grRatio.GetXaxis()->SetTitle("E_{cms}, GeV");

    grRatio.DrawClone("ap");
    grResidual.DrawClone("ap");

    TCanvas canv("canv");
    bcs->GetYaxis()->SetTitle("#sigma^{(data)}_{born}, nb");
    bcs->GetXaxis()->SetTitle("E_{cms}, GeV");
    bcs->SetTitle("Markers --- data, Line --- Fit");
    bcs->DrawClone("AP");
    f_bcs->DrawClone("same");
    canv.SetLogy();
    canv.DrawClone();
    // ToyMC: start
    std::vector<double> xsex_rand = {};
    for(int j = 0; j < 10; j++) {
        //  gRandom->SetSeed(6 + );
        // for(int i = 0; i < N_points; i++)
        // {
        //     x = bcs->GetPointX(i);
        //     y = bcs->GetPointY(i);
        //     yErr = bcs->GetErrorY(i);
        //     f_y = f_bcs->Eval(x);

        //     xs.push_back(x);
        //     ratios.push_back((y - f_y) / f_y);
        //     residuals.push_back(y - f_y);

        //     ratioErrs.push_back(yErr / f_y);
        //     residualErrs.push_back(yErr);

        //     xsex_rand.push_back(gRandom->Gaus(f_vcs->Eval(x), vcs->GetErrorY(i)));
        // }
        // std::cout << "{";
        // for(const auto& el : xsex_rand)
        // {
        //     std::cout << el << ", ";
        // }
        // xsex_rand.clear();
        // std::cout << "}," << std::endl;
    }
    // ToyMC: end

    return 0;
}