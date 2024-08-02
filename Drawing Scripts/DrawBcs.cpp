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
    // std::string filename = "bcs_kskl_track_rec_eff_sig_extr_syst.root";
    // std::string filename = "bcs_kskl_sig_extr_syst.root";
    std::string filename = "bcs_kskl_track_rec_eff_sig_extr_syst.root";
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
    // grResidual.DrawClone("ap");

    TCanvas canv("canv");
    bcs->GetYaxis()->SetTitle("#sigma^{(data)}_{born}, nb");
    bcs->GetXaxis()->SetTitle("E_{cms}, GeV");
    bcs->SetTitle("Markers --- data, Line --- Fit");
    bcs->DrawClone("AP");
    f_bcs->DrawClone("same");
    canv.SetLogy();
    canv.DrawClone();

    // ToyMC: start
    // std::vector<double> xsex_rand = {};
    // gRandom->SetSeed(12);
    // for(int j = 0; j < 40; j++) {
    //     for(int i = 0; i < N_points; i++)
    //     {
    //         x = bcs->GetPointX(i);
    //         y = bcs->GetPointY(i);
    //         yErr = bcs->GetErrorY(i);
    //         f_y = f_bcs->Eval(x);

    //         xs.push_back(x);
    //         ratios.push_back((y - f_y) / f_y);
    //         residuals.push_back(y - f_y);

    //         ratioErrs.push_back(yErr / f_y);
    //         residualErrs.push_back(yErr);
    //         if(i > 5 && i < N_points - 7)
    //         { xsex_rand.push_back(gRandom->Gaus(f_vcs->Eval(x), 9.26736)); }
    //         else
    //         { xsex_rand.push_back(vcs->GetPointY(i)); }
    //     }
    //     std::cout << "{";
    //     for(const auto& el : xsex_rand)
    //     {
    //         std::cout << el << ", ";
    //     }
    //     xsex_rand.clear();
    //     std::cout << "}," << std::endl;
    // }
    // // ToyMC: end

    // std::vector<double> ind = {0.4, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0};
    // std::vector<double> m = {1019.442, 1019.451, 1019.453, 1019.451, 1019.451, 1019.456, 1019.455, 1019.459, 1019.443, 1019.446, 1019.450, 1019.45, 1019.459, 1019.448, 1019.449, 1019.454, 1019.448, 1019.440, 1019.448, 1019.445, 1019.455, 1019.460, 1019.455, 1019.441, 1019.447, 1019.443, 1019.461, 1019.442, 1019.445, 1019.447, 1019.443, 1019.446, 1019.452, 1019.462, 1019.455, 1019.451, 1019.446, 1019.444, 1019.442};
    // std::vector<double> m_err = {0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005};
    // TGraphErrors toy_mc(ind.size(), ind.data(), m.data(), zeroes.data(), m_err.data());
    // toy_mc.GetXaxis()->SetTitle("sample");
    // toy_mc.GetYaxis()->SetTitle("M^{(fit)}_{#phi}, MeV");
    // toy_mc.DrawClone();
    // std::cout << grResidual.GetRMS(2) << std::endl;
    // std::cout << grResidual.GetMean(2) << std::endl;
    return 0;
}