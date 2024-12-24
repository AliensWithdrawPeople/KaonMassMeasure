#include <vector>
#include <iostream>
#include <utility>
#include <string>
#include <iomanip>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TRandom.h>
#include <TCanvas.h>

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooPolynomial.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooChi2Var.h"
#include "RooTFnBinding.h"
#include "RooCFunction1Binding.h"
#include "RooCFunction3Binding.h"
#include "RooWrapperPdf.h"

using namespace RooFit;

// Run exclusively from the machine with the ROOT properly tuned to work with FFTW (in my case ubuntu via WSL).
int eval_sigma()
{
    const std::vector<std::string> energyPoints = {"505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514"};
    // std::string energyPoint = "510"; 
    int rebin_factor = 80;
    std::vector<double> sigmas = {}; 
    std::vector<double> sigma_errs = {}; 
    // for(const auto& energyPoint : energyPoints)
    {
        std::string energyPoint = "508";
        // auto filename_exp = "/mnt/c/work/Science/BINP/Kaon Mass Measure/hists/Exp/tests/Hists_Exp" + energyPoint + "_std.root";
        // auto filename_exp = "/mnt/c/work/Science/BINP/Kaon Mass Measure/hists/MC/tests/Hists_MC" + energyPoint + "_std_smear_1.16sigma.root";
        auto filename_exp = "/mnt/c/work/Science/BINP/Kaon Mass Measure/hists/MC/tests/Hists_MC" + energyPoint + "_std.root";
        // auto filename_MC = "/mnt/c/work/Science/BINP/Kaon Mass Measure/hists/MC/tests/Hists_MC" + energyPoint + "_std_no_E_smear.root";
        auto filename_MC = "/mnt/c/work/Science/BINP/Kaon Mass Measure/hists/MC/tests/Hists_MC" + energyPoint + "_std_no_E_smear.root";
        
        auto exp = TFile::Open(filename_exp.c_str());
        auto exp_distr = exp->Get<TH1D>("hMass");
        exp_distr->SetName("exp_distr");
        exp_distr->Rebin(rebin_factor);

        auto MC = TFile::Open(filename_MC.c_str());
        auto MC_distr = MC->Get<TH1D>("hMass");
        MC_distr->SetName("MC_distr");
        MC_distr->Rebin(rebin_factor);

        // auto exp_2d_distr = MC->Get<TH2D>("hPiPlusThetaVsKsTheta");
        // exp_2d_distr->RebinY(8);
        // auto exp_distr = exp_2d_distr->ProjectionY("exp_distr");

        RooRealVar Mks("Mks", "Mks", 490, 505, "M_{K_{S}}, MeV");

        RooDataHist exp_hist("exp_hist", "exp_hist", Mks, exp_distr);
        const auto hist = &exp_hist;
        double hist_events = hist->sumEntries();

        auto frame = Mks.frame();
        Mks.setRange("full", 490, 505);
        Mks.setRange("short", 496, 499);
        Mks.setRange("super_short", 496.5, 498.2);
        // MC truth distribution: start
        RooDataHist sig_hist("sig_hist", "sig_hist", Mks, MC_distr);
        RooHistPdf sig_pdf("sig_pdf", "sig_pdf", Mks, sig_hist, 3); 
        // MC truth distribution: end

        // Resolution (gaussian distribution): start
        // RooRealVar mean_corr("mean_corr", "mean_corr", -0.047, -1, 1);
        RooRealVar mean_corr("mean_corr", "mean_corr", -0.0, -1, 1);
        RooRealVar sigma_corr("sigma_corr", "sigma_corr", 1e-1, 1e-6, 5);
        RooGaussian gauss_pdf("gauss", "gauss", Mks, mean_corr, sigma_corr);
        // mean_corr.setConstant(true);
        // sigma_corr.setConstant(true);

        // Resolution (gaussian distribution): end

        RooFFTConvPdf sig_conv_pdf("sig_conv_pdf", "sig_conv_pdf", Mks, sig_pdf, gauss_pdf);
    // **************** Signal + gaussian smearing: End **************************

        RooRealVar n_sig("n_sig", "n_sig", 0.9 * hist_events, 0, 1.2 * hist_events);

        RooAddPdf model("model", "MC convolved with gaussian distr", RooArgSet(sig_conv_pdf), RooArgSet(n_sig));
        auto res = model.fitTo(*hist, Save(), Range("short"), PrintLevel(0));
        std::cout << "Nll = " << res->minNll() << std::endl;

        hist->plotOn(frame);
        model.plotOn(frame);
        model.plotOn(frame, Components(sig_conv_pdf), LineStyle(kDashDotted));
        hist->plotOn(frame);
        frame->Draw();

        std::cout << "Resolution: Mean = " << mean_corr.getVal() << std::endl;
        std::cout << "Resolution: Mean error = " << mean_corr.getError() << std::endl;
        std::cout << "Resolution: Sigma = " << sigma_corr.getVal() << std::endl;
        std::cout << "Resolution: Sigma error = " << sigma_corr.getError() << std::endl;
        std::cout << "hist_events = " << hist_events << std::endl;
        std::cout << "chi2/ndf = " << frame->chiSquare(2) << std::endl;

        sigmas.push_back(sigma_corr.getVal());
        sigma_errs.push_back(sigma_corr.getError());
    }

    for(const auto& sigma : sigmas)
    {
        std::cout << sigma << ", ";
    }

    std::cout << std::endl;
    for(const auto& sigma : sigma_errs)
    {
        std::cout << sigma << ", ";
    }
    std::cout << std::endl;
    return 0;
}
