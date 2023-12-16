#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooChi2Var.h"

#include <vector>
#include <TCanvas.h>

using namespace RooFit;

// Run exclusively from the machine with the ROOT properly tuned to work with FFTW (in my case ubuntu via WSL).
int Xsec_FitVertexRho()
{
    std::vector<double> ev_sig = {};
    std::vector<double> ev_sig_err = {};
    std::vector<double> ev_bkg = {};
    std::vector<double> ev_bkg_err = {};

    std::vector<std::string> ens = {"501", "503", "505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514", "517", "520", "525", "530"};
    std::string energy = "509.5";
    auto c = new TCanvas("canv","canv",800,400);  
    // c->Divide(4, 4);
    // for(int i = 0; i < ens.size(); i++)
    // {
        auto pad = c->cd(7 + 1 - 8);
        pad->SetLogy();
        // std::string energy = ens[i];
        auto exp = TFile::Open(("./tr_ph/PhiXSection/Kn" + energy + ".root").c_str());
        auto kstlen_exp = new TH1D("kstlen_exp", "kstlen_exp", 100, 0, 6);
        auto tr_exp = exp->Get<TTree>("Kn");
        tr_exp->Draw("kstlen >> kstlen_exp", "", "goff");

        auto MC = TFile::Open(("./tr_ph/PhiXSection/MC/MC_Kn" + energy + ".root").c_str());
        auto kstlen_MC = new TH1D("kstlen_MC", "kstlen_MC", 100, 0, 6);
        auto tr_MC = MC->Get<TTree>("Kn_MC");
        tr_MC->Draw("kstlen >> kstlen_MC", "", "goff");
        
        
        RooRealVar kstlen("kstlen", "kstlen", 0, 6, "MeV/c^{2}");
        kstlen.setRange("left_slope", 430, 505);
        kstlen.setRange("right_slope", 490, 570);
        kstlen.setRange("max", 0, 6);
        RooDataHist exp_hist("kstlen_exp", "kstlen_exp", kstlen, kstlen_exp);
        RooDataHist mc_hist("kstlen_MC", "kstlen_MC", kstlen, kstlen_MC);

    // **************** Polynomial background: Start **************************
        RooRealVar a0("a0", "a0", 0.0, -1e7, 1e7);
        RooRealVar a1("a1", "a1", 2, -1000, 1000);
        RooRealVar a2("a2", "a2", 0, -1000, 1000);
        RooRealVar a3("a3", "a3", 0, -1000, 1000);
        RooPolynomial bkg("bkg", "Background", kstlen, RooArgSet());
        // RooPolynomial bkg("bkg", "Background", mass, RooArgSet(a0, a1, a2, a3));
    // **************** Polynomial background: End **************************

        // auto frame = new RooPlot("frame", ("Energy point = " + energy).c_str(), mass, 420, 580, 160);
        auto frame = kstlen.frame();
        
    // **************** Signal + gaussian smearing: Start **************************
        RooDataHist sig_hist("sig_hist", "sig_hist", kstlen, kstlen_MC);
        RooHistPdf sig_pdf("sig_pdf", "sig_pdf", kstlen, sig_hist, 3);

        RooRealVar mean_corr("mean_corr", "mean_corr", -0.7, -5, 5);
        RooRealVar sigma_corr("sigma_corr", "sigma_corr", 1e-6, 1e-6, 5);
        RooGaussian gauss_pdf("gauss", "gauss", kstlen, mean_corr, sigma_corr);
        // sigma_corr.setConstant(true);
        RooFFTConvPdf sig_conv_pdf("sig_conv_pdf", "sig_conv_pdf", kstlen, sig_pdf, gauss_pdf);
    // **************** Signal + gaussian smearing: End **************************

        RooRealVar n_sig("n_sig", "n_sig", kstlen_exp->GetEntries(), 0, 1.2 * kstlen_exp->GetEntries());
        RooRealVar n_bkg("n_bkg", "n_bkg", 0, 0, 1.2 * kstlen_exp->GetEntries());
        n_bkg.setConstant(true);        
        RooAddPdf model("model", "MC + bkg(pol2)", RooArgList(sig_pdf, bkg), RooArgList(n_sig, n_bkg));
        auto res = model.fitTo(mc_hist, Save(), PrintLevel(0), Range("max"));
        // std::cout << "Nll = " << res->minNll() << std::endl;


        model.plotOn(frame, Normalization(n_sig.getVal() + n_bkg.getVal()));
        model.plotOn(frame, Components(bkg), LineStyle(kDashed), LineColor(kRed), Normalization(n_bkg.getVal()));
        model.plotOn(frame, Components(sig_conv_pdf), LineStyle(kDashDotted), Normalization(n_sig.getVal()));
        mc_hist.plotOn(frame);
        frame->GetXaxis()->SetTitle("K_{S} Mass, MeV/c^{2}");
        frame->GetYaxis()->SetRangeUser(0.1, kstlen_exp->GetEntries() / 2.);
        frame->DrawClone();
        
        std::cout << "integral = " << n_sig.getVal() + n_bkg.getVal() << std::endl;
        std::cout << "mass_exp events = " << kstlen_exp->GetEntries() << std::endl;
        std::cout << "chi2/ndf = " << frame->chiSquare(6) << std::endl;
        ev_sig.push_back(n_sig.getVal());
        ev_sig_err.push_back(n_sig.getError());
        ev_bkg.push_back(n_bkg.getVal());
        ev_bkg_err.push_back(n_bkg.getError());
        exp->Close();  
    // }

    // TFile output("res_1.root", "recreate");
    // c->Write();
    // output.Save();
    // output.Close();

    auto print = [](const std::vector<double>& vec, std::string name) {
        std::cout << name << " = [";
        for(const auto& el:vec)
        { std::cout << el << ", "; } 
        std::cout << "]" << std::endl;
    };

    print(ev_sig, "ev_sig");
    print(ev_sig_err, "ev_sig_err");
    print(ev_bkg, "ev_bkg");
    print(ev_bkg_err, "ev_bkg_err");
    return 0;
}
