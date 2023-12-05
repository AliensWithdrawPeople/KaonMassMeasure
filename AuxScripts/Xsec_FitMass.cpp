#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TRandom.h>

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
int Xsec_FitMass()
{
    std::vector<double> status = {};
    std::vector<double> chi2 = {};
    std::vector<double> ev_sig = {};
    std::vector<double> ev_sig_err = {};
    std::vector<double> ev_bkg = {};
    std::vector<double> ev_bkg_err = {};

    std::vector<std::string> ens = {"501", "503", "505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514", "517", "520", "525", "530"};
    // std::string energy = "503";
    auto c = new TCanvas("canv","canv",800,400);  
    c->Divide(4, 2);
    for(int i = 0; i < 16; i++)
    {
        auto pad = c->cd(i + 1);
        pad->SetLogy();
        std::string energy = ens[i];
        auto exp = TFile::Open(("./tr_ph/PhiXSection/Kn" + energy + ".root").c_str());
        auto mass_exp = new TH1D("mass_exp", "mass_exp", 80, 420, 580);
        auto tr_exp = exp->Get<TTree>("Kn");
        tr_exp->Draw("mass >> mass_exp", "", "goff");

        auto MC = TFile::Open(("./tr_ph/PhiXSection/MC/MC_Kn" + energy + ".root").c_str());
        auto mass_MC = new TH1D("mass_MC", "mass_MC", 160, 420, 580);
        auto mass_MC_conv = new TH1D("mass_MC_conv", "mass_MC_conv", 160, 420, 580);
        auto tr_MC = MC->Get<TTree>("Kn_MC");
        tr_MC->Draw("mass >> mass_MC", "", "goff");
        for(int i = 0; i < mass_MC->GetNbinsX(); i++)
        {
            for(int j = 0; j < mass_MC->GetBinContent(i); j++)
            { mass_MC_conv->Fill(gRandom->Gaus(mass_MC->GetBinCenter(i) - 1., 0.3)); }
        }
        
        
        RooRealVar mass("mass", "mass", 420, 580, "MeV/c^{2}");
        mass.setRange("left_slope", 430, 505);
        mass.setRange("right_slope", 490, 570);
        mass.setRange("short", 430, 570);
        mass.setRange("max", 420, 580);
        RooDataHist exp_hist("exp_hist", "exp_hist", mass, mass_exp);
        RooDataHist mc_conv_hist("mc_conv_hist", "mc_conv_hist", mass, mass_MC_conv);
        const auto hist = &exp_hist;
        double hist_events = hist->sumEntries();
    // **************** Polynomial background: Start **************************
        RooRealVar a0("a0", "a0", 2.7, -1e7, 1e7);
        RooRealVar a1("a1", "a1", 2.9, -1000, 1000);
        RooRealVar a2("a2", "a2", 0, -1000, 1000);
        RooRealVar a3("a3", "a3", 0, -1000, 1000);
        // RooPolynomial bkg("bkg", "Background", mass, RooArgSet());
        RooPolynomial bkg("bkg", "Background", mass, RooArgSet(a0, a1));
    // **************** Polynomial background: End **************************

        // auto frame = new RooPlot("frame", ("Energy point = " + energy).c_str(), mass, 420, 580, 160);
        auto frame = mass.frame();
        
    // **************** Signal + gaussian smearing: Start **************************
        RooDataHist sig_hist("sig_hist", "sig_hist", mass, mass_MC);
        RooHistPdf sig_pdf("sig_pdf", "sig_pdf", mass, sig_hist, 3);

        RooRealVar mean_corr("mean_corr", "mean_corr", -0.7, -5, 5);
        RooRealVar sigma_corr("sigma_corr", "sigma_corr", 1e-1, 1e-6, 5);
        RooGaussian gauss_pdf("gauss", "gauss", mass, mean_corr, sigma_corr);
        // sigma_corr.setConstant(true);
        RooFFTConvPdf sig_conv_pdf("sig_conv_pdf", "sig_conv_pdf", mass, sig_pdf, gauss_pdf);
    // **************** Signal + gaussian smearing: End **************************

        RooRealVar n_sig("n_sig", "n_sig", hist_events, 0, 1.2 * hist_events);
        RooRealVar n_bkg("n_bkg", "n_bkg", 1.e2, 0, 1.2 * hist_events);        
        RooAddPdf model("model", "MC + bkg(pol2)", RooArgList(sig_conv_pdf, bkg), RooArgList(n_sig, n_bkg));
        auto res = model.fitTo(*hist, Save(), PrintLevel(0), Range("max"), MaxCalls(1000));
        // std::cout << "Nll = " << res->minNll() << std::endl;


        model.plotOn(frame, Normalization(n_sig.getVal() + n_bkg.getVal()));
        model.plotOn(frame, Components(bkg), LineStyle(kDashed), LineColor(kRed), Normalization(n_bkg.getVal()));
        model.plotOn(frame, Components(sig_conv_pdf), LineStyle(kDashDotted), Normalization(n_sig.getVal()));
        hist->plotOn(frame);

        frame->SetTitle(("Energy = " + energy + " MeV").c_str());
        frame->GetYaxis()->SetTitle("Events / 2 MeV/c^{2}");
        frame->GetXaxis()->SetTitle("K_{S} Mass, MeV/c^{2}");
        frame->GetYaxis()->SetRangeUser(0.1, mass_exp->GetEntries() / 2.);
        frame->DrawClone();
        
        std::cout << "integral = " << n_sig.getVal() + n_bkg.getVal() << std::endl;
        std::cout << "hist_events = " << hist_events << std::endl;
        std::cout << "chi2/ndf = " << frame->chiSquare(6) << std::endl;
        status.push_back(res->status());
        chi2.push_back(frame->chiSquare(6));
        ev_sig.push_back(n_sig.getVal());
        ev_sig_err.push_back(n_sig.getError());
        ev_bkg.push_back(n_bkg.getVal());
        ev_bkg_err.push_back(n_bkg.getError());
        exp->Close();  
    }

    TFile output("res_1.root", "recreate");
    c->Write();
    output.Save();
    output.Close();

    auto print = [](const std::vector<double>& vec, std::string name) {
        std::cout << name << " = [";
        for(const auto& el:vec)
        { std::cout << el << ", "; } 
        std::cout << "]" << std::endl;
    };

    print(status, "status");
    print(chi2, "chi2_ndf");
    print(ev_sig, "ev_sig");
    print(ev_sig_err, "ev_sig_err");
    print(ev_bkg, "ev_bkg");
    print(ev_bkg_err, "ev_bkg_err");
    return 0;
}
