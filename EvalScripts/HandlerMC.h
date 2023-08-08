#ifndef HandlerMC_
#define HandlerMC_h

#include "Evaluator.h"
#include "Tree.h"
#include "HistContainer.h"
#include "Spline.h"
#include "CalcSigmas.h"

#include <memory>

#include "TMath.h"


class HandlerMC: Evaluator
{
private:
    std::string energyPoint;
    std::optional<double> energy;
    bool verbose;
    std::unique_ptr<Tree> tree;
    HistContainer container;

    // Sigma_psi(lnY, kstheta) from gaus fit.
    std::vector<std::vector<double>> vSigmaMatrixFit;
    
    /**
    * Correction to psi_reco (dpsi = psi_reco - psi_gen).
    * psi_corrected = psi_reco - deltaPsi_RecGenDiff->Eval(kstheta_reco).
    */
    Spline deltaPsi_RecGenDiff;

    // Resolution of theta as a function of ksTheta.
    Spline thetaSigma;

    /**
    * Creates deltaPsi_RecGenDiff -- correction to psi_reco (dpsi = psi_reco - psi_gen).
    * @param tree reference to Tree object storing data about events,
    * @param spline_filename name of .root file where spline will be saved.
    */
    static Spline CreateDeltaPsiSpline(std::unique_ptr<Tree> const &tree, std::string spline_filename = "", bool save = false);

public:
    HandlerMC(std::string fKsKl, std::string energyPoint, std::optional<double> meanEnergy, 
                bool saveSplines = false, bool isVerbose = true);
    HandlerMC(std::string fKsKl, std::string energyPoint, bool saveSplines = false, bool isVerbose = true): 
            HandlerMC(fKsKl, energyPoint, std::nullopt, saveSplines, isVerbose) {}

};

HandlerMC::HandlerMC(std::string fKsKl, std::string energyPoint, std::optional<double> meanEnergy, 
                    bool saveSplines, bool isVerbose = true): energyPoint{energyPoint}, energy{meanEnergy}, verbose{isVerbose} 
{
    tree = std::unique_ptr<Tree>(new Tree(fKsKl));

    container.Add("hMlnY", new TH2D("hMlnY", "M(lnY)", 300, -0.3, 0.3, 600, 480, 520));
    container.Add("hDeltaM", new TProfile("hDeltaM", "DeltaM(lnY)", 40, -1, 1, -1, 1));
    container.Add("hMlnYpfx", new TProfile("hMlnYpfx","Profile of M versus lnY", 30, -1, 1, 490, 505));
    container.Add("MPsi", new TH2D("MPsi", "M(Psi)", 200, 2, TMath::Pi(), 200, 480, 520));
    container.Add("hM_CrAnglelnY", new TH2D("hM_CrAnglelnY", "M_CrAngle(lnY)", 30, -0.4, 0.4, 40000, 490, 515));
    container.Add("hPsilnY", new TH2D("hPsilnY", "Psi(lnY)", 200, -0.4, 0.4, 10000, 2.4, 3.3));
    container.Add("hPsi1", new TH1D("hPsi1", "Psi distr", 20000, -1, 6.3));
    container.Add("hEnergySpectrum", new TH1D("hEnergySpectrum", "hEnergySpectrum", 6000, 480, 540));
    // After profile cut
    container.Add("hEnergySpectrumCut", new TH1D("hEnergySpectrumCut", "hEnergySpectrumCut", 6000, 480, 540));

    container.Add("hMassVsKsTheta", new TH2D("hMassVsKsTheta", "M vs KsTheta", 600, -1.57, 1.57, 40000, 480, 520));
    container.Add("hKsThetaVsLnY", new TH2D("hKsThetaVsLnY", "KsTheta vs lnY", 300, -0.6, 0.6, 315, 0, 3.15));

    container.Add("hThetaDiffRecGen", new TH2D("hThetaDiffRecGen", "hThetaDiffRecGen", 1200, -6., 6., 20000, -1, 1));
    container.Add("hDiffRecGen_cowboy", new TH2D("hDiffRecGen_cowboy", "hDiffRecGen_cowboy", 1200, -6., 6., 20000, -1, 1));
    container.Add("hDiffRecGen_sailor", new TH2D("hDiffRecGen_sailor", "hDiffRecGen_sailor", 1200, -6., 6., 20000, -1, 1));

    container.Add("hMlnYpfx_cowboy", new TProfile("hMlnYpfx_cowboy", "Profile of M versus lnY, cowboy", 30, -1, 1, 490, 505));
    container.Add("hMlnYpfx_sailor", new TProfile("hMlnYpfx_sailor", "Profile of M versus lnY, sailor", 30, -1, 1, 490, 505));


    container.Get("hMlnY").value()->GetYaxis()->SetTitle("M_{K^{0}_{S}}, #frac{MeV}{c^{2}}");
    container.Get("hMlnY").value()->GetXaxis()->SetTitle("ln(Y)");
    container.Get("hMlnYpfx").value()->GetYaxis()->SetTitle("M_{K^{0}_{S}}, #frac{MeV}{c^{2}}");
    container.Get("hMlnYpfx").value()->GetXaxis()->SetTitle("ln(Y)");

    container.Get("hM_CrAnglelnY").value()->GetXaxis()->SetTitle("lnY");
    container.Get("hM_CrAnglelnY").value()->GetYaxis()->SetTitle("M_{K^{0}_{S}}, #frac{MeV}{c^{2}}");

    deltaPsi_RecGenDiff = CreateDeltaPsiSpline(tree, 
                        "C:/work/Science/BINP/Kaon Mass Measure/ksdpsi_splines/spline_" + energyPoint + ".root", saveSplines);

    //TODO: Create and fill vector of good entries.
    vSigmaMatrixFit = Sigmas::GetSigmaMatrix(tree, goodEntries);
    auto [ksThetaBinCenters, deltaTheta, deltaThetaError] = Sigmas::GetThetaSigmas(tree, goodEntries);
    thetaSigma = Spline(ksThetaBinCenters, deltaTheta);
}

Spline HandlerMC::CreateDeltaPsiSpline(std::unique_ptr<Tree> const &tree, std::string spline_filename, bool save)
{
    std::vector<Double_t> psi_delta = {};
    std::vector<Double_t> ksTheta_bins = {};
    TProfile psi("psi", "deltaPsi", 50, 0, 3.14, -0.1, 0.1);

    for(const auto &entry : *tree)
    {
        if(abs(tree->reco.Y - 1) > 1e-9 && tree->reco.piPos.nhit > 10 && tree->reco.piNeg.nhit > 10 && fabs(tree->reco.ks.theta - TMath::Pi() / 2) < 0.5)
        { psi.Fill(tree->reco.ks.theta, tree->reco.ksdpsi - tree->gen.ksdpsi); }
    }


    for(int i = 0; i < psi.GetNbinsX(); i++)
    {
        ksTheta_bins.push_back(psi.GetBinCenter(i));
        psi_delta.push_back(psi.GetBinContent(i));
    }

    return Spline(ksTheta_bins, psi_delta);
}

#endif