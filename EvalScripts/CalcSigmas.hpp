#ifndef CalcSigmas_hpp
#define CalcSigmas_hpp

#include "Tree.hpp"

#include <memory>
#include <tuple>
#include <iomanip>

#include "TMath.h"
#include "TH1D.h"
#include "TVector3.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

namespace Sigmas {

    /// @returns Bin numbers in the sigma^{(M_K)}_{psi} [lnY][ksTheta] matrix.
    std::pair<int, int> GetSigmaMatrix_bin(double lnY, double ksTheta)
    {
        auto binY = lnY > -0.3? int(floor((lnY + 0.3) / 0.075)) : -1;
        auto binKsTheta = ksTheta > 0.57? int(floor((ksTheta - 0.57) / 0.25)) : -1;
        binY = binY < 8? binY : -1;
        binKsTheta = binKsTheta < 8? binKsTheta : -1;
        return {binY, binKsTheta};
    }

    /// @brief Calc sigma^{(M_K)}_{psi} [lnY][ksTheta].
    /// @param tree -- pointer to a tree with data;
    /// @param goodEntries -- numbers (ids) of entries that satisfy all selection cuts.
    /// @param left_window_width -- number of sigmas from the mean to the left: mean - left_window_width * sigma.
    /// @param right_window_width -- number of sigmas from the mean to the right: mean + left_window_width * sigma.
    /// @param isVerbose -- prints out the matrix if isVerbose == true.
    /// @returns std::vector<std::vector<double>> vSigmaFitMatrix with shape (8, 12).
    std::vector<std::vector<double>> GetSigmaMatrix(std::unique_ptr<Tree> const &tree, const std::vector<int> &goodEntries, double left_window_width, double right_window_width, bool isVerbose)
    {
        // Distributions of ksdpsi (angle between pions' momentums) in different bins (lnY, kstheta).
        std::vector<std::vector<std::unique_ptr<TH1D>>> psiDistrs;
        std::vector<std::vector<double>> vSigmaMatrixFit;
        for(int i = 0; i < 8; i++)
        {
            vSigmaMatrixFit.push_back({});
            psiDistrs.push_back(std::vector<std::unique_ptr<TH1D>>());
            for(int j = 0; j < 8; j++)
            {
                std::string histName = "py_bin_lnY" + std::to_string(i) + "_kstheta" + std::to_string(j); 
                psiDistrs[i].push_back(std::make_unique<TH1D>(TH1D(histName.c_str(), histName.c_str(), 2000, -3.1415, 3.1415))); 
            }
        }

        for(const auto &entry : goodEntries)
        {   
            tree->GetEntry(entry);
            auto [binY, binKsTheta] = GetSigmaMatrix_bin(log(tree->reco.Y), tree->reco.ks.theta);
            if(binY != -1 && binKsTheta != -1)
            { psiDistrs[binY][binKsTheta]->Fill(tree->reco.ksdpsi); } 
        }

        TFitResultPtr r;
        for(int i = 0; i < 8; i++)
        {
            for(int j = 0; j < 8; j++)
            {
                if(psiDistrs[i][j]->GetEntries() < 150)
                { 
                    vSigmaMatrixFit[i].push_back(psiDistrs[i][j]->GetStdDev()); 
                    continue;
                }

                if(psiDistrs[i][j]->GetEntries() < 600)
                { psiDistrs[i][j]->Rebin(2 * int(3 - psiDistrs[i][j]->GetEntries() / 300)); }

                auto leftBorder = psiDistrs[i][j]->GetBinCenter(psiDistrs[i][j]->GetMaximumBin()) - left_window_width * psiDistrs[i][j]->GetStdDev();
                auto rightBorder = psiDistrs[i][j]->GetBinCenter(psiDistrs[i][j]->GetMaximumBin()) + right_window_width * psiDistrs[i][j]->GetStdDev();
                r = psiDistrs[i][j]->Fit("gaus", "SEQ", "", leftBorder, rightBorder);
                r = psiDistrs[i][j]->Fit("gaus", "SEQ", "", r->Parameter(1) - left_window_width * r->Parameter(2), r->Parameter(1) + right_window_width * r->Parameter(2));
                vSigmaMatrixFit[i].push_back(r->Parameter(2));
            }
        }

        if(isVerbose)
        {
            auto default_precision = std::cout.precision();
            std::cout << "vSigmaMatrixFit:" << std::endl;
            for(int i = 0; i < 8; i++)
            {
                for(int j = 0; j < 8; j++)
                { std::cout << std::left << std::setprecision(4) << std::setw(8) << vSigmaMatrixFit[i][j] << " "; }
                std::cout << std::endl;
            }
            std::cout << std::setprecision(default_precision) << std::endl;
        }
        return vSigmaMatrixFit;
    }


    int GetThetaBin(double theta)
    {
        auto bin = int(floor((theta - 1.1) / 0.094));
        return bin > -1 && bin < 10? bin : -1;
    }

    /// @brief Calc sigma^{(psi)}_{theta} and its error for different ksTheta.
    /// @param tree -- pointer to a tree with data;
    /// @param goodEntries -- numbers of entries that satisfy all selection cuts.
    /// @returns tuple<ksTheta bin centers, sigmas, errors of sigma>.
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> GetThetaSigmas(std::unique_ptr<Tree> const &tree, const std::vector<int> &goodEntries)
    {
        // Theta distributions for different ksTheta.
        std::vector<std::unique_ptr<TH1D>> piThetaProjs;

        for(int i = 0; i < 10; i++)
        { piThetaProjs.push_back(std::make_unique<TH1D>(TH1D(("piThetaProjs" + std::to_string(i + 1)).c_str(), ("piThetaProjs" + std::to_string(i + 1)).c_str(), 1000, -0.1, 0.1))); }

        for(const auto &entry : goodEntries)
        {
            tree->GetEntry(entry);

            if(auto bin = GetThetaBin(tree->reco.piPos.theta); bin != -1)
            { piThetaProjs[bin]->Fill(tree->reco.piPos.theta - tree->gen.piPos.theta); }

            if(auto bin = GetThetaBin(tree->reco.piNeg.theta); bin != -1)
            { piThetaProjs[bin]->Fill(tree->reco.piNeg.theta - tree->gen.piNeg.theta); }
        }

        std::vector<double> binCenters;
        std::vector<double> sigmas;
        std::vector<double> sigmaErrors;
        for(int i = 0; i < piThetaProjs.size(); i++)
        {
            binCenters.push_back(1.1 + (i + i + 1) / 2 * 0.094 - TMath::Pi() / 2);
            sigmas.push_back(piThetaProjs[i]->GetRMS());
            sigmaErrors.push_back(piThetaProjs[i]->GetRMSError());
        }

        return {binCenters, sigmas, sigmaErrors};
    }

    /// @brief Calc sigma^{(psi)}_{phi}.
    /// @return sigmaPhi.
    double GetPhiSigma()
    {
        return 7.12804e-03;
    }

    /// @brief Calc sigma^_{Y}.
    /// @return sigmaY.
    double GetYSigma()
    {
        return 0.018;
    }
}

#endif