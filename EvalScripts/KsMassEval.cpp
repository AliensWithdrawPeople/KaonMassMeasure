#include <chrono>
#include <ctime> 
#include <iostream> 

#include "HandlerMC.hpp"
#include "HandlerExp.hpp"


int KsMassEval()
{
    auto start = std::chrono::system_clock::now();

    const std::vector<double> meanEnergies_vec = {504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 511.444, 513.864};
    // Means of Ks energy spectrums
    const std::vector<double> meanEnergiesSpectrum_vec = {504.683, 507.762, 508.323, 508.885, 509.445, 509.841, 510.263, 510.694,  510.694, 512.297};
    const std::vector<double> meanEnergiesErr = {0.007, 0.007, 0.008, 0.009, 0.004, 0.005, 0.007, 0.009, 0.009, 0.009};


    const std::vector<double> deltaE_RC_Smeared = {0.105504, 0.0760349, 0.069423, 0.0598971, 0.0769491, 0.118673, 0.197142, 0.336072, 0.467132, 1.54063};
    const std::vector<std::string> energyPoints = {"505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514"};

    std::map<std::string, std::pair<double, double>> meanEnergies;
    std::map<std::string, std::pair<double, double>> meanEnergiesSpectrum;
    std::map<std::string, double> radiativeCorrections;
    for(int i = 0; i < energyPoints.size(); i++)
    {   
        meanEnergies[energyPoints[i]] = std::make_pair(meanEnergies_vec[i], meanEnergiesErr[i]); 
        meanEnergiesSpectrum[energyPoints[i]] = std::make_pair(meanEnergiesSpectrum_vec[i], meanEnergiesErr[i]); 
        radiativeCorrections[energyPoints[i]] = deltaE_RC_Smeared[i]; 
    }

    std::string energyPoint = "509.5";
    
    std::string fileNameMC = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/MC/KsKl_Smeared/New formfactor/MC" + energyPoint + ".root";
    auto handlerMC = new HandlerMC(fileNameMC, energyPoint, 0.27, meanEnergies[energyPoint].first, true, true);
    auto [mass, massErr] = handlerMC->Eval();
    std::cout << mass << " + " << massErr << std::endl;
    handlerMC->Draw("hMlnYpfx", {-0.4, 0.4});
    handlerMC->SaveHists("C:/work/Science/BINP/Kaon Mass Measure/hists/MC/Hists_MC" + energyPoint + ".root");
    // handlerMC->SaveSplines("C:/work/Science/BINP/Kaon Mass Measure/splines/spline_" + energyPoint + ".root");
    delete handlerMC;

    // std::string fileNameExp = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/expKsKl/exp" + energyPoint + ".root";
    // auto handlerExp = new HandlerExp(fileNameExp, energyPoint, 0.27, meanEnergies[energyPoint].first, radiativeCorrections[energyPoint], true);
    // auto [mass, massErr] = handlerExp->Eval();
    // std::cout << mass << " + " << massErr << std::endl;
    // handlerExp->Draw("hMlnYpfx", {-0.4, 0.4});
    // handlerExp->SaveHists("C:/work/Science/BINP/Kaon Mass Measure/hists/Exp/Hists_Exp" + energyPoint + ".root");
    // delete handlerExp;

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start; 
    std::cout << "exec time = " << diff.count() << std::endl;
    return 0;
}