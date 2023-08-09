#include <chrono>
#include <ctime> 
#include <iostream> 

#include "HandlerMC.h"


int KsMassEval()
{
    auto start = std::chrono::system_clock::now();

    const std::vector<double> meanEnergies_vec = {504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 513.864};
    // Means of Ks energy spectrums
    const std::vector<double> meanEnergiesSpectrum_vec = {504.683, 507.762, 508.323, 508.885, 509.445, 509.841, 510.263, 510.694, 512.297};
    const std::vector<double> meanEnergiesErr = {0.007, 0.007, 0.008, 0.009, 0.004, 0.005, 0.007, 0.009, 0.009};

    // RC without energy shift.
    // const std::vector<double> deltaM_RC_Smeared = {0.112, 0.081, 0.077, 0.068, 0.076, 0.108, 0.185, 0.324, 1.462};

    // RC with energy shifted by 134 keV
    const std::vector<double> deltaM_RC_Smeared = {0.1021, 0.0966, 0.087, 0.0759, 0.0808, 0.1336, 0.2191, 0.3726, 1.5223};
    const std::vector<std::string> energyPoints = {"505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "514"};

    std::map<std::string, std::pair<double, double>> meanEnergies;
    std::map<std::string, std::pair<double, double>> meanEnergiesSpectrum;
    std::map<std::string, double> radiativeCorrections;
    for(int i = 0; i < energyPoints.size(); i++)
    {   
        meanEnergies[energyPoints[i]] = std::make_pair(meanEnergies_vec[i], meanEnergiesErr[i]); 
        meanEnergiesSpectrum[energyPoints[i]] = std::make_pair(meanEnergiesSpectrum_vec[i], meanEnergiesErr[i]); 
        radiativeCorrections[energyPoints[i]] = deltaM_RC_Smeared[i]; 
    }

    std::string energyPoint = "509";
    std::string fileName = "C:/work/Science/BINP/Kaon Mass Measure/tr_ph/MC/KsKl_Smeared/MC" + energyPoint + ".root";
    // fileName = "tr_ph/MC/KsKl_Smeared/Field/MC" + energyPoint + "_Field.root";

    auto handler = new HandlerMC(fileName, energyPoint, 0.27, std::nullopt, true);
    auto [mass, massErr] = handler->Eval();
    std::cout << mass << " + " << massErr << std::endl;
    handler->Draw("hMlnYpfx", {-0.4, 0.4});
    // handler->SaveHists("hists" + energyPoint + ".root");
    delete handler;

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start; 
    std::cout << "exec time = " << diff.count() << std::endl;
    return 0;
}