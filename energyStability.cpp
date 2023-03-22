#include "C:\work\Science\BINP\Kaon Mass Measure\Energy.h"

int energyStability()
{
    gROOT->Reset();
    auto start = std::chrono::system_clock::now();
    const std::vector<double> meanEnergies = {504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 513.864};
    const std::vector<double> meanEnergiesErr = {0.007, 0.007, 0.008, 0.009, 0.004, 0.005, 0.007, 0.009, 0.009};
    const std::vector<std::string> energyPoints = {"505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "514"};

    std::map<std::string, Energy*> enDict;
    // for(int i = 0; i < energyPoints.size(); i++)
    // { 
    //     enDict[energyPoints[i]] = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKpKm/kchExp" + energyPoints[i] + ".root", 
    //                                                                                             meanEnergies[i], meanEnergiesErr[i], 30, 2.9);        
    // }

    enDict["510"] = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKpKm/KpKm_NoVertex/kchExp510_noVertex.root", meanEnergies[5], meanEnergiesErr[5], 30, 0);

    enDict["510"]->DrawGraph();

    // std::cout << "{";
    // for(auto &handler : enDict)
    // {
    //     auto shift = handler.second->GetEnergyDiff();
    //     for (auto k : shift)
    //     { std::cout << -k.second << ", "; }
    // }
    // std::cout << "}" << std::endl;

    for(auto &handler : enDict)
    { delete handler.second; }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start; 
    std::cout << "exec time = " << diff.count() << std::endl;
    return 0;
}