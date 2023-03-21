#include "C:\work\Science\BINP\Kaon Mass Measure\Energy.h"

int energyStability()
{
    gROOT->Reset();
    auto start = std::chrono::system_clock::now();
    const std::vector<double> meanEnergies = {504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 513.864};
    const std::vector<double> meanEnergiesErr = {0.007, 0.007, 0.008, 0.009, 0.004, 0.005, 0.007, 0.009, 0.009};
    const std::vector<std::string> energyPoints = {"505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "514"};
    std::vector<Energy*> handlers;

    for(int i = 0; i < energyPoints.size(); i++)
    { handlers.push_back(new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKpKm/kchExp" + energyPoints[i] + ".root", meanEnergies[i], meanEnergiesErr[i], 30, 2.9)); }
    // { handlers.push_back(new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/MC/Kch/MC" + energyPoints[i] + "_FSR.root", meanEnergies[i], 0.007, 1, 0, false)); }

    handlers[8]->DrawGraph();

    // std::cout << "{";
    // for(auto &handler : handlers)
    // {
    //     auto shift = handler->GetEnergyDiff();
    //     for (auto k : shift)
    //     { std::cout << -k.second << ", "; }
    // }
    // std::cout << "}" << std::endl;

    for(auto &handler : handlers)
    { delete handler; }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start; 
    std::cout << "exec time = " << diff.count() << std::endl;
    return 0;
}