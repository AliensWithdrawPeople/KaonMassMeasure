#include <vector>
#include <unordered_set>
#include <map>
#include <chrono>
#include <ctime>
#include <iostream>

#include "TROOT.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"


class Energy
{
private:
    TTree *kTr; // Charged kaons tree.
    std::map<int, TH1D*> enHists; // map<runnum, hEnergy> == Energy distribution for each run.

    std::map<int, Float_t> runAvgEmeas; // map<runnum, emeas> == emeas for each run.
    std::map<int, Float_t> runAvgDEmeas; // map<runnum, demeas> == demeas for each run.
    std::map<int, Float_t> runAvgEval; // map<runnum, enVal> == Kch energy for each run.
    std::map<int, Float_t> runAvgEerr; // map<runnum, enErr> == Kch energy error for each run.
    std::vector<std::vector<int>> emeasRunGroups; // Runs grouped by emeas.
    std::vector<std::vector<int>> enRunGroups; // Runs grouped by Kch avg energy.

    std::vector<double> runs;
    Int_t runnum;
    Float_t emeas;
    Float_t demeas;
    Float_t ebeam;
    Float_t tdedx[2];
    Float_t tptot[2];
    const double kchMass = 493.677;

    int FillHists(); // This function fills enHists.
    void DivideIntoGroops();
    std::map<int, Float_t> AverageKchEnergy();
    void CalcEnergies(); // Calculates energies for all runs and fills runAvgEval and runAvgEerr
    std::pair<double, double> Eval(int run); // Evaluates mean energy and energy error for a run

public:
    int GetGroupsNum();
    void PrintGroups();
    void DrawGraph(int graphNum = -1);
    Energy(std::string fChargedK);
    ~Energy();
};

Energy::Energy(std::string fChargedK)
{
    TFile *file = TFile::Open(fChargedK.c_str());
    kTr = (TTree *)file->Get("kChargedTree");

    kTr->SetBranchAddress("ebeam", &ebeam);
    kTr->SetBranchAddress("emeas", &emeas);
    kTr->SetBranchAddress("demeas", &demeas);
    kTr->SetBranchAddress("runnum", &runnum);
    kTr->SetBranchAddress("tdedx", tdedx);
    kTr->SetBranchAddress("tptot", tptot);

    FillHists();
    DivideIntoGroops();
    CalcEnergies();
}

Energy::~Energy()
{
    for (auto& [run, hist]: enHists) 
    { delete hist; }
    enHists.clear();
    delete kTr;
}

int Energy::FillHists()
{
    std::unordered_set<double> runs_;
    for(int i = 0; i < kTr->GetEntriesFast(); i++)
    {
        kTr->GetEntry(i);
        runs_.insert(double(runnum));
        if(enHists.count(runnum) <= 0)
        { 
            auto tmpStr = ("hEn_" + std::to_string(runnum)).c_str();
            enHists[runnum] = new TH1D(tmpStr, tmpStr, 3000, 490, 520);
        }
        enHists[runnum]->Fill(sqrt(tptot[0] * tptot[0] + kchMass * kchMass));
        enHists[runnum]->Fill(sqrt(tptot[1] * tptot[1] + kchMass * kchMass));
        runAvgEmeas[runnum] = emeas;
        runAvgDEmeas[runnum] = demeas;
    }
    runs.insert(runs.end(), runs_.begin(), runs_.end());
    return 0;
}

std::pair<double, double> Energy::Eval(int run)
{
    if(enHists.count(run) > 0)
    { return std::make_pair(enHists[run]->GetMean(), enHists[run]->GetMeanError()); }
    return std::make_pair(0., 0.);
}

void Energy::CalcEnergies()
{
    std::pair<double, double> tmpPair;
    for (auto& [run, hist]: enHists) 
    {
        tmpPair = Eval(run);
        runAvgEval[run] = tmpPair.first;
        runAvgEerr[run] = tmpPair.second;
    }
}

void Energy::DivideIntoGroops()
{
    emeasRunGroups = {{int(runs[0])}};
    for(int i = 1; i < runs.size(); i++)
    {
        if(fabs(runAvgEmeas[int(runs[i])] - runAvgEmeas[int(runs[i - 1])]) < 1e-6)
        { emeasRunGroups.back().push_back(int(runs[i])); }
        else
        { emeasRunGroups.push_back({int(runs[i])}); }
    }

    for(int i = 0; i < emeasRunGroups.size(); i++)
    {
        if(emeasRunGroups[i].size() < 4)
        {
            enRunGroups.push_back(emeasRunGroups[i]);
            continue; 
        }
        auto it = emeasRunGroups[i].begin();
        while(it + 3 < emeasRunGroups[i].end())
        { 
            enRunGroups.push_back(std::vector<int>(it, it + 3)); 
            it = it + 3;
        }
        enRunGroups.push_back(std::vector<int>(it - 3, emeasRunGroups[i].end())); 
    }
}

std::map<int, Float_t> Energy::AverageKchEnergy()
{
    std::pair<double, double> tmpPair;
    std::map<int, Float_t> averaged;
    for(auto &group : enRunGroups)
    {
        for(auto run : group)
        {
            tmpPair = Eval(run);
           // Add averaging Kch energy by group of runs. 
        }
    }
}

void Energy::DrawGraph(int graphNum = -1)
{
    std::vector<double> enVals;
    std::vector<double> enErrs;
    std::vector<double> emeasVals;
    std::vector<double> emeasErrs;
    std::pair<double, double> tmpPair;
    for (auto& [run, hist]: enHists) 
    {
        // Change to std::transform(...).
        tmpPair = Eval(run);
        enVals.push_back(tmpPair.first + 4);
        emeasVals.push_back(runAvgEmeas[run]);
        emeasErrs.push_back(0);
        enErrs.push_back(tmpPair.second);
    }

    std::vector<double> zeroes(runs.size(), 0.0);
    TGraphErrors grEnVsRun(runs.size(), runs.data(), enVals.data(), zeroes.data(), enErrs.data());
    TGraphErrors grEmeasVsRun(runs.size(), runs.data(), emeasVals.data(), zeroes.data(), emeasErrs.data());
    grEnVsRun.DrawClone("AP");
    grEmeasVsRun.SetMarkerColor(kRed);
    grEmeasVsRun.DrawClone("P same");
}

int Energy::GetGroupsNum()
{ return emeasRunGroups.size(); }

void Energy::PrintGroups()
{
    for(auto &group : emeasRunGroups)
    {
        for(auto rnum : group)
        { std::cout << rnum << ", "; } 
        std::cout << std::endl;
    }
}

int energyStability()
{
    gROOT->Reset();
    auto start = std::chrono::system_clock::now();

    auto enHandler = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKpKm/expKch_509.root");
    enHandler->DrawGraph();
    delete enHandler;

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start; 
    std::cout << "exec time = " << diff.count() << std::endl;
    return 0;
}