#include <vector>
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
    std::map<int, Float_t> runAvgEmeas; // map<runnum, hEnergy> == Energy distribution for each run.
    std::map<int, Float_t> runAvgDEmeas; // map<runnum, hEnergy> == Energy distribution for each run.
    Int_t runnum;
    Float_t emeas;
    Float_t demeas;
    Float_t ebeam;
    Float_t tdedx[2];
    Float_t tptot[2];
    const double kchMass = 493.677;

    int FillHists(); // This function fills enHists.
    std::pair<double, double> Eval(int run); // Evaluate mean energy and energy error for a run

public:
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
    for(int i = 0; i < kTr->GetEntriesFast(); i++)
    {
        kTr->GetEntry(i);
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

    return 0;
}

std::pair<double, double> Energy::Eval(int run)
{
    if(enHists.count(run) > 0)
    { return std::make_pair(enHists[run]->GetMean(), enHists[run]->GetMeanError()); }
    return std::make_pair(0., 0.);
}

void Energy::DrawGraph(int graphNum = -1)
{
    
    std::vector<double> runs;
    std::vector<double> enVals;
    std::vector<double> emeasVals;
    std::vector<double> enErrs;
    std::vector<double> emeasErrs;
    std::pair<double, double> tmpPair;
    for (auto& [run, hist]: enHists) 
    {
        tmpPair = Eval(run);
        runs.push_back(run);
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

int energyStability()
{
    gROOT->Reset();
    auto start = std::chrono::system_clock::now();

    auto enHandler = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKch_509.root");
    enHandler->DrawGraph();
    delete enHandler;

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start; 
    std::cout << "exec time = " << diff.count() << std::endl;
    return 0;
}