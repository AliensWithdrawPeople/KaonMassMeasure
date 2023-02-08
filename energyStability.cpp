#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <iostream>

#include "TROOT.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"


class Energy
{
private:
    // Charged kaons tree.
    TTree *kTr;
    // map<runnum, hEnergy> == Energy distribution for each run.
    std::map<int, TH1D*> enHists; 

    // map<runnum, pair<emeas, demeas>> == emeas for each run.
    std::map<int, std::pair<Float_t, Float_t>> runAvgEmeas; 
    // map<runnum, pair<enVal, enErr>> == Kch energy for each run.
    std::map<int, std::pair<Float_t, Float_t>> runAvgEval;
    std::vector<std::vector<int>> emeasRunGroups; // Runs grouped by emeas.
    std::vector<std::vector<int>> enRunGroups; // Runs grouped by Kch avg energy.
    double energyShift;
    double comptonEnergyMean;
    double comptonEnergyError;

    std::vector<double> runs;
    Int_t runnum;
    Float_t emeas;
    Float_t demeas;
    Float_t ebeam;
    Float_t tdedx[2];
    Float_t tptot[2];
    const double kchMass = 493.677;

    /*
    * This function fills enHists, calculates energies 
    * and errors for all runs and fills runAvgEval and runAvgEerr with pair<enVal, enErr>.
    */
    int FillHists();
    // Divide according to emeas (laser system) groups.
    int DivideIntoGroups(int maxGroupSize = 4);
    std::map<int, std::pair<Float_t, Float_t>> AverageKchEnergy();
    /* 
    * Evaluates mean energy and energy error for a run. 
    * Return pair<energy, energy error>.
    */
    std::pair<double, double> Eval(int run); 
    std::pair<double, double> Eval(std::vector<int> &runGroup); 

public:
    int GetGroupsNum();
    void PrintGroups();
    void DrawGraph(int graphNum = -1);

    Energy(std::string fChargedK, double comptonEnergyMean, double comptonEnergyError, int maxGroupSize = 4, double shiftToKchEnergy = 4.);
    ~Energy();
};

Energy::Energy(std::string fChargedK, double comptonEnergyMean, double comptonEnergyError, int maxGroupSize = 4, double shiftToKchEnergy = 4.)
{
    TFile *file = TFile::Open(fChargedK.c_str());
    kTr = (TTree *)file->Get("kChargedTree");
    energyShift = shiftToKchEnergy;
    this->comptonEnergyMean = comptonEnergyMean;
    this->comptonEnergyError = comptonEnergyError;
    kTr->SetBranchAddress("ebeam", &ebeam);
    kTr->SetBranchAddress("emeas", &emeas);
    kTr->SetBranchAddress("demeas", &demeas);
    kTr->SetBranchAddress("runnum", &runnum);
    kTr->SetBranchAddress("tdedx", tdedx);
    kTr->SetBranchAddress("tptot", tptot);

    FillHists();
    DivideIntoGroups(maxGroupSize);
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
    std::set<double> runs_;
    for(int i = 0; i < kTr->GetEntriesFast(); i++)
    {
        kTr->GetEntry(i);
        if(fabs(demeas) < 1e-8 || emeas < 100)
        { continue; }
        runs_.insert(double(runnum));
        if(enHists.count(runnum) <= 0)
        { 
            auto tmpStr = ("hEn_" + std::to_string(runnum)).c_str();
            enHists[runnum] = new TH1D(tmpStr, tmpStr, 3000, 490, 520);
        }
        enHists[runnum]->Fill(sqrt(tptot[0] * tptot[0] + kchMass * kchMass));
        enHists[runnum]->Fill(sqrt(tptot[1] * tptot[1] + kchMass * kchMass));
        runAvgEmeas[runnum] = std::make_pair(emeas, demeas);
    }
    runs.insert(runs.end(), runs_.begin(), runs_.end());
    // Calculate energies and errors and put them to a dictionary. 
    for (auto& [run, hist]: enHists) 
    { runAvgEval[run] = Eval(run); }
    return 0;
}

std::pair<double, double> Energy::Eval(int run)
{
    if(enHists.count(run) > 0)
    { 
        auto tmpHist = new TH1D(*enHists[run]);
        tmpHist->Rebin(8);
        auto res = tmpHist->Fit("gaus", "SEQM", "", tmpHist->GetBinCenter(tmpHist->GetMaximumBin()) - tmpHist->GetRMS(), 
                                                    tmpHist->GetBinCenter(tmpHist->GetMaximumBin()) + 2 * tmpHist->GetRMS());
        res = tmpHist->Fit("gaus", "SEQM", "", res->Parameter(1) - res->Parameter(2), res->Parameter(1) + 2 * res->Parameter(2));
        if(res->Chi2()/res->Ndf() > 1.5)
        { std::cout << run << " run: chi2 / ndf = " << res->Chi2()/res->Ndf() << std::endl; }
        delete tmpHist;
        return std::make_pair(res->Parameter(1), res->ParError(1)); 
        // return std::make_pair(enHists[run]->GetMean(), enHists[run]->GetMeanError()); 
    }
    return std::make_pair(-1, 0.);
}

std::pair<double, double> Energy::Eval(std::vector<int> &runGroup)
{
    if(runGroup.size() == 0)
    { return std::make_pair(-1, 0.); }
    auto tmpHist = new TH1D(("tmpHist" + std::to_string(runGroup[0])).c_str(), 
                            (std::to_string(runGroup[0]) + "-" + std::to_string(runGroup.back())).c_str(), 3000, 490, 520);
    for(int i = 0; i < tmpHist->GetNbinsX(); i++)
    {
        for(auto run : runGroup)
        { tmpHist->SetBinContent(i, tmpHist->GetBinContent(i) + enHists[run]->GetBinContent(i)); }
    }
    tmpHist->Rebin(8);
    auto res = tmpHist->Fit("gaus", "SEQM", "", tmpHist->GetBinCenter(tmpHist->GetMaximumBin()) - tmpHist->GetRMS(), 
                                                tmpHist->GetBinCenter(tmpHist->GetMaximumBin()) + 1.5 * tmpHist->GetRMS());
    res = tmpHist->Fit("gaus", "SEQM", "", res->Parameter(1) - res->Parameter(2), res->Parameter(1) + 1.5 * res->Parameter(2));
    if(res->Chi2()/res->Ndf() > 1.5)
    { std::cout << runGroup[0] << "-" << runGroup.back() << " runs: chi2 / ndf = " << res->Chi2()/res->Ndf() << std::endl; }
    // delete tmpHist;
    return std::make_pair(res->Parameter(1), res->ParError(1)); 
}

int Energy::DivideIntoGroups(int maxGroupSize)
{
    emeasRunGroups = {{int(runs[0])}};
    for(int i = 1; i < runs.size(); i++)
    {
        if(fabs(runAvgEmeas[int(runs[i])].first - runAvgEmeas[int(runs[i - 1])].first) < 1e-8)
        { emeasRunGroups.back().push_back(int(runs[i])); }
        else
        { emeasRunGroups.push_back({int(runs[i])}); }
    }
    
    for(int i = 0; i < emeasRunGroups.size(); i++)
    {
        if(emeasRunGroups[i].size() <= maxGroupSize)
        {
            enRunGroups.push_back(emeasRunGroups[i]);
            continue; 
        }
        auto it = emeasRunGroups[i].begin();
        while(it + maxGroupSize < emeasRunGroups[i].end())
        { 
            enRunGroups.push_back(std::vector<int>(it, it + maxGroupSize)); 
            it = it + maxGroupSize;
        }
        enRunGroups.push_back(std::vector<int>(it - maxGroupSize, emeasRunGroups[i].end())); 
    }
    return 0;
}

std::map<int, std::pair<Float_t, Float_t>> Energy::AverageKchEnergy()
{
    std::map<int, std::pair<Float_t, Float_t>> averaged;
    Float_t tmpEnergy = 0;
    Float_t tmpEnergyErr = 0;
    for(auto &group : enRunGroups)
    {
        // tmpEnergy = 0;
        // tmpEnergyErr = 0;
        // for(auto run : group)
        // {
        //     tmpEnergy += runAvgEval[run].first;
        //     tmpEnergyErr += 1 / runAvgEval[run].second / runAvgEval[run].second;
        // }
        // tmpEnergy = tmpEnergy / group.size();
        // tmpEnergyErr = sqrt(1 / tmpEnergyErr);
        auto tmpPair = Eval(group);
        tmpEnergy = tmpPair.first;
        tmpEnergyErr = tmpPair.second;
        for(auto run : group)
        { averaged[run] = std::make_pair(tmpEnergy, tmpEnergyErr); }
    }
    return averaged;
}

void Energy::DrawGraph(int graphNum = -1)
{
    std::vector<double> runnums;
    std::vector<double> enVals;
    std::vector<double> enErrs;
    std::vector<double> emeasVals;
    std::vector<double> emeasErrs;
    auto enValsAveraged = AverageKchEnergy();

    std::transform(runAvgEmeas.begin(), runAvgEmeas.end(), std::back_inserter(runnums),
                        [](std::pair<int, std::pair<Float_t, Float_t>> val) { return double(val.first); });

    std::transform(runAvgEmeas.begin(), runAvgEmeas.end(), std::back_inserter(emeasVals),
                        [](std::pair<int, std::pair<Float_t, Float_t>> val) { return double(val.second.first); });
    std::transform(runAvgEmeas.begin(), runAvgEmeas.end(), std::back_inserter(emeasErrs),
                        [](std::pair<int, std::pair<Float_t, Float_t>> val) { return double(val.second.second); });
    
    std::transform(enValsAveraged.begin(), enValsAveraged.end(), std::back_inserter(enVals),
                        [&](std::pair<int, std::pair<Float_t, Float_t>> val) { return double(val.second.first + energyShift); });
    std::transform(enValsAveraged.begin(), enValsAveraged.end(), std::back_inserter(enErrs),
                        [](std::pair<int, std::pair<Float_t, Float_t>> val) { return double(val.second.second); });
    
    std::vector<double> zeroes(runnums.size(), 0.0);
    TGraphErrors grEnVsRun(runnums.size(), runnums.data(), enVals.data(), zeroes.data(), enErrs.data());
    TGraphErrors grEmeasVsRun(runnums.size(), runnums.data(), emeasVals.data(), zeroes.data(), emeasErrs.data());

    TGraphErrors grComptonEnVsRun;
    grComptonEnVsRun.AddPoint(runs[0] - 10, comptonEnergyMean);
    grComptonEnVsRun.SetPointError(0, 0, comptonEnergyError);
    grComptonEnVsRun.SetLineColor(kBlue);
    grComptonEnVsRun.SetMarkerColor(kBlue);
    TLine comptonMean(runs[0], comptonEnergyMean, runs.back(), comptonEnergyMean);
    comptonMean.SetLineColor(kYellow);
    
    grEmeasVsRun.GetXaxis()->SetTitle("Run");
    grEmeasVsRun.GetYaxis()->SetTitle("Energy, MeV");
    grEmeasVsRun.SetTitle("Red -- emeas, black -- Kch E_inv, yellow(blue) -- compton mean for energy point");
    grEmeasVsRun.SetMarkerColor(kRed);
    grEmeasVsRun.SetLineColor(kRed);
    grEmeasVsRun.DrawClone("AP");
    grEnVsRun.DrawClone("P same");
    grComptonEnVsRun.DrawClone("P same");
    comptonMean.DrawClone("same");
}

int Energy::GetGroupsNum()
{ return emeasRunGroups.size(); }

void Energy::PrintGroups()
{
    for(auto &group : emeasRunGroups)
    {
        std::cout << "size = " << group.size() << "; ";
        for(auto rnum : group)
        { std::cout << rnum << ", "; } 
        std::cout << std::endl;
    }
}

int energyStability()
{
    gROOT->Reset();
    auto start = std::chrono::system_clock::now();

    auto enHandler = new Energy("C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKpKm/expKch_509.root", 508.957, 0.009, 3);
    enHandler->DrawGraph();
    delete enHandler;

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start; 
    std::cout << "exec time = " << diff.count() << std::endl;
    return 0;
}