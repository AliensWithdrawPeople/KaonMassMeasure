#pragma once
#ifndef Energy_h
#define Energy_h

#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <iostream>
#include <fstream>

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
    std::vector<int> badRuns;

    std::vector<double> runs;
    Int_t runnum;
    Float_t emeas;
    Float_t demeas;
    Float_t ebeam;
    Float_t tdedx[2];
    Float_t tptot[2];

    const double kchMass = 493.677;
    bool verbose;
    bool isExp;

    /*
    * This function fills enHists, calculates energies 
    * and errors for all runs and fills runAvgEval and runAvgEerr with pair<enVal, enErr>.
    */
    int FillHists();
    // Divide according to emeas (laser system) groups. Return number of Kch energy groups.
    int DivideIntoGroups(int maxGroupSize = 4);
    std::map<int, std::pair<Float_t, Float_t>> AverageKchEnergy(bool isValForEachRun = false);
    /* 
    * Evaluates mean energy and energy error for a run. 
    * Return pair<energy, energy error>.
    */
    std::pair<double, double> Eval(int run, bool isVerbose = true); 
    std::pair<double, double> Eval(std::vector<int> &runGroup, bool isVerbose = true); 
    int ReadBadRuns(std::string filename);

public:
    int GetGroupsNum();
    void PrintGroups();
    // Draws a Energy vs Run. Returns a pair (KchMeanEnergy, KchMeanEnergyErr).
    std::pair<double, double> DrawGraph(int graphNum = -1);
    // Return map of (groupNum, KchEnergy(groupNum) + energyShift - comptonEnergyMean) pairs.
    std::map<int, Float_t> GetEnergyDiff();

    Energy(std::string fChargedK, std::string fBadRunsList, double comptonEnergyMean, double comptonEnergyError, int maxGroupSize = 4, double shiftToKchEnergy = 4., bool isExp = true, bool isVerbose = false);
    ~Energy();
};
#endif

#define Energy_cpp
Energy::Energy(std::string fChargedK, std::string fBadRunsList, double comptonEnergyMean, double comptonEnergyError, int maxGroupSize, double shiftToKchEnergy, bool isExp, bool isVerbose)
{
    TFile *file = TFile::Open(fChargedK.c_str());
    kTr = (TTree *)file->Get("kChargedTree");
    
    energyShift = shiftToKchEnergy;
    verbose = isVerbose;
    this->isExp = isExp;
    this->comptonEnergyMean = comptonEnergyMean;
    this->comptonEnergyError = comptonEnergyError;
    kTr->SetBranchAddress("ebeam", &ebeam);
    kTr->SetBranchAddress("emeas", &emeas);
    kTr->SetBranchAddress("demeas", &demeas);
    kTr->SetBranchAddress("runnum", &runnum);
    kTr->SetBranchAddress("tdedx", tdedx);
    kTr->SetBranchAddress("tptot", tptot);
    if(fBadRunsList != "")
    { ReadBadRuns(fBadRunsList); }
    FillHists();
    DivideIntoGroups(maxGroupSize);
}

int Energy::ReadBadRuns(std::string filename)
{
    std::ifstream input(filename);
    if(input.is_open() && input.good())
    {
        std::istream_iterator<double> start(input), end;
        badRuns.insert(badRuns.end(), start, end);
    }
    input.close();
    return badRuns.size();
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

        if(std::find(badRuns.begin(), badRuns.end(), runnum) != badRuns.end())
        { continue; }
        if(isExp && (fabs(demeas) < 1e-8 || emeas < 100))
        { continue; }
        runs_.insert(double(runnum));
        if(enHists.count(runnum) <= 0)
        { 
            std::string tmpStr = "hEn_" + std::to_string(runnum);
            enHists[runnum] = new TH1D(tmpStr.c_str(), tmpStr.c_str(), 2000, 480, 520);
        }
        enHists[runnum]->Fill(sqrt(tptot[0] * tptot[0] + kchMass * kchMass));
        enHists[runnum]->Fill(sqrt(tptot[1] * tptot[1] + kchMass * kchMass));
        runAvgEmeas[runnum] = std::make_pair(emeas, demeas);
    }
    runs.insert(runs.end(), runs_.begin(), runs_.end());
    // Calculate energies and errors and put them to a dictionary. 
    for (auto& [run, hist]: enHists) 
    { runAvgEval[run] = Eval(run, verbose); }
    return 0;
}

std::pair<double, double> Energy::Eval(int run, bool isVerbose)
{
    if(enHists.count(run) > 0)
    { 
        auto tmpHist = new TH1D(*enHists[run]);
        tmpHist->Rebin(8);

        // int bin = 0;
        // double leftBorder = tmpHist->GetBinWithContent(tmpHist->GetMaximum() * 0.4, bin, 0, tmpHist->GetMaximumBin(), 1e-2);
        // leftBorder = tmpHist->GetBinCenter(bin);

        // double rightBorder = tmpHist->GetBinWithContent(tmpHist->GetMaximum() * 0.9, bin, tmpHist->GetMaximumBin(), tmpHist->GetNbinsX()-1, 1e-2);
        // rightBorder = tmpHist->GetBinCenter(bin);

        // auto res = tmpHist->Fit("gaus", "SEQM", "", leftBorder, rightBorder);

        auto res = tmpHist->Fit("gaus", "SEQM", "", tmpHist->GetBinCenter(tmpHist->GetMaximumBin()) - tmpHist->GetRMS(), 
                                                    tmpHist->GetBinCenter(tmpHist->GetMaximumBin()) + 2 * tmpHist->GetRMS());                                            
        if(res == -1)
        { return std::make_pair(0, 0); }

        // leftBorder = tmpHist->GetBinWithContent(res->Parameter(0) * 0.4, bin, 0, tmpHist->FindBin(res->Parameter(1)), 1e-2);
        // leftBorder = tmpHist->GetBinCenter(bin);

        // rightBorder = tmpHist->GetBinWithContent(res->Parameter(0) * 0.9, bin, tmpHist->FindBin(res->Parameter(1)), tmpHist->GetNbinsX()-1, 1e-2);
        // rightBorder = tmpHist->GetBinCenter(bin);
        // res = tmpHist->Fit("gaus", "SEQM", "", leftBorder, rightBorder);


        res = tmpHist->Fit("gaus", "SEQM", "", res->Parameter(1) - res->Parameter(2), res->Parameter(1) + 2 * res->Parameter(2));
        if(res == -1)
        { return std::make_pair(0, 0); }

        if(verbose && res->Chi2()/res->Ndf() > 1.5)
        { std::cout << run << " run: chi2 / ndf = " << res->Chi2()/res->Ndf() << std::endl; }
        delete tmpHist;
        return std::make_pair(res->Parameter(1), res->ParError(1)); 
        // return std::make_pair(enHists[run]->GetMean(), enHists[run]->GetMeanError()); 
    }
    return std::make_pair(-1, 0.);
}

std::pair<double, double> Energy::Eval(std::vector<int> &runGroup, bool isVerbose)
{
    if(runGroup.size() == 0)
    { return std::make_pair(-1, 0.); }
    auto tmpHist = new TH1D(("tmpHist" + std::to_string(runGroup[0]) + "to" + std::to_string(runGroup.back())).c_str(), 
                            ("tmpHist" + std::to_string(runGroup[0]) + "to" + std::to_string(runGroup.back())).c_str(), 2000, 480, 520);
    for(int i = 0; i < tmpHist->GetNbinsX(); i++)
    {
        for(auto run : runGroup)
        { tmpHist->SetBinContent(i, tmpHist->GetBinContent(i) + enHists[run]->GetBinContent(i)); }
    }

    tmpHist->Rebin(4);
    auto res = tmpHist->Fit("gaus", "SEQM", "", tmpHist->GetBinCenter(tmpHist->GetMaximumBin()) - tmpHist->GetRMS(), 
                                                tmpHist->GetBinCenter(tmpHist->GetMaximumBin()) + 2 * tmpHist->GetRMS());
    if(res == -1) { return std::make_pair(0, 0); }                                            
    res = tmpHist->Fit("gaus", "SEQM", "", res->Parameter(1) - res->Parameter(2), res->Parameter(1) + 2 * res->Parameter(2));
    if(res == -1) { return std::make_pair(0, 0); }                                            
    res = tmpHist->Fit("gaus", "SEQM", "", res->Parameter(1) - res->Parameter(2), res->Parameter(1) + 2 * res->Parameter(2));
    if(res == -1) { return std::make_pair(0, 0); }       
    
    if(verbose && res->Chi2()/res->Ndf() > 1.5)
    { std::cout << runGroup[0] << "-" << runGroup.back() << " runs: chi2 / ndf = " << res->Chi2()/res->Ndf() << std::endl; }

    delete tmpHist;
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
        auto it = emeasRunGroups[i].begin();
        while(it + maxGroupSize < emeasRunGroups[i].end())
        { 
            enRunGroups.push_back(std::vector<int>(it, it + maxGroupSize)); 
            it = it + maxGroupSize;
        }
        enRunGroups.push_back(std::vector<int>(it, emeasRunGroups[i].end()));
    }
    return enRunGroups.size();
}

std::map<int, std::pair<Float_t, Float_t>> Energy::AverageKchEnergy(bool isValForEachRun)
{
    std::map<int, std::pair<Float_t, Float_t>> averaged;
    Float_t tmpEnergy = 0;
    Float_t tmpEnergyErr = 0;
    for(auto &group : enRunGroups)
    {
        auto tmpPair = Eval(group, verbose);
        tmpEnergy = tmpPair.first;
        tmpEnergyErr = tmpPair.second;
        if(isValForEachRun)
        {
            for(auto run : group)
            { averaged[run] = std::make_pair(tmpEnergy, tmpEnergyErr); }
            continue;
        }
        else
        { averaged[group[0]] = std::make_pair(tmpEnergy, tmpEnergyErr); }
    }
    return averaged;
}

std::pair<double, double> Energy::DrawGraph(int graphNum)
{
    std::vector<double> groupNums;
    std::vector<double> groupNumsErr;
    std::vector<double> enVals;
    std::vector<double> enErrs;
    auto enValsAveraged = AverageKchEnergy();
    for(auto &group : enRunGroups)
    {
        groupNums.push_back((group.back() + group[0]) / 2);
        groupNumsErr.push_back(fabs(group.back() - group[0]) / 2);
    }

    std::vector<double> emeasGroupNums;
    std::vector<double> emeasGroupNumsErr;
    std::vector<double> emeasVals;
    std::vector<double> emeasErrs;
    for(auto &group : emeasRunGroups)
    {
        emeasGroupNums.push_back((group.back() + group[0]) / 2);
        emeasGroupNumsErr.push_back(fabs(group.back() - group[0]) / 2);
        emeasVals.push_back(runAvgEmeas[group[0]].first);
        emeasErrs.push_back(runAvgEmeas[group[0]].second);
    }
    // energyShift = emeasVals[0] - enValsAveraged[emeasRunGroups[0][0]].first;
    std::transform(enValsAveraged.begin(), enValsAveraged.end(), std::back_inserter(enVals),
                        [&](std::pair<int, std::pair<Float_t, Float_t>> val) { return double(val.second.first + energyShift); });
    std::transform(enValsAveraged.begin(), enValsAveraged.end(), std::back_inserter(enErrs),
                        [](std::pair<int, std::pair<Float_t, Float_t>> val) { return double(val.second.second); });

    TGraphErrors grEnVsRun(groupNums.size(), groupNums.data(), enVals.data(), groupNumsErr.data(), enErrs.data());
    TGraphErrors grEmeasVsRun(emeasGroupNums.size(), emeasGroupNums.data(), emeasVals.data(), emeasGroupNumsErr.data(), emeasErrs.data());

    TGraphErrors grComptonEnVsRun;
    grComptonEnVsRun.AddPoint(runs[0] - 10, comptonEnergyMean);
    grComptonEnVsRun.SetPointError(0, 0, comptonEnergyError);
    grComptonEnVsRun.AddPoint(runs.back() + 10, comptonEnergyMean);
    grComptonEnVsRun.SetPointError(1, 0, comptonEnergyError);

    grComptonEnVsRun.SetLineColor(kBlue);
    grComptonEnVsRun.SetMarkerColor(kBlue);
    
    grEmeasVsRun.GetXaxis()->SetTitle("Run");
    grEmeasVsRun.GetYaxis()->SetTitle("Energy, MeV");
    grEmeasVsRun.SetTitle("Red -- emeas, black -- Kch E, blue band -- compton mean");
    grEmeasVsRun.SetMarkerColor(kRed);
    grEmeasVsRun.SetLineColor(kRed);

    grComptonEnVsRun.SetFillColor(kBlue);
    grComptonEnVsRun.SetFillStyle(3005);
    
    if(fabs(energyShift) <= 1e-6)
    { grEmeasVsRun.GetYaxis()->SetRangeUser(comptonEnergyMean - 6, comptonEnergyMean + 2); }

    grEnVsRun.SetName("grKchEnergy");
    grComptonEnVsRun.SetName("grComptonMeanEnergy");
    grEmeasVsRun.SetName("grEmeas");
    grEmeasVsRun.DrawClone("AP");
    grComptonEnVsRun.DrawClone("L3 same");
    grEnVsRun.DrawClone("P same");

    auto r = grEnVsRun.Fit("pol0", "SQME");
    return std::pair<double, double>(r->Parameter(0), r->ParError(0));
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

std::map<int, Float_t> Energy::GetEnergyDiff()
{
    auto tmp = AverageKchEnergy(true);
    double mean = 0;
    double meanErr = 0;
    std:: map<int, Float_t> diff;
    for(auto& [run, energy] : tmp)
    { 
        mean += energy.first;
        if(energy.second > 1e-7)
        { meanErr += 1 / energy.second / energy.second; }
    }
    mean = mean / tmp.size();
    meanErr = sqrt(1 / meanErr);

    for(auto& [run, energy] : tmp)
    { diff[run] = energy.first - mean; }
    // { diff[run] = energy.first + energyShift - comptonEnergyMean; }

    return diff;
}