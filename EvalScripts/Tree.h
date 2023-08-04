#ifndef Tree_h
#define Tree_h

#include <memory>
#include "TFile.h"
#include "TTree.h"

struct track {
    double theta;
    double phi;
    int nhit;
};

class Tree 
{
private:
    std::unique_ptr<TTree> tree;
    int Nentries;
    int entry;
    
public:
    int runnum; 
    double emeas; 
    // Kaon energy calculated as invariant mass of pi+pi-.
    double etrue;

    struct {
        // Momentum ratio = P1/P2, where P1 is the momentum of pi+, P2 is the momentum of pi-.
        double Y;
        double ksdpsi;
        track piPos;
        track piNeg;
        track ks;
    } reco, gen;

    Tree(std::string filename);

    int GetEntry(int entryNum)
    { return tree->GetEntry(entryNum); }

    /* 
    * Gets next entry if the tree has one and returns true. 
    * Otherwise returns false and set next entry to entryNum = 0.
    */
    bool NextEntry();
    bool SetEntry(int entryNum);
    int GetCurrentEntry(); 
};

Tree::Tree(std::string filename) 
{
    tree = std::make_unique<TTree>(TFile::Open(filename.c_str())->Get<TTree>("ksTree"));
    Nentries = tree->GetEntries();
    entry = -1;

    tree->SetBranchAddress("emeas", &emeas);
    tree->SetBranchAddress("etrue", &etrue);
    tree->SetBranchAddress("runnum", &runnum);
    
    // reco
    tree->SetBranchAddress("Y", &reco.Y);
    tree->SetBranchAddress("ksdpsi", &reco.ksdpsi);
    tree->SetBranchAddress("kstheta", &reco.ks.theta);
    tree->SetBranchAddress("ksphi", &reco.ks.phi);

    tree->SetBranchAddress("nhitPos", &reco.piPos.nhit);
    tree->SetBranchAddress("piPhiPos", &reco.piPos.phi);
    tree->SetBranchAddress("piThetaPos", &reco.piPos.theta);

    tree->SetBranchAddress("nhitNeg", &reco.piNeg.nhit);
    tree->SetBranchAddress("piPhiNeg", &reco.piNeg.phi);
    tree->SetBranchAddress("piThetaNeg", &reco.piNeg.theta);

    // gen
    tree->SetBranchAddress("Y_gen", &gen.Y);
    tree->SetBranchAddress("ksdpsi_gen", &gen.ksdpsi);
    tree->SetBranchAddress("kstheta_gen", &gen.ks.theta);
    tree->SetBranchAddress("ksphi_gen", &gen.ks.phi);

    tree->SetBranchAddress("nhitPos_gen", &gen.piPos.nhit);
    tree->SetBranchAddress("piPhiPos_gen", &gen.piPos.phi);
    tree->SetBranchAddress("piThetaPos_gen", &gen.piPos.theta);

    tree->SetBranchAddress("nhitNeg_gen", &gen.piNeg.nhit);
    tree->SetBranchAddress("piPhiNeg_gen", &gen.piNeg.phi);
    tree->SetBranchAddress("piThetaNeg_gen", &gen.piNeg.theta);
}

bool Tree::NextEntry()
{ 
    entry++;
    if(entry < Nentries)
    {   
        GetEntry(entry);
        return true;
    }

    entry = -1; 
    return false; 
}

bool Tree::SetEntry(int entryNum)
{
    if(entryNum < Nentries)
    {
        entry = entryNum;
        return true;
    }
    return false;
}

int Tree::GetCurrentEntry() 
{ return entry; } 

#endif