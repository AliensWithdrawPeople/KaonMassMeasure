#ifndef Tree_hpp
#define Tree_hpp

#include <memory>
#include <optional>
#include <stdexcept>

#include "TFile.h"
#include "TTree.h"


struct track {
    float theta;
    float phi;
    int nhit;
};

class Tree 
{
public:
    /// @brief iterator for the Tree class. Mostly usable in range-based for loop. 
    class iterator {
    private:
        int entryNum = 0;
        Tree* tr;

    public:
        iterator(int entryNum, Tree* tr): entryNum(entryNum), tr(tr) 
        {}

        iterator operator++() 
        { 
            ++entryNum; 
            tr->GetEntry(entryNum);
            return *this; 
        }

        bool operator!=(const iterator & other) const 
        { return entryNum != other.entryNum;  }

        const int& operator*() const 
        { return entryNum; }
    };

private:
    std::unique_ptr<TTree> tree;

public:
    int Nentries;

    int runnum; 
    float emeas; 
    // Kaon energy calculated as invariant mass of pi+pi-.
    float etrue;

    struct Data {
        // Momentum ratio = P1/P2, where P1 is the momentum of pi+, P2 is the momentum of pi-.
        float Y;
        float ks_len;
        float ksdpsi;
        track piPos;
        track piNeg;
        track ks;
    } reco, gen;

    Data* data = &reco;

    Tree(std::string filename, bool isExp = false);

    iterator begin() { return iterator(0, this); }
    iterator end() { return iterator(Nentries, this); }

    int GetEntry(int entryNum)
    { return tree->GetEntry(entryNum); }
    
    void SetReco(bool yes) 
    { data = yes? &reco : &gen; }
};

Tree::Tree(std::string filename, bool isExp): tree(std::unique_ptr<TTree>(TFile::Open(filename.c_str())->Get<TTree>("ksTree"))), Nentries{int(tree->GetEntries())}
{
    tree->SetBranchAddress("emeas", &emeas);
    tree->SetBranchAddress("etrue", &etrue);
    tree->SetBranchAddress("runnum", &runnum);
    
    // reco
    tree->SetBranchAddress("Y", &reco.Y);
    tree->SetBranchAddress("ks_len", &reco.ks_len);
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
    if(!isExp)
    {
        tree->SetBranchAddress("Y_gen", &gen.Y);
        tree->SetBranchAddress("ksdpsi_gen", &gen.ksdpsi);
        tree->SetBranchAddress("kstheta_gen", &gen.ks.theta);
        tree->SetBranchAddress("ksphi_gen", &gen.ks.phi);

        // tree->SetBranchAddress("nhitPos", &gen.piPos.nhit);
        gen.piPos.nhit = reco.piPos.nhit;
        tree->SetBranchAddress("piPhiPos_gen", &gen.piPos.phi);
        tree->SetBranchAddress("piThetaPos_gen", &gen.piPos.theta);

        // tree->SetBranchAddress("nhitNeg", &gen.piNeg.nhit);
        gen.piNeg.nhit = reco.piNeg.nhit;
        tree->SetBranchAddress("piPhiNeg_gen", &gen.piNeg.phi);
        tree->SetBranchAddress("piThetaNeg_gen", &gen.piNeg.theta);
    }
}

#endif