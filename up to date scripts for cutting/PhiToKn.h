//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jul  2 19:27:21 2023 by ROOT version 6.28/04
// from TTree ksPrelim/Preliminary cut tree
// found on file: ./Prelim501.root
//////////////////////////////////////////////////////////

#ifndef PhiToKn_h
#define PhiToKn_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class PhiToKn {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         emeas0;
   Float_t         demeas0;
   Float_t         emeas;
   Float_t         demeas;
   Int_t           runnum;
   Int_t           nks;
   Int_t           tcharge[15];
   Int_t           ksvind[15][2];
   Float_t         kspiphi[15][2];
   Float_t         kspith[15][2];
   Float_t         kspipt[15][2];
   Float_t         ksminv[15];
   Float_t         kstlen[15];
   Float_t         ksalign[15];

   // List of branches
   TBranch        *b_emeas0;   //!
   TBranch        *b_demeas0;   //!
   TBranch        *b_emeas;   //!
   TBranch        *b_demeas;   //!
   TBranch        *b_runnum;   //!
   TBranch        *b_nks;   //!
   TBranch        *b_tcharge;   //!
   TBranch        *b_ksvind;   //!
   TBranch        *b_kspiphi;   //!
   TBranch        *b_kspith;   //!
   TBranch        *b_kspipt;   //!
   TBranch        *b_ksminv;   //!
   TBranch        *b_kstlen;   //!
   TBranch        *b_ksalign;   //!

   PhiToKn(TTree *tree=0);
   virtual ~PhiToKn();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(std::string output_fname);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PhiToKn_cxx
PhiToKn::PhiToKn(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("./Prelim501.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("./Prelim501.root");
      }
      f->GetObject("ksPrelim",tree);

   }
   Init(tree);
}

PhiToKn::~PhiToKn()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PhiToKn::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PhiToKn::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void PhiToKn::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("emeas0", &emeas0, &b_emeas0);
   fChain->SetBranchAddress("demeas0", &demeas0, &b_demeas0);
   fChain->SetBranchAddress("emeas", &emeas, &b_emeas);
   fChain->SetBranchAddress("demeas", &demeas, &b_demeas);

   fChain->SetBranchAddress("runnum", &runnum, &b_runnum);
   fChain->SetBranchAddress("nks", &nks, &b_nks);
   fChain->SetBranchAddress("tcharge", tcharge, &b_tcharge);
   fChain->SetBranchAddress("ksvind", ksvind, &b_ksvind);
   fChain->SetBranchAddress("kspiphi", kspiphi, &b_kspiphi);
   fChain->SetBranchAddress("kspith", kspith, &b_kspith);
   fChain->SetBranchAddress("kspipt", kspipt, &b_kspipt);
   fChain->SetBranchAddress("ksminv", ksminv, &b_ksminv);
   fChain->SetBranchAddress("kstlen", kstlen, &b_kstlen);
   fChain->SetBranchAddress("ksalign", ksalign, &b_ksalign);
   Notify();
}

Bool_t PhiToKn::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PhiToKn::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PhiToKn::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PhiToKn_cxx
