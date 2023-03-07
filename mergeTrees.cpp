#include <iostream>
#include <TTree.h>
#include <TList.h>
#include <TFile.h>

int mergeTrees()
{
    TList *list = new TList();
    TFile *file1 = TFile::Open("tr_ph/mcgpj/tr_ph v9/EnergySmearing/510.5/tr_ph_run000006.root");
    list->Add((TTree *)file1->Get("tr_ph"));

    TFile *file2 = TFile::Open("tr_ph/mcgpj/tr_ph v9/EnergySmearing/510.5/tr_ph_run000007.root");
    list->Add((TTree *)file2->Get("tr_ph"));

    TFile *file3 = TFile::Open("tr_ph/mcgpj/tr_ph v9/EnergySmearing/510.5/tr_ph_run000008.root");
    list->Add((TTree *)file3->Get("tr_ph"));

    TFile *file4 = TFile::Open("tr_ph/mcgpj/tr_ph v9/EnergySmearing/510.5/tr_ph_run000009.root");
    list->Add((TTree *)file4->Get("tr_ph"));

    TFile *file5 = TFile::Open("tr_ph/mcgpj/tr_ph v9/EnergySmearing/510.5/tr_ph_run000010.root");
    list->Add((TTree *)file5->Get("tr_ph"));

    TFile *file6 = TFile::Open("tr_ph/mcgpj/tr_ph v9/EnergySmearing/510.5/tr_ph_run000001.root");
    list->Add((TTree *)file5->Get("tr_ph"));

    TFile *file7 = TFile::Open("tr_ph/mcgpj/tr_ph v9/EnergySmearing/510.5/tr_ph_run000002.root");
    list->Add((TTree *)file5->Get("tr_ph"));

    TFile *file8 = TFile::Open("tr_ph/mcgpj/tr_ph v9/EnergySmearing/510.5/tr_ph_run000003.root");
    list->Add((TTree *)file5->Get("tr_ph"));

    TFile *file9 = TFile::Open("tr_ph/mcgpj/tr_ph v9/EnergySmearing/510.5/tr_ph_run000004.root");
    list->Add((TTree *)file5->Get("tr_ph"));

    TFile *top = new TFile("tr_ph/mcgpj/tr_ph v9/EnergySmearing/MCGPJ_kskl510.5_Merged_9points.root", "recreate");
    TTree *newtree = TTree::MergeTrees(list);
    newtree->SetName("tr_ph_merged");
    top->Write();
    top->Save();
    return 0;
}