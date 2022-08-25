#define kskl2bGen_cxx
#include "kskl2bGen.h"
#include <TH2.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TVector3.h>

void kskl2bGen::Loop(std::string histFileName)
{
    //   In a ROOT session, you can do:
    //      root> .L ksklCut_mcgpj.C
    //      root> ksklCut_mcgpj t
    //      root> t.GetEntry(12); // Fill t data members with entry number 12
    //      root> t.Show();       // Show values of entry 12
    //      root> t.Show(16);     // Read and show values of entry 16
    //      root> t.Loop();       // Loop on all entries
    //

    //     This is the loop skeleton where:
    //    jentry is the global entry number in the chain
    //    ientry is the entry number in the current Tree
    //  Note that the argument to GetEntry must be:
    //    jentry for TChain::GetEntry
    //    ientry for TTree::GetEntry and TBranch::GetEntry
    //
    //       To read only selected branches, Insert statements like:
    // METHOD1:
    //    fChain->SetBranchStatus("*",0);  // disable all branches
    //    fChain->SetBranchStatus("branchname",1);  // activate branchname
    // METHOD2: replace line
    //    fChain->GetEntry(jentry);       //read all branches
    // by  b_branchname->GetEntry(ientry); //read only this branch
    if (fChain == 0)
        return;

    TFile *top = new TFile(histFileName.c_str(), "recreate");
    auto tNew = new TTree("ksTree", "Cutted tr_ph (Ks mass meas important data)");

    Float_t dpsi;
    Float_t Y;
    Float_t eff;
    Float_t p1; Float_t p2;
    tNew->Branch("emeas", &emeas, "emeas/F");
    tNew->Branch("demeas", &demeas, "demeas/F");
    tNew->Branch("runnum", &runnum, "runnum/I");
    tNew->Branch("ksdpsi", &dpsi, "dpsi/F");
    tNew->Branch("Y", &Y, "Y/F");

    Double_t halfPi = TMath::Pi() / 2;
    int NgoodTr = 0;
    int NgoodTrS = 0;
    double cutChi2r = 15.;
    double cutChi2z = 10.;
    int cutNhitMin = 10;
    int cutNhitMax = 30;
    double cutRmin = 0.1;
    double cutRmax = 6.0;
    double cutZtrack = 12.;
    double cutPtot = 40;

    double cosKlKs = 0;
    TVector3 ks;
    TVector3 kl;
    // dPhi and dTheta - phi and theta angle between Ks and Kl respectively
    double dPhi = 0;
    double dTheta = 0;
    int tmpCounter = 0;
    std::vector<Int_t> ksCand = {};

    auto hist = new TH2D("hist", "", 1000, 0, 600, 1000, 0, 600);
    auto histKlCands = new TH1D("KlCands", "number of Kl candidates for one Ks candidate", 5, 0, 5);
    auto histKsCands = new TH1D("KsCands", "number of Ks candidates", 5, 0, 5);
    // Theta = kl.Theta()
    auto hdPhiTheta = new TH2D("hdPhiTheta", "", 600, 0, TMath::Pi(), 600, -TMath::Pi(), TMath::Pi());
    auto hdThetadPhi = new TH2D("hdThetadPhi", "", 600, -TMath::Pi(), TMath::Pi(), 600, -TMath::Pi(), TMath::Pi());
    auto hClEdPhi = new TH2D("hClEdPhi", "", 600, -TMath::Pi(), TMath::Pi(), 600, 0, 600);
    auto hPsi = new TH1D("hPsi", "", 628, 0, 6.28);
    auto hPhi = new TH1D("hPhi", "", 1000, -0.2, 0.2);
    auto hDeltaMom = new TH1D("hDeltaMom", "Lorentz delta mom", 1000, -0.2, 0.2);
    auto hDeltaMomVsDeltaPhi1 = new TH2D("hDeltaMomVsDeltaPhi1", "delta mom vs delta phi", 1000, -0.1, 0.1, 1000, -10, 10);
    auto hDeltaMomVsDeltaPhi2 = new TH2D("hDeltaMomVsDeltaPhi2", "delta mom vs delta phi Cowboy", 1000, -0.1, 0.1, 1000, -10, 10);
    auto hDeltaMomVsDeltaPhi3 = new TH2D("hDeltaMomVsDeltaPhi3", "delta mom vs delta phi Sailor", 1000, -0.1, 0.1, 1000, -10, 10);

    auto hE = new TH1D("hE", "Corrected energy", 50000, 480, 515);

    hdThetadPhi->GetXaxis()->SetTitle("#Delta#phi, rad");
    hdThetadPhi->GetYaxis()->SetTitle("#Delta#theta, rad");

    hdPhiTheta->GetXaxis()->SetTitle("#theta of Kl, rad");
    hdPhiTheta->GetYaxis()->SetTitle("#Delta#phi, rad");

    hPsi->GetXaxis()->SetTitle("Angle between motion vectors Ks and Kl, rad");

    hClEdPhi->GetXaxis()->SetTitle("#Delta#phi, rad");
    hClEdPhi->GetYaxis()->SetTitle("Cluster Energy,  MeV");

    TVector3 piPos;
    TVector3 piNeg;
    TVector3 piPosRec;
    TVector3 piNegRec;
    TVector3 field(0., 0., 1.);
    int counter1 = 0;
    int counter2 = 0;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0)
            break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        for (int i = 0; i < nt; i++)
        {
            if (tptot[i] > cutPtot && fabs(trho[i]) < cutRmax && fabs(tz[i]) < cutZtrack &&
                tchi2r[i] < cutChi2r && tchi2z[i] < cutChi2z && tnhit[i] > cutNhitMin && tnhit[i] < cutNhitMax)
            { NgoodTr++; }
        }

        if (NgoodTr == 2)
        { NgoodTrS++; }
        
    
        if (NgoodTr == 2 && is_coll != 1 )
        {
            for(int k = 0; k < nks; k++)
            {
                if(ksalign[k] > 0.85 && (tdedx[ksvind[k][0]] + tdedx[ksvind[k][1]]) / 2 < 5000 &&
                abs(kspith[k][0] - TMath::Pi() / 2) <= 0.7 && 
                abs(kspith[k][1] - TMath::Pi() / 2) <= 0.7 &&
                //20 - half of the linear size of Drift Chamber
                //(20 - ksz0[0]) * fabs(TMath::Tan(kspith[0][0])) > 15 && (20 - ksz0[0]) * fabs(TMath::Tan(kspith[0][1])) > 15 &&
                //kspipt[k][0] > 120 && kspipt[k][1] > 120 && 
                //kspipt[k][0] < 350 && kspipt[k][1] < 350 &&
                tcharge[ksvind[k][0]] * tcharge[ksvind[k][1]] < 0 && kstype[k] == 0) // Added kstype[k] == 0.
                {
                    ks.SetMagThetaPhi(1, ksth[k], ksphi[k]);

                    for(int j = 0; j < nph; j++)
                    {
                        // phth0 and phphi0 are theta and phi respectively angles of KL candidate.
                        kl.SetMagThetaPhi(phrho[j], phth0[j], phphi0[j]);
                        // Shift to the center of phi decay. 
                        kl.SetX(kl.X() - xbeam);
                        kl.SetY(kl.Y() - ybeam);
                        kl.SetZ(kl.Z() - ksz0[k]);

                        hPsi->Fill(ks.Angle(kl));

                        dPhi = ks.DeltaPhi(kl);
                        dTheta = ks.Theta() + kl.Theta() - TMath::Pi();
                        hdThetadPhi->Fill(dPhi, dTheta);
                        // phen0 - cluster energy. So this is the energy deposition of Kl candidate.
                        hClEdPhi->Fill(dPhi, phen0[j]); 
                        hdPhiTheta->Fill(kl.Theta(), dPhi);

                        if((dPhi < -TMath::Pi() + 1 || dPhi > TMath::Pi() - 1) && fabs(dTheta) < 1 && phen0[j] > 40)
                        { tmpCounter++; } 
                    }        
        
                    histKlCands->Fill(tmpCounter);
                    if(tmpCounter > 0 || 1) // '|| 1' if there is no Kl cut
                    { ksCand.push_back(k); }
                    tmpCounter = 0;
                }
            }

            histKsCands->Fill(ksCand.size());
            if(ksCand.size() > 0)
            {
                if (tcharge[ksvind[0][0]] > 0)
                {
                    Y = kspipt[0][0] / kspipt[0][1];
                    hist->Fill(kspipt[0][0], kspipt[0][1]);
                    p1 = kspipt[0][0]; p2 = kspipt[0][1];
                    piPosRec.SetMagThetaPhi(kspipt[0][0], kspith[0][0], kspiphi[0][0]);
                    piNegRec.SetMagThetaPhi(kspipt[0][1], kspith[0][1], kspiphi[0][1]);
                }
                else
                {
                    Y = kspipt[0][1] / kspipt[0][0];
                    hist->Fill(kspipt[0][1], kspipt[0][0]);
                    p1 = kspipt[0][1]; p2 = kspipt[0][0];
                    piPosRec.SetMagThetaPhi(kspipt[0][1], kspith[0][1], kspiphi[0][1]);
                    piNegRec.SetMagThetaPhi(kspipt[0][0], kspith[0][0], kspiphi[0][0]);
                }

                for(int i = 0; i < nsim; i++)
                {
                    if(simtype[i] == 211 && simorig[i] == 310)
                    { 
                        piPos.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); 
                        counter1++;
                    }

                    if(simtype[i] == -211 && simorig[i] == 310)
                    { 
                        piNeg.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); 
                        counter2++;
                    }

                }
                dpsi = ksdpsi[0];
                tNew->Fill();
                hE->Fill(sqrt(139.57 * 139.57 + piNeg.Mag2()) + sqrt(139.57 * 139.57 + piPos.Mag2()));
                // emeas = sqrt(139.57 * 139.57 + piNeg.Mag2()) + sqrt(139.57 * 139.57 + piPos.Mag2());  
                if(piPos.Cross(field).XYvector().DeltaPhi(piNeg.XYvector()) < TMath::Pi() / 2)
                { 
                    // hE->Fill(sqrt(139.57 * 139.57 + piNeg.Mag2()) + sqrt(139.57 * 139.57 + piPos.Mag2()));
                    // dpsi = piNeg.Angle(piPos);
                    // tNew->Fill(); 
                    hPhi->Fill(piPos.Phi() - piPosRec.Phi()); 
                    hPhi->Fill(piNeg.Phi() - piNegRec.Phi()); 

                    hDeltaMomVsDeltaPhi2->Fill(piPos.XYvector().DeltaPhi(piNeg.XYvector()) - piPosRec.XYvector().DeltaPhi(piNegRec.XYvector()), 
                                            (piPos.Mag() - piNeg.Mag()) - (piPosRec.Mag() - piNegRec.Mag()));
                }

                if(piPos.Cross(field).XYvector().DeltaPhi(piNeg.XYvector()) > TMath::Pi() / 2)
                {  
                    hDeltaMomVsDeltaPhi3->Fill(piPos.XYvector().DeltaPhi(piNeg.XYvector()) - piPosRec.XYvector().DeltaPhi(piNegRec.XYvector()), 
                                            (piPos.Mag() - piNeg.Mag()) - (piPosRec.Mag() - piNegRec.Mag())); 
                }

                // hDeltaMomVsDeltaPhi2->Fill(piPos.Phi() - piPosRec.Phi(), piPos.Mag() - piPosRec.Mag());
                // hDeltaMomVsDeltaPhi3->Fill(piNeg.Phi() - piNegRec.Phi(), piNeg.Mag() - piNegRec.Mag());
                hDeltaMomVsDeltaPhi1->Fill(piPos.XYvector().DeltaPhi(piNeg.XYvector()) - piPosRec.XYvector().DeltaPhi(piNegRec.XYvector()), 
                                            (piPos.Mag() - piNeg.Mag()) - (piPosRec.Mag() - piNegRec.Mag()));

                hDeltaMom->Fill( (piPos.Mag() - piPosRec.Mag()) / piPos.Mag());
            }
            ksCand.clear();
            ksCand.shrink_to_fit();
        }

        NgoodTr = 0;
    }

    eff = (double) tNew->GetEntriesFast() / nentries;
    std::cout << "nentries " << nentries << std::endl;
    std::cout << "e_mc = " <<  eff << std::endl;
    std::cout << "NgoodTrS = " << NgoodTrS << std::endl;

    std::cout << "Energy point: " << emeas << "; Mean Energy = " << hE->GetMean() << " +/- " << hE->GetMeanError() << std::endl;

    hist->GetYaxis()->SetTitle("P_{#pi^{+}} [MeV/c]");
    hist->GetXaxis()->SetTitle("P_{#pi^{-}} [MeV/c]");
    hist->Draw("COL");
    hdThetadPhi->Draw("COL");

    // hDeltaMomVsDeltaPhi2->GetYaxis()->SetTitle("#DeltaP_{#pi^{+}} [MeV/c]");
    // hDeltaMomVsDeltaPhi1->Draw("COL");
    // auto c = new TCanvas("deltaP vs deltaPhi", "#pi^{-}", 200, 10, 600, 400);
    // c->Divide(2, 1);
    // c->cd(1);
    // hDeltaMomVsDeltaPhi2->Draw("COL");
    // c->cd(2);
    // hDeltaMomVsDeltaPhi3->Draw("COL");


    top->Write();
    top->Save();
}