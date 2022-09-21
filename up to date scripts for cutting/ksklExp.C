#define ksklExp_cxx
#include "ksklExp.h"
#include <TH2.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TVector3.h>

void ksklExp::Loop(std::string histFileName)
{
//   In a ROOT session, you can do:
//      root> .L ksklExp.C
//      root> ksklExp t
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
//by  b_branchname->GetEntry(ientry); //read only this branch
    if (fChain == 0)
        return;

    TFile *top = new TFile(histFileName.c_str(), "recreate");
    auto tNew = new TTree("ksTree", "Cutted tr_ph (Ks mass meas important data)");

    Float_t dpsi;
    Float_t Y;
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

    TVector3 ks;
    TVector3 kl;
    // dPhi and dTheta - phi and theta angle between Ks and Kl respectively
    double dPhi = 0;
    double dTheta = 0;
    int tmpCounter = 0;
    std::vector<Int_t> ksCand = {};

    auto hist = new TH2D("hist", "", 1000, 0, 600, 1000, 0, 600);
    auto hSumMomTr = new TH1D("hSumMomTr", "Sum of moms", 4000, 0, 1200);
    auto histKlCands = new TH1D("histKlCands", "number of Kl candidates for one Ks candidate", 5, 0, 5);
    auto histKsCands = new TH1D("histKsCands", "number of Ks candidates", 5, 0, 5);
    // Theta = kl.Theta()
    auto hdPhiTheta = new TH2D("hdPhiTheta", "", 600, 0, TMath::Pi(), 600, -TMath::Pi(), TMath::Pi());
    auto hdThetadPhi = new TH2D("hdThetadPhi", "", 1200, 0, 2*TMath::Pi(), 1200, -TMath::Pi(), TMath::Pi());
    auto hClEdPhi = new TH2D("hClEdPhi", "", 600, 0, 2 * TMath::Pi(), 600, 0, 600);
    auto hPsiUncutted = new TH1D("hPsiUncutted", "", 628, 0, 3.15);
    auto hPsiCutted = new TH1D("hPsiCutted", "", 628, 0, 3.15);
    // Phi angle between pions

    auto hKsKlPhi = new TH1D("hKsKlPhi", "", 1000, -1.5, 1.5);

    auto hKsKlTheta = new TH1D("hKsKlTheta", "", 1000, -1.5, 1.5);

    auto hHit = new TH1D("hHit", "nhits", 40, 0, 40);
    auto hTrackTheta = new TH1D("hTrackTheta", "hTrackTheta", 250, -TMath::Pi() / 2, TMath::Pi() / 2);

    hdThetadPhi->GetXaxis()->SetTitle("#Delta#phi, rad");
    hdThetadPhi->GetYaxis()->SetTitle("#Delta#theta, rad");

    hdPhiTheta->GetXaxis()->SetTitle("#theta of Kl, rad");
    hdPhiTheta->GetYaxis()->SetTitle("#Delta#phi, rad");

    hPsiUncutted->GetXaxis()->SetTitle("Angle between motion vectors Ks and Kl, rad");

    hClEdPhi->GetXaxis()->SetTitle("#Delta#phi, rad");
    hClEdPhi->GetYaxis()->SetTitle("Cluster Energy, MeV");

    TVector3 piPos;
    TVector3 piNeg;
    TVector3 field(0., 0., 1.);

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
                // kspipt[k][0] > 120 && kspipt[k][1] > 120 && 
                // kspipt[k][0] < 350 && kspipt[k][1] < 350 &&
		(kspipt[k][0]+kspipt[k][1])<500 &&
                tcharge[ksvind[k][0]] * tcharge[ksvind[k][1]] < 0 && kstype[k] == 0) // Added kstype[k] == 0.
                {
                    ks.SetMagThetaPhi(1, ksth[k], ksphi[k]);

                    for(int j = 0; j < nph; j++)
                    {
                        // if(phen0[j] < 40)
                        // { continue; }
                        // phth0 and phphi0 are theta and phi respectively angles of KL candidate.
                        kl.SetMagThetaPhi(phrho[j], phth0[j], phphi0[j]);
                        // Shift to the center of phi decay. 
                        kl.SetX(kl.X() - xbeam);
                        kl.SetY(kl.Y() - ybeam);
                        kl.SetZ(kl.Z() - ksz0[k]);

                        hPsiUncutted->Fill(ks.Angle(kl));

                        dPhi = ks.DeltaPhi(kl);
                        dTheta = ks.Theta() + kl.Theta() - TMath::Pi();
                        // hdThetadPhi->Fill(dPhi, dTheta);
                        // phen0 - cluster energy. So this is the energy deposition of Kl candidate.
                        // hClEdPhi->Fill(dPhi, phen0[j]); 
                        hdPhiTheta->Fill(kl.Theta(), dPhi);

                        for(int i = 0; i < nsim; i++)
                        {
                            if(simtype[i] == 211 && simorig[i] == 310)
                            { piPos.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); }
                            if(simtype[i] == -211 && simorig[i] == 310)
                            { piNeg.SetMagThetaPhi(simmom[i], simtheta[i], simphi[i]); }
                        }
                        

                        if(dPhi < 0)
                        { 
                            hdThetadPhi->Fill(dPhi + 2*TMath::Pi(), dTheta); 
                            hClEdPhi->Fill(dPhi + 2*TMath::Pi(), phen0[j]);
                        }
                        else
                        { 
                            hdThetadPhi->Fill(dPhi, dTheta); 
                            hClEdPhi->Fill(dPhi, phen0[j]);
                        }

                        if((dPhi < -TMath::Pi() + 1 || dPhi > TMath::Pi() - 1) && fabs(dTheta) < 0.1 && phen0[j] > 40)
                        { 
                            tmpCounter++; 
                            hPsiCutted->Fill(ks.Angle(kl));
                            // hdThetadPhi->Fill(dPhi, dTheta);
                            if(dPhi < 0)
                            { hKsKlPhi->Fill(dPhi + TMath::Pi()); }
                            else
                            { hKsKlPhi->Fill(dPhi - TMath::Pi()); }
                            hKsKlTheta->Fill(dTheta);
                        } 
                    }        
        
                    histKlCands->Fill(tmpCounter);
                    if(tmpCounter > 0) // '|| 1' if there is no Kl cut
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
                    
                    piPos.SetMagThetaPhi(kspipt[0][0], kspith[0][0], kspiphi[0][0]);
                    piNeg.SetMagThetaPhi(kspipt[0][1], kspith[0][1], kspiphi[0][1]);
                }
                else
                {
                    Y = kspipt[0][1] / kspipt[0][0];
                    hist->Fill(kspipt[0][1], kspipt[0][0]);
                    p1 = kspipt[0][1]; p2 = kspipt[0][0];
                    piPos.SetMagThetaPhi(kspipt[0][1], kspith[0][1], kspiphi[0][1]);
                    piNeg.SetMagThetaPhi(kspipt[0][0], kspith[0][0], kspiphi[0][0]);
                }
                hSumMomTr->Fill(kspipt[0][1] + kspipt[0][0]);

                dpsi = ksdpsi[0];
                tNew->Fill();
                if(piPos.Cross(field).XYvector().DeltaPhi(piNeg.XYvector()) < TMath::Pi() / 2)
                { 
                    // cowboy
                    // dpsi = piNeg.Angle(piPos);
                    // tNew->Fill(); 
                }

                if(piPos.Cross(field).XYvector().DeltaPhi(piNeg.XYvector()) > TMath::Pi() / 2)
                {  
                    // sailor
                }
            }
            ksCand.clear();
            ksCand.shrink_to_fit();
        }
        NgoodTr = 0;
    }

    // Drawing hits with cuts
    fChain->Draw("tnhit>>hHit", "", "goff");

    auto hitCutLine1 = new TLine(10, 0, 10, 10000);
    hitCutLine1->SetLineColor(kBlue);
    hitCutLine1->SetLineWidth(4);
    auto hitCutLine2 = new TLine(30, 0, 30, 10000);
    hitCutLine2->SetLineColor(kBlue);
    hitCutLine2->SetLineWidth(4);

    hHit->GetXaxis()->SetTitle("Number of hits");
    // hHit->Draw();
    // hitCutLine1->Draw("same");
    // hitCutLine2->Draw("same");
    // Drawing hits with cuts 

    // ****************************************************

    // Drawing track theta with cuts' lines 
    fChain->Draw("tth - TMath::Pi() / 2>>hTrackTheta", "", "goff");

    auto thetaCutLine1 = new TLine(-0.7, 0, -0.7, 1600);
    thetaCutLine1->SetLineColor(kBlue);
    thetaCutLine1->SetLineWidth(4);
    auto thetaCutLine2 = new TLine(0.7, 0, 0.7, 1600);
    thetaCutLine2->SetLineColor(kBlue);
    thetaCutLine2->SetLineWidth(4);

    hTrackTheta->GetXaxis()->SetTitle("#theta - #frac{#pi}{2}, rad");
    // hTrackTheta->Draw();
    // thetaCutLine1->Draw("same");
    // thetaCutLine2->Draw("same");
    // Drawing track theta with cuts' lines 

    // ****************************************************

    // Drawing KsKl dThetadPhi with cuts' lines
    auto dThetadPhiCutLine1 = new TLine(TMath::Pi() - 1, -0.3, TMath::Pi() + 1, -0.3);
    auto dThetadPhiCutLine2 = new TLine(TMath::Pi() - 1, -0.3, TMath::Pi() - 1, 0.3);
    auto dThetadPhiCutLine3 = new TLine(TMath::Pi() + 1, -0.3, TMath::Pi() + 1, 0.3);
    auto dThetadPhiCutLine4 = new TLine(TMath::Pi() - 1, 0.3, TMath::Pi() + 1, 0.3);
    dThetadPhiCutLine1->SetLineWidth(4);
    dThetadPhiCutLine2->SetLineWidth(4);
    dThetadPhiCutLine3->SetLineWidth(4);
    dThetadPhiCutLine4->SetLineWidth(4);
    // hdThetadPhi->Draw("col");
    // dThetadPhiCutLine1->Draw("same");
    // dThetadPhiCutLine2->Draw("same");
    // dThetadPhiCutLine3->Draw("same");
    // dThetadPhiCutLine4->Draw("same");
    // Drawing KsKl dThetadPhi with cuts' lines

    // ****************************************************

    // Drawing KsKl Cluster Energy vs dPhi with cuts' lines
    auto hClEdPhiCutLine1 = new TLine(TMath::Pi() - 1, 40, TMath::Pi() + 1, 40);
    auto hClEdPhiCutLine2 = new TLine(TMath::Pi() - 1, 40, TMath::Pi() - 1, 600);
    auto hClEdPhiCutLine3 = new TLine(TMath::Pi() + 1, 40, TMath::Pi() + 1, 600);
    hClEdPhiCutLine1->SetLineWidth(4);
    hClEdPhiCutLine2->SetLineWidth(4);
    hClEdPhiCutLine3->SetLineWidth(4);
    // hClEdPhi->Draw("col");
    // hClEdPhiCutLine1->Draw("same");
    // hClEdPhiCutLine2->Draw("same");
    // hClEdPhiCutLine3->Draw("same");
    // Drawing KsKl Cluster Energy vs dPhi with cuts' lines

    std::cout << "nentries " << nentries << std::endl;
    std::cout << "NgoodTrS = " << NgoodTrS << std::endl;

    hist->GetYaxis()->SetTitle("P_{#pi^{+}} [MeV/c]");
    hist->GetXaxis()->SetTitle("P_{#pi^{-}} [MeV/c]");
    hist->Draw("COL");
    // hdThetadPhi->Draw("COL");
    // hPsiUncutted->Draw();
    // hPsiCutted->SetLineColor(kGreen);
    // hPsiCutted->Draw("same");


    // hKsKlTheta->Draw();

    top->Write();
    top->Save();
}
