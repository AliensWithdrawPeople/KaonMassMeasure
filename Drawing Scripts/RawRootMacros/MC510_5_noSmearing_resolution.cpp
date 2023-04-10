#ifdef __CLING__
#pragma cling optimize(0)
#endif

#include "TCanvas.h"
#include "TH1D.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TF1.h"
#include "TStyle.h"

void MC510_5_noSmearing_resolution()
{
//=========Macro generated from canvas: MlnY/Mass(lnY)
//=========  (Mon Apr 10 15:13:35 2023) by ROOT version 6.26/06
   TCanvas *MlnY = new TCanvas("MlnY", "Mass(lnY)",200,10,600,400);
   gStyle->SetOptFit(1);
   MlnY->Range(2.449402,-12.12152,2.787716,63.63797);
   MlnY->SetFillColor(0);
   MlnY->SetBorderMode(0);
   MlnY->SetBorderSize(2);
   MlnY->SetGridx();
   MlnY->SetGridy();
   MlnY->SetTickx(1);
   MlnY->SetTicky(1);
   MlnY->SetLeftMargin(0.12);
   MlnY->SetRightMargin(0.05);
   MlnY->SetTopMargin(0.05);
   MlnY->SetBottomMargin(0.16);
   MlnY->SetFrameBorderMode(0);
   MlnY->SetFrameBorderMode(0);
   
   TH1D *py6__2 = new TH1D("py6__2","Psi6(lnY)",250,2.4,3.3);
   py6__2->SetBinContent(39,2);
   py6__2->SetBinContent(43,1);
   py6__2->SetBinContent(45,2);
   py6__2->SetBinContent(46,1);
   py6__2->SetBinContent(47,1);
   py6__2->SetBinContent(49,4);
   py6__2->SetBinContent(50,3);
   py6__2->SetBinContent(51,5);
   py6__2->SetBinContent(52,10);
   py6__2->SetBinContent(53,15);
   py6__2->SetBinContent(54,20);
   py6__2->SetBinContent(55,20);
   py6__2->SetBinContent(56,33);
   py6__2->SetBinContent(57,26);
   py6__2->SetBinContent(58,39);
   py6__2->SetBinContent(59,51);
   py6__2->SetBinContent(60,46);
   py6__2->SetBinContent(61,42);
   py6__2->SetBinContent(62,57);
   py6__2->SetBinContent(63,45);
   py6__2->SetBinContent(64,46);
   py6__2->SetBinContent(65,39);
   py6__2->SetBinContent(66,29);
   py6__2->SetBinContent(67,23);
   py6__2->SetBinContent(68,17);
   py6__2->SetBinContent(69,13);
   py6__2->SetBinContent(70,14);
   py6__2->SetBinContent(71,6);
   py6__2->SetBinContent(72,7);
   py6__2->SetBinContent(73,7);
   py6__2->SetBinContent(74,3);
   py6__2->SetBinContent(75,4);
   py6__2->SetBinContent(76,4);
   py6__2->SetBinContent(77,1);
   py6__2->SetBinContent(78,1);
   py6__2->SetBinContent(79,1);
   py6__2->SetBinContent(82,1);
   py6__2->SetBinContent(83,1);
   py6__2->SetBinContent(90,1);
   py6__2->SetBinContent(95,1);
   py6__2->SetEntries(642);
   
   TF1 *gaus3 = new TF1("gaus","gaus",2.581322,2.656672, TF1::EAddToList::kNo);
   gaus3->SetFillColor(19);
   gaus3->SetFillStyle(0);
   gaus3->SetMarkerStyle(20);
   gaus3->SetLineWidth(3);
   gaus3->SetChisquare(10.48726);
   gaus3->SetNDF(18);
   gaus3->GetXaxis()->SetNdivisions(515);
   gaus3->GetXaxis()->SetLabelFont(42);
   gaus3->GetXaxis()->SetLabelOffset(0.0125);
   gaus3->GetXaxis()->SetLabelSize(0.04);
   gaus3->GetXaxis()->SetTitleSize(0.05);
   gaus3->GetXaxis()->SetTitleOffset(1);
   gaus3->GetXaxis()->SetTitleFont(42);
   gaus3->GetYaxis()->SetNdivisions(515);
   gaus3->GetYaxis()->SetLabelFont(42);
   gaus3->GetYaxis()->SetLabelOffset(0.0125);
   gaus3->GetYaxis()->SetLabelSize(0.04);
   gaus3->GetYaxis()->SetTitleSize(0.05);
   gaus3->GetYaxis()->SetTitleOffset(1);
   gaus3->GetYaxis()->SetTitleFont(42);
   gaus3->SetParameter(0,48.71202);
   gaus3->SetParError(0,2.637446);
   gaus3->SetParLimits(0,0,0);
   gaus3->SetParameter(1,2.618461);
   gaus3->SetParError(1,0.0008191932);
   gaus3->SetParLimits(1,0,0);
   gaus3->SetParameter(2,0.01787566);
   gaus3->SetParError(2,0.000768446);
   gaus3->SetParLimits(2,0,0.1617282);
   gaus3->SetParent(py6__2);
   py6__2->GetListOfFunctions()->Add(gaus3);
   
   TPaveStats *ptstats = new TPaveStats(0.62,0.595,0.98,0.995,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(2);
   ptstats->SetFillColor(0);
   ptstats->SetLineWidth(2);
   ptstats->SetTextAlign(12);
   TText *ptstats_LaTex = ptstats->AddText("py6");
   ptstats_LaTex->SetTextSize(0.04088889);
   ptstats_LaTex = ptstats->AddText("Entries = 642    ");
   ptstats_LaTex = ptstats->AddText("Mean  =   2.62");
   ptstats_LaTex = ptstats->AddText("Std Dev   = 0.02148");
   ptstats_LaTex = ptstats->AddText("#chi^{2} / ndf = 10.4873 / 18");
   ptstats_LaTex = ptstats->AddText("Prob  = 0.914836");
   ptstats_LaTex = ptstats->AddText("Constant = 48.712 #pm 2.637 ");
   ptstats_LaTex = ptstats->AddText("Mean     = 2.61846 #pm 0.00082 ");
   ptstats_LaTex = ptstats->AddText("Sigma    = 0.0178757 #pm 0.0007684 ");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(111111);
   ptstats->Draw();
   py6__2->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(py6__2);
   py6__2->SetLineWidth(2);
   py6__2->SetMarkerStyle(20);
   py6__2->GetXaxis()->SetRange(26,103);
   py6__2->GetXaxis()->SetNdivisions(515);
   py6__2->GetXaxis()->SetLabelFont(42);
   py6__2->GetXaxis()->SetLabelOffset(0.0125);
   py6__2->GetXaxis()->SetLabelSize(0.04);
   py6__2->GetXaxis()->SetTitleSize(0.05);
   py6__2->GetXaxis()->SetTitleOffset(1);
   py6__2->GetXaxis()->SetTitleFont(42);
   py6__2->GetYaxis()->SetNdivisions(515);
   py6__2->GetYaxis()->SetLabelFont(42);
   py6__2->GetYaxis()->SetLabelOffset(0.0125);
   py6__2->GetYaxis()->SetLabelSize(0.04);
   py6__2->GetYaxis()->SetTitleSize(0.05);
   py6__2->GetYaxis()->SetTitleOffset(1);
   py6__2->GetYaxis()->SetTitleFont(42);
   py6__2->GetZaxis()->SetNdivisions(6);
   py6__2->GetZaxis()->SetLabelFont(42);
   py6__2->GetZaxis()->SetLabelSize(0.04);
   py6__2->GetZaxis()->SetTitleSize(0.05);
   py6__2->GetZaxis()->SetTitleOffset(1);
   py6__2->GetZaxis()->SetTitleFont(42);
   py6__2->Draw("");
   
   TPaveText *pt = new TPaveText(0.1517997,0.8362761,0.2856025,0.8908507,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(2);
   pt->SetLineWidth(2);
   TText *pt_LaTex = pt->AddText("Psi6(lnY)");
   pt->Draw();
   MlnY->Modified();
   MlnY->cd();
   MlnY->SetSelected(MlnY);
}
