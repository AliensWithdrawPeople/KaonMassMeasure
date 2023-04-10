#ifdef __CLING__
#pragma cling optimize(0)
#endif

#include "TCanvas.h"
#include "TH1D.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TF1.h"
#include "TStyle.h"

void MC510_5_Smeared5_resolution()
{
//=========Macro generated from canvas: MlnY/Mass(lnY)
//=========  (Mon Apr 10 15:15:44 2023) by ROOT version 6.26/06
   TCanvas *MlnY = new TCanvas("MlnY", "Mass(lnY)",0,23,1280,649);
   gStyle->SetOptFit(1);
   MlnY->Range(2.404424,-110.3696,2.907557,579.4405);
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
   
   TH1D *py6__1 = new TH1D("py6__1","Psi6(lnY)",250,2.4,3.3);
   py6__1->SetBinContent(24,1);
   py6__1->SetBinContent(30,1);
   py6__1->SetBinContent(31,1);
   py6__1->SetBinContent(33,1);
   py6__1->SetBinContent(36,1);
   py6__1->SetBinContent(37,1);
   py6__1->SetBinContent(38,1);
   py6__1->SetBinContent(39,3);
   py6__1->SetBinContent(41,1);
   py6__1->SetBinContent(42,8);
   py6__1->SetBinContent(43,1);
   py6__1->SetBinContent(44,2);
   py6__1->SetBinContent(45,8);
   py6__1->SetBinContent(46,8);
   py6__1->SetBinContent(47,19);
   py6__1->SetBinContent(48,18);
   py6__1->SetBinContent(49,27);
   py6__1->SetBinContent(50,30);
   py6__1->SetBinContent(51,52);
   py6__1->SetBinContent(52,78);
   py6__1->SetBinContent(53,120);
   py6__1->SetBinContent(54,136);
   py6__1->SetBinContent(55,164);
   py6__1->SetBinContent(56,243);
   py6__1->SetBinContent(57,318);
   py6__1->SetBinContent(58,426);
   py6__1->SetBinContent(59,419);
   py6__1->SetBinContent(60,513);
   py6__1->SetBinContent(61,497);
   py6__1->SetBinContent(62,519);
   py6__1->SetBinContent(63,441);
   py6__1->SetBinContent(64,433);
   py6__1->SetBinContent(65,352);
   py6__1->SetBinContent(66,293);
   py6__1->SetBinContent(67,246);
   py6__1->SetBinContent(68,174);
   py6__1->SetBinContent(69,134);
   py6__1->SetBinContent(70,102);
   py6__1->SetBinContent(71,83);
   py6__1->SetBinContent(72,90);
   py6__1->SetBinContent(73,54);
   py6__1->SetBinContent(74,23);
   py6__1->SetBinContent(75,33);
   py6__1->SetBinContent(76,22);
   py6__1->SetBinContent(77,26);
   py6__1->SetBinContent(78,13);
   py6__1->SetBinContent(79,12);
   py6__1->SetBinContent(80,14);
   py6__1->SetBinContent(81,8);
   py6__1->SetBinContent(82,9);
   py6__1->SetBinContent(83,9);
   py6__1->SetBinContent(84,4);
   py6__1->SetBinContent(85,2);
   py6__1->SetBinContent(87,6);
   py6__1->SetBinContent(88,3);
   py6__1->SetBinContent(89,2);
   py6__1->SetBinContent(90,4);
   py6__1->SetBinContent(91,1);
   py6__1->SetBinContent(92,3);
   py6__1->SetBinContent(93,3);
   py6__1->SetBinContent(94,3);
   py6__1->SetBinContent(95,2);
   py6__1->SetBinContent(96,4);
   py6__1->SetBinContent(97,1);
   py6__1->SetBinContent(98,1);
   py6__1->SetBinContent(100,2);
   py6__1->SetBinContent(102,1);
   py6__1->SetBinContent(103,1);
   py6__1->SetBinContent(104,1);
   py6__1->SetBinContent(106,1);
   py6__1->SetBinContent(107,1);
   py6__1->SetBinContent(114,1);
   py6__1->SetEntries(6235);
   
   TF1 *gaus1 = new TF1("gaus","gaus",2.582889,2.656513, TF1::EAddToList::kNo);
   gaus1->SetFillColor(19);
   gaus1->SetFillStyle(0);
   gaus1->SetMarkerStyle(20);
   gaus1->SetLineWidth(3);
   gaus1->SetChisquare(31.87409);
   gaus1->SetNDF(17);
   gaus1->GetXaxis()->SetNdivisions(515);
   gaus1->GetXaxis()->SetLabelFont(42);
   gaus1->GetXaxis()->SetLabelOffset(0.0125);
   gaus1->GetXaxis()->SetLabelSize(0.04);
   gaus1->GetXaxis()->SetTitleSize(0.05);
   gaus1->GetXaxis()->SetTitleOffset(1);
   gaus1->GetXaxis()->SetTitleFont(42);
   gaus1->GetYaxis()->SetNdivisions(515);
   gaus1->GetYaxis()->SetLabelFont(42);
   gaus1->GetYaxis()->SetLabelOffset(0.0125);
   gaus1->GetYaxis()->SetLabelSize(0.04);
   gaus1->GetYaxis()->SetTitleSize(0.05);
   gaus1->GetYaxis()->SetTitleOffset(1);
   gaus1->GetYaxis()->SetTitleFont(42);
   gaus1->SetParameter(0,493.5905);
   gaus1->SetParError(0,8.780403);
   gaus1->SetParLimits(0,0,0);
   gaus1->SetParameter(1,2.619266);
   gaus1->SetParError(1,0.0002515676);
   gaus1->SetParLimits(1,0,0);
   gaus1->SetParameter(2,0.01705428);
   gaus1->SetParError(2,0.0002480881);
   gaus1->SetParLimits(2,0,0.1542052);
   gaus1->SetParent(py6__1);
   py6__1->GetListOfFunctions()->Add(gaus1);
   
   TPaveStats *ptstats = new TPaveStats(0.62,0.595,0.98,0.995,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(2);
   ptstats->SetFillColor(0);
   ptstats->SetLineWidth(2);
   ptstats->SetTextAlign(12);
   TText *ptstats_LaTex = ptstats->AddText("py6");
   ptstats_LaTex->SetTextSize(0.04088889);
   ptstats_LaTex = ptstats->AddText("Entries = 6235   ");
   ptstats_LaTex = ptstats->AddText("Mean  =  2.621");
   ptstats_LaTex = ptstats->AddText("Std Dev   = 0.02276");
   ptstats_LaTex = ptstats->AddText("#chi^{2} / ndf = 31.8741 / 17");
   ptstats_LaTex = ptstats->AddText("Prob  = 0.0155991");
   ptstats_LaTex = ptstats->AddText("Constant = 493.59 #pm 8.78 ");
   ptstats_LaTex = ptstats->AddText("Mean     = 2.61927 #pm 0.00025 ");
   ptstats_LaTex = ptstats->AddText("Sigma    = 0.0170543 #pm 0.0002481 ");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(111111);
   ptstats->Draw();
   py6__1->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(py6__1);
   py6__1->SetLineColor(2);
   py6__1->SetLineWidth(2);
   py6__1->SetMarkerStyle(20);
   py6__1->GetXaxis()->SetRange(19,134);
   py6__1->GetXaxis()->SetNdivisions(515);
   py6__1->GetXaxis()->SetLabelFont(42);
   py6__1->GetXaxis()->SetLabelOffset(0.0125);
   py6__1->GetXaxis()->SetLabelSize(0.04);
   py6__1->GetXaxis()->SetTitleSize(0.05);
   py6__1->GetXaxis()->SetTitleOffset(1);
   py6__1->GetXaxis()->SetTitleFont(42);
   py6__1->GetYaxis()->SetNdivisions(515);
   py6__1->GetYaxis()->SetLabelFont(42);
   py6__1->GetYaxis()->SetLabelOffset(0.0125);
   py6__1->GetYaxis()->SetLabelSize(0.04);
   py6__1->GetYaxis()->SetTitleSize(0.05);
   py6__1->GetYaxis()->SetTitleOffset(1);
   py6__1->GetYaxis()->SetTitleFont(42);
   py6__1->GetZaxis()->SetNdivisions(6);
   py6__1->GetZaxis()->SetLabelFont(42);
   py6__1->GetZaxis()->SetLabelSize(0.04);
   py6__1->GetZaxis()->SetTitleSize(0.05);
   py6__1->GetZaxis()->SetTitleOffset(1);
   py6__1->GetZaxis()->SetTitleFont(42);
   py6__1->Draw("");
   
   TPaveText *pt = new TPaveText(0.01,0.9404546,0.1437124,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(2);
   pt->SetLineWidth(2);
   TText *pt_LaTex = pt->AddText("Psi6(lnY)");
   pt->Draw();
   MlnY->Modified();
   MlnY->cd();
   MlnY->SetSelected(MlnY);
}
