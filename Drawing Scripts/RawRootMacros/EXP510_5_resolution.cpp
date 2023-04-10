#ifdef __CLING__
#pragma cling optimize(0)
#endif

#include "TCanvas.h"
#include "TH1D.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TF1.h"
#include "TStyle.h"


void EXP510_5_resolution()
{
//=========Macro generated from canvas: MlnY/Mass(lnY)
//=========  (Mon Apr 10 15:17:55 2023) by ROOT version 6.26/06
   TCanvas *MlnY = new TCanvas("MlnY", "Mass(lnY)",0,23,1280,649);
   gStyle->SetOptFit(1);
   MlnY->Range(2.405552,-104.1245,2.839287,546.6539);
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
   py6__1->SetBinContent(25,1);
   py6__1->SetBinContent(27,1);
   py6__1->SetBinContent(28,2);
   py6__1->SetBinContent(29,3);
   py6__1->SetBinContent(32,2);
   py6__1->SetBinContent(33,2);
   py6__1->SetBinContent(34,3);
   py6__1->SetBinContent(35,3);
   py6__1->SetBinContent(37,1);
   py6__1->SetBinContent(38,1);
   py6__1->SetBinContent(39,1);
   py6__1->SetBinContent(40,3);
   py6__1->SetBinContent(41,4);
   py6__1->SetBinContent(42,4);
   py6__1->SetBinContent(43,6);
   py6__1->SetBinContent(44,4);
   py6__1->SetBinContent(45,6);
   py6__1->SetBinContent(46,7);
   py6__1->SetBinContent(47,7);
   py6__1->SetBinContent(48,22);
   py6__1->SetBinContent(49,28);
   py6__1->SetBinContent(50,40);
   py6__1->SetBinContent(51,40);
   py6__1->SetBinContent(52,70);
   py6__1->SetBinContent(53,93);
   py6__1->SetBinContent(54,134);
   py6__1->SetBinContent(55,147);
   py6__1->SetBinContent(56,240);
   py6__1->SetBinContent(57,287);
   py6__1->SetBinContent(58,365);
   py6__1->SetBinContent(59,375);
   py6__1->SetBinContent(60,468);
   py6__1->SetBinContent(61,443);
   py6__1->SetBinContent(62,447);
   py6__1->SetBinContent(63,384);
   py6__1->SetBinContent(64,350);
   py6__1->SetBinContent(65,310);
   py6__1->SetBinContent(66,266);
   py6__1->SetBinContent(67,188);
   py6__1->SetBinContent(68,158);
   py6__1->SetBinContent(69,129);
   py6__1->SetBinContent(70,97);
   py6__1->SetBinContent(71,87);
   py6__1->SetBinContent(72,50);
   py6__1->SetBinContent(73,40);
   py6__1->SetBinContent(74,32);
   py6__1->SetBinContent(75,19);
   py6__1->SetBinContent(76,28);
   py6__1->SetBinContent(77,14);
   py6__1->SetBinContent(78,15);
   py6__1->SetBinContent(79,10);
   py6__1->SetBinContent(80,15);
   py6__1->SetBinContent(81,11);
   py6__1->SetBinContent(82,6);
   py6__1->SetBinContent(83,8);
   py6__1->SetBinContent(84,4);
   py6__1->SetBinContent(85,3);
   py6__1->SetBinContent(86,4);
   py6__1->SetBinContent(87,4);
   py6__1->SetBinContent(89,3);
   py6__1->SetBinContent(90,3);
   py6__1->SetBinContent(91,1);
   py6__1->SetBinContent(93,2);
   py6__1->SetBinContent(94,1);
   py6__1->SetBinContent(96,1);
   py6__1->SetBinContent(97,1);
   py6__1->SetBinContent(102,1);
   py6__1->SetBinContent(103,1);
   py6__1->SetBinContent(104,1);
   py6__1->SetBinContent(106,1);
   py6__1->SetBinContent(107,1);
   py6__1->SetBinContent(110,1);
   py6__1->SetBinContent(111,1);
   py6__1->SetBinContent(114,1);
   py6__1->SetEntries(5513);
   
   TF1 *gaus1 = new TF1("gaus","gaus",2.582296,2.656529, TF1::EAddToList::kNo);
   gaus1->SetFillColor(19);
   gaus1->SetFillStyle(0);
   gaus1->SetMarkerStyle(20);
   gaus1->SetLineWidth(3);
   gaus1->SetChisquare(37.14289);
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
   gaus1->SetParameter(0,429.965);
   gaus1->SetParError(0,8.192908);
   gaus1->SetParLimits(0,0,0);
   gaus1->SetParameter(1,2.619064);
   gaus1->SetParError(1,0.0002753365);
   gaus1->SetParLimits(1,0,0);
   gaus1->SetParameter(2,0.0173616);
   gaus1->SetParError(2,0.0002769931);
   gaus1->SetParLimits(2,0,0.1560027);
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
   ptstats_LaTex = ptstats->AddText("Entries = 5513   ");
   ptstats_LaTex = ptstats->AddText("Mean  =   2.62");
   ptstats_LaTex = ptstats->AddText("Std Dev   = 0.02333");
   ptstats_LaTex = ptstats->AddText("#chi^{2} / ndf = 37.1429 / 17");
   ptstats_LaTex = ptstats->AddText("Prob  = 0.00321829");
   ptstats_LaTex = ptstats->AddText("Constant = 429.965 #pm 8.193 ");
   ptstats_LaTex = ptstats->AddText("Mean     = 2.61906 #pm 0.00028 ");
   ptstats_LaTex = ptstats->AddText("Sigma    = 0.0173616 #pm 0.0002770 ");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(111111);
   ptstats->Draw();
   py6__1->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(py6__1);
   py6__1->SetLineColor(4);
   py6__1->SetLineWidth(2);
   py6__1->SetMarkerColor(4);
   py6__1->SetMarkerStyle(20);
   py6__1->SetMarkerSize(1.4);
   py6__1->GetXaxis()->SetRange(17,116);
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
   py6__1->Draw("E");
   
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
