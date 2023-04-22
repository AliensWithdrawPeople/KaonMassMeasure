#ifdef __CLING__
#pragma cling optimize(0)
#endif

#include "TCanvas.h"
#include "TH1D.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TStyle.h"

void EnergyControl510()
{
//=========Macro generated from canvas: c1/c1
//=========  (Sat Apr 22 11:18:22 2023) by ROOT version 6.26/06
   TCanvas *c1 = new TCanvas("c1", "c1",0,23,2048,1081);
   gStyle->SetOptFit(1);
   c1->Range(61249.82,509.6962,61938.01,510.1176);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetGridx();
   c1->SetGridy();
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.12);
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.16);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   Double_t grEmeas_fx1001[38] = {
   61385,
   61397,
   61410,
   61420,
   61425,
   61434,
   61447,
   61458,
   61468,
   61479,
   61487,
   61497,
   61511,
   61525,
   61540,
   61551,
   61568,
   61577,
   61589,
   61605,
   61621,
   61636,
   61650,
   61666,
   61682,
   61697,
   61712,
   61726,
   61740,
   61754,
   61767,
   61779,
   61791,
   61805,
   61819,
   61832,
   61842,
   61851};
   Double_t grEmeas_fy1001[38] = {
   510.049,
   509.997,
   509.985,
   510.031,
   510.029,
   509.988,
   509.914,
   509.936,
   509.925,
   509.916,
   509.944,
   509.978,
   509.914,
   509.906,
   509.927,
   509.924,
   509.992,
   510,
   509.944,
   509.94,
   509.951,
   509.967,
   509.974,
   509.968,
   509.92,
   509.94,
   510.007,
   509.952,
   509.993,
   509.991,
   509.994,
   509.981,
   509.947,
   509.912,
   509.975,
   509.887,
   509.958,
   509.916};
   Double_t grEmeas_fex1001[38] = {
   5,
   6.5,
   5.5,
   3,
   1,
   6.5,
   5,
   4,
   5.5,
   4,
   3,
   6,
   7,
   6.5,
   7,
   3.5,
   4,
   4,
   7.5,
   7.5,
   7,
   7,
   6.5,
   7,
   8,
   6.5,
   7,
   6.5,
   6.5,
   6.5,
   5.5,
   5,
   6.5,
   6.5,
   6.5,
   5,
   3,
   5};
   Double_t grEmeas_fey1001[38] = {
   0.031,
   0.044,
   0.044,
   0.026,
   0.029,
   0.027,
   0.035,
   0.031,
   0.025,
   0.027,
   0.028,
   0.029,
   0.024,
   0.02,
   0.03,
   0.032,
   0.031,
   0.028,
   0.029,
   0.038,
   0.034,
   0.029,
   0.026,
   0.029,
   0.032,
   0.032,
   0.027,
   0.029,
   0.031,
   0.031,
   0.031,
   0.031,
   0.03,
   0.03,
   0.031,
   0.031,
   0.035,
   0.035};
   TGraphErrors *gre = new TGraphErrors(38,grEmeas_fx1001,grEmeas_fy1001,grEmeas_fex1001,grEmeas_fey1001);
   gre->SetName("grEmeas");
   gre->SetTitle("Red -- emeas, black -- Kch E, blue band -- compton mean");
   gre->SetFillStyle(1000);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ff0000");
   gre->SetLineColor(ci);
   gre->SetLineWidth(3);

   ci = TColor::GetColor("#ff0000");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(1.5);
   
   TH1F *Graph_grEmeas1001 = new TH1F("Graph_grEmeas1001","Red -- emeas, black -- Kch E, blue band -- compton mean",100,61332.4,61903.6);
   Graph_grEmeas1001->SetMinimum(509.7637);
   Graph_grEmeas1001->SetMaximum(510.0965);
   Graph_grEmeas1001->SetDirectory(0);
   Graph_grEmeas1001->SetStats(0);
   Graph_grEmeas1001->SetLineWidth(3);
   Graph_grEmeas1001->SetMarkerStyle(20);
   Graph_grEmeas1001->SetMarkerSize(1.5);
   Graph_grEmeas1001->GetXaxis()->SetTitle("Run");
   Graph_grEmeas1001->GetXaxis()->SetNdivisions(515);
   Graph_grEmeas1001->GetXaxis()->SetLabelFont(42);
   Graph_grEmeas1001->GetXaxis()->SetLabelOffset(0.0125);
   Graph_grEmeas1001->GetXaxis()->SetLabelSize(0.04);
   Graph_grEmeas1001->GetXaxis()->SetTitleSize(0.05);
   Graph_grEmeas1001->GetXaxis()->SetTitleOffset(1);
   Graph_grEmeas1001->GetXaxis()->SetTitleFont(42);
   Graph_grEmeas1001->GetYaxis()->SetTitle("Energy, MeV");
   Graph_grEmeas1001->GetYaxis()->SetNdivisions(515);
   Graph_grEmeas1001->GetYaxis()->SetLabelFont(42);
   Graph_grEmeas1001->GetYaxis()->SetLabelOffset(0.0125);
   Graph_grEmeas1001->GetYaxis()->SetLabelSize(0.04);
   Graph_grEmeas1001->GetYaxis()->SetTitleSize(0.05);
   Graph_grEmeas1001->GetYaxis()->SetTitleOffset(1);
   Graph_grEmeas1001->GetYaxis()->SetTitleFont(42);
   Graph_grEmeas1001->GetZaxis()->SetNdivisions(6);
   Graph_grEmeas1001->GetZaxis()->SetLabelFont(42);
   Graph_grEmeas1001->GetZaxis()->SetLabelSize(0.04);
   Graph_grEmeas1001->GetZaxis()->SetTitleSize(0.05);
   Graph_grEmeas1001->GetZaxis()->SetTitleOffset(1);
   Graph_grEmeas1001->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_grEmeas1001);
   
   gre->Draw("ap");
   
   Double_t grComptonMeanEnergy_fx1002[2] = {
   61370,
   61866};
   Double_t grComptonMeanEnergy_fy1002[2] = {
   509.956,
   509.956};
   Double_t grComptonMeanEnergy_fex1002[2] = {
   0,
   0};
   Double_t grComptonMeanEnergy_fey1002[2] = {
   0.005,
   0.005};
   gre = new TGraphErrors(2,grComptonMeanEnergy_fx1002,grComptonMeanEnergy_fy1002,grComptonMeanEnergy_fex1002,grComptonMeanEnergy_fey1002);
   gre->SetName("grComptonMeanEnergy");
   gre->SetTitle("");

   ci = TColor::GetColor("#0000ff");
   gre->SetFillColor(ci);
   gre->SetFillStyle(3005);

   ci = TColor::GetColor("#0000ff");
   gre->SetLineColor(ci);
   gre->SetLineWidth(4);

   ci = TColor::GetColor("#0000ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(1.5);
   
   TH1F *Graph_grComptonMeanEnergy1002 = new TH1F("Graph_grComptonMeanEnergy1002","",100,61320.4,61915.6);
   Graph_grComptonMeanEnergy1002->SetMinimum(509.95);
   Graph_grComptonMeanEnergy1002->SetMaximum(509.962);
   Graph_grComptonMeanEnergy1002->SetDirectory(0);
   Graph_grComptonMeanEnergy1002->SetStats(0);
   Graph_grComptonMeanEnergy1002->SetLineWidth(3);
   Graph_grComptonMeanEnergy1002->SetMarkerStyle(20);
   Graph_grComptonMeanEnergy1002->SetMarkerSize(1.5);
   Graph_grComptonMeanEnergy1002->GetXaxis()->SetNdivisions(515);
   Graph_grComptonMeanEnergy1002->GetXaxis()->SetLabelFont(42);
   Graph_grComptonMeanEnergy1002->GetXaxis()->SetLabelOffset(0.0125);
   Graph_grComptonMeanEnergy1002->GetXaxis()->SetLabelSize(0.04);
   Graph_grComptonMeanEnergy1002->GetXaxis()->SetTitleSize(0.05);
   Graph_grComptonMeanEnergy1002->GetXaxis()->SetTitleOffset(1);
   Graph_grComptonMeanEnergy1002->GetXaxis()->SetTitleFont(42);
   Graph_grComptonMeanEnergy1002->GetYaxis()->SetNdivisions(515);
   Graph_grComptonMeanEnergy1002->GetYaxis()->SetLabelFont(42);
   Graph_grComptonMeanEnergy1002->GetYaxis()->SetLabelOffset(0.0125);
   Graph_grComptonMeanEnergy1002->GetYaxis()->SetLabelSize(0.04);
   Graph_grComptonMeanEnergy1002->GetYaxis()->SetTitleSize(0.05);
   Graph_grComptonMeanEnergy1002->GetYaxis()->SetTitleOffset(1);
   Graph_grComptonMeanEnergy1002->GetYaxis()->SetTitleFont(42);
   Graph_grComptonMeanEnergy1002->GetZaxis()->SetNdivisions(6);
   Graph_grComptonMeanEnergy1002->GetZaxis()->SetLabelFont(42);
   Graph_grComptonMeanEnergy1002->GetZaxis()->SetLabelSize(0.04);
   Graph_grComptonMeanEnergy1002->GetZaxis()->SetTitleSize(0.05);
   Graph_grComptonMeanEnergy1002->GetZaxis()->SetTitleOffset(1);
   Graph_grComptonMeanEnergy1002->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_grComptonMeanEnergy1002);
   
   gre->Draw("l3 ");
   
   Double_t grKchEnergy_fx1003[38] = {
   61385,
   61397,
   61410,
   61420,
   61425,
   61434,
   61447,
   61458,
   61468,
   61479,
   61487,
   61497,
   61511,
   61525,
   61540,
   61551,
   61568,
   61577,
   61589,
   61605,
   61621,
   61636,
   61650,
   61666,
   61682,
   61697,
   61712,
   61726,
   61740,
   61754,
   61767,
   61779,
   61791,
   61805,
   61819,
   61832,
   61842,
   61851};
   Double_t grKchEnergy_fy1003[38] = {
   509.9471,
   509.8676,
   509.8506,
   509.8612,
   509.85,
   509.8426,
   509.8315,
   509.8272,
   509.8254,
   509.8312,
   509.8107,
   509.8217,
   509.8168,
   509.809,
   509.8021,
   509.814,
   509.8857,
   509.8914,
   509.8639,
   509.8493,
   509.8536,
   509.842,
   509.835,
   509.8281,
   509.8371,
   509.8256,
   509.8347,
   509.8407,
   509.8437,
   509.8401,
   509.8357,
   509.8267,
   509.8338,
   509.826,
   509.8208,
   509.8196,
   509.8316,
   509.8235};
   Double_t grKchEnergy_fex1003[38] = {
   5,
   6.5,
   5.5,
   3,
   1,
   6.5,
   5,
   4,
   5.5,
   4,
   3,
   6,
   7,
   6.5,
   7,
   3.5,
   4,
   4,
   7.5,
   7.5,
   7,
   7,
   6.5,
   7,
   8,
   6.5,
   7,
   6.5,
   6.5,
   6.5,
   5.5,
   5,
   6.5,
   6.5,
   6.5,
   5,
   3,
   5};
   Double_t grKchEnergy_fey1003[38] = {
   0.007598007,
   0.006406982,
   0.006726283,
   0.01095321,
   0.0157649,
   0.00702742,
   0.009230571,
   0.008969479,
   0.006959426,
   0.00857276,
   0.01019226,
   0.006889328,
   0.006739956,
   0.007094312,
   0.007486487,
   0.009902036,
   0.009232044,
   0.009892941,
   0.006648914,
   0.006292443,
   0.007011882,
   0.006823654,
   0.008104713,
   0.006529504,
   0.006965077,
   0.006663254,
   0.006169513,
   0.006605253,
   0.006369451,
   0.006073671,
   0.007057939,
   0.007747112,
   0.006597415,
   0.006600739,
   0.006784084,
   0.007532346,
   0.008937997,
   0.007294001};
   gre = new TGraphErrors(38,grKchEnergy_fx1003,grKchEnergy_fy1003,grKchEnergy_fex1003,grKchEnergy_fey1003);
   gre->SetName("grKchEnergy");
   gre->SetTitle("Graph");
   gre->SetFillStyle(1000);
   gre->SetLineWidth(3);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(1.5);
   
   TH1F *Graph_grKchEnergy1003 = new TH1F("Graph_grKchEnergy1003","Graph",100,61332.4,61903.6);
   Graph_grKchEnergy1003->SetMinimum(509.7786);
   Graph_grKchEnergy1003->SetMaximum(509.9707);
   Graph_grKchEnergy1003->SetDirectory(0);
   Graph_grKchEnergy1003->SetStats(0);
   Graph_grKchEnergy1003->SetLineWidth(3);
   Graph_grKchEnergy1003->SetMarkerStyle(20);
   Graph_grKchEnergy1003->SetMarkerSize(1.5);
   Graph_grKchEnergy1003->GetXaxis()->SetNdivisions(515);
   Graph_grKchEnergy1003->GetXaxis()->SetLabelFont(42);
   Graph_grKchEnergy1003->GetXaxis()->SetLabelOffset(0.0125);
   Graph_grKchEnergy1003->GetXaxis()->SetLabelSize(0.04);
   Graph_grKchEnergy1003->GetXaxis()->SetTitleSize(0.05);
   Graph_grKchEnergy1003->GetXaxis()->SetTitleOffset(1);
   Graph_grKchEnergy1003->GetXaxis()->SetTitleFont(42);
   Graph_grKchEnergy1003->GetYaxis()->SetNdivisions(515);
   Graph_grKchEnergy1003->GetYaxis()->SetLabelFont(42);
   Graph_grKchEnergy1003->GetYaxis()->SetLabelOffset(0.0125);
   Graph_grKchEnergy1003->GetYaxis()->SetLabelSize(0.04);
   Graph_grKchEnergy1003->GetYaxis()->SetTitleSize(0.05);
   Graph_grKchEnergy1003->GetYaxis()->SetTitleOffset(1);
   Graph_grKchEnergy1003->GetYaxis()->SetTitleFont(42);
   Graph_grKchEnergy1003->GetZaxis()->SetNdivisions(6);
   Graph_grKchEnergy1003->GetZaxis()->SetLabelFont(42);
   Graph_grKchEnergy1003->GetZaxis()->SetLabelSize(0.04);
   Graph_grKchEnergy1003->GetZaxis()->SetTitleSize(0.05);
   Graph_grKchEnergy1003->GetZaxis()->SetTitleOffset(1);
   Graph_grKchEnergy1003->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_grKchEnergy1003);
   
   gre->Draw("p ");
   
   TPaveText *pt = new TPaveText(0.2150538,0.8966825,0.914956,0.9526066,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(2);
   pt->SetLineWidth(3);
   TText *pt_LaTex = pt->AddText("Red -- emeas, black -- Kch E, blue band -- compton mean");
   pt->Draw();
   
   TF1 *ke11 = new TF1("ke1","364.312 / (x-61332.6) / (x-61332.6) + 506.203",61379,61560, TF1::EAddToList::kDefault);
   ke11->SetFillColor(19);
   ke11->SetFillStyle(0);
   ke11->SetMarkerStyle(20);
   ke11->SetMarkerSize(1.5);
   ke11->SetLineWidth(3);
   ke11->GetXaxis()->SetNdivisions(515);
   ke11->GetXaxis()->SetLabelFont(42);
   ke11->GetXaxis()->SetLabelOffset(0.0125);
   ke11->GetXaxis()->SetLabelSize(0.04);
   ke11->GetXaxis()->SetTitleSize(0.05);
   ke11->GetXaxis()->SetTitleOffset(1);
   ke11->GetXaxis()->SetTitleFont(42);
   ke11->GetYaxis()->SetNdivisions(515);
   ke11->GetYaxis()->SetLabelFont(42);
   ke11->GetYaxis()->SetLabelOffset(0.0125);
   ke11->GetYaxis()->SetLabelSize(0.04);
   ke11->GetYaxis()->SetTitleSize(0.05);
   ke11->GetYaxis()->SetTitleOffset(1);
   ke11->GetYaxis()->SetTitleFont(42);
   ke11->Draw("same");
   
   TF1 *ke12 = new TF1("ke1","364.312 / (x-61332.6) / (x-61332.6) + 506.203 + 3.6",61379,61560, TF1::EAddToList::kDefault);
   ke12->SetFillColor(19);
   ke12->SetFillStyle(0);
   ke12->SetMarkerStyle(20);
   ke12->SetMarkerSize(1.5);
   ke12->SetLineWidth(3);
   ke12->GetXaxis()->SetNdivisions(515);
   ke12->GetXaxis()->SetLabelFont(42);
   ke12->GetXaxis()->SetLabelOffset(0.0125);
   ke12->GetXaxis()->SetLabelSize(0.04);
   ke12->GetXaxis()->SetTitleSize(0.05);
   ke12->GetXaxis()->SetTitleOffset(1);
   ke12->GetXaxis()->SetTitleFont(42);
   ke12->GetYaxis()->SetNdivisions(515);
   ke12->GetYaxis()->SetLabelFont(42);
   ke12->GetYaxis()->SetLabelOffset(0.0125);
   ke12->GetYaxis()->SetLabelSize(0.04);
   ke12->GetYaxis()->SetTitleSize(0.05);
   ke12->GetYaxis()->SetTitleOffset(1);
   ke12->GetYaxis()->SetTitleFont(42);
   ke12->Draw("same");
   
   TF1 *ke13 = new TF1("ke1","364.312 / (x-61332.6) / (x-61332.6) + 506.203 + 3.6",61379,61560, TF1::EAddToList::kDefault);
   ke13->SetFillColor(19);
   ke13->SetFillStyle(0);
   ke13->SetMarkerStyle(20);
   ke13->SetMarkerSize(1.5);
   ke13->GetXaxis()->SetNdivisions(515);
   ke13->GetXaxis()->SetLabelFont(42);
   ke13->GetXaxis()->SetLabelOffset(0.0125);
   ke13->GetXaxis()->SetLabelSize(0.04);
   ke13->GetXaxis()->SetTitleSize(0.05);
   ke13->GetXaxis()->SetTitleOffset(1);
   ke13->GetXaxis()->SetTitleFont(42);
   ke13->GetYaxis()->SetNdivisions(515);
   ke13->GetYaxis()->SetLabelFont(42);
   ke13->GetYaxis()->SetLabelOffset(0.0125);
   ke13->GetYaxis()->SetLabelSize(0.04);
   ke13->GetYaxis()->SetTitleSize(0.05);
   ke13->GetYaxis()->SetTitleOffset(1);
   ke13->GetYaxis()->SetTitleFont(42);
   ke13->Draw("same");
   
   TF1 *ce14 = new TF1("ce1","364.312 / (x-61332.6) / (x-61332.6) + 509.925",61379,61560, TF1::EAddToList::kDefault);
   ce14->SetFillColor(19);
   ce14->SetFillStyle(0);
   ce14->SetMarkerStyle(20);
   ce14->SetMarkerSize(1.5);
   ce14->SetLineColor(2);
   ce14->GetXaxis()->SetNdivisions(515);
   ce14->GetXaxis()->SetLabelFont(42);
   ce14->GetXaxis()->SetLabelOffset(0.0125);
   ce14->GetXaxis()->SetLabelSize(0.04);
   ce14->GetXaxis()->SetTitleSize(0.05);
   ce14->GetXaxis()->SetTitleOffset(1);
   ce14->GetXaxis()->SetTitleFont(42);
   ce14->GetYaxis()->SetNdivisions(515);
   ce14->GetYaxis()->SetLabelFont(42);
   ce14->GetYaxis()->SetLabelOffset(0.0125);
   ce14->GetYaxis()->SetLabelSize(0.04);
   ce14->GetYaxis()->SetTitleSize(0.05);
   ce14->GetYaxis()->SetTitleOffset(1);
   ce14->GetYaxis()->SetTitleFont(42);
   ce14->Draw("same");
   
   TF1 *ke25 = new TF1("ke2","598.053 / (x-61479.6) / (x-61479.6) + 506.217 + 3.6",61561,61685, TF1::EAddToList::kDefault);
   ke25->SetFillColor(19);
   ke25->SetFillStyle(0);
   ke25->SetMarkerStyle(20);
   ke25->SetMarkerSize(1.5);
   ke25->GetXaxis()->SetNdivisions(515);
   ke25->GetXaxis()->SetLabelFont(42);
   ke25->GetXaxis()->SetLabelOffset(0.0125);
   ke25->GetXaxis()->SetLabelSize(0.04);
   ke25->GetXaxis()->SetTitleSize(0.05);
   ke25->GetXaxis()->SetTitleOffset(1);
   ke25->GetXaxis()->SetTitleFont(42);
   ke25->GetYaxis()->SetNdivisions(515);
   ke25->GetYaxis()->SetLabelFont(42);
   ke25->GetYaxis()->SetLabelOffset(0.0125);
   ke25->GetYaxis()->SetLabelSize(0.04);
   ke25->GetYaxis()->SetTitleSize(0.05);
   ke25->GetYaxis()->SetTitleOffset(1);
   ke25->GetYaxis()->SetTitleFont(42);
   ke25->Draw("same");
   
   TF1 *ce26 = new TF1("ce2","598.053 / (x-61479.6) / (x-61479.6) + 509.928",61561,61685, TF1::EAddToList::kDefault);
   ce26->SetFillColor(19);
   ce26->SetFillStyle(0);
   ce26->SetMarkerStyle(20);
   ce26->SetMarkerSize(1.5);
   ce26->SetLineColor(2);
   ce26->GetXaxis()->SetNdivisions(515);
   ce26->GetXaxis()->SetLabelFont(42);
   ce26->GetXaxis()->SetLabelOffset(0.0125);
   ce26->GetXaxis()->SetLabelSize(0.04);
   ce26->GetXaxis()->SetTitleSize(0.05);
   ce26->GetXaxis()->SetTitleOffset(1);
   ce26->GetXaxis()->SetTitleFont(42);
   ce26->GetYaxis()->SetNdivisions(515);
   ce26->GetYaxis()->SetLabelFont(42);
   ce26->GetYaxis()->SetLabelOffset(0.0125);
   ce26->GetYaxis()->SetLabelSize(0.04);
   ce26->GetYaxis()->SetTitleSize(0.05);
   ce26->GetYaxis()->SetTitleOffset(1);
   ce26->GetYaxis()->SetTitleFont(42);
   ce26->Draw("same");
   
   TF1 *ke37 = new TF1("ke3","0.0101 * sin(0.0397 * (x-61704.6)) + 506.231 + 3.6",61690,61856, TF1::EAddToList::kDefault);
   ke37->SetFillColor(19);
   ke37->SetFillStyle(0);
   ke37->SetMarkerStyle(20);
   ke37->SetMarkerSize(1.5);
   ke37->GetXaxis()->SetNdivisions(515);
   ke37->GetXaxis()->SetLabelFont(42);
   ke37->GetXaxis()->SetLabelOffset(0.0125);
   ke37->GetXaxis()->SetLabelSize(0.04);
   ke37->GetXaxis()->SetTitleSize(0.05);
   ke37->GetXaxis()->SetTitleOffset(1);
   ke37->GetXaxis()->SetTitleFont(42);
   ke37->GetYaxis()->SetNdivisions(515);
   ke37->GetYaxis()->SetLabelFont(42);
   ke37->GetYaxis()->SetLabelOffset(0.0125);
   ke37->GetYaxis()->SetLabelSize(0.04);
   ke37->GetYaxis()->SetTitleSize(0.05);
   ke37->GetYaxis()->SetTitleOffset(1);
   ke37->GetYaxis()->SetTitleFont(42);
   ke37->Draw("same");
   
   TF1 *ce38 = new TF1("ce3","0.0101 * sin(0.0397 * (x-61704.6)) + 509.960",61690,61856, TF1::EAddToList::kDefault);
   ce38->SetFillColor(19);
   ce38->SetFillStyle(0);
   ce38->SetMarkerStyle(20);
   ce38->SetMarkerSize(1.5);
   ce38->SetLineColor(2);
   ce38->GetXaxis()->SetNdivisions(515);
   ce38->GetXaxis()->SetLabelFont(42);
   ce38->GetXaxis()->SetLabelOffset(0.0125);
   ce38->GetXaxis()->SetLabelSize(0.04);
   ce38->GetXaxis()->SetTitleSize(0.05);
   ce38->GetXaxis()->SetTitleOffset(1);
   ce38->GetXaxis()->SetTitleFont(42);
   ce38->GetYaxis()->SetNdivisions(515);
   ce38->GetYaxis()->SetLabelFont(42);
   ce38->GetYaxis()->SetLabelOffset(0.0125);
   ce38->GetYaxis()->SetLabelSize(0.04);
   ce38->GetYaxis()->SetTitleSize(0.05);
   ce38->GetYaxis()->SetTitleOffset(1);
   ce38->GetYaxis()->SetTitleFont(42);
   ce38->Draw("same");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
