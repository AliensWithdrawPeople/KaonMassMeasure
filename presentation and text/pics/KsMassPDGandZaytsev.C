#ifdef __CLING__
#pragma cling optimize(0)
#endif
void KsMassPDGandZaytsev()
{
//=========Macro generated from canvas: c1_n4/c1
//=========  (Thu Sep  1 21:22:17 2022) by ROOT version 6.26/04
   TCanvas *c1_n4 = new TCanvas("c1_n4", "c1",123,97,1291,581);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1_n4->Range(503.8109,497.5077,508.0733,497.868);
   c1_n4->SetFillColor(0);
   c1_n4->SetBorderMode(0);
   c1_n4->SetBorderSize(2);
   c1_n4->SetGridy();
   c1_n4->SetTickx(1);
   c1_n4->SetTicky(1);
   c1_n4->SetLeftMargin(0.1616688);
   c1_n4->SetRightMargin(0.007822686);
   c1_n4->SetTopMargin(0.05);
   c1_n4->SetBottomMargin(0.1173403);
   c1_n4->SetFrameBorderMode(0);
   c1_n4->SetFrameBorderMode(0);
   
   Double_t CMD_fx1001[6] = {
   505,
   505.5,
   506,
   506.5,
   507,
   507.5};
   Double_t CMD_fy1001[6] = {
   497.742,
   497.661,
   497.625,
   497.634,
   497.583,
   497.607};
   Double_t CMD_fex1001[6] = {
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t CMD_fey1001[6] = {
   0.085,
   0.033,
   0.031,
   0.024,
   0.021,
   0.017};
   TGraphErrors *gre = new TGraphErrors(6,CMD_fx1001,CMD_fy1001,CMD_fex1001,CMD_fey1001);
   gre->SetName("CMD");
   gre->SetTitle("CMD");
   gre->SetFillStyle(1000);
   gre->SetLineWidth(2);
   gre->SetMarkerStyle(20);
   
   TH1F *Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001 = new TH1F("Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001","CMD (1985)",100,504.5,510.5);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->SetMinimum(497.55);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->SetMaximum(497.85);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->SetDirectory(0);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->SetStats(0);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->SetLineWidth(2);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->SetMarkerStyle(20);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetXaxis()->SetRange(1,59);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetXaxis()->SetNdivisions(0);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetXaxis()->SetLabelOffset(0.0125);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetXaxis()->SetLabelSize(0.04);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetXaxis()->SetTitleOffset(1);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetYaxis()->SetTitle("M_{K^{0}}, #frac{MeV}{c^{2}}");
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetYaxis()->CenterTitle(true);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetYaxis()->SetNdivisions(512);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetYaxis()->SetLabelOffset(0.0125);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetYaxis()->SetLabelSize(0.06);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetZaxis()->SetNdivisions(6);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetZaxis()->SetLabelSize(0.04);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetZaxis()->SetTitleSize(0.05);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetZaxis()->SetTitleOffset(1);
   Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph_Graph_Graph_Graph_CMDsPoP1985cP10011001100310051001);
   
   
   TF1 *PrevFitTMP1002 = new TF1("PrevFitTMP","pol0",504.5,510.5, TF1::EAddToList::kNo);
   PrevFitTMP1002->SetFillColor(19);
   PrevFitTMP1002->SetFillStyle(0);
   PrevFitTMP1002->SetMarkerStyle(20);
   PrevFitTMP1002->SetLineWidth(3);
   PrevFitTMP1002->SetChisquare(7.449172);
   PrevFitTMP1002->SetNDF(5);
   PrevFitTMP1002->GetXaxis()->SetNdivisions(515);
   PrevFitTMP1002->GetXaxis()->SetLabelFont(42);
   PrevFitTMP1002->GetXaxis()->SetLabelOffset(0.0125);
   PrevFitTMP1002->GetXaxis()->SetLabelSize(0.04);
   PrevFitTMP1002->GetXaxis()->SetTitleSize(0.05);
   PrevFitTMP1002->GetXaxis()->SetTitleOffset(1);
   PrevFitTMP1002->GetXaxis()->SetTitleFont(42);
   PrevFitTMP1002->GetYaxis()->SetNdivisions(6);
   PrevFitTMP1002->GetYaxis()->SetLabelFont(42);
   PrevFitTMP1002->GetYaxis()->SetLabelOffset(0.0125);
   PrevFitTMP1002->GetYaxis()->SetLabelSize(0.04);
   PrevFitTMP1002->GetYaxis()->SetTitleSize(0.05);
   PrevFitTMP1002->GetYaxis()->SetTitleOffset(1);
   PrevFitTMP1002->GetYaxis()->SetTitleFont(42);
   PrevFitTMP1002->SetParameter(0,497.6153);
   PrevFitTMP1002->SetParError(0,0.01022692);
   PrevFitTMP1002->SetParLimits(0,0,0);
   PrevFitTMP1002->SetParent(gre);
   gre->GetListOfFunctions()->Add(PrevFitTMP1002);
   
   TPaveStats *ptstats = new TPaveStats(0.62,0.835,0.98,0.995,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(2);
   ptstats->SetFillColor(0);
   ptstats->SetLineWidth(2);
   ptstats->SetTextAlign(12);
   TText *ptstats_LaTex = ptstats->AddText("#chi^{2} / ndf = 7.449 / 5");
   ptstats_LaTex = ptstats->AddText("Prob  = 0.1893");
   ptstats_LaTex = ptstats->AddText("p0       = 497.6 #pm 0.01023 ");
   ptstats->SetOptStat(0);
   ptstats->SetOptFit(111111);
   ptstats->Draw();
   gre->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(gre->GetListOfFunctions());
   gre->Draw(" ap");
   TLine *line = new TLine(504.53,497.611,508,497.611);
   line->SetLineColor(4);
   line->SetLineWidth(2);
   line->Draw();
   TText *text = new TText(505.2,497.65,"CMD '85");
   text->SetTextFont(42);
   text->SetTextAngle(90);
   text->Draw();
   text = new TText(505.7,497.65,"CMD '87");
   text->SetTextFont(42);
   text->SetTextAngle(90);
   text->Draw();
   text = new TText(506.2,497.62,"NA48 '02");
   text->SetTextFont(42);
   text->SetTextAngle(90);
   text->Draw();
   text = new TText(506.7,497.625,"CMD-2 '03");
   text->SetTextFont(42);
   text->SetTextAngle(90);
   text->Draw();
   text = new TText(507.2,497.63,"KLOE '07");
   text->SetTextFont(42);
   text->SetTextAngle(90);
   text->Draw();
   text = new TText(507.7,497.63,"CLEO-c data '14");
   text->SetTextFont(42);
   text->SetTextAngle(90);
   text->Draw();
   c1_n4->Modified();
   c1_n4->cd();
   c1_n4->SetSelected(c1_n4);
}
