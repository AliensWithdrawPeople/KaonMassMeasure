#ifdef __CLING__
#pragma cling optimize(0)
#endif
void RawTrackHits()
{
//=========Macro generated from canvas: c1/c1
//=========  (Wed Sep  7 12:09:45 2022) by ROOT version 6.26/04
   TCanvas *c1 = new TCanvas("c1", "c1",0,23,1536,793);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(0.6626506,-5735.605,36.80723,30111.93);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.12);
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.16);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   TH1D *hHit__1 = new TH1D("hHit__1","nhits",40,0,40);
   hHit__1->SetBinContent(6,8027);
   hHit__1->SetBinContent(7,7777);
   hHit__1->SetBinContent(8,7152);
   hHit__1->SetBinContent(9,6888);
   hHit__1->SetBinContent(10,6589);
   hHit__1->SetBinContent(11,6488);
   hHit__1->SetBinContent(12,6594);
   hHit__1->SetBinContent(13,7116);
   hHit__1->SetBinContent(14,9676);
   hHit__1->SetBinContent(15,14926);
   hHit__1->SetBinContent(16,21965);
   hHit__1->SetBinContent(17,26971);
   hHit__1->SetBinContent(18,25725);
   hHit__1->SetBinContent(19,18565);
   hHit__1->SetBinContent(20,9868);
   hHit__1->SetBinContent(21,3883);
   hHit__1->SetBinContent(22,1355);
   hHit__1->SetBinContent(23,515);
   hHit__1->SetBinContent(24,261);
   hHit__1->SetBinContent(25,154);
   hHit__1->SetBinContent(26,120);
   hHit__1->SetBinContent(27,104);
   hHit__1->SetBinContent(28,100);
   hHit__1->SetBinContent(29,83);
   hHit__1->SetBinContent(30,77);
   hHit__1->SetBinContent(31,63);
   hHit__1->SetBinContent(32,79);
   hHit__1->SetBinContent(33,73);
   hHit__1->SetBinContent(34,68);
   hHit__1->SetBinContent(35,55);
   hHit__1->SetBinContent(36,55);
   hHit__1->SetBinContent(37,52);
   hHit__1->SetBinContent(38,39);
   hHit__1->SetBinContent(39,37);
   hHit__1->SetBinContent(40,33);
   hHit__1->SetBinContent(41,120);
   hHit__1->SetEntries(191653);
   hHit__1->SetLineWidth(2);
   hHit__1->SetMarkerStyle(20);
   hHit__1->GetXaxis()->SetTitle("Number of hits");
   hHit__1->GetXaxis()->SetRange(6,35);
   hHit__1->GetXaxis()->SetNdivisions(515);
   hHit__1->GetXaxis()->SetLabelFont(42);
   hHit__1->GetXaxis()->SetLabelOffset(0.0125);
   hHit__1->GetXaxis()->SetLabelSize(0.04);
   hHit__1->GetXaxis()->SetTitleSize(0.05);
   hHit__1->GetXaxis()->SetTitleOffset(1);
   hHit__1->GetXaxis()->SetTitleFont(42);
   hHit__1->GetYaxis()->SetNdivisions(6);
   hHit__1->GetYaxis()->SetLabelFont(42);
   hHit__1->GetYaxis()->SetLabelOffset(0.0125);
   hHit__1->GetYaxis()->SetLabelSize(0.04);
   hHit__1->GetYaxis()->SetTitleSize(0.05);
   hHit__1->GetYaxis()->SetTitleOffset(1);
   hHit__1->GetYaxis()->SetTitleFont(42);
   hHit__1->GetZaxis()->SetNdivisions(6);
   hHit__1->GetZaxis()->SetLabelFont(42);
   hHit__1->GetZaxis()->SetLabelSize(0.04);
   hHit__1->GetZaxis()->SetTitleSize(0.05);
   hHit__1->GetZaxis()->SetTitleOffset(1);
   hHit__1->GetZaxis()->SetTitleFont(42);
   hHit__1->Draw("");
   TLine *line = new TLine(10,0,10,27000);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#0000ff");
   line->SetLineColor(ci);
   line->SetLineWidth(4);
   line->Draw();
   line = new TLine(30,0,30,27027.26);

   ci = TColor::GetColor("#0000ff");
   line->SetLineColor(ci);
   line->SetLineWidth(4);
   line->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
