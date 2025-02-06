#include <TFile.h>
#include <TGraphErrors.h>
#include <TRatioPlot.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TRoot.h> 
#include <TLegend.h>
#include <TAxis.h>
#include <TMath.h>
#include <TLatex.h>
#include <TLine.h>

#include <vector>
#include <string>
#include <iomanip>
#include <sstream>
#include <numeric>
#include <iostream>

int draw_born()
{
    auto bottom_pad_scale = 0.4;
    auto top_pad_scale = 1 - bottom_pad_scale;
    auto min_x = 1.000; // Minimum E_{cms} on plots
    auto max_x = 1.063; // Maximum E_{cms} on plots

    // Setup: start
    double n_params = 3;

    std::string filename = "bcs_kskl_sig_extr_syst.root";
    std::string output_filename = "C:/work/Science/KaonMassMeasure/presentation and text/paper/figs_v2/Born_xsec.pdf";
    // Setup: end

    auto canvas = new TCanvas("c1","", 800, 800);
    auto bottom_pad = new TPad("bottom_pad","bottom",0., 0., 1., bottom_pad_scale); 
    bottom_pad->Draw();

    bottom_pad->SetBottomMargin(0.4);

    bottom_pad->SetBorderMode(0);

    auto top_pad = new TPad("top_pad","top",0., bottom_pad_scale, 1., 1.);  
    top_pad->Draw();
    top_pad->SetBottomMargin(0.00001);
    top_pad->SetBorderMode(0);

    top_pad->cd();

    std::string fullFilename = "C:/work/Science/KaonMassMeasure/presentation and text/paper/figs_v2/scripts/" + filename; 
    auto file = TFile::Open(fullFilename.c_str());
    auto vcs = file->Get<TGraphErrors>("vcs");
    auto bcs = file->Get<TGraphErrors>("bcs");


    auto f_bcs = file->Get<TF1>("f_bcs");
    f_bcs->SetRange(1.002, 1.060);
    int N_points = bcs->GetN();
    double x = 0;
    double x_err = 0;
    double y = 0;
    double y_err = 0;

    std::vector<double> energies = {};
    std::vector<double> energy_spreads = {};

    std::vector<double> chi2 = {};
    std::vector<double> xsec = {};
    std::vector<double> xsec_err = {};
    std::vector<double> pull = {};
    std::vector<double> pull_err = {};

    for(int i = 0; i < N_points; i++)
    {
        x = bcs->GetPointX(i);
        y = bcs->GetPointY(i);
        x_err = bcs->GetErrorX(i);
        y_err = bcs->GetErrorY(i);

        energies.push_back(x);
        energy_spreads.push_back(x_err);
        xsec.push_back(y);
        xsec_err.push_back(y_err);
        pull.push_back((y - f_bcs->Eval(x))/y_err);
        pull_err.push_back(0);

        chi2.push_back((y - f_bcs->Eval(x))/y_err * (y - f_bcs->Eval(x))/y_err);
    }

    std::vector<double> zeroes(100, 0.0);

    auto gr_bcs = new TGraphErrors(energies.size(), energies.data(), xsec.data(), zeroes.data(), xsec_err.data());
    auto gr_pull = new TGraphErrors(energies.size(), energies.data(), pull.data(), zeroes.data(), pull_err.data());

    auto chi2_val = std::accumulate(chi2.begin(), chi2.end(), 0.);
    auto ndf_val = chi2.size() - n_params;
    std::stringstream chi2_stream;
    chi2_stream << "#chi^{2}/ndf = " << std::fixed << std::setprecision(1) << chi2_val  << std::setprecision(0)<< "/" << ndf_val;
    std::string chi2_str = chi2_stream.str();

    std::cout << chi2_str << std::endl; 
    gr_bcs->SetName("data");
    gr_bcs->GetXaxis()->SetTitle("E_{cms}, GeV");
    gr_bcs->GetYaxis()->SetTitle("#sigma_{Born}, nb");
    
    gr_bcs->GetYaxis()->SetRangeUser(2, 7e3);
    gr_bcs->GetXaxis()->SetRangeUser(min_x, max_x);
    
    gr_bcs->GetYaxis()->SetTitleSize(gStyle->GetTitleSize("y") / top_pad_scale);
    gr_bcs->GetYaxis()->SetLabelSize(gStyle->GetLabelSize("y") / top_pad_scale);
    gr_bcs->GetXaxis()->SetLabelSize(0.0);

    gr_bcs->GetYaxis()->SetNdivisions(505);
    gr_bcs->GetYaxis()->SetTitleOffset(0.9);


    TLatex lat; 
    double x_label = 0.2;
    double y_label = 0.85;
    lat.SetNDC();
    lat.SetTextFont(42);
    lat.SetTextSize(gStyle->GetTitleSize("y") / top_pad_scale);

    gr_bcs->DrawClone("AP");
    f_bcs->SetLineColor(kBlack);
    f_bcs->DrawClone("same");
    std::string description = "#font[72]{CMD-3} #chi^{2}/ndf = 9.2/13";
    lat.DrawLatex(x_label, y_label, description.c_str());
    top_pad->SetLogy();

    bottom_pad->cd();
    gr_pull->SetName("delta");
    gr_pull->GetXaxis()->SetRangeUser(min_x, max_x);

    gr_pull->GetXaxis()->SetTitle("E_{cms}, GeV");

    gr_pull->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("x") / bottom_pad_scale);
    gr_pull->GetYaxis()->SetLabelSize(gStyle->GetLabelSize("y") / bottom_pad_scale);

    gr_pull->GetXaxis()->SetLabelOffset(gStyle->GetLabelOffset("x") / bottom_pad_scale);

    gr_pull->GetXaxis()->SetTitleSize(gStyle->GetTitleSize("x") / bottom_pad_scale);
    gr_pull->GetYaxis()->SetTitleSize(gStyle->GetTitleSize("y") / bottom_pad_scale);

    gr_pull->GetXaxis()->SetTitleOffset(1.2);
    gr_pull->GetYaxis()->SetTitleOffset(0.6);

    gr_pull->SetFillColor(kGray);
    TLine zero_pull(min_x, 0, max_x, 0);
    zero_pull.SetLineWidth(2);

    gr_pull->GetYaxis()->SetTitle("Pull");
    gr_pull->GetYaxis()->SetNdivisions(505);
    gr_pull->DrawClone("AB");
    zero_pull.DrawClone("same");
    
    // canvas->SetCanvasSize(1000,1000);
    canvas->SetLeftMargin(0.3);
    canvas->Draw();
    canvas->SaveAs(output_filename.c_str());
    return 0;
}