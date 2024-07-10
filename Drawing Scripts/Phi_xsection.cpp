#include "TF1.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TFile.h"
#include "TCanvas.h"
#include <TLegend.h>
#include <TRandom.h>

int Phi_xsection()
{
    TCanvas canv("canv");
/************ phi -> KsKl cross section vs energy ************/
    std::vector<double> zeroes(100, 0.0);

    std::vector<double> energy = {1.001859, 1.00596, 1.0096, 1.015724, 1.016808, 1.017914, 1.019056, 1.019912, 1.020916, 1.02207, 1.022896, 1.027728, 1.0338, 1.03979, 1.0498, 1.06};
    std::vector<double> energySpread = {0.0003492, 0.0003625, 0.0003626, 0.0003937, 0.0003473, 0.0003556, 0.0003819, 0.0004158, 0.0004071, 0.0003791, 0.0003914, 0.0003672, 0.00036, 0.00041, 0.000417, 0.000378};

    // With Bkg
    std::vector<double> xsec_vis_with_bkg = {4.72124, 10.16127, 24.95628, 225.61187, 387.14261, 667.53118, 998.52905, 989.26049, 788.19274, 530.25236, 407.61139, 151.73518, 81.20724, 53.87495, 34.92748, 25.41842};
    std::vector<double> xsec_vis_with_bkg_err = {0.197209, 0.18273, 0.576638, 1.509155, 1.439604, 2.117351, 1.900085, 2.073628, 2.662238, 2.300914, 2.195585, 1.226125, 0.870177, 0.734421, 0.538556, 0.484794};
    auto grVisXsec_with_bkg = new TGraphErrors(energy.size(), energy.data(), xsec_vis_with_bkg.data(), energySpread.data(), xsec_vis_with_bkg_err.data());

    // New FF: Const bkg fit
    // std::vector<double> xsec_vis2 = {2.55493, 8.57642, 22.91284, 220.21836, 378.98208, 657.72612, 984.25605, 974.97838, 778.57844, 518.9548, 397.77018, 145.34911, 75.96123, 49.65799, 32.40645, 23.18225};
    // std::vector<double> xsec_vis_err2 = {0.144931, 0.167673, 0.552109, 1.487758, 1.418073, 2.096495, 1.879297, 2.051185, 2.637251, 2.267091, 2.161583, 1.197029, 0.838923, 0.703349, 0.517876, 0.462314};

    // std::vector<double> xsec_vis2 = {2.2188, 8.42657, 23.47782, 223.86866, 384.72266, 667.19337, 999.73643, 992.35409, 790.57739, 529.32899, 407.17639, 150.40073, 78.14517, 51.40992, 31.92956, 22.7146};
    // std::vector<double> xsec_vis_err2 = {0.135943, 0.167415, 0.56196, 1.524775, 1.458682, 2.182089, 1.963416, 2.15065, 2.740329, 2.353383, 2.244232, 1.236502, 0.860266, 0.724327, 0.519018, 0.46168};

    std::vector<double> xsec_vis2 = {3.41088, 8.74448, 23.34787, 223.87477, 385.92215, 667.67413, 1000.62134, 993.56095, 791.73284, 530.64294, 408.19207, 150.85117, 80.21025, 53.24817, 36.41072, 27.77034};
    std::vector<double> xsec_vis_err2 = {0.168654, 0.170589, 0.560372, 1.5248, 1.461933, 2.183167, 1.964766, 2.152644, 2.743444, 2.357467, 2.247874, 1.238595, 0.872728, 0.73802, 0.55603, 0.512258};
    auto grVisXsec2 = new TGraphErrors(energy.size(), energy.data(), xsec_vis2.data(), energySpread.data(), xsec_vis_err2.data());

    // New FF: MC + pol2 fit
    std::vector<double> xsec_vis = {2.50786, 8.68968, 23.62778, 226.54159, 389.94455, 675.32426, 1010.04235, 1001.09009, 800.55697, 533.56961, 411.19664, 150.52218, 78.77594, 51.6739, 33.70135, 23.99367};
    std::vector<double> xsec_vis_err = {0.143587, 0.168791, 0.560805, 1.512828, 1.446984, 2.133892, 1.916838, 2.092197, 2.694341, 2.310825, 2.20793, 1.220632, 0.855789, 0.718334, 0.528582, 0.470581};
    

    std::vector<std::vector<double>> xsec_vis_rand = {
        {2.57247, 8.26814, 23.7833, 214.22, 371.415, 661.6, 1002.1, 1005.63, 767.857, 514.892, 400.077, 143.449, 74.9358, 49.0609, 31.3562, 22.4868},
        {2.53891, 8.21763, 23.6718, 213.696, 373.168, 659.13, 1005.98, 1004.85, 767.48, 514.691, 396.675, 146.836, 74.7731, 51.4805, 31.9079, 22.9335},
        {2.65096, 8.35693, 24.1935, 216.075, 370.053, 661.301, 1003.24, 1005.54, 764.776, 515.852, 397.117, 144.543, 73.8681, 51.3029, 31.7333, 23.0143},
        {2.25664, 8.14144, 22.2421, 216.255, 370.938, 662.34, 1004.25, 1003.54, 766.256, 514.037, 391.492, 143.172, 73.9338, 49.9087, 32.8378, 22.9119},
        {2.10354, 8.29098, 24.335, 212.946, 370.61, 660.946, 1002.45, 1004.99, 767.042, 515.535, 398.914, 145.173, 74.9507, 50.9194, 32.0915, 22.8704},
        {2.56283, 8.5742, 23.892, 212.933, 370.64, 660.141, 1002.63, 1003.12, 770.437, 516.511, 393.859, 146.238, 75.1617, 49.9506, 31.7336, 23.2989},
        {2.36254, 8.59378, 23.3287, 214.386, 371.825, 663.649, 1001.15, 1006.41, 767.016, 514.436, 401.07, 144.669, 74.899, 49.7429, 32.3396, 22.3556},
        {2.1385, 8.41783, 22.2859, 212.741, 369.358, 661.476, 1002.67, 1007.54, 761.996, 515.402, 395.479, 144.256, 75.656, 48.9775, 32.5585, 22.9618},
        {2.58969, 8.40348, 23.341, 212.913, 368.226, 664.675, 1001.84, 1005.31, 763.226, 513.022, 397.493, 144.108, 75.9234, 49.852, 31.6514, 21.9755},
        {2.38908, 8.71652, 23.1769, 214.035, 371.545, 661.188, 1002.23, 1005.13, 769.987, 514.326, 397.18, 145.398, 74.6333, 49.823, 31.5205, 23.5013},
    };

    // New FF + track eff accounted for (514 > E > 505): Const bkg fit
    std::vector<double> xsec_vis2_tr_eff = {2.44829, 9.03045, 23.34802, 224.52051, 386.89303, 676.26567, 1023.60277, 1018.18031, 814.36251, 534.28396, 418.82626, 150.77984, 85.33161, 55.83057, 38.1369, 28.22991};
    std::vector<double> xsec_vis_err2_tr_eff = {0.158988, 0.209823, 0.793483, 2.400864, 2.154328, 3.484861, 2.86498, 3.368048, 4.563764, 4.240477, 4.733099, 3.64503, 3.408002, 3.018571, 2.267439, 2.017582};
    auto grVisXsec2_tr_eff = new TGraphErrors(energy.size(), energy.data(), xsec_vis2_tr_eff.data(), energySpread.data(), xsec_vis_err2_tr_eff.data());
    TFile *top = new TFile("C:/work/Science/BINP/Kaon Mass Measure/PhiMesonFit/vcs/vcs_kskl_track_rec_eff.root", "recreate");
    grVisXsec2_tr_eff->GetYaxis()->SetTitle("#sigma_{vis}, nb");
    grVisXsec2_tr_eff->GetXaxis()->SetTitle("E_{avg}, MeV");
    grVisXsec2_tr_eff->SetName("vcs");
    grVisXsec2_tr_eff->SetMarkerColor(kRed);
    grVisXsec2_tr_eff->Write();
    top->Write();
    top->Save();

    // auto grVisXsec = new TGraphErrors(energy.size(), energy.data(), xsec_vis.data(), energySpread.data(), xsec_vis_err.data());
    // for(int i = 0; i < xsec_vis_rand.size(); i++)
    // {
    //     auto sample = xsec_vis_rand[i];
    //     TFile *top = new TFile(("C:/work/Science/BINP/Kaon Mass Measure/PhiMesonFit/vcs/ToyMC/vcs_kskl_rand_sample" + std::to_string(i) + ".root").c_str(), "recreate");
    //     auto grVisXsec = new TGraphErrors(energy.size(), energy.data(), sample.data(), energySpread.data(), xsec_vis_err.data());
    //     grVisXsec->GetYaxis()->SetTitle("#sigma_{vis}, nb");
    //     grVisXsec->GetXaxis()->SetTitle("E_{avg}, MeV");
    //     grVisXsec->SetName("vcs");
    //     grVisXsec->SetMarkerColor(kRed);
    //     grVisXsec->DrawClone("AP");
        
    //     grVisXsec->Write();
    //     top->Write();
    //     top->Save();
    // }
    // auto grVisXsec = new TGraphErrors(energy.size(), energy.data(), xsec_vis_rand.data(), energySpread.data(), xsec_vis_err.data());
    // grVisXsec2_tr_eff->GetYaxis()->SetTitle("#sigma_{vis}, nb");
    // grVisXsec2_tr_eff->GetXaxis()->SetTitle("E_{avg}, MeV");
    // grVisXsec2_tr_eff->SetName("vcs");
    // grVisXsec2_tr_eff->SetMarkerColor(kRed);
    // grVisXsec2->SetMarkerColor(kBlue);
    // grVisXsec2_tr_eff->DrawClone("AP");
    // grVisXsec2->DrawClone("P same");
    
    // grVisXsec2_tr_eff->Write();
    // top->Write();
    // top->Save();


    std::vector<double> sample = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, -1, 10};
    std::vector<double> mass = {1019.442,  1019.437, 1019.437, 1019.424, 1019.437, 1019.440, 1019.439, 1019.432, 1019.435, 1019.440};
    std::vector<double> mass_err = {0.005, 0.006, 0.005, 0.005, 0.007, 0.005, 0.005, 0.005, 0.005, 0.005};
    std::vector<double> width = {4.29, 4.29, 4.29, 4.27, 4.29, 4.3, 4.27, 4.28, 4.30, 4.32};
    std::vector<double> width_err = {0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010};
    auto gr_toyMC_mass = new TGraphErrors(sample.size(), sample.data(), mass.data(), zeroes.data(), mass_err.data());
    auto gr_toyMC_width = new TGraphErrors(sample.size(), sample.data(), width.data(), zeroes.data(), width_err.data());

    // gr_toyMC_mass->DrawClone("AP");

/************ Efficiency vs energy ************/
    // Outdated (old Kaon FF in MCGPJ) efficiency.
    std::vector<double> efficiency_old = {0.25922, 0.25138, 0.24396, 0.23382, 0.23255, 0.23099, 0.23027, 0.22983, 0.22779, 0.22737, 0.22637, 0.22331, 0.21986, 0.21707, 0.21214, 0.20391};
    std::vector<double> efficiency_err_old = {0.00044, 0.000436, 0.000431, 0.00045, 0.000424, 0.000423, 0.000423, 0.000451, 0.000298, 0.000421, 0.00042, 0.000437, 0.000416, 0.000414, 0.000531, 0.000524};
    
    auto grEff_old = new TGraphErrors(energy.size(), energy.data(), efficiency_old.data(), zeroes.data(), efficiency_err_old.data());
    grEff_old->SetMarkerColor(kRed);
    grEff_old->SetName("grEff_old");

    // Up-to-date efficiency.
    std::vector<double> efficiency = {0.25945, 0.25114, 0.2444, 0.2339, 0.23235, 0.23133, 0.23111, 0.22996, 0.22806, 0.22772, 0.22709, 0.22395, 0.22032, 0.21754, 0.21277, 0.2047};
    std::vector<double> efficiency_err = {0.00048, 0.00047, 0.00047, 0.00049, 0.00049, 0.00028, 0.00025, 0.00025, 0.00047, 0.00047, 0.00047, 0.00044, 0.00061, 0.0006, 0.00059, 0.00058};
    
    auto grEff = new TGraphErrors(energy.size(), energy.data(), efficiency.data(), zeroes.data(), efficiency_err.data());
    grEff->SetName("grEff");

    grEff->GetYaxis()->SetTitleOffset(1.2);
    grEff->GetYaxis()->SetTitle("#varepsilon_{MC}");
    grEff->GetXaxis()->SetTitle("E_{cms}, GeV");
    grEff->SetName("eff");
    grEff->SetMarkerColor(kBlue);
    grEff->SetTitle("Selection efficiency of e^{+}e^{-}#rightarrowK_{S}K_{L} events");

    // auto legend = new TLegend(0.1,0.7,0.48,0.9);
    // legend->AddEntry(grEff_old,"With old Kaon FF in MCGPJ","ep");
    // legend->AddEntry(grEff,"With new Kaon FF in MCGPJ","ep");

    grEff->DrawClone("AP");
    // grEff_old->DrawClone("P same");
    // legend->DrawClone();

/************ Background vs energy ************/
    // Const bkg; chi2 fit
    // std::vector<double> bckg = {0.5448, 0.32112, 0.61351, 0.71244, 1.23553, 1.67868, 2.71171, 2.57525, 2.08286, 1.31758, 1.41573, 0.88321, 0.8125,  0.76288, 0.71542, 0.84304};
    // std::vector<double> bckg_err = {0.03429, 0.01627, 0.04449, 0.0373, 0.03013, 0.04185, 0.03296, 0.03633, 0.04488, 0.04353, 0.0528,  0.04186, 0.03907, 0.03983, 0.03503, 0.03965};
    
    // Const bkg, only right sideband (range = (540, 576)); NLL fit
    std::vector<double> bckg2 = {0.69951, 0.47594, 0.42957, 0.55441, 0.86295, 0.77104, 1.04199, 1.10664, 1.00505, 0.88549, 0.81114, 1.00279, 1.04139, 1.02398, 1.2634, 1.2652};
    std::vector<double> bckg_err2 = {0.03887, 0.01981, 0.03722, 0.0329, 0.02517, 0.02833, 0.02039, 0.02377, 0.03114, 0.03567, 0.03994, 0.04461, 0.04424, 0.04616, 0.04659, 0.0486};
    auto grBckg_const = new TGraphErrors(energy.size(), energy.data(), bckg2.data(), zeroes.data(), bckg_err2.data());
    grBckg_const->SetMarkerColor(kRed);

    std::vector<double> bckg_const_4pi = {0.642, 0.443, 0.35795, 0.46776, 0.81084, 0.75258, 0.97121, 1.03651, 0.93662, 0.88549, 0.7588, 0.5482, 0.7449, 0.64325, 0.70352, 0.59548};
    std::vector<double> bckg_const_4pi_err = {0.03723, 0.01911, 0.03397, 0.03021, 0.0244, 0.02799, 0.01969, 0.023, 0.03006, 0.03567, 0.03863, 0.03296, 0.0374, 0.03657, 0.03474, 0.03331};
    auto grBckg_const_4pi = new TGraphErrors(energy.size(), energy.data(), bckg_const_4pi.data(), zeroes.data(), bckg_const_4pi_err.data());
    grBckg_const_4pi->SetMarkerColor(kBlack);

    // MC + bkg fit
    std::vector<double> bckg = {0.60188, 0.40148, 0.39382, 0.23414, 0.23994, 0.21656, 0.29378, 0.20309, 0.1284, 0.35398, 0.41952, 0.56358, 0.7418, 0.77178, 0.96882, 1.05519};
    std::vector<double> bckg_err = {0.03891, 0.02062, 0.04516, 0.05381, 0.02114, 0.09759, 0.00026, 0.00369, 0.04822, 0.04436, 0.00083, 0.04753, 0.05373, 0.05456, 0.05113, 0.05295};

    auto grBckg = new TGraphErrors(energy.size(), energy.data(), bckg.data(), zeroes.data(), bckg_err.data());

    std::vector<double> energy_bkg_calculated = {1.0105, 1.013, 1.0151, 1.0161, 1.0172, 1.0172, 1.018, 1.0191, 1.0192, 1.0194, 1.0199, 1.0212, 1.0213, 1.0221, 1.0227, 1.0233, 1.0253, 1.028, 1.0291, 1.0339, 1.04, 1.0499, 1.0509, 1.06};
    std::vector<double> bkg_calculated = {0.35, 0.365, 0.399, 0.442, 0.516, 0.515, 0.623, 0.748, 0.75, 0.754, 0.746, 0.71, 0.708, 0.64, 0.597, 0.575, 0.521, 0.487, 0.48, 0.461, 0.556, 0.548, 0.548, 0.545};
    auto grBkg_calculated = new TGraphErrors(energy_bkg_calculated.size(), energy_bkg_calculated.data(), bkg_calculated.data(), zeroes.data(), zeroes.data());
    grBkg_calculated->SetName("grBkg_calculated");

    grBckg_const->GetYaxis()->SetTitle("#sigma_{bkg}, nb");
    grBckg_const->GetXaxis()->SetTitle("E_{cms}, GeV");
    grBckg->SetName("bckg");
    grBckg->SetMarkerColor(kBlue);

    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->AddEntry(grBckg_const, "pol0 bkg from right sideband","ep");
    legend->AddEntry(grBckg_const_4pi, "pol0 bkg from right sideband + 4pi bkg","ep");
    legend->AddEntry(grBckg,"MC + pol2 bkg fit","ep");
    legend->AddEntry(grBkg_calculated,"Estimated Bkg: K^{+}K^{-}, #pi^{+}#pi^{-}#pi^{+}#pi^{-}, #pi^{+}#pi^{-}#pi^{0}#pi^{0}","l");
    legend->SetBorderSize(4);
    legend->SetLineColor(kBlack);
    gStyle->SetLegendBorderSize(1);

    grBckg_const->DrawClone("AP");
    grBckg->DrawClone("P same");
    grBckg_const_4pi->DrawClone("P same");
    grBkg_calculated->DrawClone("same");

    legend->DrawClone("same l");


    canv.DrawClone();
    return 0;
}