#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TCanvas.h>
#include <TRatioPlot.h>

#include <vector>


int CountEvents()
{
    std::string energy = "509.5";

    TCanvas canv("canv", "canv", 800, 600);
    auto exp = TFile::Open(("C:/work/Science/BINP/Kaon Mass Measure/tr_ph/PhiXSection/Kn" + energy + ".root").c_str());
    auto mass_exp = new TH1D("mass_exp", "mass_exp", 50, 460, 540);
    auto tr_exp = exp->Get<TTree>("Kn");
    tr_exp->Draw("mass >> mass_exp", "", "goff");

    auto MC = TFile::Open(("C:/work/Science/BINP/Kaon Mass Measure/tr_ph/PhiXSection/MC/MC_Kn" + energy + ".root").c_str());
    auto mass_MC = new TH1D("mass_MC", "mass_MC", 100, 460, 540);
    auto tr_MC = MC->Get<TTree>("Kn_MC");
    tr_MC->Draw("mass >> mass_MC", "", "goff");

    auto sig_func = new TF1("sig_func", "[0] / x[0] / sqrt(2 * 3.14 * ([1]*[1])) * exp(-(log(x[0]) - [2]) * (log(x[0]) - [2]) / 2 / ([1]*[1])) + gaus(3) + pol3(6)", 460, 520);

    auto sig_bkg_func = new TF1("sig_bkg_func", "[0] / ([1] + [2] + [3] + [4]) * (\
        [1] / sqrt(2 * 3.14 * ([5]*[5] + [14]*[14])) * exp(-(x[0] - [6] - [13]) * (x[0] - [6] - [13])/2/([5]*[5] + [14]*[14])) + \
        [2] / sqrt(2 * 3.14 * ([7]*[7] + [14]*[14])) * exp(-(x[0] - [8] - [13]) * (x[0] - [8] - [13])/2/([7]*[7] + [14]*[14])) + \
        [3] / sqrt(2 * 3.14 * ([9]*[9] + [14]*[14])) * exp(-(x[0] - [10] - [13]) * (x[0] - [10] - [13])/2/([9]*[9] + [14]*[14])) + \
        [4] / sqrt(2 * 3.14 * ([11]*[11] + [14]*[14])) * exp(-(x[0] - [12] - [13]) * (x[0] - [12] - [13])/2/([11]*[11] + [14]*[14]))\
    ) + [15] + [16] * (x[0] - 496.) + [17] * (x[0] - 496.) * (x[0] - 496.)", 460, 540);

    // auto sig_bkg_func = new TF1("sig_bkg_func", "[0] / sqrt(2 * 3.14 * ([1]*[1] + [3]*[3])) * exp(-(x[0] - [2] - [4]) * (x[0] - [2] - [4])/2/([1]*[1] + [3]*[3])) + [5]", 460, 540);

    double pars_alt[11] = {1000, 0.005, 6, 10, 500, 1, 0, 0, 0, 0, 0};
    double pars_alt2[6] = {1000, 2, 497.6, -1, 0.5};

    double pars[18] = {1000, 100, 50, 10, 1, 1, 500, 1, 490, 1, 510, 1, 520, 0, 0, 0, 0, 0};
    sig_bkg_func->SetParameters(pars);
    sig_bkg_func->FixParameter(13, 0.0);
    sig_bkg_func->FixParameter(14, 0.0);

    sig_bkg_func->FixParameter(15, 0.0);
    sig_bkg_func->FixParameter(16, 0.0);
    sig_bkg_func->FixParameter(17, 0.0);

    auto res_MC = mass_MC->Fit(sig_bkg_func, "SQLME0");
    res_MC = mass_MC->Fit(sig_bkg_func, "SQLME0");
    res_MC = mass_MC->Fit(sig_bkg_func, "SQLME");
    std::cout << "MC chi2/ndf = " << res_MC->Chi2() << "/" << res_MC->Ndf() << std::endl;
    std::cout << "mass_MC signal integral = " << sig_bkg_func->Integral(460, 540) << std::endl;
    mass_MC->DrawClone("P E");
    // mass_exp->DrawNormalized("same", 187133);


    for(int i = 1; i < 13; i++)
    { sig_bkg_func->FixParameter(i, res_MC->Parameter(i)); }

    sig_bkg_func->ReleaseParameter(13);
    sig_bkg_func->ReleaseParameter(14);

    // sig_bkg_func->ReleaseParameter(15);
    // sig_bkg_func->ReleaseParameter(16);
    // sig_bkg_func->ReleaseParameter(17);

    // auto res_exp = mass_exp->Fit(sig_bkg_func, "SQLME");
    // mass_exp->DrawClone("PE");
    // mass_exp->GetXaxis()->SetTitle("x");
    // // auto ratio_plot = new TRatioPlot(mass_exp);
    // // ratio_plot->Draw();
    // // ratio_plot->GetLowerRefYaxis()->SetTitle("ratio");
    // // ratio_plot->GetUpperRefYaxis()->SetTitle("entries");
    // // ratio_plot->GetUpperPad()->SetLogy();

    // sig_bkg_func->FixParameter(15, 0.0);
    // sig_bkg_func->FixParameter(16, 0.0);
    // sig_bkg_func->FixParameter(17, 0.0);

    // sig_bkg_func->SetLineColor(kBlue);
    // sig_bkg_func->SetLineStyle(kDashed);
    // sig_bkg_func->DrawClone("same");

    // std::cout << "mass_exp signal integral = " << sig_bkg_func->Integral(460, 540) << std::endl;
    // std::cout << "mass_exp signal integral = " << sig_bkg_func->IntegralError(460, 540) << std::endl;


// **************************** Pars_alt2 ***************************
    // sig_bkg_func->SetParameters(pars_alt2);
    // sig_bkg_func->FixParameter(3, 0.0);
    // sig_bkg_func->FixParameter(4, 0.0);
    // sig_bkg_func->FixParameter(5, 0.0);

    // auto res_MC = mass_MC->Fit(sig_bkg_func, "SQLME");
    // res_MC = mass_MC->Fit(sig_bkg_func, "SQLME");
    // res_MC = mass_MC->Fit(sig_bkg_func, "SQLME");
    // std::cout << "mass_MC signal integral = " << sig_bkg_func->Integral(460, 540) << std::endl;
    // mass_MC->DrawClone();

    // for(int i = 1; i < 13; i++)
    // { sig_bkg_func->FixParameter(i, res_MC->Parameter(i)); }

    // sig_bkg_func->ReleaseParameter(3);
    // sig_bkg_func->ReleaseParameter(4);
    // sig_bkg_func->ReleaseParameter(5);


    // auto res_exp = mass_exp->Fit(sig_bkg_func, "SQLME");
    // mass_exp->DrawClone();

    // sig_bkg_func->FixParameter(5, 0.0);
    // std::cout << "mass_exp signal integral = " << sig_bkg_func->Integral(460, 540) << std::endl;
    // std::cout << "mass_exp signal integral = " << sig_bkg_func->IntegralError(460, 540) << std::endl;

    canv.SetLogy();
    canv.DrawClone();
    return 0;
}