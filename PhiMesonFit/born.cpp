#include <complex>
#include <functional>
#include <vector>
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

double FAS_ASPO(double TwoE)
{
    double Polinom(int, double, double*);

    const int K_max = 4;
    int I_mark = 100;
    int Mode = 1;

    bool Calc = true;

    double A[K_max] = {-6.12394622, 25.0341405, -34.1311022, 15.5413717};
    double B[K_max] = {5.29354148, -7.90990714, -2.26007613, 5.21453902};
    double C[K_max] = {-0.56436115415016, 2.69953793, -4.32966739, 2.33116866};
    double D[K_max] = {-0.0548334238,  0.31600391, -0.609523718, 0.393667808};

    double A1[K_max] = {-4.91624401 , 19.8655606, -26.9136128 , 12.2412286};
    double D1[K_max] = {-0.00794774189,0.0522269164,-0.114526409 , 0.0838126536};

    double LM[6] = {1.1,0.875,0.75,0.62,0.52,0.46};

    double POL, TEMP1, TEMP2; 

    if ( Calc ) {
        TEMP2 = D[1]*LM[4]+D[2]*LM[4]*LM[4]+D[3]*LM[4]*LM[4]*LM[4];
        TEMP1 = D1[0]+D1[1]*LM[4]+D1[2]*LM[4]*LM[4]+D1[3]*LM[4]*LM[4]*LM[4];
        D[0]  = TEMP1 - TEMP2;

        TEMP2 = C[1]*LM[3]+C[2]*LM[3]*LM[3]+C[3]*LM[3]*LM[3]*LM[3];
        TEMP1 = D[0]+D[1]*LM[3]+D[2]*LM[3]*LM[3]+D[3]*LM[3]*LM[3]*LM[3];
        C[0]  = TEMP1 - TEMP2;

        TEMP2 = A1[1]*LM[2]+A1[2]*LM[2]*LM[2]+A1[3]*LM[2]*LM[2]*LM[2];
        TEMP1 = C[0]+C[1]*LM[2]+C[2]*LM[2]*LM[2]+C[3]*LM[2]*LM[2]*LM[2];
        A1[0] = TEMP1 - TEMP2;

        TEMP1 = A1[0]+A1[1]*LM[1]+A1[2]*LM[1]*LM[1]+A1[3]*LM[1]*LM[1]*LM[1];
        TEMP2 = A[1]*LM[1]+A[2]*LM[1]*LM[1]+A[3]*LM[1]*LM[1]*LM[1];
        A[0]  = TEMP1 - TEMP2;
        
        TEMP1 = A[0] + A[1]*LM[0]+A[2]*LM[0]*LM[0]+A[3]*LM[0]*LM[0]*LM[0];
        TEMP2 = B[1]*LM[0]+B[2]*LM[0]*LM[0]+B[3]*LM[0]*LM[0]*LM[0];
        B[0]  = TEMP1 - TEMP2;

        Calc = false;
    }
    POL = 0.0;
    if ( TwoE >= LM[0] ) {
        for( int i = 0; i < K_max; i++ ) {
        if ( Mode == 1 ) POL += std::pow(TwoE,(double)i)*B[i];
        if ( Mode == 2 ) POL += std::pow(TwoE,(double)(i-1))*i*B[i];
        }
    }

    if ( TwoE >= LM[1] && TwoE < LM[0] ) {
        for( int i = 0; i < K_max; i++ ) {
        if ( Mode == 1 ) POL += std::pow(TwoE,(double)i)*A[i];
        if ( Mode == 2 ) POL += std::pow(TwoE,(double)(i-1))*i*A[i];
        }
    }

    if ( TwoE >= LM[2] && TwoE < LM[1] ) {
        for( int i = 0; i < K_max; i++ ) {
        if ( Mode == 1 ) POL += std::pow(TwoE,(double)i)*A1[i];
        if ( Mode == 2 ) POL += std::pow(TwoE,(double)(i-1))*i*A1[i];
        }
    }

    if ( TwoE >= LM[3] && TwoE < LM[2] ) {
        for( int i = 0; i < K_max; i++ ) {
        if ( Mode == 1 ) POL += std::pow(TwoE,(double)i)*C[i];
        if ( Mode == 2 ) POL += std::pow(TwoE,(double)(i-1))*i*C[i];
        }
    }

    if ( TwoE >= LM[4] && TwoE < LM[3] ) {
        for( int i = 0; i < K_max; i++ ) {
        if ( Mode == 1 ) POL += std::pow(TwoE,(double)i)*D[i];
        if ( Mode == 2 ) POL += std::pow(TwoE,(double)(i-1))*i*D[i];
        }
    }

    if ( TwoE >= LM[5] && TwoE < LM[4] ) {
        for( int i = 0; i < K_max; i++ ) {
        if ( Mode == 1 ) POL += std::pow(TwoE,(double)i)*D1[i];
        if ( Mode == 2 ) POL += std::pow(TwoE,(double)(i-1))*i*D1[i];
        }
    }
    if ( TwoE < LM[5] ) {
        if ( Mode == 1 ) POL = D1[0]+D1[1]*LM[5]+D1[2]*LM[5]*LM[5]+D1[3]*LM[5]*LM[5]*LM[5];
        if ( Mode == 2 ) POL = 3.0*LM[5]*LM[5]*D1[3]+2.0*LM[5]*D1[2]+D1[1]; 
    }


    double Fval = (POL/0.393728*(0.00749/0.0361478));

    return Fval;
}

double SigmaBorn(double *x, double *par) {
    double s = x[0] * x[0];
    double m_phi = par[0];
    double gamma_phi = par[1];
    double B_phi_gamma_phi_ee = par[2];

    auto pi = 3.14159265359;
    auto alpha = 1 / 137.035999;
    
    double m_K0 = 497.614;
    double m_KC = 493.677;
    double m_pi = 139.570;
    
    const double m_rho = 775.26;
    const double m_omega = 782.66;
    const double gamma_rho = 8.68; // todo
    const double gamma_omega = 149.10; // todo

    // K^0 momentum
    auto p_K = [&m_K0](double s){return std::sqrt(std::complex<double>(s - 4 * m_K0 * m_K0, 0)); };

    // Propagator 
    auto D = [](double s, double m_V, std::function<double(double)> Gamma_V) { return m_V * m_V - s + std::complex<double>(0.0, -std::sqrt(s) * Gamma_V(s)); };

    auto Gamma_rho = [&gamma_rho, &m_rho, &m_pi](double s) {
        return gamma_rho * (m_rho * m_rho / s) * std::pow((s / 4 - m_pi * m_pi) / (m_rho * m_rho / 4 - m_pi * m_pi), 3./2.);
    };

    auto Gamma_omega = [&gamma_omega, &m_omega, &m_pi, &m_K0](double s) { 
        auto B_2pi = 1.53e-2;   
        auto B_pi0gamma = 8.35e-2;   
        auto B_3pi = 89.2e-2;   

        auto m_pi0 = 134.9768;
        
        auto part_2pi = B_2pi * (m_omega * m_omega / s) * std::pow((s / 4 - m_pi * m_pi) / (m_omega * m_omega / 4 - m_pi * m_pi), 3./2.);
        auto part_pi0gamma = B_pi0gamma * std::pow((std::sqrt(s) * (1 - m_pi0 * m_pi0 / s) / 2) / (m_omega * (1 - m_pi0 * m_pi0 / m_omega / m_omega) / 2), 3);
        auto part_3pi = B_3pi * std::sqrt(s) / m_omega * FAS_ASPO(std::sqrt(s) / 1e3) / FAS_ASPO(m_omega / 1e3); // Omega_3pi todo!
        return gamma_omega * (part_2pi + part_pi0gamma + part_3pi);
    };

    auto Gamma_phi = [&gamma_phi, &m_phi, &m_pi, &m_K0, &m_KC](double s) {
        auto B_2pi = 1.53e-2;   
        auto B_eta_gamma = 1.3e-2;   
        auto B_3pi = 89.2e-2;   
        auto B_KpKm = 49.1e-2;
        auto B_KsKl = 33.9e-2;
        auto m_pi0 = 134.9768;
        auto m_eta = 547.862;
        
        auto part_KpKm = B_KpKm * (m_phi * m_phi / s) * std::pow((s / 4 - m_KC * m_KC) / (m_phi * m_phi / 4 - m_KC * m_KC), 3./2.); 
        auto part_KsKl = B_KpKm * (m_phi * m_phi / s) * std::pow((s / 4 - m_K0 * m_K0) / (m_phi * m_phi / 4 - m_K0 * m_K0), 3./2.); 
        auto part_eta_gamma = B_eta_gamma * std::pow((std::sqrt(s) * (1 - m_eta * m_eta / s) / 2) / (m_phi * (1 - m_eta * m_eta / m_phi / m_phi) / 2), 3);
        auto part_3pi = B_3pi * std::sqrt(s) / m_phi * FAS_ASPO(std::sqrt(s) / 1e3) / FAS_ASPO(m_phi / 1e3); // Phi_3pi todo!
        return gamma_phi * (part_KpKm + part_KsKl + part_eta_gamma + part_3pi);
    };

    auto g_phi_gamma = std::sqrt(3 * std::pow(m_phi, 3) / 4 / pi / alpha);
    auto g_phi_KK = std::sqrt(6 * pi * m_phi * m_phi * gamma_phi * B_phi_gamma_phi_ee / std::pow(p_K(m_phi * m_phi), 3) );
    auto g_phi = g_phi_gamma * g_phi_KK;
    
    auto g_phi_KK_0 = std::sqrt(6 * pi * m_phi * m_phi * gamma_phi * 0.339 / std::pow(p_K(m_phi * m_phi), 3) );
    auto g_rho = std::sqrt(3 * std::pow(m_rho, 3) * 7.04 / 4. / pi / alpha) * (-g_phi_KK_0 / std::sqrt(2));
    auto g_omega = std::sqrt(3 * std::pow(m_omega, 3) * 0.60 / 4. / pi / alpha) * (g_phi_KK_0 / std::sqrt(2));
    
    auto rho_part = g_rho / D(s, m_rho, Gamma_rho);
    auto omega_part = g_omega / D(s, m_omega, Gamma_omega);
    auto phi_part = g_phi / D(s, m_phi, Gamma_phi);

    auto tot = rho_part + omega_part + phi_part;
    auto val = 8 * pi * alpha / (3 * std::pow(s, 5./2.)) * std::pow(p_K(m_phi * m_phi).real(), 3.) 
            * (tot.real() * tot.real() + tot.imag() * tot.imag());
    // std::cout << std::sqrt(s) << " : " << val << std::endl;
    return val * 1e9;
}

int born()
{
    std::vector<double> energies = {1004.066, 1010.466, 1012.955, 1015.068, 1016.105, 1017.155, 1017.156, 1018.046, 1019.118,  1019.214, 1019.421, 1019.902, 
                1021.222, 1021.309, 1022.078,  1022.744, 1023.264,  1025.320, 1027.956, 1029.090, 1033.907, 1040.028,  1049.864,  1050.862,  1059.947};
    std::vector<double> cross_sections = {6.87, 42.16, 96.74, 219.53, 366.33, 628.15, 624.76, 996.62, 1413.65,  1433.05, 1434.84, 1341.91,  833.20, 
                        807.54, 582.93,  443.71, 377.77,  199.26, 115.93, 96.96, 50.12, 31.27,  16.93,  17.47,  12.09};
    std::vector<double> xsec_err = {0.42, 0.47, 1., 5.02, 3.33, 2.95, 9.89, 4.28,  6.02, 15.03, 18.40, 4.74, 4.89, 10.36, 4.03, 4.38, 
                                    5.31, 4.97, 1.7, 3., 1.26, 1.01, 0.5, 0.94, 0.71};
    std::vector<double> zeroes(100, 0.0);
    
    TCanvas c("canv", "canv");
    TGraphErrors a(energies.size(), energies.data(), cross_sections.data(), zeroes.data(), xsec_err.data());
    auto bb = new TF1("born", SigmaBorn, 1000, 1100, 3);
    bb->SetParameter(0, 1019.457);
    bb->SetParameter(1, 4.240);
    bb->SetParameter(2, 0.428);

    a.Fit(bb, "ME");
    a.DrawClone("AP");
    // bb->DrawClone("same");

    c.DrawClone();
    return 0;
}