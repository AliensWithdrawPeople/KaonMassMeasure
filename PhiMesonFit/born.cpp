#include <complex>
#include <functional>

double SigmaBorn(double s, double m_phi, double gamma_phi, double B_phi_gamma_phi_ee) {
    auto pi = 3.14159265359;
    auto alpha = 1 / 137.035999;
    
    double m_K = 497.614;
    double m_pi = 139.570;
    
    const double m_rho = 775.26;
    const double m_omega = 782.66;
    const double gamma_rho = 1; // todo
    const double gamma_omega = 1; // todo

    // K^0 momentum
    auto p_K = [&m_K](double s){return std::sqrt(std::complex<double>(s - 4 * m_K * m_K, 0)); };

    // Propagator 
    auto D = [](double s, double m_V, std::function<double(double)> Gamma_V) { return m_V * m_V - s - std::complex<double>(0.0, Gamma_V(s)); };

    auto Gamma_rho = [&gamma_rho, &m_rho, &m_pi](double s) {
        return gamma_rho * (m_rho * m_rho / s) * std::pow((s / 4 - m_pi * m_pi) / (m_rho * m_rho / 4 - m_pi * m_pi), 3./2.);
    };

    auto Gamma_omega = [&gamma_omega, &m_omega, &m_pi, &m_K](double s) { 
        auto B_2pi = 1.53e-2;   
        auto B_pi0gamma = 8.35e-2;   
        auto B_3pi = 89.2e-2;   

        auto m_pi0 = 134.9768;
        
        auto part_2pi = B_2pi * (m_omega * m_omega / s) * std::pow((s / 4 - m_pi * m_pi) / (m_omega * m_omega / 4 - m_pi * m_pi), 3./2.);
        auto part_pi0gamma = B_pi0gamma * std::pow((std::sqrt(s) * (1 - m_pi0 * m_pi0 / s) / 2) / (m_omega * (1 - m_pi0 * m_pi0 / m_omega / m_omega) / 2), 3);
        auto part_3pi = B_3pi * std::sqrt(s) / m_omega * Phi_3pi; // Phi_3pi todo!
        return gamma_omega * (part_2pi + part_pi0gamma + part_3pi);
    };

    auto Gamma_phi = [&gamma_phi, &m_phi, &m_pi, &m_K](double s) {
        auto B_2pi = 1.53e-2;   
        auto B_eta_gamma = 1.3e-2;   
        auto B_3pi = 89.2e-2;   
        auto B_KpKm = 49.1e-2;
        auto B_KsKl = 33.9e-2;
        auto m_pi0 = 134.9768;
        
        auto part_2pi = B_2pi * (m_omega * m_omega / s) * std::pow((s / 4 - m_pi * m_pi) / (m_omega * m_omega / 4 - m_pi * m_pi), 3./2.);
        auto part_pi0gamma = B_pi0gamma * std::pow((std::sqrt(s) * (1 - m_pi0 * m_pi0 / s) / 2) / (m_omega * (1 - m_pi0 * m_pi0 / m_omega / m_omega) / 2), 3);
        auto part_3pi = B_3pi * std::sqrt(s) / m_omega * Phi_3pi; // Phi_3pi todo!
        return gamma_omega * (part_2pi + part_pi0gamma + part_3pi);
    };

    auto g_phi_gamma = std::sqrt(3 * std::pow(m_phi, 3) * B_phi_gamma_phi_ee / 4 / pi / alpha);
    auto g_phi_KK = std::sqrt(6 * pi * m_phi * m_phi * gamma_phi / std::pow(p_K(m_omega * m_omega), 3) );
    auto g_phi = g_phi_gamma * g_phi_KK;
    
    auto g_rho = std::sqrt(3 * std::pow(m_rho, 3) * 7.04 / 4 / pi / alpha) * (g_phi_KK / std::sqrt(2));
    auto g_omega = std::sqrt(3 * std::pow(m_omega, 3) * 0.60 / 4 / pi / alpha) * (-g_phi_KK / std::sqrt(2));
    
    auto rho_part = g_rho / D(s, m_rho, Gamma_rho);
    auto omega_part = g_omega / D(s, m_omega, Gamma_omega);
    auto phi_part = g_phi / D(s, m_phi, Gamma_phi);


    return 0;
}