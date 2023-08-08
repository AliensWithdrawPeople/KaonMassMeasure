#ifndef MassFunc_h
#define MassFunc_h

#include "misc.h"

#include <math.h>

namespace FullRecMassFunc
{
    // Mass of pi meson in MeV/c^2.
    constexpr double mass_pi = 139.57039; 

    /// @brief Variables of Ks mass function (full reconstruction method).
    enum Var {psi, energy, Y};

    /// @brief Evaluate Ks mass function
    /// @param dpsi Angle between two pions tracks in radians;
    /// @param energy Beam energy in MeV
    /// @param Y Ration of pions' momentums (Y = p_{pi^+} / p_{pi^-});
    /// @return Value of Ks mass function (full reconstruction method).
    static double Eval(double dpsi, double energy, double Y)
    {
        auto E_2 = energy * energy; 
        auto eta_2 = std::pow((1 - Y * Y) / (1 + Y * Y), 2); 
        
        auto part1 = 1 / eta_2 * (1 + cos(dpsi) * sqrt(1 - eta_2));
        auto part2 = 1 - sqrt(1 - eta_2 * (1 - 4 * mass_pi * mass_pi / E_2));
        auto tot = E_2 * (1 - part1 * part2);
        return sqrt(tot);
    }

    /// @brief Evaluate numerical partial derivative of Ks Mass function with respect to var variable
    /// in point (dpsi, energy, Y) of corresponding order.
    /// @param var Evaluate derivative with respect to this variable;
    /// @param dpsi Angle between two pions tracks in radians;
    /// @param energy Beam energy in MeV
    /// @param Y Ration of pions' momentums (Y = p_{pi^+} / p_{pi^-});
    /// @param order Order of partial derivative.
    /// @return Numerical partial derivative of Ks Mass function (full reconstruction method).
    double Derivative(Var var, double dpsi, double energy, double Y, int order)
    {
        double res = 0;
        double dx = 1e-6;
        
        switch (var)
        {
        case Var::psi:
            res = misc::derivative( [&energy, &Y](double psi){ 
                return Eval(psi, energy, Y); 
                }, dpsi, dx, order);
            break;

        case Var::Y:
            res = misc::derivative( [&energy, &dpsi](double Y_val){ 
                return Eval(dpsi, energy, Y_val); 
                }, Y, dx, order);
            break;
        
        case Var::energy:
            res = misc::derivative( [&dpsi, &Y](double E){ 
                return Eval(dpsi, E, Y); 
                }, energy, dx, order);
            break;
        
        default:
            res = misc::derivative( [&energy, &Y](double psi){ 
                    return Eval(psi, energy, Y); 
                    }, dpsi, dx, order);
            break;
        }   

        return res; 
    }
}

#endif