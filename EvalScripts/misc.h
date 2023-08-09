#ifndef misc_h
#define misc_h

#include <functional>
#include <vector>
#include <numeric>

namespace misc {
    std::vector<int> pascalsTriangle(int n) 
    {
        std::vector<int> triangle(n + 1, 1);
        for (int i = 1; i < n; i++) 
        {
            for (int j = i; j >= 1; j--) 
            { triangle[j] += triangle[j - 1]; }
        }
        return triangle;
    }

    // Derivative with Richardson extrapolation.
    double derivative(std::function<double(double)> func, double x, double dx, int n) 
    {
        std::vector<int> binomial = pascalsTriangle(n);
        auto deriv = [&binomial, &func, &x, &n](double dx) {
            double result = 0.0;
            for (int i = 0; i <= n; i++) 
            { result += pow(-1, i) * binomial[i] * func(x + (n / 2. - i) * dx); }
            return result / pow(dx, n);
        };
        
        return (4 * deriv(dx / 2) - deriv(dx)) / 3;
    }

    double derivative2(std::function<double(double)> func, double x, double dx, int n) 
    {
        return (func(x + dx) - 2 * func(x) + func(x - dx)) / dx / dx;
    }
}

#endif