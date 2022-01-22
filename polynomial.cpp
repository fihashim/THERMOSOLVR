#include "cluster2.h"
#include "atomic_data.h"
#include "polynomial.h"
#include<complex>
#include<iterator>
#include<vector>
#include<limits>
#include<numeric>
#include<iostream>
#include<iomanip>
#include<tuple>

std::tuple<bool, double>booldoubletuple(bool success, double x) {
    return  std::make_tuple(success, x);
}

bool definitelygreaterthan(double a, double b, const double &epsilon)
{
    return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool definitelylessthan(double a, double b, const double &epsilon)
{
    return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

double horner_derivative(int n, const std::vector<double>poly, double x)
{
    double result = poly[n];
    double pdiff = 0.00;

    for (std::vector<double>::size_type t = poly.size()- 2; t != (std::vector<double>::size_type) -1; t--)
    {
        pdiff = pdiff * x + result;
        result = result * x + poly[t];
    }
    return pdiff; 
}

double horner_polynomial(int n, const std::vector<double>&poly, double x)
{
    double result = poly[n];
    double pdiff = 0.00;
    for (std::vector<int>::size_type t = poly.size()- 2;t != (std::vector<int>::size_type) -1; t--)
    {
        pdiff = pdiff * x + result;
        result = result * x + poly[t];
    };
    return result; 
}

std::tuple<bool, double>newton(int n, const std::vector<double>&poly, double x, int newton_iterations, bool success)
{
    double result;
    double pdiff; 
    double dx;
    double new_pop;
    double xnew;
    std::tuple<bool, double>final_newton;
    double newton_convergence = global_eps;
    for (int i = 0; i < newton_iterations; i++)
    {
        pdiff = horner_derivative(n-1, poly, x);
        result = horner_polynomial(n-1, poly, x);

        if (std::abs(result) <= newton_convergence)
        {
            success = true; 
            final_newton = booldoubletuple(success,x);
            return final_newton;
        }

        dx = -result/pdiff;
        xnew = x + dx;

        while(definitelylessthan(xnew, 0.00,global_eps) or definitelygreaterthan(xnew, 1.00, global_eps))
        {
            dx = 0.5 * dx;
            xnew = x + dx;
        }

        x = xnew;
        final_newton = booldoubletuple(success,x);
    }
    success = false; 
    return final_newton; 
}

typedef std::numeric_limits<double>dbl;
std::vector<std::complex<double>>solve_polynomial3(std::vector<double>coeffs)
{
    std::vector<std::complex<double>>roots; 
    std::vector<std::complex<double>>cof;
    std::complex<double>im(0.00,1.00);
    std::complex<double>one_third(1.00/3.00, 0.00);
    double two_power_one_third;
    two_power_one_third = pow(2.00,1.00/3.00);
    double sqrt3 = sqrt(3.00);
    const std::complex<double>temp(0.00,0.00);
    std::cout.precision(dbl::max_digits10 - 2);

    double cof1 = -coeffs[2]/3.00/coeffs[3];
    std::complex<double>cof_complex1(cof1);
    double cof2 = -pow(coeffs[2],2.00) + 3.00 * coeffs[3]*coeffs[1];
    std::complex<double>cof_complex2(cof2);
    double cof3 = (-2.00 * pow(coeffs[2],3.00) + 9.00 * coeffs[3] * coeffs[2] * coeffs[1] - 27.00 * coeffs[0] * pow(coeffs[3],2.00));
    std::complex<double>cof_complex3(cof3);
    double cof4 = 4.00 * pow(cof2,3.00) + pow(cof3,2.00);
    std::complex<double>cof_complex4(cof4);
    double cof5 = 3.00 * coeffs[3];
    std::complex<double>cof_complex5(cof5);
    std::complex<double>cof_complex6(pow(cof_complex3 + sqrt(cof_complex4), one_third));
    
    cof.push_back(cof_complex1);
    cof.push_back(cof_complex2);
    cof.push_back(cof_complex3);
    cof.push_back(cof_complex4);
    cof.push_back(cof_complex5);
    cof.push_back(cof_complex6);
  
    std::complex<double>root1(cof_complex1 - two_power_one_third * cof_complex2/(cof_complex5*cof_complex6) + cof_complex6/two_power_one_third/cof_complex5);
    std::complex<double>root2(cof_complex1 + (1.00 + im * sqrt3) * cof_complex2 / (pow(two_power_one_third,2.00) * cof_complex5 * cof_complex6) - (1.00 - im * sqrt3) * cof_complex6 / (2.0 * cof_complex5 * two_power_one_third));
    std::complex<double>root3(cof_complex1 + (1.00 - im * sqrt3) * cof_complex2/(pow(two_power_one_third,2.00) * cof_complex5 * cof_complex6) - (1.00 + im * sqrt3) * cof_complex6/(2.00 * cof_complex5 * two_power_one_third));

    roots.push_back(root1);
    roots.push_back(root2);
    roots.push_back(root3);
    
    return roots; 
}