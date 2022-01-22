#include<tuple>
#include<vector>
#include<complex>

std::tuple<bool, double>booldoubletuple(bool success, double x);
bool definitelygreaterthan(double a, double b, const double &epsilon);
bool definitelylessthan(double a, double b, const double &epsilon);
double horner_derivative(int n, const std::vector<double>poly, double x);
double horner_polynomial(int n, const std::vector<double>&poly, double x);
std::tuple<bool, double>newton(int n, const std::vector<double>&poly, double x, int newton_iterations, bool success);
std::vector<std::complex<double>>solve_polynomial3(std::vector<double>coefficients);






