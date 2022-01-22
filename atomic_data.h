#ifndef ATOMIC_H
#define ATOMIC_H
#include<map>
#include<cmath>
#include<string>


const std::map<std::string, double>periodic_table= 
{
    { "H", 1.0080000162124634 },
    { "He", 4.003 },
    { "Li", 6.94 },
    { "Be", 9.012 },
    { "B", 10.81 },
    { "C", 12.01 },
    { "N", 14.01 },
    { "O", 16.00 },
    { "F", 19.00 },
    { "Ne", 20.18 },
    { "Na", 22.99 },
    { "Mg", 24.31 },
    { "Al", 26.98 },
    { "Si", 28.09 },
    { "P", 30.97 },
    { "S", 32.06 },
    { "Cl", 35.45 },
    { "Ar", 39.95 }
};

const double pi = 4.0*atan(1.0);
const double planck = 6.62606957e-34; //J*s
const double kb = 1.3806488e-23; //J/K
const double speed_of_light = 299792458.0; // m/s
const double avogadro = 6.0221413e23; // 1/mol
const double amu = 1.660538921e-27; //kg
const double pressure = 1.01325; 
const double hbar = planck / (2 * pi);
const double factor = planck * 100.0 * speed_of_light/kb;
const double global_eps = 1.0e-10; 

#endif