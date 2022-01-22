#ifndef GLOBAL_H
#define GLOBAL_H
#include<vector>
#include "auxiliary2.h"

class range_c;

class global_c
{
    public:
    int qce_iterations, newton_iterations, nconverged;
    std::vector<double>degree;
    double press, mtot, vexcl;
    std::vector<double>ntot;
    double max_deviation, vdamp;
    range_c amf, bxv, temp;
    std::vector<double>amf_pure;
    std::vector<double>bxv_pure;
    std::vector<double>monomer_amounts; 
    bool progress_bar;     
};

extern class global_c global_data; 

#endif
