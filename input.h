#ifndef INPUT_H
#define INPUT_H
#include "auxiliary2.h"
#include<iostream>
#include<algorithm>
#include<vector>

class range_c; 


class input_c
{
    // input read from the [system] section
    int components;

    // input read from the [ensemble] section
    range_c temperature_c; 
    double pressure_c;
    std::vector<double>monomer_amounts_c; 

    // input read from the [qce] section
    range_c amf_c, bxv_c;
    std::vector<double>amf_pure;
    std::vector<double>bxv_pure;
    double max_deviation, volume_damping_factor;
    int qce_iterations;
    int newton_iterations;

    // input read from the [reference] section
    bool compare, compare_density, compare_isobar, compare_phase_transition;
    double ref_phase_transition_c, ref_density_c, ref_density_temperature_c; 
    double ref_phase_transition_weight_c, ref_density_weight_c, ref_isobar_weight_c; 
    std::vector<double>ref_isobar_temperature_c, ref_isobar_volume_c;
    std::vector<std::string>ref_isobar_file; 

    // input read from the [output] section
    bool contrib, helmholtz_contrib, internal_contrib, entropy_contrib, cv_contrib;
    bool progress_bar_c; 

    public:
    input_c()
    {
        // defaults for the [system] section
        components = 1; 

        // defaults for the [ensemble] section
        temperature_c.set_range(temperature_c, 298.15, 298.15, 1); // in K
        pressure_c = 1.01325; // in bar 
        monomer_amounts_c.resize(components);
        monomer_amounts_c.push_back(1.0 / components); // in mol 

        // defaults for [qce] section
        amf_c.set_range(amf_c, 0.0, 0.0, 1);
        bxv_c.set_range(bxv_c, 1.0, 1.0, 1);
        amf_pure.resize(components);
        bxv_pure.resize(components);
        amf_pure.push_back(0.0);
        bxv_pure.push_back(0.0);
        max_deviation = 1.0e-9;
        volume_damping_factor = 0.01;
        qce_iterations = 100;
        newton_iterations = 500; 

        // defaults for [reference] section
        compare = false;
        compare_isobar = false;
        compare_density = false; 
        compare_phase_transition = false;
        ref_isobar_weight_c = 1.0;
        ref_density_weight_c = 1.0;
        ref_phase_transition_weight_c = 1.0;

        // defaults for [output] section
        contrib = false;
        helmholtz_contrib = false;
        internal_contrib = false;
        entropy_contrib = false;
        cv_contrib = false;
        progress_bar_c = true; 
    }
    void set_components(const int component);
    void set_amf(const double amf_start, const double amf_end, const int amf_step);
    void set_bxv(const double bxv_start, const double bxv_end, const int bxv_step);
    void set_temperature(const double temp_low, const double temp_high, const int temp_step);
    void set_pressure(const double pressure);
    void set_monomer_amounts(const double monomer_amounts);
    // void set_max_deviation(const double &max_deviation);
    void set_contributions();
    void set_compare_density(const double ref_density_temperature, const double ref_density, const std::string ref_density_weight);
    void set_phase_transition(const double phase_transition_temp, const std::string phase_transition_weight);
    void set_isobar(const std::string isobar_weight);
    void set_ref_isobar_volume(const std::vector<double>&ref_isobar_volume);
    void set_ref_isobar_temp(const std::vector<double>&ref_isobar_temp);
    
    // getter functions
    std::vector<double>get_ref_isobar_volume_c();
    int get_components();
    double get_pressure(); 
    range_c get_amf(); 
    range_c get_bxv(); 
    range_c get_temperature();
    std::vector<double>get_monomer_amounts();
    double get_max_deviations(); 
    double get_volume_damping_factor(); 
    int get_qce_iterations();
    int get_newton_iterations();
    bool get_progress_bar(); 
    std::vector<double>get_amf_pure(); 
    std::vector<double>get_bxv_pure();
    bool get_compare(); 
    bool get_compare_isobar(); 
    double get_compare_isobar_weight();
    std::vector<double>get_ref_isobar_temperature_c();
    bool get_compare_density(); 
    bool get_contributions(); 
    bool get_helmholtz_contrib();
    bool get_internal_contrib(); 
    bool get_entropy_contrib(); 
    bool get_cv_contrib();
    double get_ref_density_weight(); 
    double get_ref_density(); 
    double get_ref_density_temperature(); 
    bool get_compare_phase_transition(); 
    double get_ref_phase_transition_weight_c();
    double get_ref_phase_transition(); 
    double get_ref_isobar_weight_c(); 


    

};

    input_c process_input(std::string inputfile);
    void check_input(input_c &input_object); 



#endif


