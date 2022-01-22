#include "input.h"
#include "config2.h"
#include "auxiliary2.h"
#include<iterator>

// generic template to print any type of vector. 
template<typename T>
void printVector(const T& t) {
    std::copy(t.cbegin(), t.cend(), std::ostream_iterator<typename T::value_type>(std::cout, "\n "));
};

void input_c::set_components(int component)
{
    components = component; 
};

void input_c::set_amf(double amf_start, double amf_end, int amf_step)
{
    amf_c.set_range(amf_c, amf_start, amf_end, amf_step);
};

void input_c::set_bxv(double bxv_start, double bxv_end, int bxv_step)
{
    bxv_c.set_range(bxv_c, bxv_start, bxv_end, bxv_step);
};

void input_c::set_temperature(double temp_low, double temp_high, int temp_step)
{
    temperature_c.set_range(temperature_c, temp_low, temp_high, temp_step);
};

void input_c::set_pressure(double pressure)
{
    pressure_c = pressure; 
}

void input_c::set_monomer_amounts(double monomer_amounts)
{
    std::replace(monomer_amounts_c.begin(), monomer_amounts_c.end(), monomer_amounts_c[0], monomer_amounts);
}

void input_c::set_contributions()
{
    contrib = true;
    cv_contrib = true;
    internal_contrib = true;
    helmholtz_contrib = true;
    entropy_contrib = true; 
}

void input_c::set_compare_density(double ref_density_temperature, double ref_density, std::string ref_density_weight)
{
    compare_density = true;
    compare = true; 
    ref_density_temperature_c = ref_density_temperature;
    ref_density_c = ref_density;
    if (ref_density_weight == "#")
    {
        ref_density_weight_c = 1.0;
    }
    else
    {
        ref_density_weight_c = stod(ref_density_weight); 
    }
}

void input_c::set_phase_transition(double phase_transition_temp, std::string phase_transition_weight)
{
    compare_phase_transition = true; 
    compare = true;
    ref_phase_transition_c = phase_transition_temp;
    if (phase_transition_weight == "#")
    {
        ref_phase_transition_weight_c = 1.0; 
    }
    else
    {
        ref_phase_transition_weight_c = stod(phase_transition_weight);
    }

}

void input_c::set_isobar(std::string isobar_weight)
{
    compare_isobar = true; 
    if (isobar_weight == "#")
    {
        ref_isobar_weight_c = 1.0;
    }
    else
    {
        ref_isobar_weight_c = stod(isobar_weight);
    }
}

void input_c::set_ref_isobar_volume(const std::vector<double>&ref_isobar_volume)
{
    ref_isobar_volume_c = ref_isobar_volume; 
}

void input_c::set_ref_isobar_temp(const std::vector<double>&ref_isobar_temp)
{
    ref_isobar_temperature_c = ref_isobar_temp;
}

int input_c::get_components()
{
    return components; 
}

double input_c::get_pressure()
{
    return pressure_c; 
}

range_c input_c::get_amf()
{
    return amf_c;
}

std::vector<double>input_c::get_ref_isobar_volume_c()
{
    return ref_isobar_volume_c;
}

std::vector<double>input_c::get_ref_isobar_temperature_c()
{
    return ref_isobar_temperature_c;
}

range_c input_c::get_bxv()
{
    return bxv_c;
}

range_c input_c::get_temperature()
{
    return temperature_c; 
}

std::vector<double> input_c::get_monomer_amounts()
{
    return monomer_amounts_c;
}

bool input_c::get_contributions()
{
    return contrib; 
}

bool input_c::get_helmholtz_contrib()
{
    return helmholtz_contrib;
}

bool input_c::get_internal_contrib()
{
    return internal_contrib;
}

bool input_c::get_entropy_contrib()
{
    return entropy_contrib;
}

bool input_c::get_cv_contrib()
{
    return cv_contrib; 
}

double input_c::get_max_deviations()
{
    return max_deviation; 
}

double input_c::get_volume_damping_factor()
{
    return volume_damping_factor;
}

double input_c::get_ref_isobar_weight_c()
{
    return ref_isobar_weight_c; 
}

int input_c::get_qce_iterations()
{
    return qce_iterations; 
}

int input_c::get_newton_iterations()
{
    return newton_iterations;
}

bool input_c::get_progress_bar()
{
    return progress_bar_c; 
}

std::vector<double>input_c::get_amf_pure()
{
    return amf_pure;
}

std::vector<double>input_c::get_bxv_pure()
{
    return bxv_pure;
}

bool input_c::get_compare()
{
    return compare; 
};

bool input_c::get_compare_isobar()
{
    return compare_isobar; 
}

double input_c::get_compare_isobar_weight()
{
    return ref_isobar_weight_c; 
}

bool input_c::get_compare_density()
{
    return compare_density; 
}

double input_c::get_ref_density_weight()
{
    return ref_density_weight_c;
}

double input_c::get_ref_density()
{
    return ref_density_c; 
}

double input_c::get_ref_density_temperature()
{
    return ref_density_temperature_c; 
}

bool input_c::get_compare_phase_transition()
{
    return compare_phase_transition; 
}

double input_c::get_ref_phase_transition_weight_c()
{
    return ref_phase_transition_weight_c; 
}

double input_c::get_ref_phase_transition()
{
    return ref_phase_transition_c; 
}

input_c process_input(std::string inputfile)
{
    std::string input_contents;
    input_c input_obj; 
    input_contents = readFile(inputfile);
    input_obj = get_config_arguments(input_contents);
    check_input(input_obj);
return input_obj; 
}

void check_input(input_c &input_object)
{
    if (input_object.get_amf().first < 0.00)
    {
        std::cout << "unphysical value for keyword amf" << std::endl;
        exit(1);  
    }
    range_c amf = input_object.get_amf(); 
    if (not check_range(amf))
    {
        std::cout << "illegal range specification for keyword amf" << std::endl; 
    }

    range_c bxv = input_object.get_bxv();
    if (not check_range(bxv))
    {
        std::cout << "illegal range specification for keyword bxv" << std::endl;
        exit(1);
    }

    if (input_object.get_bxv().first <= 0.00)
    {
        std::cout << "unphysical value for keyword bxv" << std::endl; 
        exit(1);
    }

    if (input_object.get_newton_iterations() <= 0)
    {
        std::cout << "newton iterations must be > 0" << std::endl; 
    } 

    if (input_object.get_temperature().first < 0.00)
    {
        std::cout << "minimum temperature must be > 0" << std::endl;
    }

    range_c temp = input_object.get_temperature(); 
    if (not check_range(temp))
    {
        std::cout << "illegel range specification for keyword temperature" << std::endl; 
    }

    if (input_object.get_pressure() < 0.0)
    {
        std::cout << "unphysical value for keyword pressure" << std::endl; 
    }

    if (input_object.get_components() <= 0)
    {
        std::cout << "unphysical value for keyword components" << std::endl;
    }

    if (input_object.get_monomer_amounts()[0] < 0.00)
    {
        std::cout << "monomer amounts must be > 0" << std::endl; 
    }

    if (input_object.get_max_deviations() <= 0.00)
    {
        std::cout << "max deviations must be > 0" << std::endl; 
    }

    if (input_object.get_volume_damping_factor() <= 0.00 or input_object.get_volume_damping_factor() >= 1.0)
    {
        std::cout << "volume damping factor must be between 0.01 to 0.99" << std::endl; 
    }

    if (input_object.get_amf().num * input_object.get_bxv().num > 1 and not input_object.get_compare())
    {
        std::cout << "amf and/or bxv interval specified, but no reference section" << std::endl;
    }

    // check density
    if (input_object.get_compare_density())
    {
        if (input_object.get_ref_density() <= 0.00)
        {
            std::cout << "unphysical value for reference density" << std::endl; 
        }

        if (input_object.get_ref_density_weight() < 0.0)
        {
            std::cout << "unphysical value for reference density weight" << std::endl; 
        }

        if (input_object.get_ref_density_temperature() <= 0.00)
        {
            std::cout << "unphysical value for reference density temperature" << std::endl;
        }

        if (input_object.get_ref_density_temperature() <= input_object.get_temperature().first or input_object.get_ref_density_temperature() >= input_object.get_temperature().last)
        {
            std::cout << "density reference temperature must be within the investigated temperature range" << std::endl; 
        }
    }

    // check phase transition
    if (input_object.get_compare_phase_transition())
    {
        if (input_object.get_ref_phase_transition() < 0) 
        {
            std::cout << "unphysical value for reference phase transition" << std::endl;
        }

        if (input_object.get_ref_phase_transition_weight_c() < 0.0)
        {
            std::cout << "unphysical value for phase transition weight" << std::endl; 
        }

        if (input_object.get_ref_phase_transition() <= input_object.get_temperature().first or input_object.get_ref_phase_transition() >= input_object.get_temperature().last)
        {
            std::cout << "reference temperature of phase transition must be within the investigated temperature range" << std::endl; 
        }

        if (input_object.get_temperature().num < 2)
        {
            std::cout << "need two temperature points for phase transition determination" << std::endl; 
        }
    }
    
    // check isobar
    if (input_object.get_compare_isobar())
    {
        if (input_object.get_ref_isobar_weight_c() < 0.0)
        {
            std::cout << "unphysical value for isobar weight in reference section" << std::endl;
        }

        std::vector<double>ref_temp = input_object.get_ref_isobar_temperature_c();
        bool zeros = std::any_of(ref_temp.begin(), ref_temp.end(), [](int i) { return i<0; });
        if (zeros)
        {
            std::cout << "unphysical value for isobar reference temperatures in reference section" << std::endl; 
        }

        std::vector<double>ref_volume = input_object.get_ref_isobar_volume_c();
        bool zeros_vol = std::any_of(ref_volume.begin(), ref_volume.end(), [](int i) { return i<0; });
        if (zeros_vol)
        {
            std::cout << "unphysical value for isobar reference volumes in reference section" << std::endl; 
        }

        if (ref_temp[0] < input_object.get_temperature().first or ref_temp[ref_temp.size() - 1] > input_object.get_temperature().last)
        {
            std::cout << "reference isobar temperatures must be within the investigated temperature range" << std::endl; 
        }
    }
}