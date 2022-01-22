#include "qce.h"
#include "global.h"
#include "input.h"
#include "/opt/homebrew/include/omp.h"
#include "thermo.h"
#include "cluster2.h"
#include "atomic_data.h"
#include "auxiliary2.h" 
#include "polynomial.h"
#include "partfunctions.h"
#include<numeric>
#include<filesystem>
#include<map>
#include<functional>
#include<iomanip>
#include<iterator>
#include<set>
#include<vector>
#include<tuple>
#include<algorithm>
#include<fstream>


class global_c global_data; 

std::vector<std::filesystem::path>getdatfiles(std::filesystem::path const &root, std::string const &ext)
{
    std::vector<std::filesystem::path>datfiles;

    if (std::filesystem::exists(root) && std::filesystem::is_directory(root))
    {
        for (auto const & entry : std::filesystem::recursive_directory_iterator(root))
        {
            if (std::filesystem::is_regular_file(entry) && entry.path().extension() == ext)
            {
                datfiles.emplace_back(entry.path().filename());
            }
        }
    }
    return datfiles; 
}

std::tuple<bool, std::vector<double>>boolvecttuple(bool success, std::vector<double>population) {
    return  std::make_tuple(success, population);
};

// template<class T>
// void printVector(std::vector<T> const&input)
// {
//     for (auto i = input.begin(); i != input.end(); ++i)
//         std::cout << *i << std::endl;
// }

void isobar_c::set_error(const double &error)
{
    error_c = error;
}

void isobar_c::set_amf(const double &amf)
{
    amf_c = amf; 
}

void isobar_c::set_bxv(const double &bxv)
{
    bxv_c = bxv; 
}

void isobar_c::set_temperature(double temperature, int i)
{
    temp_c[i] = temperature; 
}

void isobar_c::set_gibbs(double gibbs, int i)
{
    gibbs_c[i] = gibbs; 
}

void isobar_c::set_volume(double vol, int i)
{
    vol_c[i] = vol; 
}

void isobar_c::set_solution(int solution, int& i)
{
    solution_c[i] = solution; 
}

void isobar_c::set_populations(std::vector<double>populations, int& i)
{
    populations_c[i] = populations; 
}

void isobar_c::set_lnq(std::vector<pf_c>lnq, int& i)
{
    lnq_c[i] = lnq; 
}

void isobar_c::set_converged(bool converged, int i)
{
    converged_c[i] = converged;
}
void isobar_c::resize_temperature(int &num)
{
    temp_c.resize(num);
};

void isobar_c::resize_volume(int &num)
{
    vol_c.resize(num);
}

void isobar_c::resize_gibbs(int &num)
{
    gibbs_c.resize(num);
}

void isobar_c::resize_converged(int &num)
{
    converged_c.resize(num);
}

void isobar_c::resize_solution(int &num)
{
    solution_c.resize(num);
}

void isobar_c::resize_lnq(int &num, int clustersize)
{
    lnq_c.resize(num);
    for (int i = 0; i < num; ++i)
    {
        lnq_c[i].resize(clustersize);
    }
}

void isobar_c::resize_populations(int &num, int clustersize)
{
    populations_c.resize(num);
    for (int i = 0; i < num; i++)
    {
        populations_c[i].resize(clustersize);
    }
}

std::vector<double>isobar_c::get_temperature()
{
    return temp_c; 
}

double& isobar_c::get_temperature_element(int i)
{
    return temp_c[i];
}

double& isobar_c::get_amf()
{
    return amf_c; 
}

double& isobar_c::get_bxv()
{
    return bxv_c; 
}

double& isobar_c::get_volume(int i)
{
    return vol_c[i];
}

int isobar_c::get_solution(int i)
{
    return solution_c[i];
}

int isobar_c::get_converged_element(int i)
{
    return converged_c[i];
}

double& isobar_c::get_gibbs(int i)
{
    return gibbs_c[i];
}

std::vector<bool>& isobar_c::get_converged()
{
    return converged_c;
}

int count_converge(std::vector<bool>&converged_c)
{
    int converge_count = 0; 
    for (int i = 0; i < converged_c.size(); i++)
    {
        if (converged_c[i] == true)
        {
            converge_count += 1; 
        }
    }
    return converge_count; 
}
std::vector<double>& isobar_c::get_populations(int i)
{
    return populations_c[i];
}

std::vector<pf_c>& isobar_c::get_lnq(int i )
{
    return lnq_c[i];
}

std::vector<std::vector<pf_c>>& isobar_c::get_lnq_vector()
{
    return lnq_c;
}

std::vector<double>& isobar_c::get_volume_vec()
{
    return vol_c; 
}

std::vector<std::vector<double>>& isobar_c::get_populations_vec()
{
    return populations_c; 
}

double& isobar_c::get_error()
{
    return error_c; 
}


void initialize_conserved_quantities(cluster_c &monomer)
{ 
    // particle number
    global_data.ntot.resize(monomer.get_composition());
    global_data.ntot[0] = global_data.monomer_amounts[0] * avogadro;

    // mass
    global_data.mtot = global_data.mtot + global_data.monomer_amounts[0] * monomer.get_molmass(); 
    global_data.mtot = global_data.mtot / 1000.00; 

    // excluded volume
    global_data.vexcl = global_data.monomer_amounts[0] * monomer.get_volume();
    global_data.vexcl = global_data.vexcl * avogadro * 1.0e-30; 
}

void initialize_degree(std::vector<cluster_c>&clusterset)
{
    global_data.degree.resize(1);
    for (int i = 0; i < clusterset.size(); i++)
    {
        global_data.degree[0] = clusterset[i].get_composition();
    }
}

void initialize_populations(std::vector<double>&populations, int &no_clusters)
{
    populations.resize(no_clusters);
    populations.insert(populations.begin(), global_data.monomer_amounts[0] * avogadro);
}



double compare_density(isobar_c& ib, double &temp, double &density)
{
    double vol, vcalc, error;
    vol = global_data.mtot / (1.00e3 * density);
    vcalc = determine_volume(ib, temp);
    error = pow(((vcalc - vol) / vol), 2.00);
    return error; 
}

double determine_volume(isobar_c &ib, double &temp)
{
    double x, predicted_volume;
    int index;

    for (int i = 0; i < ib.get_temperature().size(); i++)
    {
        if (ib.get_temperature_element(i) >= temp)
        {
            index = i;
            break; 
        }
    }
    int index_1 = index + 1;
    x = (temp - ib.get_temperature_element(index)) / (ib.get_temperature_element(index_1) - ib.get_temperature_element(index));
    predicted_volume = (1.00 - x) * ib.get_volume(index) + x * ib.get_volume(index_1);
    return predicted_volume; 
}

double compare_isobar(isobar_c& ib, std::vector<double>&temp, std::vector<double>&vol)
{
    std::vector<double>temp_vec = temp;
    std::vector<double>vol_vec = vol; 
    double error, vcalc;
    error = 0.00; 
    int no_temp = temp_vec.size();
    for (int i = 0; i < temp_vec.size(); i++)
    {
        vcalc = determine_volume(ib, temp_vec[i]);
        double expt_diff = (vcalc - vol_vec[i]) / vol_vec[i]; 
        double expt_diff_sqrd = pow(expt_diff, 2.00);
        error += expt_diff_sqrd; 
    } 
    error = error / no_temp; 
    return error; 
}

double determine_phase_transition(isobar_c& ib)
{
    double dv_max, dv, predicted_phase_transition;
    int itemp_max;
    itemp_max = 1; 
    dv_max = 0.00;
    for (int i = 1; i < ib.get_temperature().size(); i++)
    {
        bool first_converged = ib.get_converged_element(i);
        bool second_converged = ib.get_converged_element(i-1);
        if (first_converged and second_converged)
        {
            dv = ib.get_volume(i) - ib.get_volume(i-1);
            if (dv > dv_max)
            {
                dv_max = dv;
                itemp_max = i;
            }
        }
    }

    // peacemaker equation
    predicted_phase_transition = 0.50 * (ib.get_temperature_element(itemp_max) + ib.get_temperature_element(itemp_max-1));

    // old qce_program
    // predicted_phase_transition =ib.get_temperature_element(itemp_max);

    return predicted_phase_transition; 
}

double compare_phase_transition(isobar_c& ib, double& pt)
{
    double error, ptcalc;
    ptcalc = determine_phase_transition(ib);
    ib.predicted_bp = ptcalc; 
    error = pow((ptcalc - pt) / pt, 2.00);

    return error; 
}

double ideal_gas_volume(double &temp)
{
    double ideal_gas_volume;
    double ntot_sum;
    for (auto ntot : global_data.ntot)
    { 
        ntot_sum += ntot;    
    }

    ideal_gas_volume = ntot_sum * kb * temp / global_data.press; 
    return ideal_gas_volume; 
}

std::vector<double>calculate_remaining_populations(std::vector<double>&monomer_populations, std::vector<cluster_c>&clusterset, std::vector<pf_c>&lnq)
{
    std::vector<double>populations;
    populations = monomer_populations; 
    double tmp; 
    for (int i = 1; i < clusterset.size(); i++)
    {
        tmp = tmp + clusterset[i].get_composition() * (log(monomer_populations[0]) - lnq[0].get_qtot());
        tmp = tmp + lnq[i].get_qtot(); 

        populations.push_back(exp(tmp));
        tmp = 0.00; 
    }
    return populations; 
}

std::tuple<bool, std::vector<double> >calculate_populations(std::vector<double>&populations, std::vector<cluster_c>&clusterset, std::vector<pf_c>&lnq, int &cyclus, cluster_c &monomer)
{
    std::vector<double>pop_coeffs; 
    std::vector<double>coeffs1_pop; 
    std::vector<double>mass_coeffs; 
    std::vector<double>monomer_populations; 
    std::tuple<bool, std::vector<double>>final_pop; 
    bool success; 
    std::cout.precision(12);

    for (int i = 0; i < clusterset.size(); i++)
    {
        double coeff = 0.00; 
        coeff = lnq[i].get_qtot();
        coeff = coeff - (clusterset[i].get_composition() * lnq[0].get_qtot()); 
        coeff = coeff + (clusterset[i].get_composition() * log(global_data.ntot[0]));
        pop_coeffs.push_back(exp(coeff) * clusterset[i].get_composition() / global_data.ntot[0]); 
        // std::cout << std::scientific << std::uppercase << pop_coeffs[i] << std::endl;
    }

    if (cyclus == 1)
    {
        monomer_populations.push_back(global_data.monomer_amounts[0] * avogadro / global_data.ntot[0]);
    }
    else
    {
        monomer_populations.push_back(1.0e-2 * global_data.monomer_amounts[0] * avogadro / global_data.ntot[0]);
    }
   
    int n = global_data.degree[0] + 1;
    coeffs1_pop.resize(n);
    std::fill(coeffs1_pop.begin(), coeffs1_pop.end(), 0.00);
    coeffs1_pop[0] = -1.00;
    std::map<int, double>coeffs_composition_map; 

    for (int i = 0; i < clusterset.size(); i++)
    {
        coeffs_composition_map.insert(std::pair<int,double>(clusterset[i].get_composition(), pop_coeffs[i]));
    }

    std::map<int,double>::iterator itr; 
    for (auto &comp : coeffs_composition_map)
    {
        coeffs1_pop[comp.first] = comp.second; 
    }
   
    std::tie(success, monomer_populations[0]) = newton(n, coeffs1_pop, monomer_populations[0], global_data.newton_iterations, success);
 
    if (monomer_populations[0] < 0.00)
    {
        success = false;
    }
    
    if (success == true)
    {
        monomer_populations[0] = monomer_populations[0] * global_data.ntot[0];
        populations = calculate_remaining_populations(monomer_populations, clusterset, lnq);        
    }

    final_pop = boolvecttuple(success, populations);
    
    return final_pop; 
}

std::vector<int>isobar_c::get_solution_vec()
{
    return solution_c; 
}

typedef std::numeric_limits<double>dbl;
std::tuple<bool, double>calculate_volume(double &vol, double &vdamp, double &amf, double &bxv, double &temp, std::vector<double>&populations)
{
    std::cout.precision(dbl::max_digits10 - 2);
    bool success;
    double sum_of_populations; 
    std::vector<double>coeffs;
    std::vector<std::complex<double>>roots;
    bool valid_root; 
    std::vector<bool>valid_roots; 
    std::multimap<bool, double>valid_roots_map; 
    std::vector<double>amf_pure_local;
    std::vector<bool>valid_roots_vec;     

    for (auto amf_value : global_data.amf_pure)
    {
        if (amf_value <= 0.00)
        {
            amf_value = amf; 
        }
    }
    
    coeffs.push_back(amf * bxv * global_data.vexcl * pow(global_data.ntot[0], 2));
    coeffs.push_back(-amf * pow(global_data.ntot[0], 2));

    for (int i = populations.size() - 1; i > -1; i--)
    {
        sum_of_populations += populations[i];
    }
    
    // calculate the coefficients
    coeffs.push_back(kb * temp * sum_of_populations + global_data.press * bxv * global_data.vexcl);
    coeffs.push_back(-global_data.press);
 
    roots = solve_polynomial3(coeffs);

    for (int i = 0; i < roots.size(); i++)
    {
        if (abs(imag(roots[i])) <= global_eps and real(roots[i]) - bxv * global_data.vexcl >= global_eps)
        {
            valid_root = true;
            valid_roots_vec.push_back(valid_root);
            valid_roots_map.insert(std::pair<bool, double>(valid_root, real(roots[i])));
        }    
        else
        {
            valid_root = false;
            valid_roots_vec.push_back(valid_root);
            valid_roots_map.insert(std::pair<bool, double>(valid_root, real(roots[i])));
        }
    }

    typedef std::function<bool(std::pair<bool, double>, std::pair<bool, double>)> Comparator; 
    if (std::any_of(valid_roots_vec.begin(), valid_roots_vec.end(), [](bool valid_roots){return valid_roots ==true;}))
    {
        if (vdamp <= 1.00)
        {
            Comparator compfunctor = [](std::pair<bool, double> elem1, std::pair<bool, double> elem2)
            {
                return elem1.second < elem2.second; 
            };

            std::set<std::pair<bool, double>, Comparator>valid_roots_sorted(valid_roots_map.begin(), valid_roots_map.end(), compfunctor);
            for (auto element : valid_roots_sorted)
            {
                if (element.first == true)
                {
                    vol = element.second; 
                }
            }
        }

        else
        {
            Comparator complessthan = [](std::pair<bool, double>elem1, std::pair<bool, double>elem2)
            {
                return elem1.second > elem2.second; 
            };
            std::set<std::pair<bool, double>, Comparator>valid_roots_sorted(valid_roots_map.begin(), valid_roots_map.end(), complessthan);
            for (auto &element : valid_roots_sorted)
            {
                if (element.first == true)
                {
                    vol = element.second; 
                }
            }
        }
        success = true;
    }
    else
    {
        success = false; 
    }

    return booldoubletuple(success, vol);
}

std::tuple<bool, double>check_convergence(double &gibbs, double &vol, double &temp, std::vector<double>&populations, std::vector<pf_c>&lnq)
{
    bool converged = false;
    double new_gibbs = 0.00;
    double deviation = 0.00; 
    for (int i = 0; i < populations.size(); i++)
    {
        new_gibbs = new_gibbs - ln_factorial(populations[i]) + populations[i] * lnq[i].get_qtot(); 
    }
    new_gibbs = -kb * temp * new_gibbs + global_data.press * vol; 
    deviation = new_gibbs / gibbs - 1.00;
    gibbs = new_gibbs; 

    if (abs(deviation) <= global_data.max_deviation)
    {
        converged = true; 
    }

    return std::make_tuple(converged, gibbs);
}



void qce_prepare(input_c input_object, cluster_c &monomer, std::vector<cluster_c>&clusterset, reference_c &ref_object)
{
    global_data.press = input_object.get_pressure() * 1.0e5; 
    global_data.max_deviation = input_object.get_max_deviations(); 
    global_data.vdamp = 0.01; 
    global_data.qce_iterations = input_object.get_qce_iterations();
    global_data.newton_iterations = input_object.get_newton_iterations();
    global_data.amf = input_object.get_amf();
    global_data.bxv = input_object.get_bxv();
    global_data.temp = input_object.get_temperature();
    global_data.progress_bar = input_object.get_progress_bar(); 

    global_data.amf_pure.resize(monomer.get_composition());
    double amf_pure = input_object.get_amf_pure()[0];
    global_data.amf_pure[0] = amf_pure / pow(avogadro, 2);

    global_data.bxv_pure.resize(monomer.get_composition());
    double bxv_pure = input_object.get_bxv_pure()[0]; 
    global_data.bxv_pure[0] = bxv_pure; 
    
    global_data.monomer_amounts = input_object.get_monomer_amounts(); 

    initialize_conserved_quantities(monomer);
    global_data.degree.resize(monomer.get_composition());
    initialize_degree(clusterset);
    global_data.nconverged = 0; 

    // assign reference data
    ref_object.compare = input_object.get_compare(); 
    ref_object.compare_isobar = input_object.get_compare_isobar();

    if (ref_object.compare_isobar)
    {
        ref_object.isobar_weight = input_object.get_compare_isobar_weight();
        ref_object.isobar_temperature = input_object.get_ref_isobar_temperature_c(); 
        ref_object.isobar_volume = input_object.get_ref_isobar_volume_c(); 
        for (int i = 0; i < ref_object.isobar_volume.size(); i++)
        {
            ref_object.isobar_volume[i] = 1.0e-3 * ref_object.isobar_volume[i]; 
        }
    }

    ref_object.compare_density = input_object.get_compare_density(); 
    if (ref_object.compare_density)
    {
        ref_object.density_weight = input_object.get_ref_density_weight();
        ref_object.density = input_object.get_ref_density();
        ref_object.density_temperature = input_object.get_ref_density_temperature(); 
    }

    ref_object.compare_phase_transition = input_object.get_compare_phase_transition(); 
    if (ref_object.compare_phase_transition)
    {
        ref_object.phase_transition_weight = input_object.get_ref_phase_transition_weight_c();
        ref_object.phase_transition = input_object.get_ref_phase_transition(); 
    }
}

isobar_c qce_iteration(double &temp, double &bxv, double &amf, double &vdamp, double& v0, isobar_c &ib, std::vector<cluster_c>&clusterset,  std::vector<pf_c>&lnq, cluster_c &monomer, int cyclus)
{
    double vol, gibbs; 
    std::vector<double>populations;
    bool converged;
    int iteration;
    bool success;
    double vdamp_local;
    double old_vol; 

    // initialize volume, populations, gibbs free enthalpy. 
    vol = v0;
    vdamp_local = vdamp;
    gibbs = std::numeric_limits<float>::max();
    int no_clusters = clusterset.size() - 1; 
    initialize_populations(populations, no_clusters);

    converged = false;
    for (int i = 0; i < global_data.qce_iterations; i++)
    {
        if (i == 0)
        {
            calculate_lnq(lnq, vol, temp, amf, bxv, clusterset);
        }
        else
        {
            update_lnq(lnq, amf, bxv, temp, vol, clusterset);
        }

        std::tie(success, populations) = calculate_populations(populations, clusterset, lnq, cyclus, monomer);
        
         if (not success)
        {
            vol = v0 * vdamp_local;
            vdamp_local = vdamp_local * vdamp;
            initialize_populations(populations, no_clusters);
            continue; 
        }
        
        old_vol = vol;

        std::tie(success, vol) = calculate_volume(vol, vdamp, amf, bxv, temp, populations);

        if (not success)
        {
            vol = v0 * vdamp_local;
            vdamp_local = vdamp_local * vdamp; 
            initialize_populations(populations, no_clusters);
            continue;
        }
        
        std::tie(converged, gibbs) = check_convergence(gibbs, vol, temp, populations, lnq);
        
        if (converged == true)
        {
            ib.iter_vol = vol;
            ib.iter_gibbs = gibbs;
            ib.iter_populations = populations; 
            ib.iter_lnq = lnq; 
            ib.iter_converged = converged; 
            populations.shrink_to_fit();
            lnq.shrink_to_fit();
            break; 
        }
    }
    return ib; 
}

typedef std::numeric_limits<double>dbl;
void qce_main(isobar_c &ib, std::vector<cluster_c>&clusterset, cluster_c &monomer, reference_c &reference)
{
    double v0, vdamp, vol, gibbs, error;
    bool converged, copy;
    std::vector<pf_c>lnq;
    lnq.resize(clusterset.size());
    std::vector<double>populations[clusterset.size()];
    int itemp = ib.get_temperature().size();
    std::cout.precision(dbl::max_digits10);

    // Cycle One - high T to low T (gas phase cycle)
    vdamp = 1.00 - global_data.vdamp; 
    ib.iter_converged = false;
    for (int i = itemp - 1; i > -1; i--)
    {
        if (not ib.iter_converged)
        {
            v0 = ideal_gas_volume(ib.get_temperature_element(i));
        }
        else
        {
            v0 = ib.iter_vol; 
        }

        ib = qce_iteration(ib.get_temperature_element(i), ib.get_bxv(), ib.get_amf(), vdamp, v0, ib, clusterset, lnq, monomer, 1);
        
        //copy results
        ib.set_gibbs(ib.iter_gibbs, i);
        ib.set_volume(ib.iter_vol, i);
        ib.set_lnq(ib.iter_lnq, i); 
        ib.set_converged(ib.iter_converged, i);
        ib.set_populations(ib.iter_populations, i); 
        ib.set_solution(1, i);
        if (ib.iter_converged)
        {
            int current_solution = ib.get_solution(i); 
            ib.set_solution(ib.get_solution(i) + 100, i); 
        }
    }

    // Cycle Two - low T to high T (liquid phase cycle)
    vdamp = 1.00 + global_data.vdamp; 
    ib.iter_converged = false; // this used to be converged = false. 

    for (int i = 0; i < itemp; i++)
    {
        if (not ib.iter_converged)
        {
            v0 = 1.0e-2 * ideal_gas_volume(ib.get_temperature_element(i));
        }
        else
        {
            v0 = fmin(ib.iter_vol, 1.00e-1f * ideal_gas_volume(ib.get_temperature_element(i)));
        }

        ib = qce_iteration(ib.get_temperature_element(i), ib.get_bxv(), ib.get_amf(), vdamp, v0, ib, clusterset, lnq, monomer, 2);

        if (ib.iter_converged)
        {
            int current_solution = ib.get_solution(i); 
            ib.set_solution(current_solution + 10, i);
            if (not ib.get_converged_element(i))
            { 
                copy = true; 
            }
            else if (ib.iter_gibbs < ib.get_gibbs(i))
            {
                copy = true;
            }
            else
            {
                copy = false; 
            }

            if (copy == true)
            {
                ib.set_gibbs(ib.iter_gibbs, i);
                ib.set_volume(ib.iter_vol, i);
                ib.set_lnq(ib.iter_lnq, i);
                ib.set_converged(ib.iter_converged, i);
                ib.set_populations(ib.iter_populations, i); 
                ib.set_solution(ib.get_solution(i) + 1, i);
            }
        }
    }

    // determine isobar quality, if necessary. 
    ib.set_error(0.00);

    if (reference.compare_density)
    {
        error = compare_density(ib, reference.density_temperature, reference.density);
        ib.density_error = error; 
        ib.set_error(ib.get_error() + reference.density_weight * error);
    }

    if (reference.compare_isobar)
    {
        error = compare_isobar(ib, reference.isobar_temperature, reference.isobar_volume);
        ib.isobar_error = error; 
        ib.set_error(ib.get_error() + reference.isobar_weight * error);
    }

    if (reference.compare_phase_transition)
    {
        error = compare_phase_transition(ib, reference.phase_transition);
        ib.pt_error = error; 
        ib.set_error(ib.get_error() + reference.phase_transition_weight * error); 
    }
}

isobar_c qce_start(std::vector<cluster_c>&clusterset, cluster_c &monomer, reference_c &reference)
{
    std::vector<cluster_c>clusterset_temp = clusterset; 
    int iamf, ibxv, itemp, nr_isobars_computed, nr_isobars_total; 
    isobar_c best_ib; 
    isobar_c ib;

    best_ib.set_error(std::numeric_limits<float>::max());
    nr_isobars_computed = 0;
    nr_isobars_total = global_data.amf.num * global_data.bxv.num;
    

    //write amf + bxv errors to a file;
    std::ofstream erroroutput;
    erroroutput.open("amf_bxverrors.dat");
    erroroutput << "amf" << std::setw(22) << "bxv" << std::setw(25) << "dens_error" << std::setw(22) << "isobar_error" << std::setw(16) << "pt_error" << std::setw(22) << "final error" << std::endl; 
    #pragma omp parallel for default(none) private(ib, iamf, ibxv, itemp) shared(global_data, avogadro, clusterset, monomer, best_ib, reference, nr_isobars_computed, nr_isobars_total, erroroutput) collapse(2), schedule (guided)
    for (int iamf = 0; iamf < global_data.amf.num; iamf++)
    {
        for (int ibxv = 0; ibxv < global_data.bxv.num; ibxv++)
        {
            // allocate sizes. 
            ib.resize_temperature(global_data.temp.num);
            ib.resize_volume(global_data.temp.num);
            ib.resize_gibbs(global_data.temp.num);
            ib.resize_converged(global_data.temp.num);
            ib.resize_solution(global_data.temp.num);
            ib.resize_lnq(global_data.temp.num, clusterset.size());
            ib.resize_populations(global_data.temp.num, clusterset.size());

            ib.set_amf(((global_data.amf.first + ((iamf) * global_data.amf.delta)) / pow(avogadro,2)));
            ib.set_bxv(global_data.bxv.first + ((ibxv) * global_data.bxv.delta));
            for (int itemp = 0; itemp < global_data.temp.num; itemp++)
            {
                double temp = global_data.temp.first + ((itemp) * global_data.temp.delta);
                ib.set_temperature(temp, itemp);
            }        
            qce_main(ib, clusterset, monomer, reference);
            erroroutput << std::scientific << std::uppercase << std::setprecision(3) << ib.get_amf() * pow(avogadro,2) << std::setw(20) << ib.get_bxv() << std::setw(20) << ib.density_error;
            erroroutput << std::setw(20) << std::setprecision(3) << ib.isobar_error << std::setw(20) << std::setprecision(3) << std::scientific << std::uppercase << ib.pt_error << std::setw(20) << std::scientific << std::uppercase << std::setprecision(3) << ib.get_error() << std::endl; 
            #pragma omp critical 
            global_data.nconverged = global_data.nconverged + count_converge(ib.get_converged());
            if (ib.get_error() < best_ib.get_error())
            {
                best_ib = ib; 
            }
           
            #pragma omp atomic
            nr_isobars_computed = nr_isobars_computed + 1; 

            // #if defined(_OPENMP)
            // if (omp_get_thread_num() == 0)
            // {
            //     progress_bar(nr_isobars_computed, nr_isobars_total, global_data.progress_bar, false);
            // }
            // #else
            // {
            //     progress_bar(nr_isobars_computed, nr_isobars_total, global_data.progress_bar, false);
            // }
            // #endif
        }
    }
    erroroutput.close(); 
    progress_bar(nr_isobars_computed, nr_isobars_total, global_data.progress_bar, true);
    return best_ib; 
}

void post_processing(int& temp_num, std::vector<cluster_c>&clusterset, isobar_c &best_ib, input_c &input)
{
    std::vector<pf_c>d; 
    std::vector<pf_c>dd; 
    std::vector<std::vector<pf_c>>dlnq_clust; 
    std::vector<std::vector<pf_c>>ddlnq_clust; 
    std::vector<double>dvol; 

    dvol.resize(temp_num);
    d.resize(clusterset.size());
    dd.resize(clusterset.size());
    for (int i = 0; i < temp_num; i++)
    {
        calculate_dlnq(d, clusterset, best_ib.get_amf(), best_ib.get_temperature_element(i), best_ib.get_volume(i));
        calculate_ddlnq(dd, clusterset, best_ib.get_amf(), best_ib.get_temperature_element(i), best_ib.get_volume(i));
        dlnq_clust.push_back(d);
        ddlnq_clust.push_back(dd);
    }
    dvol = derivative(best_ib.get_volume_vec(), best_ib.get_temperature());

    // calculate and write thermodynamic quantities
    int clusterset_size = clusterset.size();
    std::vector<double>temp = best_ib.get_temperature(); 
    double best_ib_bxv = best_ib.get_bxv() * global_data.vexcl; 
    std::vector<int>best_ib_solution = best_ib.get_solution_vec(); 
    std::vector<std::vector<double>>best_ib_populations = best_ib.get_populations_vec(); 
    std::vector<double>volume_vec = best_ib.get_volume_vec(); 
    std::vector<std::vector<pf_c>>lnq_vec = best_ib.get_lnq_vector();
    std::vector<bool>best_ib_converged = best_ib.get_converged(); 

    calculate_thermo(temp_num, clusterset_size, best_ib.get_populations_vec(), global_data.press, temp, best_ib.get_volume_vec(), best_ib_bxv, lnq_vec, dlnq_clust, ddlnq_clust, dvol, best_ib_solution, best_ib.get_converged(), input);
    calculate_cluster_entropy("Stot_contrib.dat", best_ib, lnq_vec, dlnq_clust, clusterset, "qtot");
    calculate_cluster_entropy("Strans_contrib.dat", best_ib, lnq_vec, dlnq_clust, clusterset, "qtrans");
    calculate_cluster_entropy("Svib_contrib.dat", best_ib, lnq_vec, dlnq_clust, clusterset, "qvib");
    calculate_cluster_entropy("Srot_contrib.dat", best_ib, lnq_vec, dlnq_clust, clusterset, "qrot");
    calculate_cluster_entropy("Sint_contrib.dat", best_ib, lnq_vec, dlnq_clust, clusterset, "qint");

    // normalise populations
    std::vector<std::vector<double>>pop_norm; 
    double sum_of_composition = 0.00; 
    pop_norm.resize(temp_num, std::vector<double>(clusterset_size));

    for (int i = 0; i < clusterset.size(); i++)
    {
        sum_of_composition += clusterset[i].get_composition(); 
    }

    for (int i = 0; i < temp_num; i++)
    {
        for (int j = 0; j < clusterset_size; j++)
        {
            pop_norm[i][j] = best_ib_populations[i][j] * clusterset[j].get_composition() / global_data.ntot[0];
        }
    }

    // Write populations
    std::ofstream population_output;
    population_output.open("populations.dat");
    population_output << std::left << "#" << std::setw(20) << "T/K" << std::setw(20);
    for (int i = 0; i < clusterset_size; i++)
    {
        population_output << std::setw(20) << clusterset[i].get_label();
        if (i == clusterset_size - 1)
        {
            population_output << std::endl; 
        }
    }

    for (int j = 0; j < temp_num; j++)
    {
        population_output << std::left << std::setw(20) << std::setprecision(6) << std::scientific << std::uppercase << temp[j];
        for (int k = 0; k < clusterset_size; k++)
        {
            population_output << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << pop_norm[j][k];
            if (k == clusterset_size - 1)
            {
                population_output << '\n';
            }
        }
    }
    population_output.close();

    // calculate concentrations
    std::vector<std::vector<double>>conc; 
    conc.resize(temp_num, std::vector<double>(clusterset_size));
    for (int i = 0; i < temp_num; i++)
    {
        for (int j = 0; j < clusterset_size; j++)
        {
            conc[i][j] = best_ib_populations[i][j] / (avogadro * best_ib.get_volume(i) * 1.0e3);
        }
    }

    // write concentrations 
    std::ofstream concentration_output; 
    concentration_output.open("concentrations.dat");
    concentration_output << std::left << "#" << std::setw(20) << "T/K" << std::setw(20);
    for (int i = 0; i < clusterset_size; i++)
    {
        concentration_output << std::setw(20) << clusterset[i].get_label(); 
        if (i == clusterset_size - 1)
        {
            concentration_output << std::endl; 
        }
    }
    for (int j = 0; j < temp_num; j++)
    {
        concentration_output << std::left << std::setw(20) << std::scientific << std::setprecision(6) << std::uppercase << temp[j];
        for (int k = 0; k < clusterset_size; k++)
        {
            concentration_output << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << conc[j][k];
            if (k == clusterset_size - 1)
            {
                concentration_output << '\n';
            }
        }
    }
    concentration_output.close();
}  

void print_input(input_c &input_object, isobar_c &best_ib, reference_c &reference)
{
    std::ofstream outputfile;
    outputfile.open("output.dat");
    outputfile << "Using the following input: " << std::endl; 
    outputfile << std::endl; 
    outputfile << "[system]" << std::endl;
    outputfile << '\t' << "components: " << input_object.get_components() << std::endl; 
    outputfile << std::endl;
    outputfile << "[ensemble]" << std::endl;
    outputfile << '\t' << "pressure: " << input_object.get_pressure() << " [bar]" << std::endl;
    outputfile << '\t' << "temperature: " << input_object.get_temperature().first << " ... " << input_object.get_temperature().last << " [K]" << std::endl; 
    outputfile << '\t' << "monomer amounts: " << std::setprecision(4) << input_object.get_monomer_amounts()[0] << " [mol]" << std::endl; 
    outputfile << std::endl;
    outputfile << "[thermoinput]" << std::endl;
    outputfile << '\t' << "amf: " <<  best_ib.get_amf()  * pow(avogadro,2) << " [J*m^3/mol^2]" << std::endl;
    outputfile << '\t' << "bxv: " << best_ib.get_bxv() << std::endl; 
    outputfile << '\t' << "error: " << std::setprecision(6) << best_ib.get_error() << std::endl; 
    outputfile << '\t' << "calculated phase transition: " << std::setprecision(5) << determine_phase_transition(best_ib) << " K" << std::endl;
    outputfile << '\t' << "calculated density: " << std::setprecision(6) << 1.0e3 * global_data.mtot / (1.0e6 * determine_volume(best_ib, reference.density_temperature)) << "g/cm^3" << std::endl; 
    outputfile << '\t' << "Number of converged iterations: " << count_converge(best_ib.get_converged()) << std::endl; 
    outputfile << '\t' << "Number of converged iterations (total): " << global_data.nconverged << "/" << global_data.amf.num * global_data.bxv.num * global_data.temp.num << std::endl; 
    outputfile << '\t' << "maximum number of THERMOSOLVR iterations: " << global_data.qce_iterations << std::endl;
    outputfile << '\t' << "maximum number of Newton-Raphson iterations: " << global_data.newton_iterations << std::endl; 
    outputfile << '\t' << "maximum relative deviation: " << std::setprecision(6) << std::scientific << std::uppercase << global_data.max_deviation << std::endl; 
    outputfile << '\t' << "volume damping factor: " << std::setprecision(6) << std::scientific << std::uppercase << global_data.vdamp << std::endl; 
    outputfile << std::endl; 
    outputfile.close(); 
}

void print_clusterset(std::vector<cluster_c>&clusterset, cluster_c &monomer)
{
    std::ofstream finalclusterset;
    finalclusterset.open("clusterset.dat");
    finalclusterset << "Using the following clusterset: " << std::endl;
    finalclusterset << std::endl; 
    for (int i = 0; i < clusterset.size(); i++)
    {
        finalclusterset << "[" << clusterset[i].get_label() << "]" << std::endl;
        finalclusterset << '\t' << "composition: " << clusterset[i].get_composition()<< std::endl; 
        finalclusterset << '\t' << "sigma: " << clusterset[i].get_sigma() << std::endl;
        finalclusterset << '\t' << "energy: " << clusterset[i].get_energy() << " [kJ/mol]" << std::endl;
        if (i == 0)
        {
            finalclusterset << '\t' << "volume: " << monomer.get_volume() << " [A^3]" << std::endl;
        }
        else
        {
            finalclusterset << '\t' << "volume: " << monomer.get_volume() * clusterset[i].get_composition() << " [A^3]" << std::endl; 
        }
        finalclusterset << '\t' << "mass: " << std::fixed << std::setprecision(4) << clusterset[i].get_molmass() << " [amu]" << std::endl; 

        finalclusterset << '\t' << "inertia: " << std::fixed << std::setprecision(4) << clusterset[i].get_momofinertia(2) << " ";
        finalclusterset << clusterset[i].get_momofinertia(1) << " ";
        finalclusterset << clusterset[i].get_momofinertia(0) << " [amu*Angstrom^2]" << std::endl;
        
        finalclusterset << '\t' << "frequencies: " << std::fixed << std::setprecision(4) << clusterset[i].get_frequency_element(0) << " ... " << clusterset[i].get_frequency_element(clusterset[i].get_frequencies().size()-1) << " [1/cm]" << std::endl;
        
    
        finalclusterset << std::endl; 
    }
    finalclusterset.close(); 
}

void qce_finalize(isobar_c &best_ib, std::vector<cluster_c>&clusterset, reference_c &reference, input_c &input, cluster_c &monomer, int option, int clusterset_no)
{
    if (option == 1)
    {
        namespace fs = std::filesystem;
        std::cout << "clusterset " << std::to_string(clusterset_no) << ": [ "; 
        for (int i = 0; i < clusterset.size(); i++)
        {
            std::cout << clusterset[i].get_label() << " ";
            if (i == clusterset.size() - 1)
            {
                std::cout << "]" << std::endl; 
            }
        }
         // perform post processing
        post_processing(global_data.temp.num, clusterset, best_ib, input);
        print_clusterset(clusterset, monomer);
        print_input(input, best_ib, reference);
        fs::path p = fs::current_path();

        // convert clusterset directory from string to path
        std::string clustersetpath_string = "test_clustersets/clusterset" + std::to_string(clusterset_no);
        const std::filesystem::path clustersetpath = clustersetpath_string; 

        // move *.dat files for each clusterset into their corresponding directories. 
        std::vector<std::filesystem::path>outputfiles;
        outputfiles = getdatfiles(p, ".dat");
        for (auto const & entry : outputfiles)
        {
            if (std::filesystem::exists(p/entry))
            {
                fs::rename(p/entry, p/clustersetpath/entry);
            }
        } 
    }

    if (option == 2)
    {
        std::cout << "Number of converged iterations: " << count_converge(best_ib.get_converged()) << std::endl; 

        if (reference.compare)
        {
            std::cout << "Number of converged iterations (total): " << global_data.nconverged << "/" << global_data.amf.num * global_data.bxv.num * global_data.temp.num << std::endl; 
        }

        std::cout << "Best isobar found for: amf =" << std::setprecision(6) << best_ib.get_amf() * pow(avogadro, 2.00) << " J*m^3/mol^2, bxv = " << best_ib.get_bxv() << std::endl;
        std::cout << "Error: " << std::setprecision(6) << std::scientific << std::uppercase << best_ib.get_error() << std::endl;
        if (reference.compare_phase_transition)
        {
            std::cout << "Calculated phase transition: " << determine_phase_transition(best_ib) << " K" << std::endl; 
        }
        
        if (reference.compare_density)
        {
            std::cout << "Calculated density: " << 1.0e3 * global_data.mtot / (1.0e6 * determine_volume(best_ib, reference.density_temperature)) << " g/cm3" << std::endl; 
        }
        std::cout << std::endl;

        // perform post processing
        post_processing(global_data.temp.num, clusterset, best_ib, input);
        print_clusterset(clusterset, monomer);
        print_input(input, best_ib, reference);

        // create folder for dat files
        std::filesystem::create_directories("outputfiles"); 
        
        // get outputfiles folder path
        std::filesystem::path p = std::filesystem::current_path();
        std::filesystem::path outputfilepath = p/"outputfiles";
        
        // get all *dat files and move into outputfiles folder. 
        std::vector<std::filesystem::path>outputfiles;
        outputfiles = getdatfiles(p, ".dat");

        for (auto const & entry : outputfiles)
        {
            if (std::filesystem::exists(p/entry))
            {
                std::filesystem::rename(p/entry, p/outputfilepath/entry);
            }
        } 

    }
}
