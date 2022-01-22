#ifndef QCE_H
#define QCE_H
#include "global.h"
#include "partfunctions.h"
#include "input.h"
#include "cluster2.h"

class reference_c
{
    public:
    bool compare, compare_isobar, compare_density, compare_phase_transition;
    double density_weight, isobar_weight, phase_transition_weight;
    double phase_transition, density, density_temperature;
    std::vector<double>isobar_temperature, isobar_volume; 
};

class isobar_c
{
    double amf_c, bxv_c, error_c;
    std::vector<int>solution_c;
    std::vector<bool>converged_c;
    std::vector<double>vol_c;
    std::vector<double>temp_c;
    std::vector<double>gibbs_c;
    std::vector<std::vector<double>>populations_c;
    std::vector<std::vector<pf_c>>lnq_c; 
    
    public:
    double density_error, pt_error, isobar_error; 
    double iter_vol, iter_gibbs, iter_converged;
    double predicted_bp;
    std::vector<double>iter_populations; 
    std::vector<pf_c>iter_lnq; 
    void set_error(const double &error); 
    void set_amf(const double &amf);
    void set_bxv(const double &bxv);
    void set_temperature(double temperature, int i);
    void set_volume(double vol, int i);
    void set_gibbs(double gibbs, int i);
    void set_converged(bool converged, int i);
    void set_populations(std::vector<double>populations, int& i);
    void set_solution(int solution, int &i);
    void set_lnq(std::vector<pf_c>lnq, int& i);
    void resize_temperature(int &num);
    void resize_volume(int &num);
    void resize_gibbs(int &num);
    void resize_converged(int &num);
    void resize_solution(int &num);
    void resize_lnq(int &num, int clustersize);
    void resize_populations(int &num, int clustersize);

    double& get_amf();
    double& get_bxv();
    std::vector<double>get_temperature();
    double& get_temperature_element(int i);
    double& get_gibbs(int i);
    double& get_volume(int i);
    int get_solution(int i);
    double &get_error(); 
    int get_converged_element(int i);
    std::vector<bool>& get_converged();
    std::vector<int>get_solution_vec();
    std::vector<std::vector<double>>& get_populations_vec();
    std::vector<double>& get_volume_vec();
    std::vector<double>& get_populations(int i);
    std::vector<pf_c>& get_lnq(int i);
    std::vector<std::vector<pf_c>>& get_lnq_vector();
};

    std::tuple<bool, std::vector<double> >boolvecttuple(bool success, std::vector<double>population);
    void initialize_conserved_quantities(cluster_c &monomer);
    void qce_prepare(input_c input_object, cluster_c &monomer, std::vector<cluster_c>&clusterset, reference_c &ref_object);
    void initialize_degree(std::vector<cluster_c>&clusterset);
    double ideal_gas_volume(double &temp);
    std::tuple<bool, std::vector<double> >calculate_populations(std::vector<double>&populations, std::vector<cluster_c>&clusterset, std::vector<pf_c>&lnq, int &cyclus, cluster_c &monomer);
    std::vector<double>calculate_remaining_populations(std::vector<double>&monomer_populations, std::vector<cluster_c>&clusterset, std::vector<pf_c>&lnq);
    isobar_c qce_start(std::vector<cluster_c>&clusterset,cluster_c &monomer, reference_c &reference);
    void qce_main(isobar_c &ib, std::vector<cluster_c>&clusterset, cluster_c &monomer, reference_c &reference);
    isobar_c qce_iteration(double &temp, double &bxv, double &amf, double &vdamp, double& v0, isobar_c &ib, std::vector<cluster_c>&clusterset, std::vector<pf_c>&lnq, cluster_c &monomer, int cyclus);
    std::tuple<bool, double>calculate_volume(double &vol, double &vdamp, double &amf, double &bxv, double &temp, std::vector<double>&populations);
    std::tuple<bool, double>check_convergence(double &gibbs, double &vol, double &temp, std::vector<double>&populations, std::vector<pf_c>&lnq);
    double determine_volume(isobar_c &ib, double &temp);
    double compare_density(isobar_c& ib, double &temp, double &density);
    double compare_isobar(isobar_c& ib, std::vector<double>&temp, std::vector<double>&vol);
    double determine_phase_transition(isobar_c& ib);
    double compare_phase_transition(isobar_c& ib, double &pt);
    void post_processing(int& temp_num, std::vector<cluster_c>&clusterset, isobar_c &best_ib);
    void print_clusterset(std::vector<cluster_c>&clusterset, cluster_c &monomer);
    void print_input(input_c &input_object, isobar_c &best_ib, reference_c &reference);
    void qce_finalize(isobar_c &best_ib, std::vector<cluster_c>&clusterset, reference_c &reference, input_c &input, cluster_c &monomer, int option, int clusterset_no);

#endif
