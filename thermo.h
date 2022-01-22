#include<vector>
#include<string>
#include "partfunctions.h"
#include "input.h"


std::vector<double>calculate_internal_energy(int &temp_num, std::vector<double>&temp, std::vector<double>&dlnq, std::vector<double>&u);
std::vector<double>calculate_enthalpy(int &temp_num, std::vector<double>&temp, std::vector<double>&dlnq, std::vector<double>&vol, double &press, std::vector<double>&h);
std::vector<double>calculate_entropy(int &temp_num, std::vector<double>&temp, std::vector<double>&lnq, std::vector<double>&dlnq, std::vector<double>&s);
std::vector<double>calculate_cv(int &temp_num, std::vector<double>&temp, std::vector<double>&dlnq, std::vector<double>&ddlnq, std::vector<double>&cv);
std::vector<double>calculate_cp(int &temp_num, std::vector<double>&temp, std::vector<double>&dlnq, std::vector<double>&ddlnq, std::vector<double>&dvol, double &press, std::vector<double>&cp);
std::vector<double>calculate_gibbs_enthalpy(int &temp_num, std::vector<double>&temp, std::vector<double>&lnq, std::vector<double>&vol, double &press, std::vector<double>&g);
std::vector<double>calculate_helmholtz_energy(int &temp_num, std::vector<double>&temp, std::vector<double>&lnq, std::vector<double>&a);
void calculate_cluster_entropy(std::string filename, isobar_c &ib, std::vector<std::vector<pf_c>>&lnq_clust, std::vector<std::vector<pf_c>>&dlnq_clust, std::vector<cluster_c>&clusterset, std::string pftype);
void write_contributions(int &temp_num, std::vector<bool>&conevrged, std::vector<double>&temp, std::vector<double>&a, std::vector<double>&vib, std::vector<double>&rot, std::vector<double>&trans, std::vector<double>&elec, std::vector<double>&mf, std::vector<double>&indi, std::string header, double factor, std::string filename);
void write_volume(int &temp_num, std::vector<bool>&converged, std::vector<double>&temp, std::vector<double>&vol, double &vexcl, std::vector<double>&alpha, std::vector<int>&solution);
void write_thermo0(int &temp_num, std::vector<bool>&converged, std::vector<double>&temp, std::vector<double>&a, std::vector<double>&g);
void write_thermo1(int &temp_num, std::vector<bool>&converged, std::vector<double>&temp, std::vector<double>&u, std::vector<double>&h, std::vector<double>&s);
void write_thermo2(int &temp_num, std::vector<bool>&converged, std::vector<double>&temp, std::vector<double>&cp, std::vector<double>&cv);
std::vector<double>calculate_expansion_coefficient(int& temp_num, std::vector<double>&vol, std::vector<double>&dvol, std::vector<double>&alpha);
std::vector<double>add_lnq_indi(int& temp_num, int &clusterset_size, std::vector<std::vector<double>>&populations, std::vector<double>&lnq_sys);
std::vector<double>calculate_lnq_sys(int& temp_num, int &clusterset_size, std::vector<std::vector<double>>&populations, std::vector<std::vector<pf_c>>&lnq_clust, std::string lnq_type);
void calculate_thermo(int &temp_num, int &clusterset_size, std::vector<std::vector<double>>&populations, double &press, std::vector<double>&temp, std::vector<double>&vol, double &vexcl, std::vector<std::vector<pf_c>>&lnq_clust, std::vector<std::vector<pf_c>>&dlnq_clust, std::vector<std::vector<pf_c>>&ddlnq_clust, std::vector<double>&dvol,std::vector<int>&solution, std::vector<bool>&converged, input_c &input);














