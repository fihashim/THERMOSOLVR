#ifndef PARTFUNCTIONS_H
#define PARTFUNCTIONS_H
#include "cluster2.h"
// class global_c; 
// class cluster_c; 

class pf_c
{
    double qtrans_c, qvib_c, qrot_c, qelec_c, qint_c, qtot_c; 
    
    public:
    void set_qtrans(double qtrans);
    void set_qvib(double qvib);
    void set_qrot(double qrot);
    void set_qelec(double qelec);
    void set_qint(double qint);
    void set_qtot(double qtot);
    double& get_qtrans();
    double get_qvib();
    double get_qrot();
    double get_qelec();
    double get_qint(); 
    double get_qtot(); 
    void calculate_dlnq();
    void calculate_ddlnq(); 
};

void calculate_lnqtrans(pf_c &lnq, double &vol, double &temp, double &bxv, cluster_c &cluster);
void calculate_dlnqtrans(pf_c &dlnq, double &temp);
void calculate_ddlnqtrans(pf_c &ddlnq, double &temp);
void calculate_lnqvib(pf_c &lnq, double &temp, cluster_c &cluster);
void calculate_dlnqvib(pf_c &ddlnq, double &temp, cluster_c &cluster);
void calculate_ddlnqvib(pf_c &ddlnq, double &temp, cluster_c &cluster);
void calculate_lnqrot(pf_c &lnq, double &temp, cluster_c &cluster);
void calculate_dlnqrot(pf_c &dlnq, double &temp, cluster_c &cluster);
void calculate_ddlnqrot(pf_c &dlnq, double &temp, cluster_c &cluster);
void calculate_lnqelec(pf_c &lnq, double &temp, cluster_c &cluster);
void calculate_dlnqelec(pf_c &dlnq, double &temp, cluster_c &cluster);
void calculate_ddlnqelec(pf_c &dlnq, double &temp, cluster_c &cluster);
void calculate_lnqint(pf_c &lnq, double &temp, double &amf, double &vol, cluster_c &cluster);
void calculate_dlnqint(pf_c &dlnq, double &temp, double &amf, double &vol, cluster_c &cluster);
void calculate_ddlnqint(pf_c &ddlnq, double &temp, double &amf, double &vol, cluster_c &cluster);
double calculate_emf(int &composition, double &amf, double &vol,cluster_c &cluster);
void calculate_lnq(std::vector<pf_c>&lnq, double &vol, double &temp, double &amf, double &bxv, std::vector<cluster_c>&clusterset);
void update_lnq(std::vector<pf_c>&lnq, double &amf, double &bxv, double &temp, double &vol, std::vector<cluster_c>&clusterset);
void calculate_dlnq(std::vector<pf_c>&dlnq, std::vector<cluster_c>&clusterset, double &amf, double &temp, double &vol);
void calculate_ddlnq(std::vector<pf_c>&ddlnq, std::vector<cluster_c>&clusterset, double &amf, double &temp, double &vol);


#endif