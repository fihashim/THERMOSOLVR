#include "partfunctions.h"
#include "cluster2.h"
#include "atomic_data.h"
#include "global.h"
#include<vector>
#include<iostream>
#include<iomanip>


void pf_c::set_qtrans(double qtrans)
{
    qtrans_c = qtrans; 
}

void pf_c::set_qvib(double qvib)
{
    qvib_c = qvib; 
}

void pf_c::set_qrot(double qrot)
{
    qrot_c = qrot; 
}

void pf_c::set_qelec(double qelec)
{
    qelec_c = qelec; 
}

void pf_c::set_qint(double qint)
{
    qint_c = qint; 
}

void pf_c::set_qtot(double qtot)
{
    qtot_c = qtot; 
}

double& pf_c::get_qtrans()
{
    return qtrans_c;
}

double pf_c::get_qelec()
{
    return qelec_c;
}

double pf_c::get_qint()
{
    return qint_c; 
}

double pf_c::get_qvib()
{
    return qvib_c;
}

double pf_c::get_qrot()
{
    return qrot_c;
}

double pf_c::get_qtot()
{
    return qtot_c; 
}

void calculate_lnqtrans(pf_c &lnq, double &vol, double& temp, double &bxv, cluster_c &cluster)
{
    double lnqtrans; 
    double lambda;
    double mass;
    // calculate de Broglie wavelength
    mass = cluster.get_molmass() * amu; 
    lambda = planck / sqrt(2.0 * pi * mass * kb * temp);
    lnqtrans = log(vol - (bxv * global_data.vexcl)) - 3.0 * log(lambda);
    lnq.set_qtrans(lnqtrans);
}

void calculate_lnqvib(pf_c &lnq, double &temp, cluster_c &cluster)
{
    std::vector<double>frequencies = cluster.get_frequencies();
    double lnqvib = 0.00; 
    double anharmonicity = cluster.get_anharmonicity(); 
    for (int j = 0; j < frequencies.size(); j++)
    {
        double t_vib = 0.00; 
        t_vib = frequencies[j] * factor;
        if (anharmonicity > 0.00)
        {
            // morse oscillator
            lnqvib = lnqvib - log(2.0 * sinh(0.5 * t_vib / temp)) + anharmonicity * (t_vib / temp)* (0.25 + 0.5 / pow(sinh(0.5 * t_vib/temp), 2));  
        }
        else
        {
            // harmonic oscillator
            lnqvib = lnqvib - 0.5 * t_vib/temp - log(1.0 - exp(-t_vib/temp));
        }
    }
    lnq.set_qvib(lnqvib);
}

void calculate_lnqrot(pf_c& lnq, double &temp, cluster_c &cluster)
{
    double lnqrot; 
    double inertia_temp[3];
   
    if (cluster.get_atom())
    {
        lnqrot = 0.00; 
    }
    else
    {
        lnqrot = 1.00; 
        inertia_temp[0] = cluster.get_momofinertia(0); 
        inertia_temp[1] = cluster.get_momofinertia(1);
        inertia_temp[2] = cluster.get_momofinertia(2);

        for (int j = 0; j < 3; j++)
        {
            double t_rot;
            t_rot = pow(hbar, 2) / (2.0 * inertia_temp[j] * amu * 1.0e-20 * kb);
            lnqrot = lnqrot * temp/t_rot;
        }
        if (cluster.get_linear())
        {
            lnqrot = sqrt(lnqrot / cluster.get_sigma());
        }
        else
        {
            lnqrot = sqrt(pi * lnqrot) / cluster.get_sigma(); 
        }
        lnqrot = log(lnqrot); 
        lnq.set_qrot(lnqrot);
    }
}

void calculate_lnqelec(pf_c &lnq, double &temp, cluster_c &cluster)
{
    double lnqelec; 
    double cluster_energy;
    cluster_energy = cluster.get_energy(); 
    lnqelec = (-1000 / avogadro) * cluster_energy / (kb * temp);
    lnq.set_qelec(lnqelec);
}

double calculate_emf(int &composition, double &amf, double &vol, cluster_c &cluster)
{
    double emf; 
    emf = -1.00 / vol * global_data.ntot[0] * cluster.get_composition() * amf; 
    return emf; 
}

void calculate_lnqint(pf_c &lnq, double &temp, double &amf, double &vol, cluster_c &cluster)
{
    double lnqint; 
    double emf; 
    int composition = cluster.get_composition();
    emf = -1.00 / vol * global_data.ntot[0] * composition * amf; // calculate_emf
    lnqint = -emf / (kb * temp);
    lnq.set_qint(lnqint); 
}

void calculate_dlnqtrans(pf_c &dlnq, double &temp)
{
    dlnq.set_qtrans(1.50 / temp);
}

void calculate_dlnqvib(pf_c &dlnq, double &temp, cluster_c &cluster)
{
    std::vector<double>frequencies = cluster.get_frequencies();
    double x = cluster.get_anharmonicity();
    double factor, t_vib, fsinh, fcosh, q;
    factor = planck * 100.00 * speed_of_light / kb;
    for (int i = 0; i < frequencies.size(); i++)
    {
        t_vib = factor * frequencies[i];
        if (x > 0.00)
        {
            fsinh = sinh(0.50 * t_vib / temp);
            fcosh = cosh(0.40 * t_vib / temp);
            q = q + 0.50 * t_vib * fcosh / (pow(temp,2.00) * fsinh) - x * t_vib / (pow(temp,2.00) * (0.25 + 0.5 / pow(fsinh, 2.00))) + 0.50 * x * pow(t_vib, 2.00) * fcosh / (pow(temp, 3.00) * pow(fsinh, 3.00));
            dlnq.set_qvib(q);
        }
        else
        {
            q = q + 0.50 * t_vib / pow(temp, 2.00) + t_vib / pow(temp,2.00) / (exp(t_vib / temp) - 1.00);
            dlnq.set_qvib(q);
        }
    }
}

void calculate_dlnqrot(pf_c &dlnq, double &temp, cluster_c &cluster)
{
    double q; 
    if (cluster.get_atom())
    {
        dlnq.set_qrot(0.00);
    }
    else if (cluster.get_linear())
    {
        dlnq.set_qrot(1.0 / temp);
    }
    else
    {
        dlnq.set_qrot(1.5 / temp);
    }
}

void calculate_dlnqelec(pf_c &dlnq, double &temp, cluster_c &cluster)
{
    dlnq.set_qelec((1000 / avogadro) * cluster.get_energy() / (kb * pow(temp,2.00)));
}

void calculate_dlnqint(pf_c &dlnq, double &temp, double &amf, double &vol, cluster_c &cluster)
{
    int composition = cluster.get_composition();
    double emf = calculate_emf(composition, amf, vol, cluster);
    dlnq.set_qint(emf / (kb * pow(temp,2.00)));
}

void calculate_ddlnqtrans(pf_c &ddlnq, double &temp)
{
    ddlnq.set_qtrans(-1.50 / pow(temp,2.00));
}

void calculate_ddlnqvib(pf_c &ddlnq, double &temp, cluster_c &cluster)
{
    std::vector<double>frequencies = cluster.get_frequencies();
    double factor, t_vib, fsinh, fcosh, q;
    double x = cluster.get_anharmonicity();
    factor = planck * 100.00 * speed_of_light / kb;
    
    for (int i = 0; i < frequencies.size(); i++)
    {
        t_vib = factor * frequencies[i];
        if (x > 0.0)
        {
            // Morse oscillator
            fsinh = sinh(0.5 * t_vib / temp);
            fcosh = cosh(0.5 * t_vib / temp);
            q = q - fcosh / fsinh * t_vib / (pow(temp,3.00)) - 2.0 * x * pow(t_vib,2.00) / pow(temp,4.00) * fcosh / pow(fsinh, 3.00) + pow(t_vib / (2.0 * pow(temp,2.00) * fsinh), 2.00) + 2.0 * x * t_vib / pow(temp,3.00) * (0.25 + 0.5 / pow(fsinh,2.00)) + x * pow(t_vib,3.00) / (4.0 * pow(temp, 5.00)) * (2.0 * pow(fcosh,2.00) + 1.0) / pow(fsinh,4.00);
            ddlnq.set_qvib(q);
        }
        else
        {
            // Harmonic oscillator
            q = q - t_vib / pow(temp,3.00) - 2.0 * t_vib / (pow(temp,3.00) * (exp(t_vib / temp) - 1.00)) + 
            pow(t_vib,2.00) * exp(t_vib/temp) / (pow(temp,4.00) * pow((exp(t_vib/temp) - 1),2.00));
            ddlnq.set_qvib(q);
        }
    }
}

void calculate_ddlnqrot(pf_c &ddlnq, double &temp, cluster_c &cluster)
{
    if (cluster.get_atom())
    {
        ddlnq.set_qrot(0.00);
    }
    else if (cluster.get_linear())
    {
        ddlnq.set_qrot(-1.0 / pow(temp,2.00));
    }
    else
    {
        ddlnq.set_qrot(-1.5 / pow(temp, 2.00));
    }
}

void calculate_ddlnqelec(pf_c &ddlnq, double &temp, cluster_c &cluster)
{
    ddlnq.set_qelec(-(2000 / avogadro) * cluster.get_energy() / (kb * pow(temp,3.00)));
}

void calculate_ddlnqint(pf_c &ddlnq, double &temp, double &amf, double &vol, cluster_c &cluster)
{
    int composition = cluster.get_composition();
    double emf;
    emf = calculate_emf(composition,amf, vol, cluster);
    ddlnq.set_qint(-2.0 * emf / (kb * pow(temp,3.00)));
}

void calculate_lnq(std::vector<pf_c>&lnq, double &vol, double &temp, double &amf, double &bxv, std::vector<cluster_c>&clusterset)
{
    lnq.resize(clusterset.size());
    for (int i = 0; i < clusterset.size(); i++)
    {
        double lnqtot; 
        double lnqtrans;
        double lnqvib;
        double lnqrot;
        double lnqelec;
        double lnqint;
        calculate_lnqtrans(lnq[i], vol, temp, bxv, clusterset[i]); 
        calculate_lnqvib(lnq[i], temp, clusterset[i]);
        calculate_lnqrot(lnq[i], temp, clusterset[i]);
        calculate_lnqelec(lnq[i], temp, clusterset[i]);
        calculate_lnqint(lnq[i], temp, amf, vol, clusterset[i]);

        lnqtrans = lnq[i].get_qtrans();
        lnqvib = lnq[i].get_qvib(); 
        lnqrot = lnq[i].get_qrot(); 
        lnqelec = lnq[i].get_qelec();
        lnqint = lnq[i].get_qint();
        lnqtot = lnqtrans + lnqvib + lnqrot + lnqelec + lnqint; 
        lnq[i].set_qtot(lnqtot);
    }

}

void update_lnq(std::vector<pf_c>&lnq, double &amf, double &bxv, double &temp, double &vol, std::vector<cluster_c>&clusterset)
{
    for (int i = 0; i < lnq.size(); i++)
    {
        double lnqtot; 
        double lnqtrans;
        double lnqvib;
        double lnqrot;
        double lnqelec;
        double lnqint;
        calculate_lnqtrans(lnq[i], vol, temp, bxv, clusterset[i]);
        calculate_lnqint(lnq[i], temp, amf, vol, clusterset[i]);
        lnqtrans = lnq[i].get_qtrans();
        lnqvib = lnq[i].get_qvib(); 
        lnqrot = lnq[i].get_qrot(); 
        lnqelec = lnq[i].get_qelec();
        lnqint = lnq[i].get_qint();
        lnqtot = lnq[i].get_qtot(); 
        lnqtot = lnqtrans + lnqvib + lnqrot + lnqelec + lnqint; 
        lnq[i].set_qtot(lnqtot);
    }
}

void calculate_dlnq(std::vector<pf_c>&dlnq, std::vector<cluster_c>&clusterset, double &amf, double &temp, double &vol)
{
    dlnq.resize(clusterset.size());
    for (int i = 0; i < dlnq.size(); i++)
    {
        double dlnqtrans;
        double dlnqvib;
        double dlnqrot;
        double dlnqelec;
        double dlnqint; 
        double dlnqtot;
        calculate_dlnqtrans(dlnq[i],temp);
        calculate_dlnqvib(dlnq[i], temp, clusterset[i]);
        calculate_dlnqrot(dlnq[i], temp, clusterset[i]);
        calculate_dlnqelec(dlnq[i], temp, clusterset[i]);
        calculate_dlnqint(dlnq[i], temp, amf, vol, clusterset[i]);

        dlnqtrans = dlnq[i].get_qtrans();        
        dlnqvib = dlnq[i].get_qvib();
        dlnqrot = dlnq[i].get_qrot();
        dlnqelec = dlnq[i].get_qelec();
        dlnqint = dlnq[i].get_qint();
        dlnqtot = dlnqtrans + dlnqvib + dlnqrot + dlnqelec + dlnqint; 
        dlnq[i].set_qtot(dlnqtot);
    }
}

void calculate_ddlnq(std::vector<pf_c>&ddlnq, std::vector<cluster_c>&clusterset, double &amf, double &temp, double &vol)
{
    ddlnq.resize(clusterset.size());
    for (int i = 0; i < ddlnq.size(); i++)
    {
        double ddlnqtrans;
        double ddlnqvib;
        double ddlnqelec;
        double ddlnqrot;
        double ddlnqint;
        double ddlnqtot; 
        calculate_ddlnqtrans(ddlnq[i], temp);
        calculate_ddlnqvib(ddlnq[i], temp, clusterset[i]);
        calculate_ddlnqrot(ddlnq[i], temp, clusterset[i]);
        calculate_ddlnqelec(ddlnq[i], temp, clusterset[i]);
        calculate_ddlnqint(ddlnq[i], temp, amf, vol, clusterset[i]);
        ddlnqtrans = ddlnq[i].get_qtrans();
        ddlnqvib = ddlnq[i].get_qvib();
        ddlnqrot = ddlnq[i].get_qrot();
        ddlnqelec = ddlnq[i].get_qelec();
        ddlnqint = ddlnq[i].get_qint();
        ddlnqtot = ddlnqtrans + ddlnqvib + ddlnqrot + ddlnqelec + ddlnqint;
        ddlnq[i].set_qtot(ddlnqtot);
    }
}