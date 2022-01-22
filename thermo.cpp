#include "qce.h"
#include "partfunctions.h"
#include "auxiliary2.h"
#include "atomic_data.h"
#include "input.h"
#include<fstream>
#include<iomanip>

std::vector<double>calculate_internal_energy(int &temp_num, std::vector<double>&temp, std::vector<double>&dlnq, std::vector<double>&u)
{
    for (int i = 0; i < temp_num; i++)
    {
        u[i] = kb * pow(temp[i],2.00) * dlnq[i];
    }
    return u; 
}

std::vector<double>calculate_helmholtz_energy(int &temp_num, std::vector<double>&temp, std::vector<double>&lnq, std::vector<double>&a)
{
    for (int i = 0; i < temp_num; i++)
    {
        a[i] = -kb * temp[i] * lnq[i];
    } 
    return a; 
}

std::vector<double>calculate_enthalpy(int &temp_num, std::vector<double>&temp, std::vector<double>&dlnq, std::vector<double>&vol, double &press, std::vector<double>&h)
{
    std::vector<double>u; 
    u.resize(temp_num);
    std::fill(u.begin(), u.end(), 0.00);
    u = calculate_internal_energy(temp_num, temp, dlnq, u);
    
    for (int i = 0; i < temp_num; i++)
    {
        h[i] = u[i] + press * vol[i];
    }
    return h; 
}

std::vector<double>calculate_entropy(int &temp_num, std::vector<double>&temp, std::vector<double>&lnq, std::vector<double>&dlnq, std::vector<double>&s)
{
    std::vector<double>u;
    u.resize(temp_num);
    std::fill(u.begin(), u.end(), 0.00);
    std::vector<double>a;
    a.resize(temp_num);
    std::fill(a.begin(), a.end(), 0.00);
    u = calculate_internal_energy(temp_num, temp, dlnq, u);
    a = calculate_helmholtz_energy(temp_num, temp, lnq, a);

    for (int i = 0; i < temp_num; i++)
    {
        s[i] = (u[i] - a[i]) / temp[i];
    }
    
    return s; 
}

std::vector<double>calculate_cv(int &temp_num, std::vector<double>&temp, std::vector<double>&dlnq, std::vector<double>&ddlnq, std::vector<double>&cv)
{
    for (int i = 0; i < temp_num; i++)
    {
        cv[i] = 2.0 * kb * temp[i] * dlnq[i] + kb * pow(temp[i],2.00) * ddlnq[i];
    }
    return cv; 
}

std::vector<double>calculate_cp(int &temp_num, std::vector<double>&temp, std::vector<double>&dlnq, std::vector<double>&ddlnq, std::vector<double>&dvol, double &press, std::vector<double>&cp)
{
    std::vector<double>cv; 
    cv.resize(temp_num);
    std::fill(cv.begin(), cv.end(), 0.00);
    cv = calculate_cv(temp_num, temp, dlnq, ddlnq, cv);
    if (temp_num > 1)
    {
        for (int i = 0; i < temp_num; i++)
        {
            cp[i] = cv[i] + dvol[i] * press; 
        }
    }
    else
    {
        std::fill(cp.begin(), cp.end(), 0.00);
    }
    return cp; 
}

std::vector<double>calculate_gibbs_enthalpy(int &temp_num, std::vector<double>&temp, std::vector<double>&lnq, std::vector<double>&vol, double &press, std::vector<double>&g)
{
    std::vector<double>a;
    a.resize(temp_num);
    std::fill(a.begin(), a.end(), 0.00);
    a = calculate_helmholtz_energy(temp_num, temp, lnq, a);
    for (int i = 0; i < temp_num; i++)
    {
        g[i] = a[i] + press * vol[i];
    }
    return g; 
}


void write_contributions(int &temp_num, std::vector<bool>&converged, std::vector<double>&temp, std::vector<double>&a, std::vector<double>&vib, std::vector<double>&rot, std::vector<double>&trans, std::vector<double>&elec, std::vector<double>&mf, std::vector<double>&indi, std::string header, double factor, std::string filename)
{
    std::ofstream contrib;
    contrib.open(filename);
    contrib << std::left << header << std::endl; 
    contrib << std::left << "#" << std::setw(20) << "T/K" << std::setw(20) << "tot" << std::setw(20) << "vib" << std::setw(20) << "rot" << std::setw(20) << "trans" << std::setw(20) << "elec" << std::setw(20) << "mf" << std::setw(20) << "indi" << std::endl; 

    for (int i = 0; i < temp_num; i++)
    {
        if (converged[i])
        {
            contrib << std::left << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << temp[i]; 
            contrib << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << factor * a[i];
            contrib << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << factor * vib[i];
            contrib << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << factor * rot[i];
            contrib << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << factor * trans[i];
            contrib << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << factor * elec[i];
            contrib << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << factor * mf[i];
            contrib << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << factor * indi[i];
            contrib << std::endl; 
        }
    }
    contrib.close(); 
}

void write_volume(int &temp_num, std::vector<bool>&converged, std::vector<double>&temp, std::vector<double>&vol, double &vexcl, std::vector<double>&alpha, std::vector<int>&solution)
{
    std::ofstream volume;
    volume.open("volume.dat");
    volume << std::left << "#" << std::setw(20) << "T/K" << std::setw(20) << "V/dm3" << std::setw(20) << "V(excl)dm3" << std::setw(20) << "alpha/K^-1" << std::setw(20) << "status" << std::endl;  
    for (int i = 0; i < temp_num; i++)
    {
        if (converged[i])
        {
            volume << std::left << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << temp[i]; 
            volume << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << 1.00e3 * vol[i]; 
            volume << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << 1.00e3 * vexcl; 
            volume << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << alpha[i]; 
            volume << std::setw(20) << solution[i]; 
            volume << std::endl; 
        }
    }
    volume.close(); 
}

void write_thermo0(int &temp_num, std::vector<bool>&converged, std::vector<double>&temp, std::vector<double>&a, std::vector<double>&g)
{
    std::ofstream thermo0;
    thermo0.open("thermo0.dat");
    thermo0 << std::left << "#" << std::setw(20) << "T/K" << std::setw(20) << "A/kJ" << std::setw(20) << "G/kJ" << std::endl; 
    for (int i = 0; i < temp_num; i++)
    {
        if (converged[i])
        {
            thermo0 << std::left << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << temp[i]; 
            thermo0 << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << 1.00e-3 * a[i];
            thermo0 << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << 1.00e-3 * g[i];
            thermo0 << std::endl; 
        }
    }
    thermo0.close(); 
}

void write_thermo1(int &temp_num, std::vector<bool>&converged, std::vector<double>&temp, std::vector<double>&u, std::vector<double>&h, std::vector<double>&s)
{
    std::ofstream thermo1;
    thermo1.open("thermo1.dat");
    thermo1 << std::left << "#" << std::setw(20) << "T/K" << std::setw(20) << "U/kJ" << std::setw(20) << "H/kJ" << std::setw(20) << "S/J*K^-1" << std::endl; 
    for (int i = 0; i < temp_num; i++)
    {
        if (converged[i])
        {
            thermo1 << std::left << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << temp[i]; 
            thermo1 << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << 1.00e-3 * u[i];
            thermo1 << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << 1.00e-3 * h[i];
            thermo1 << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << s[i];
            thermo1 << std::endl; 
        }
    }
    thermo1.close();
}

void write_thermo2(int &temp_num, std::vector<bool>&converged, std::vector<double>&temp, std::vector<double>&cp, std::vector<double>&cv)
{
    std::ofstream thermo2;
    thermo2.open("thermo2.dat");
    thermo2 << std::left << "#" << std::setw(20) << "T/K" << std::setw(20) << "Cv/J*K^-1" << std::setw(20) << "Cp/J*K^-1" << std::endl;
    for (int i = 0; i < temp_num; i++)
    {
        if (converged[i])
        {
            thermo2 << std::left << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << temp[i]; 
            thermo2 << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << cv[i];
            thermo2 << std::setw(20) << std::uppercase << std::scientific << std::setprecision(6) << cp[i];
            thermo2 << std::endl; 
        }
    }
    thermo2.close(); 
}

std::vector<double>calculate_expansion_coefficient(int& temp_num, std::vector<double>&vol, std::vector<double>&dvol, std::vector<double>&alpha)
{
    if (temp_num > 1)
    {
        for (int i = 0; i < temp_num; i++)
        {
            alpha[i] = 1.0 / vol[i] * dvol[i];
        }
    }
    else
    {
        std::fill(alpha.begin(), alpha.end(), 0.00);
    }
    return alpha; 
}

std::vector<double>add_lnq_indi(int& temp_num, int &clusterset_size, std::vector<std::vector<double>>&populations, std::vector<double>&lnq_sys)
{
    for (int i = 0; i < temp_num; i++)
    {
        for (int j = 0; j < clusterset_size; j++)
        {
            lnq_sys[i] = lnq_sys[i] - ln_factorial(populations[i][j]);
        }
    }
    return lnq_sys; 
}

std::vector<double>calculate_lnq_sys(int& temp_num, int &clusterset_size, std::vector<std::vector<double>>&populations, std::vector<std::vector<pf_c>>&lnq_clust, std::string lnq_type)
{
    std::vector<double>lnq_sys; 
    lnq_sys.resize(temp_num);
    std::fill(lnq_sys.begin(), lnq_sys.end(), 0.00);
  
    if (lnq_type == "qtrans")
    {
        for (int i = 0; i < temp_num; i++)
        {
            for (int j = 0; j < clusterset_size; j++)
            {
                lnq_sys[i] = lnq_sys[i] + populations[i][j] * lnq_clust[i][j].get_qtrans();
            }
        }
    }

    if (lnq_type == "qvib")
    {
        for (int i = 0; i < temp_num; i++)
        {
            for (int j = 0; j < clusterset_size; j++)
            {
                lnq_sys[i] = lnq_sys[i] + populations[i][j] * lnq_clust[i][j].get_qvib();
            }
        }
    }

    if (lnq_type == "qrot")
    {
        for (int i = 0; i < temp_num; i++)
        {
            for (int j = 0; j < clusterset_size; j++)
            {
                lnq_sys[i] = lnq_sys[i] + populations[i][j] * lnq_clust[i][j].get_qrot();
            }
        }
    }

    if (lnq_type == "qelec")
    {
        for (int i = 0; i < temp_num; i++)
        {
            for (int j = 0; j < clusterset_size; j++)
            {
                lnq_sys[i] = lnq_sys[i] + populations[i][j] * lnq_clust[i][j].get_qelec();
            }
        }
    }

    if (lnq_type == "qint")
    {
        std::vector<std::vector<double>>lnq_int; 
        for (int i = 0; i < temp_num; i++)
        {
            for (int j = 0; j < clusterset_size; j++)
            {
                lnq_sys[i] = lnq_sys[i] + populations[i][j] * lnq_clust[i][j].get_qint();
            }
        }
    }

    if (lnq_type == "qtot")
    {
        for (int i = 0; i < temp_num; i++)
        {
            for (int j = 0; j < clusterset_size; j++)
            {
                lnq_sys[i] = lnq_sys[i] + populations[i][j] * lnq_clust[i][j].get_qtot();
            }
        }
    }
    return lnq_sys; 
}
void calculate_cluster_entropy(std::string filename, isobar_c &ib, std::vector<std::vector<pf_c>>&lnq_clust, std::vector<std::vector<pf_c>>&dlnq_clust, std::vector<cluster_c>&clusterset, std::string pftype)
{
    std::ofstream cluster_entropy;
    cluster_entropy.open(filename);
    cluster_entropy << std::left << "#" << std::setw(20) << "T/K(J*K^-1)" << std::setw(20);

    for (int i = 0; i < clusterset.size(); i++)
    {
        cluster_entropy << std::setw(20) << clusterset[i].get_label();
        if (i == clusterset.size() - 1)
        {
            cluster_entropy << std::endl; 
        }
    }

    for (int i = 0; i < ib.get_temperature().size(); i++)
    {
        cluster_entropy << std::left << std::setw(20) << std::setprecision(6) << std::scientific << std::uppercase << ib.get_temperature_element(i);
        for (int j = 0; j < lnq_clust[0].size(); j++)
        {
            if (pftype == "qtot")
            {
                double pop_qtot;
                double pop_factorial;
                double pop_qtot_dlnq; 
                double helmholtz;
                double internal_energy; 
                double entropy; 
                pop_qtot = ib.get_populations(i)[j] * lnq_clust[i][j].get_qtot(); 
                pop_qtot_dlnq = ib.get_populations(i)[j] * dlnq_clust[i][j].get_qtot(); 
                pop_factorial = ln_factorial(ib.get_populations(i)[j]);
                helmholtz = (-kb * ib.get_temperature_element(i) * (pop_qtot - pop_factorial));
                internal_energy = (kb * pow(ib.get_temperature_element(i),2) * pop_qtot_dlnq);
                entropy = (internal_energy - helmholtz) / ib.get_temperature_element(i);
                cluster_entropy << std::setw(20) << std::uppercase << std::setprecision(6) << std::scientific << entropy; 
                if (j == ib.get_lnq_vector()[0].size() - 1)
                {
                    cluster_entropy << std::endl; 
                }
            }

            if (pftype == "qtrans")
            {
                double pop_qtot;
                double pop_factorial;
                double pop_qtot_dlnq; 
                double helmholtz;
                double internal_energy; 
                double entropy; 
                pop_qtot = ib.get_populations(i)[j] * lnq_clust[i][j].get_qtrans(); 
                pop_qtot_dlnq = ib.get_populations(i)[j] * dlnq_clust[i][j].get_qtrans(); 
                pop_factorial = ln_factorial(ib.get_populations(i)[j]);
                helmholtz = (-kb * ib.get_temperature_element(i) * (pop_qtot - pop_factorial));
                internal_energy = (kb * pow(ib.get_temperature_element(i),2) * pop_qtot_dlnq);
                entropy = (internal_energy - helmholtz) / ib.get_temperature_element(i);
                cluster_entropy << std::setw(20) << std::uppercase << std::setprecision(6) << std::scientific << entropy; 
                if (j == ib.get_lnq_vector()[0].size() - 1)
                {
                    cluster_entropy << std::endl; 
                }
            }

            if (pftype == "qvib")
            {
                double pop_qtot;
                double pop_factorial;
                double pop_qtot_dlnq; 
                double helmholtz;
                double internal_energy; 
                double entropy; 
                pop_qtot = ib.get_populations(i)[j] * lnq_clust[i][j].get_qvib(); 
                pop_qtot_dlnq = ib.get_populations(i)[j] * dlnq_clust[i][j].get_qvib(); 
                pop_factorial = ln_factorial(ib.get_populations(i)[j]);
                helmholtz = (-kb * ib.get_temperature_element(i) * (pop_qtot - pop_factorial));
                internal_energy = (kb * pow(ib.get_temperature_element(i),2) * pop_qtot_dlnq);
                entropy = (internal_energy - helmholtz) / ib.get_temperature_element(i);
                cluster_entropy << std::setw(20) << std::uppercase << std::setprecision(6) << std::scientific << entropy; 
                if (j == ib.get_lnq_vector()[0].size() - 1)
                {
                    cluster_entropy << std::endl; 
                }
            }

            if (pftype == "qrot")
            {
                double pop_qtot;
                double pop_factorial;
                double pop_qtot_dlnq; 
                double helmholtz;
                double internal_energy; 
                double entropy; 
                pop_qtot = ib.get_populations(i)[j] * lnq_clust[i][j].get_qrot(); 
                pop_qtot_dlnq = ib.get_populations(i)[j] * dlnq_clust[i][j].get_qrot(); 
                pop_factorial = ln_factorial(ib.get_populations(i)[j]);
                helmholtz = (-kb * ib.get_temperature_element(i) * (pop_qtot - pop_factorial));
                internal_energy = (kb * pow(ib.get_temperature_element(i),2) * pop_qtot_dlnq);
                entropy = (internal_energy - helmholtz) / ib.get_temperature_element(i);
                cluster_entropy << std::setw(20) << std::uppercase << std::setprecision(6) << std::scientific << entropy; 
                if (j == ib.get_lnq_vector()[0].size() - 1)
                {
                    cluster_entropy << std::endl; 
                }
            }

            if (pftype == "qint")
            {
                double pop_qtot;
                double pop_factorial;
                double pop_qtot_dlnq; 
                double helmholtz;
                double internal_energy; 
                double entropy; 
                pop_qtot = ib.get_populations(i)[j] * lnq_clust[i][j].get_qint(); 
                pop_qtot_dlnq = ib.get_populations(i)[j] * dlnq_clust[i][j].get_qint(); 
                pop_factorial = ln_factorial(ib.get_populations(i)[j]);
                helmholtz = (-kb * ib.get_temperature_element(i) * (pop_qtot - pop_factorial));
                internal_energy = (kb * pow(ib.get_temperature_element(i),2) * pop_qtot_dlnq);
                entropy = (internal_energy - helmholtz) / ib.get_temperature_element(i);
                cluster_entropy << std::setw(20) << std::uppercase << std::setprecision(6) << std::scientific << entropy; 
                if (j == ib.get_lnq_vector()[0].size() - 1)
                {
                    cluster_entropy << std::endl; 
                }
            }
        }
    }
}
void calculate_thermo(int &temp_num, int &clusterset_size, std::vector<std::vector<double>>&populations, double &press, std::vector<double>&temp, std::vector<double>&vol, double &vexcl, std::vector<std::vector<pf_c>>&lnq_clust, std::vector<std::vector<pf_c>>&dlnq_clust, std::vector<std::vector<pf_c>>&ddlnq_clust, std::vector<double>&dvol,std::vector<int>&solution, std::vector<bool>&converged, input_c &input)
{
    std::vector<double>lnq, dlnq, ddlnq;
    std::vector<double>lnqvib, dlnqvib, ddlnqvib; 
    std::vector<double>lnqrot, dlnqrot, ddlnqrot;
    std::vector<double>lnqelec, dlnqelec, ddlnqelec; 
    std::vector<double>lnqint, dlnqint, ddlnqint;
    std::vector<double>lnqtrans, dlnqtrans, ddlnqtrans; 
    std::vector<double>lnqindi, dlnqindi, ddlnqindi; 
    std::vector<double>a, g, u, s, h, alpha, cv, cp; 
    std::vector<double> vib, rot, trans, elec, mf, indi; 

    vib.resize(temp_num);    
    rot.resize(temp_num);
    trans.resize(temp_num);
    elec.resize(temp_num);
    mf.resize(temp_num);
    indi.resize(temp_num);
    a.resize(temp_num);
    g.resize(temp_num);
    u.resize(temp_num);
    s.resize(temp_num);
    h.resize(temp_num);
    alpha.resize(temp_num);
    cv.resize(temp_num);
    cp.resize(temp_num);

    std::fill(vib.begin(), vib.end(), 0.00);
    std::fill(rot.begin(), rot.end(), 0.00);
    std::fill(trans.begin(), trans.end(), 0.00);
    std::fill(elec.begin(), elec.end(), 0.00);
    std::fill(mf.begin(), mf.end(), 0.00);
    std::fill(indi.begin(), indi.end(), 0.00);
    std::fill(a.begin(), a.end(), 0.00);
    std::fill(g.begin(), g.end(), 0.00);
    std::fill(u.begin(), u.end(), 0.00);
    std::fill(s.begin(), s.end(), 0.00);
    std::fill(h.begin(), h.end(), 0.00);
    std::fill(alpha.begin(), alpha.end(), 0.00);
    std::fill(cv.begin(), cv.end(), 0.00);
    std::fill(cp.begin(), cp.end(), 0.00);

    // Calculate system partition function. 
    lnq = calculate_lnq_sys(temp_num, clusterset_size, populations, lnq_clust, "qtot");
    lnq = add_lnq_indi(temp_num, clusterset_size, populations, lnq);
    dlnq = calculate_lnq_sys(temp_num, clusterset_size, populations, dlnq_clust, "qtot");
    ddlnq = calculate_lnq_sys(temp_num, clusterset_size, populations, ddlnq_clust, "qtot");
 
    // Calculate contributions to the system partition function
    lnqtrans = calculate_lnq_sys(temp_num, clusterset_size, populations, lnq_clust, "qtrans");
    dlnqtrans = calculate_lnq_sys(temp_num, clusterset_size, populations, dlnq_clust, "qtrans");
    ddlnqtrans = calculate_lnq_sys(temp_num, clusterset_size, populations, ddlnq_clust, "qtrans");
    lnqvib = calculate_lnq_sys(temp_num, clusterset_size, populations, lnq_clust, "qvib");
    dlnqvib = calculate_lnq_sys(temp_num, clusterset_size, populations, dlnq_clust, "qvib");
    ddlnqvib = calculate_lnq_sys(temp_num, clusterset_size, populations, ddlnq_clust, "qvib");
    lnqrot = calculate_lnq_sys(temp_num, clusterset_size, populations, lnq_clust, "qrot");
    dlnqrot = calculate_lnq_sys(temp_num, clusterset_size, populations, dlnq_clust, "qrot");
    ddlnqrot = calculate_lnq_sys(temp_num, clusterset_size, populations, ddlnq_clust, "qrot");
    lnqelec = calculate_lnq_sys(temp_num, clusterset_size, populations, lnq_clust, "qelec");
    dlnqelec = calculate_lnq_sys(temp_num, clusterset_size, populations, dlnq_clust, "qelec");
    ddlnqelec = calculate_lnq_sys(temp_num, clusterset_size, populations, ddlnq_clust, "qelec");
    lnqint = calculate_lnq_sys(temp_num, clusterset_size, populations, lnq_clust, "qint");
    dlnqint = calculate_lnq_sys(temp_num, clusterset_size, populations, dlnq_clust, "qint");
    ddlnqint = calculate_lnq_sys(temp_num, clusterset_size, populations, ddlnq_clust, "qint");
    
    lnqindi.resize(temp_num);
    dlnqindi.resize(temp_num);
    ddlnqindi.resize(temp_num);
    std::fill(lnqindi.begin(), lnqindi.end(), 0.00);
    std::fill(dlnqindi.begin(), dlnqindi.end(), 0.00);
    std::fill(ddlnqindi.begin(), ddlnqindi.end(), 0.00);
    lnqindi = add_lnq_indi(temp_num, clusterset_size, populations, lnqindi);
   
    // write volume, exclusion volume and status code. 
    alpha = calculate_expansion_coefficient(temp_num, vol, dvol, alpha);
    write_volume(temp_num, converged, temp, vol, vexcl, alpha, solution);

    // Calculate and write thermodynamic functions
    a = calculate_helmholtz_energy(temp_num, temp, lnq, a);
    g = calculate_gibbs_enthalpy(temp_num, temp, lnq, vol, press, g);
    u = calculate_internal_energy(temp_num, temp, dlnq, u);
    h = calculate_enthalpy(temp_num, temp, dlnq, vol, press, h);
    s = calculate_entropy(temp_num, temp, lnq, dlnq, s);
    cv = calculate_cv(temp_num, temp, dlnq, ddlnq, cv);
    cp = calculate_cp(temp_num, temp, dlnq, ddlnq, dvol, press, cp);
    write_thermo0(temp_num, converged, temp, a, g);
    write_thermo1(temp_num, converged, temp, u, h, s);
    write_thermo2(temp_num, converged, temp, cp, cv);

    // Calculate and write contributions to the Helmholtz free energy. 
    if (input.get_contributions() or input.get_helmholtz_contrib())
    {
        vib = calculate_helmholtz_energy(temp_num, temp, lnqvib, vib);
        rot = calculate_helmholtz_energy(temp_num, temp, lnqrot, rot);
        trans = calculate_helmholtz_energy(temp_num, temp, lnqtrans, trans);
        elec = calculate_helmholtz_energy(temp_num, temp, lnqelec, elec);
        mf = calculate_helmholtz_energy(temp_num, temp, lnqint, mf);
        indi = calculate_helmholtz_energy(temp_num, temp, lnqindi, indi);
        write_contributions(temp_num, converged, temp, a, vib, rot, trans, elec, mf, indi, "A/kJ", 1.00e-3, "helmholtz_contrib.dat");
    }

    // Calculate and write contributions to the internal energy. 
    if (input.get_contributions() or input.get_internal_contrib())
    {
        vib = calculate_internal_energy(temp_num, temp, dlnqvib, vib);
        rot = calculate_internal_energy(temp_num, temp, dlnqrot, rot);
        trans = calculate_internal_energy(temp_num, temp, dlnqtrans, trans);
        elec = calculate_internal_energy(temp_num, temp, dlnqelec, elec);
        mf = calculate_internal_energy(temp_num, temp, dlnqint, mf);
        indi = calculate_internal_energy(temp_num, temp, dlnqindi, indi);
        write_contributions(temp_num, converged, temp, u, vib, rot, trans, elec, mf, indi, "U/kJ", 1.0e-3, "internal_contrib.dat");
    }
    
    // Calculate and write contributions to the entropy
    if (input.get_contributions() or input.get_entropy_contrib())
    {
        vib = calculate_entropy(temp_num, temp, lnqvib, dlnqvib, vib);
        rot = calculate_entropy(temp_num, temp, lnqrot, dlnqrot, rot);
        trans = calculate_entropy(temp_num, temp, lnqtrans, dlnqtrans, trans);
        
        std::fill(elec.begin(), elec.end(), 0.00);
        std::fill(mf.begin(), mf.end(), 0.00);

        indi = calculate_entropy(temp_num, temp, lnqindi, dlnqindi, indi);
        write_contributions(temp_num, converged, temp, s, vib, rot, trans, elec, mf, indi, "S/J*K^-1", 1.0, "entropy_contrib.dat");
    }

    if (input.get_contributions() or input.get_cv_contrib())
    {
        vib = calculate_cv(temp_num, temp, dlnqvib, ddlnqvib, vib);
        rot = calculate_cv(temp_num, temp, dlnqrot, ddlnqrot, rot);
        trans = calculate_cv(temp_num, temp, dlnqtrans, ddlnqtrans, trans);

        // electronic contributions to cv are zero. 
        std::fill(elec.begin(), elec.end(), 0.00);

        //mean field contributions to cv are zero. 
        std::fill(mf.begin(), mf.end(), 0.00);

        indi = calculate_cv(temp_num, temp, dlnqindi, ddlnqindi, indi);
        write_contributions(temp_num, converged, temp, cv, vib, rot, trans, elec, mf, indi, "Cv/J*K^-1", 1.0, "cv_contrib.dat");
    }
}
