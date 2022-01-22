#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<filesystem>
#include "input.h"


std::string readFile(std::string &input)
{
    std::ifstream inputfile;
    std::string inputstring;
    std::string line;

    inputfile.open(input);
    if (inputfile.fail())
    {
        std::cout << "error with " << input << " file" << std::endl;
        exit (1);
    }
    while (getline(inputfile, line))
    {
        inputstring += line;         
    }
    return inputstring; 
}

input_c get_config_arguments(std::string &input_contents)
{
    input_c input_object; 
    std::string input_string; 
    std::stringstream inputstream(input_contents);
    while (inputstream >> input_string)
    {
        if (input_string == "components")
        {
            std::string components;
            inputstream >> components; 
            input_object.set_components(stoi(components));
        }

        if (input_string == "amf")
        {
            std::string amf_start, amf_end, amf_step;
            inputstream >> amf_start >> amf_end >> amf_step;
            input_object.set_amf(stod(amf_start), stod(amf_end), stoi(amf_step));
        }

        if (input_string == "bxv")
        {
            std::string bxv_start, bxv_end, bxv_step;
            inputstream >> bxv_start >> bxv_end >> bxv_step;
            input_object.set_bxv(stod(bxv_start), stod(bxv_end), stoi(bxv_step));
        }

        if (input_string == "temp")
        {
            std::string temp_low, temp_high, temp_step; 
            inputstream >> temp_low >> temp_high >> temp_step;
            input_object.set_temperature(stod(temp_low), stod(temp_high), stoi(temp_step));
        }

        if (input_string == "pressure")
        {
            std::string pressure;
            inputstream >> pressure; 
            input_object.set_pressure(stod(pressure));
        }

        if (input_string == "monomer_amounts")
        {
            std::string monomer_amounts;
            inputstream >> monomer_amounts;
            input_object.set_monomer_amounts(stod(monomer_amounts));
        }

        if (input_string == "contributions")
        {
            input_object.set_contributions();
        }

        if (input_string == "density")
        {
            bool compare_density; 
            compare_density = true;
            if (compare_density == true)
            {
                std::string ref_density_temp, ref_density, ref_density_weight; 
                inputstream >> ref_density_temp >> ref_density >> ref_density_weight; 
                input_object.set_compare_density(stod(ref_density_temp), stod(ref_density), ref_density_weight); 
            }
        }

        if (input_string == "phase_transition")
        {
            bool compare_phase_transition = true;
            bool compare = true; 
            if (compare_phase_transition)
            {
                std::string phase_transition_temp, phase_transition_weight;
                inputstream >> phase_transition_temp >> phase_transition_weight;
                input_object.set_phase_transition(stod(phase_transition_temp), phase_transition_weight);
            }
        }

       if (input_string == "isobar")
        {
            bool compare_isobar = true;
            std::string exptfilename; 
            if (compare_isobar)
            {                
                std::string isobar_weight;
                inputstream >> exptfilename >> isobar_weight; 
                input_object.set_isobar(isobar_weight);
            }
            std::ifstream exptfile;
            std::string exptvolinput; 
            std::string expttemp;
            std::string exptvol;
            std::vector<double>exptvol_vec; 
            std::vector<double>expttemp_vec; 
            if (std::filesystem::exists(exptfilename))
            {
                exptfile.open(exptfilename);
                while (getline(exptfile, exptvolinput))
                {
                    std::stringstream exptstream(exptvolinput);
                    while (exptstream >> expttemp >> exptvol)
                    {
                        double vol_L = stod(exptvol);
                        exptvol_vec.push_back(vol_L);
                        expttemp_vec.push_back(stod(expttemp));
                        input_object.set_ref_isobar_volume(exptvol_vec);
                        input_object.set_ref_isobar_temp(expttemp_vec);
                    }
                }
            }
            else
            {
                std::cout << "expt_vol.inp not found" << std::endl;
                exit(1);
            }
        }        
    }
    return input_object; 
}


