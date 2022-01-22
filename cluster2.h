#ifndef CLUSTER_H
#define CLUSTER_H
#include<vector>
#include<map>
#include<string>

class input_c; 

class cluster_c
{
    int sigma_c;
    std::vector<double>frequencies_c;
    double inertia_c[3]; 
    int composition_c; 
    std::vector<std::vector<double> >xyz_c;
    std::vector<double>atomicmass_c; 
    std::vector<std::string>element_c;
    double energy_c; 
    double totalmass_c; 
    double volume_c;
    double fscale_c; 
    double anharmonicity_c; 
    bool linear_c;
    bool atom_c;
    bool monomer_c; 
    std::string label_c; 
    std::string coordpath_c;
    std::string freqpath_c;

    public:
    cluster_c()
    {
        totalmass_c = 0.0; 
        energy_c = 0.0; 
        sigma_c = 1; 
        fscale_c = 1.0;
        anharmonicity_c = 0.0;
        linear_c = false;
        atom_c = false;
        monomer_c = false; 
        volume_c = 0.0; 
    }
    
    void set_composition(const int &composition);
    void set_element(const std::vector<std::string> &element_vec);
    void set_volume(const double volume);
    void set_energy(const double &energy);
    void set_sigma(const int &sigma);
    void set_label(const std::string &label);
    void set_coordpath(const std::string &coordpath);
    void set_freqpath(const std::string &freqpath);
    void set_inertia(const double &inertia, const int &index);
    void set_frequencies(const std::vector<double>&frequencies);
    void set_coordinates(const std::vector<std::vector<double> >&coordinates);
    void set_molmass(const double &molmass);
    void set_linear(const bool &linear);
    void set_atom(const bool &atom);
    void set_monomer(); 
    void set_atomicmass(const std::vector<double> &atomic_mass);

    // getter functions. 
    int get_composition();
    double get_volume(); 
    double get_energy(); 
    double get_anharmonicity();
    int get_sigma(); 
    bool get_atom();
    bool get_linear();
    double get_momofinertia(int i);
    std::string get_label(); 
    std::string get_coordpath(); 
    std::string get_freqpath(); 
    double get_frequency_element(int i);
    std::vector<std::vector<double> >get_coordinates();
    std::vector<double>get_frequencies();
    std::vector<std::string>get_elements(); 
    std::vector<double>get_atomicmass(); 
    double get_molmass(); 
    double get_fscale();    
};

int number_of_sections(std::string &cluster_contents);
std::vector<cluster_c>get_cluster_arguments(std::string &cluster_contents);
std::vector<cluster_c>process_cluster(std::string clusterfile, input_c &input_object);
std::vector<cluster_c>cluster_range(std::vector<cluster_c> &v, int m, int n);
std::vector<cluster_c>cluster_reverse_range(std::vector<cluster_c> &v, int m, int n);
std::map<int, std::vector<cluster_c>>clustersetlist(std::vector<cluster_c>&clusterset, int &combination_no);
void extract_coordinates(std::string coordpath, cluster_c &cluster_object);
void extract_frequencies(std::string freqpath, cluster_c &cluster_object);
void extract_atomicmass(cluster_c &cluster_object);
void calculate_molmass(cluster_c &cluster_object);
void calculate_momofinertia(cluster_c &cluster_object);
void check_clusterset(std::vector<cluster_c>&clusterset, input_c &input_object); 

#endif