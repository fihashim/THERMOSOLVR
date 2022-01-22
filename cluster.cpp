#include<iostream>
#include<sstream>
#include<vector>
#include<regex>
#include<fstream>
#include<algorithm>
#include<set>
#include<map>
#include<numeric>
#include "cluster2.h"
#include "config2.h"
#include "mm_malloc.h"
#include "atomic_data.h"



template<class T>
void printVectofVects(std::vector<std::vector<T>> const &mat) {
	for (std::vector<T> row: mat) {
		for (T val: row) {
			std::cout << val << " ";
		}
		std::cout << '\n';
	}
}

void cluster_c::set_element(const std::vector<std::string> &element_vec)
{
    element_c = element_vec;
}

void cluster_c::set_label(const std::string &label)
{
    label_c = label; 
}

void cluster_c::set_volume(const double volume)
{
    volume_c = volume; 
}

void cluster_c::set_composition(const int &composition)
{
    composition_c = composition; 
}

void cluster_c::set_sigma(const int &sigma)
{
    sigma_c = sigma; 
}

void cluster_c::set_energy(const double &energy)
{
    energy_c = energy; 
}

void cluster_c::set_coordpath(const std::string &coordpath)
{
    coordpath_c = coordpath; 
}

void cluster_c::set_freqpath(const std::string &freqpath)
{
    freqpath_c = freqpath; 
}

void cluster_c::set_coordinates(const std::vector<std::vector<double>> &coordinates)
{
    xyz_c = coordinates;
}

void cluster_c::set_frequencies(const std::vector<double>&frequencies)
{
    frequencies_c = frequencies; 
}

void cluster_c::set_molmass(const double &molmass)
{
    totalmass_c = molmass; 
}

void cluster_c::set_atomicmass(const std::vector<double> &atomic_mass)
{
    atomicmass_c = atomic_mass; 
}

void cluster_c::set_linear(const bool &linear)
{
    linear_c = linear; 
}

void cluster_c::set_atom(const bool &atom)
{
    atom_c = atom; 
}

void cluster_c::set_inertia(const double &inertia, const int &index)
{
    inertia_c[index] = inertia; 
}
int cluster_c::get_composition()
{
    return composition_c;
}

double cluster_c::get_volume()
{
    return volume_c; 
}

std::string cluster_c::get_coordpath()
{
    return coordpath_c; 
}

double cluster_c::get_energy()
{
    return energy_c;
}

double cluster_c::get_fscale()
{
    return fscale_c; 
}
int cluster_c::get_sigma()
{
    return sigma_c;
}

std::string cluster_c::get_label()
{
    return label_c;
}

std::string cluster_c::get_freqpath()
{
    return freqpath_c; 
}

std::vector<std::vector<double>>cluster_c::get_coordinates()
{
    return xyz_c; 
}

std::vector<double>cluster_c::get_frequencies()
{
    return frequencies_c;
}

std::vector<std::string>cluster_c::get_elements()
{
    return element_c;
}

double cluster_c::get_molmass()
{
    return totalmass_c;
}

std::vector<double>cluster_c::get_atomicmass()
{
    return atomicmass_c; 
}

double cluster_c::get_frequency_element(int i)
{
    return frequencies_c[i];
}

double cluster_c::get_anharmonicity()
{
    return anharmonicity_c; 
}

bool cluster_c::get_atom()
{
    return atom_c; 
}

bool cluster_c::get_linear()
{
    return linear_c; 
}

double cluster_c::get_momofinertia(int i)
{
    return inertia_c[i]; 
}

std::vector<cluster_c>get_cluster_arguments(std::string &cluster_contents)
{
    std::vector<cluster_c>clusterset; 
    std::string cluster_string; 
    std::stringstream clusterstream(cluster_contents);
    cluster_c cluster_object; 
    double monomer_vol; 
    int monomer_count; 

    while (clusterstream >> cluster_string)
    {   
        unsigned open_bracket = cluster_string.find('[') + 1;
        unsigned closed_bracket = cluster_string.find(']');
        
        if (open_bracket)
        {
            std::string label; 
            label = cluster_string.substr(open_bracket, closed_bracket - open_bracket);
            cluster_object.set_label(label);
        }

        if (cluster_string == "composition")
        {
            std::string composition; 
            clusterstream >> composition; 
            cluster_object.set_composition(stoi(composition));
        }

        if (cluster_string == "sigma")
        {
            std::string sigma;
            clusterstream >> sigma;
            cluster_object.set_sigma(stoi(sigma));
        }

        if (cluster_string == "volume")
        {
            std::string volume;
            clusterstream >> volume;
            monomer_vol = stod(volume);
        }

        if (cluster_string == "energy")
        {
            std::string energy;
            clusterstream >> energy;
            cluster_object.set_energy(stod(energy));
        }

        if (cluster_string == "coordinates")
        {
            std::string coordpath;
            clusterstream >> coordpath; 
            cluster_object.set_coordpath(coordpath);
        }
        
        if (cluster_string == "frequencies")
        {
            std::string freqpath;
            clusterstream >> freqpath;
            cluster_object.set_freqpath(freqpath);
            clusterset.push_back(cluster_object);
        }
    }

    clusterset[0].set_volume(monomer_vol);
    
    return clusterset;
};



void extract_coordinates(std::string coordpath, cluster_c &cluster_object)
{
    std::string coord_contents;
    std::string coord_string; 
    std::regex coord_pattern("-?[[:digit:]]{1,2}.[[:digit:]]");
    std::regex element_pattern("[[:alpha:]]{1,2}");
    coord_contents = readFile(coordpath);

    std::stringstream coordstream(coord_contents);
    std::vector<std::string>element_vec;
    std::vector<double>xyz_vec; 
    std::vector<std::vector<double>>xyz_matrix; 
    int coord_count = 0; 
    while (coordstream >> coord_string)
    {
        // get xyz elements
        if (regex_search(coord_string, element_pattern))
        {
            std::string element; 
            element += coord_string[coord_string.length() - 1];
            element_vec.push_back(element);
            coord_string = regex_replace(coord_string, element_pattern, "");
        }
        
        // get xyz coordinates. 
        if (regex_search(coord_string, coord_pattern))
        {
            coord_count += 1; 
            if (coord_count < 4)
            {
                xyz_vec.push_back(stod(coord_string));
                if (xyz_vec.size() == 3)
                {
                    xyz_matrix.push_back(xyz_vec);
                    xyz_vec.clear(); 
                    coord_count = 0; 
                }
            }
        }
    }
    cluster_object.set_element(element_vec);
    cluster_object.set_coordinates(xyz_matrix);
}

void extract_frequencies(std::string freqpath, cluster_c &cluster_object)
{
    std::ifstream freqfile; 
    std::string freq_string; 
    std::vector<double>frequencies; 
    freqfile.open(freqpath);

    if (freqfile.fail())
    {
        std::cout << cluster_object.get_label() << " frequency file does not exist." << std::endl;
        exit(1);
    }

    while (getline(freqfile, freq_string))
    {
        std::istringstream freqstream(freq_string);
        std::string modes;
        std::string strfreq;
        std::string symmetry;
        std::string redmass;
        std::string ir; 
        while (freqstream >> modes >> strfreq >> symmetry >> redmass >> ir)
        {
            if (strfreq != "FREQ(CM**-1)")
            {
                frequencies.push_back(stod(strfreq));
            }
        }
    };
    freqfile.close(); 
    cluster_object.set_frequencies(frequencies);

}

void calculate_molmass(cluster_c &cluster_object)
{
    std::vector<std::string>elements = cluster_object.get_elements(); 
    double molmass; 
    for (auto iter = elements.begin(); iter != elements.end(); iter++)
    {
        for (auto mapiter = periodic_table.begin(); mapiter != periodic_table.end(); mapiter++)
        {
            if (*iter == mapiter ->first)
            {
                molmass = molmass + mapiter ->second; 
            }
        }
    }
    cluster_object.set_molmass(molmass);
}

void extract_atomicmass(cluster_c &cluster_object)
{
    std::vector<double>atomic_mass; 
    std::vector<std::string>elements = cluster_object.get_elements(); 
    //std::cout.precision(17);
    for (auto iter = elements.begin(); iter != elements.end(); iter++)
    {
        for (auto mapiter = periodic_table.begin(); mapiter != periodic_table.end(); mapiter++)
        {
            if (*iter == mapiter->first)
            {
                atomic_mass.push_back(mapiter->second);
            }
        }
    }
    cluster_object.set_atomicmass(atomic_mass);
}



int number_of_sections(std::string &cluster_contents)
{
    std::stringstream clusterstream(cluster_contents);
    std::string cluster_string; 
    int no_clusters = 0; 
    while (clusterstream >> cluster_string)
    {
        unsigned open_bracket = cluster_string.find('[') + 1;
        unsigned closed_bracket = cluster_string.find(']');
        
        if (open_bracket)
        {
            std::string label; 
            label = cluster_string.substr(open_bracket,closed_bracket - open_bracket);
            no_clusters += 1; 
        }
    }
    
    return no_clusters; 
}

std::vector<cluster_c>process_cluster(std::string clusterfile, input_c &input_object)
{
    std::string cluster_contents;
    std::vector<cluster_c>clusterset; 
    cluster_contents = readFile(clusterfile);
    int no_clusters; 
    no_clusters = number_of_sections(cluster_contents);
    if (no_clusters == 0)
    {
        std::cerr << "empty clusterset" << std::endl;
        exit (1);
    }

    clusterset = get_cluster_arguments(cluster_contents);

    // set volumes for clusters > 1 and extract elements / coordinates / frequencies.
    for (int i = 0; i < clusterset.size(); i++)
    {
        if (clusterset[i].get_volume() == 0.0)
        {
            clusterset[i].set_volume(clusterset[0].get_volume() * clusterset[i].get_composition());
        }
        extract_coordinates(clusterset[i].get_coordpath(), clusterset[i]);
        extract_frequencies(clusterset[i].get_freqpath(), clusterset[i]);
        extract_atomicmass(clusterset[i]);
        calculate_molmass(clusterset[i]);
        calculate_momofinertia(clusterset[i]);
    }
    check_clusterset(clusterset, input_object);
    return clusterset; 
}

typedef std::numeric_limits<double>dbl; 
void calculate_momofinertia(cluster_c &cluster_object)
{
    std::vector<double>atomicmass = cluster_object.get_atomicmass(); 
    std::vector<std::vector<double>>xyz_matrix = cluster_object.get_coordinates(); 
    std::vector<double>atom_coordinates; 
    for (int i = 0; i < atomicmass.size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            double atom_coord = atomicmass[i] * xyz_matrix[i][j];
            atom_coordinates.push_back(atom_coord);
        }
    }

    // add first 3 atomic_coordinates to com_temp
    std::vector<double>com_temp; 
    for (int i = 0; i < 3; i++)
    {
        com_temp.push_back(atom_coordinates[i]);
    }

    // calculate com for second row only. 
    double second_row; 
    for (int i = 0; i < atom_coordinates.size(); i++)
    {
        if (i == 3 or 4 or 5)
        {
            second_row = atom_coordinates[i] + atom_coordinates[i+3];
            com_temp.push_back(second_row); 
        }
    }

    com_temp.erase(com_temp.end() - (com_temp.size() - 6), com_temp.end());

    // add all elements
    for (int i = 0; i < atom_coordinates.size(); i++)
    {
        if (i > 5)
        {
            second_row = atom_coordinates[i] + com_temp[i-3];
            if (i == com_temp.size())
            {
                com_temp.push_back(second_row);
            }
        }
    }

    // add elements into 2d vector according to coordinate pattern
    std::vector<std::vector<double>>com(atomicmass.size());
    for (int i = 0; i < atomicmass.size(); i++)
    {
        std::vector<double>&inner_vector = com[i];
        for (int j = 0; j < 3; j++)
        {
            inner_vector.push_back(com_temp[3 * i + j]);
        }
    }

    for (int j = 0; j < com[0].size(); j++)
    {
        for (int i = com.size() - 3; i < com.size(); i++)
        {
            com[i][j] = com[i][j] / cluster_object.get_molmass(); 
        }
    }

    // extract last com row
    std::vector<double>last_comrow = com[com.size() - 1];
    std::vector<double>xyz_shift_temp; 
    for (int i = 0; i < xyz_matrix.size(); i++)
    {
        for (int j = 0; j < xyz_matrix[0].size(); j++)
        {
            xyz_shift_temp.push_back(xyz_matrix[i][j] - last_comrow[j]);
        }
    }

    // for (int i = 0; i < xyz_shift_temp.size(); i++)
    // {
    //     std::cout << xyz_shift_temp[i] << std::endl;
    // }

    // convert 1d vector to 2d vector
    std::vector<std::vector<double>>xyz_shift(atomicmass.size());
    for (int i = 0; i < atomicmass.size(); i++)
    {
        std::vector<double>&inner_vector = xyz_shift[i];
        for (int j = 0; j < 3; j++)
        {
            inner_vector.push_back(xyz_shift_temp[3 * i + j]); 
        }
    }

    const size_t M = 3; 
    const size_t N = 3; 
    double inertia[M][N] {};

    for (int i = 0; i < atomicmass.size(); i++)
    {
        inertia[0][0] = inertia[0][0] + atomicmass[i] * (pow(xyz_shift[i][1], 2) + pow(xyz_shift[i][2], 2));
        inertia[1][1] = inertia[1][1] + atomicmass[i] * (pow(xyz_shift[i][0], 2) + pow(xyz_shift[i][2], 2));
        inertia[2][2] = inertia[2][2] + atomicmass[i] * (pow(xyz_shift[i][0],2) + pow(xyz_shift[i][1], 2)); 
        inertia[0][1] = inertia[0][1] - atomicmass[i] * (xyz_shift[i][0] * xyz_shift[i][1]);
        inertia[0][2] = inertia[0][2] - atomicmass[i] * (xyz_shift[i][0] * xyz_shift[i][2]);
        inertia[1][2] = inertia[1][2] - atomicmass[i] * (xyz_shift[i][1] * xyz_shift[i][2]);
    }

    inertia[1][0] = inertia[0][1];
    inertia[2][0] = inertia[0][2];
    inertia[2][1] = inertia[1][2];

    double identity[3][3];
    identity[0][0] = 1;
    identity[1][1] = 1;
    identity[2][2] = 1; 
    identity[0][1] = 0;
    identity[1][0] = 0;
    identity[0][2] = 0;
    identity[2][0] = 0;
    identity[1][2] = 0;
    identity[2][1] = 0;

    double p1;
    double q;
    double p2;
    double p;
    double r;
    double phi;
    double tmp;
    double eig[3];
    std::vector<double>eig_temp; 
    double B[3][3];

    p1 = pow(inertia[0][1],2) + pow(inertia[0][2], 2) + pow(inertia[1][2],2);

    if (p1 <= global_eps)
    {
        eig_temp[0] = inertia[0][0];
        eig_temp[1] = inertia[1][1];
        eig_temp[2] = inertia[2][2];
    }
    else
    {
        q = (inertia[0][0] + inertia[1][1] + inertia[2][2]) / 3.00;
        p2 = pow(inertia[0][0] - q, 2) + pow(inertia[1][1] - q, 2) + pow(inertia[2][2] - q, 2) + (2.00 *p1);
        p = sqrt(p2 / 6.00);

        for (int j = 0; j < 3; j++)
        {
            for (int i = 0; i < 3; i++)
            {
                B[i][j] = (inertia[i][j] - q * identity[i][j]) / p;
            }
        }

        r = B[0][0] * B[1][1] * B[2][2] + B[0][1] * B[1][2] * B[2][0] + B[0][2] * B[1][0] * B[2][1] - B[0][2] * B[1][1] * B[2][0] - B[0][1] * B[1][0] * B[2][2] - B[0][0] * B[1][2] * B[2][1];
        r = 0.5 * r; 

        if (r <= -1.00)
        {
            phi = pi / 3.00;
        }
        else if (r>= 1.00)
        {
            phi = 0.00; 
        }
        else
        {
            phi = acos(r) / 3.00; 
        }
        
        eig[0] = q + 2.0 * p * cos(phi);
        eig[2] = q + 2.0 * p * cos(phi + (2.0 * pi/ 3.0));
        eig[1] = 3.0 * q - eig[0] - eig[2];
    }

    if (eig[0] < eig[1])
    {
        tmp = eig[1];
        eig[1] = eig[0];
        eig[0] = tmp; 
    }

    if (eig[0] < eig[2])
    {
        tmp = eig[2];
        eig[2] = eig[0];
        eig[0] = tmp;
    }

    if (eig[1] < eig[2])
    {
        tmp = eig[2];
        eig[2] = eig[1];
        eig[1] = tmp;
    }

    int eig_count = 0; 
    if (eig[1] <= global_eps)
    {
        eig_count += 1;
    }
    else if (eig[0] <= global_eps)
    {
        eig_count += 1; 
    }
    else if (eig[2] <= global_eps)
    {
        eig_count += 1; 
    }


    int n = 3 - eig_count;
   
    if (n == 3)
    {
        cluster_object.set_inertia(eig[0], 0);
        cluster_object.set_inertia(eig[1], 1);
        cluster_object.set_inertia(eig[2], 2);
    }

    else if (n == 2)
    {
        cluster_object.set_inertia(eig[0], 0);
        // cluster_object.inertia_c[0] = eig[0];
        cluster_object.set_linear(true); 
    }

    else if (n == 0)
    {
        cluster_object.set_atom(true); 
    }
    cluster_object.set_linear(false);
    cluster_object.set_atom(false); 
}
// this function removes clusters in reverse order (from smallest to largest)
std::vector<cluster_c>cluster_range(std::vector<cluster_c> &v, int m, int n)
{
    int k = n - m + 1;
 
    auto it = v.cbegin() + m;
    while (it != v.cend() && k--) {
        it = v.erase(it);
    }
    return v; 
}

// this function removes clusters in reverse order (from largest to smallest)
std::vector<cluster_c>cluster_reverse_range(std::vector<cluster_c> &v, int m, int n)
{
    int k = n - m + 1; 
    auto it = v.cend() - n - 1; 
    while (it != v.cbegin() && k--)
    {
        it = v.erase(it);
    }
    return v; 
}

void makeCombiUtil(std::vector<std::vector<int> >& ans,
    std::vector<int>& tmp, int n, int left, int k)
{
    // Pushing this vector to a vector of vector
    if (k == 0) {
        ans.push_back(tmp);
        return;
    }
 
    // i iterates from left to n. First time
    // left will be 1
    for (int i = left; i <= n; ++i)
    {
        tmp.push_back(i);
        makeCombiUtil(ans, tmp, n, i + 1, k - 1);
 
        // Popping out last inserted element
        // from the vector
        tmp.pop_back();
    }
}
 
// Prints all combinations of size k of numbers
// from 1 to n.
std::vector<std::vector<int> > makeCombi(int n, int k)
{
    std::vector<std::vector<int> > ans;
    std::vector<int> tmp;
    makeCombiUtil(ans, tmp, n, 1, k);
    return ans;
}

// uses reverse_range function above. 
std::map<int, std::vector<cluster_c>>clustersetlist(std::vector<cluster_c>&clusterset, int &combination_no)
{
    std::vector<std::vector<int>>combi = makeCombi(clusterset.size(), combination_no);
    std::vector<std::vector<int>>new_combi; 
    std::map<int, std::vector<cluster_c>>clustersetmap; 
    new_combi.resize(combi.size(), std::vector<int>(combi[0].size()));

    for (int i = 0; i < combi.size(); i++) 
    {
        std::vector<cluster_c>temp_clusterset;
        for (int j = 0; j < combi[i].size(); j++) 
        {
            // reassign combinations according to index of clusterset vector. 
            new_combi[i][j] = (combi[i][j]) - 1;

            // assign all clusterset combinations into individual vectors. 
            // only get combinations where monomer is present. 
            if (clusterset[new_combi[i][0]].get_composition() == 1)
            {
                temp_clusterset.emplace_back(clusterset[new_combi[i][j]]);
            }

            // assign combinations to clustersetmap. 
            if (j == combi[i].size() - 1 and ! temp_clusterset.empty())
            {
                clustersetmap[i] = temp_clusterset; 
            }  
        }
    }

    return clustersetmap; 
}

void check_clusterset(std::vector<cluster_c>&clusterset, input_c &input_object)
{
    // check monomer count
    if (clusterset[0].get_composition() != input_object.get_components())
    {
        std::cout << "invalid number of monomers" << std::endl; 
    }

    for (int i = 0; i < clusterset.size(); i++)
    {
        // check composition
        if (clusterset[i].get_composition() < 0)
        {
            std::cout << "unphysical composition value for " << clusterset[i].get_label() << std::endl; 
        }

        // check sigma 
        if (clusterset[i].get_sigma() < 0)
        {
            std::cout << "unphysical sigma value for " << clusterset[i].get_label() << std::endl; 
        }

        // check anharmonicity
        if (clusterset[i].get_anharmonicity() < 0.0)
        {
            std::cout << "unphysical anharmonicity value for " << clusterset[i].get_label() << std::endl; 
        }

        //check frequency scaling factor 
        if (clusterset[i].get_fscale() <= 0.00)
        {
            std::cout << "unphysical frequency scale value for " << clusterset[i].get_label() << std::endl; 
        }

        // check frequencies 
        std::vector<double>frequencies = clusterset[i].get_frequencies(); 
        for (int j = 0; j < frequencies.size(); j++)
        {
            if (frequencies[j] < 0.0)
            {
                std::cout << "unphysical frequency for " << clusterset[i].get_label() << std::endl;  
            }
        }
        
        // check volume 
        if  (clusterset[i].get_volume() <= 0.00)
        {
            std::cout << "unphysical volume for " << clusterset[i].get_label() << std::endl; 
        }
    }
}