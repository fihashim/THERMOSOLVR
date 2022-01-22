#include "auxiliary2.h"
#include"atomic_data.h"
#include "input.h"
#include "qce.h"
#include "cluster2.h"
#include<iostream>
#include<vector>
#include<fstream>
#include<filesystem>

void range_c::set_range(range_c &r, const double& first, const double& last, const double& num)
{
    r.first = first;
    r.last = last;
    r.num = num;
    
    if (num != 1)
    {
        r.delta = (last - first) / (num - 1);
    }
    else
    {
        r.delta = 0.00;
    }
};

double ln_factorial(double x)
{
    std::vector<double>c{0.99999999999999709182,
                57.156235665862923517,
                -59.597960355475491248,
                14.136097974741747174, 
                -0.49191381609762019978,
                0.33994649984811888699e-4,
                0.46523628927048575665e-4, 
                -0.98374475304879564677e-4,
                0.15808870322491248884e-3, 
                -0.21026444172410488319e-3, 
                0.21743961811521264320e-3, 
                -0.16431810653676389022e-3, 
                0.84418223983852743293e-4, 
                -0.26190838401581408670e-4, 
                0.36899182659531622704e-5};
    double sqrt2pi = sqrt(2.0*pi);
    double g = 607.0 / 128.0;
    double tmp;
    double series;
    double denominator;
    double ln_factorialVal; 

    tmp = x + g + 0.50;
    tmp = (x + 0.50) * log(tmp) - tmp;

    denominator = x;
    series = c[0];
    for (int i = 1; i < c.size() - 1; i++)
    {
        denominator = denominator + 1.00;
        series = series + c[i]/denominator;
    }

    ln_factorialVal = tmp + log(sqrt2pi*series);
    return ln_factorialVal; 
};

std::vector<double>derivative(std::vector<double>y, std::vector<double>x)
{
    std::vector<double>d;
    std::vector<long double>ydiff;
    std::vector<long double>ylong;
    std::vector<double>xdiff; 
    //std::cout.precision(dbl::max_digits10);
    d.resize(y.size());
    ydiff.resize(y.size());
    xdiff.resize(x.size());
    int i;

    if (y.size() == 1)
    {
        d[0] = 0.00; 
    }
    else
    {
        d[0] = (y[1] - y[0])/(x[1]-x[0]);
        for (int i = 1; i < y.size() - 1; i++)
        {
            d[i] = (y[i+1] - y[i-1])/(x[i+1] - x[i-1]);
        }
        d[y.size() - 1] = (y[y.size() - 1] - y[y.size() - 2]) / (x[y.size() - 1] - x[y.size()-2]);   
    }
    return d; 
}

bool check_range(range_c &r)
{
    bool check_range = true;
    if (r.last < r.first)
    {
        check_range = false;
    }

    if (r.num < 1)
    {
        check_range = false;
    }
    return check_range; 
}



void welcomemsg(int& option)
{
    std::cout << " Hello from THERMOSOLVR. " << std::endl;
    std::cout << std::endl; 
    std::cout << " Would you like to create cluster/input files for THERMOSOLVR? " << std::endl; 
    std::cout << " Enter: (1) Run combinations clusterset search on QCE " << std::endl; 
    std::cout << "        (2) Run QCE" << std::endl; 
    std::cin >> option;
    if (option == 1 or option == 2)
    {
        return; 
    }
    else
    {
        std::cerr << "option invalid" << std::endl;
        exit(1); 
    }
}

void process_files(input_c &input_object, std::vector<cluster_c>&clusterset, cluster_c &monomer, int argc, char** argv)
{
    input_object = process_input(argv[1]);
    clusterset = process_cluster(argv[2], input_object);
    monomer = clusterset[0];
}


void clusterset_testing(std::vector<cluster_c>&clusterset, cluster_c &monomer, input_c &input, reference_c &reference, isobar_c &best_ib, int &combination_no)
{
    // get clusterset list
    std::map<int,std::vector<cluster_c>>clustersetmap; 
    clustersetmap = clustersetlist(clusterset, combination_no);

    // create folders for each clusterset tested
    namespace fs = std::filesystem;
    for (int i = 0; i < clustersetmap.size(); i++)
    {
        fs::path p = fs::current_path();
        fs::create_directories("test_clustersets/clusterset" + std::to_string(i));
    }

    int counter = 0; // use to move clustersets into designated folders in qce_finalize. 
    // write all clusterset combinations to a file. 
    std::ofstream clustersetlist; 
    clustersetlist.open("clustersetlist.txt");
    for (auto &cluster : clustersetmap)
    {
       
        std::vector<cluster_c>test_clusterset = cluster.second; 
        qce_prepare(input, monomer, test_clusterset, reference);
        best_ib = qce_start(test_clusterset,monomer, reference); 
        qce_finalize(best_ib, test_clusterset, reference, input, monomer, 1, counter);        
        counter += 1;
        clustersetlist << "clusterset" << cluster.first << " : ";
        for (int i = 0; i < test_clusterset.size(); i++)
        {
            clustersetlist << " [" << test_clusterset[i].get_label() << "]";
            if (i == test_clusterset.size() - 1)
            {
                clustersetlist << std::endl; 
            }
        }
    }
    clustersetlist.close(); 
}

void progress_bar(int &nr_isobars_computed, int &nr_isobars_total, bool &progress_bar, bool newline)
{
    int nbars;

    if (not progress_bar)
    {
        return; 
    }

    nbars = int((double(nr_isobars_computed) * 100) / nr_isobars_total);
    std::cout << "[";
    for (int i = 1; i < nbars - 1; i++)
    {
        std::cout << "-";
    }
    std::cout << ">";
    
    for (int i = nbars + 1; i < nr_isobars_total + 1; i++)
    {
        std::cout << " "; 
    }
    std::cout << "]" << nbars << "%"; 
    if (newline)
    {
        std::cout << std::endl;
    }
}