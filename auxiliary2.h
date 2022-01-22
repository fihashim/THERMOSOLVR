#ifndef AUXILIARY_H
#define AUXILIARY_H
#include<vector>

class cluster_c;
class input_c;
class reference_c; 
class isobar_c; 

class range_c
{
    public:
    double first;
    double last;
    double delta;
    int num; 
    void set_range(range_c &r, const double& first, const double& last, const double& num);
};

double ln_factorial(double x);
std::vector<double>derivative(std::vector<double>y, std::vector<double>x);
bool check_range(range_c &r);
void welcomemsg(int &option);
void process_files(input_c &input_object, std::vector<cluster_c>&clusterset, cluster_c &monomer, int argc, char** argv);
void clusterset_testing(std::vector<cluster_c>&clusterset, cluster_c &monomer, input_c &input, reference_c &reference, isobar_c &best_ib, int &combination_no);
void progress_bar(int &nr_isobars_computed, int &nr_isobars_total, bool &progress_bar, bool newline);

#endif 

