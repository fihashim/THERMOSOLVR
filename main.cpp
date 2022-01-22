#include<iostream>
#include<string>
#include<tuple>
#include<chrono>
#include<map>
#include<set>
#include<iomanip>
#include<iterator>
#include "input.h"
#include "auxiliary2.h"
#include "config2.h"
#include "cluster2.h"
#include "thermo.h"
#include "global.h"
#include "partfunctions.h"
#include "qce.h"
#include "polynomial.h"
#include "atomic_data.h"

// to compile:
// clang++ -std=c++17 auxiliary2.cpp cluster.cpp config2.cpp input.cpp thermo.cpp polynomial.cpp partfunctions.cpp qce_2.cpp main.cpp -I /opt/homebrew/include -L /opt/homebrew/lib -fopenmp -o test
int main(int argc, char** argv)
{
    // int option = 1;
    int option;  
    int combination_no; 
    input_c input;  
    cluster_c monomer; 
    
    std::vector<cluster_c>clusterset; 
    isobar_c best_ib; 
    reference_c reference; 
    

    if (argc < 3)
    {
        std::cout << "Usage: " << "./main qce.input input.clusterset" << std::endl; 
        exit(1);
    }
    // print welcome options
    welcomemsg(option);

    // enter and check input files
    process_files(input, clusterset, monomer, argc, argv);

    if (option == 1) // run clusterset testing on qce. 
    {
        std::cout << "Enter: (2) 2-cluster combinations " << std::endl;
        std::cout << "       (3) 3-cluster combinations " << std::endl;
        std::cout << "       (4) 4-cluster combinations " << std::endl;
        std::cout << "       (5) 5-cluster combinations " << std::endl;
        std::cin >> combination_no; 
        // start clock
        auto start = std::chrono::high_resolution_clock::now();
        
        // run clusterset testing 
        clusterset_testing(clusterset, monomer, input, reference, best_ib, combination_no);

        //end clock
        auto end = std::chrono::high_resolution_clock::now();

        // print time
        std::chrono::duration<double, std::milli> ms_double = end - start;         
        std::cout << "Elapsed time in seconds: "
        << std::setprecision(3) << std::scientific << ms_double.count() / 1000
        << " s " << std::endl;
        std::cout << std::endl;   
    }   
    else if (option == 2) // run qce 
    {
        // perform qce calculations and start timer
        auto start = std::chrono::high_resolution_clock::now();
        qce_prepare(input, monomer, clusterset, reference);
        best_ib = qce_start(clusterset,monomer, reference); 
        qce_finalize(best_ib, clusterset, reference, input, monomer,2, 0);

        // end timer and print time. 
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> ms_double = end - start; 

        std::cout << "Elapsed time in seconds: "
            << std::setprecision(3) << std::scientific << ms_double.count() / 1000
            << " s " << std::endl;   
    }


}
