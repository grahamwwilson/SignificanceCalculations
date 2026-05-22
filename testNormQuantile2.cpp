//
// Test 2 of MyQuantile
//
#include "CLI11.hpp"
#include <iostream> 
#include <algorithm> //std::sort
#include <cmath>
#include <TMath.h>   //TMath::Prob
#include <vector>
#include <fstream>   
#include <cstdlib>
#include <string>
#include <iomanip>
#include "MyQuantile.h"

int main(int argc, char** argv){

// Simple test program to loop through various p-values and call MyQuantile
 
    CLI::App app{"Test MyQuantile with large values of p (Test 2)"};
    
    double divisor=2.0;
    long int nmax=60;
    
    app.add_option("--d", divisor, "Divisor value"); 
    app.add_option("--n", nmax, "Max number of evaluations"); 
              
    CLI11_PARSE(app, argc, argv);
    
    std::cout << "divisor " << divisor << std::endl;
    std::cout << "nmax " << nmax << std::endl;

    double delta = 2.0;
    double pthis;

    int num = 0;
    while (num < nmax){ 
        num += 1;
        delta = delta/divisor;
        double pthis = 1.0 - delta;        
        double Z = MyQuantile(pthis);
    }
    std::cout << "Ending program with pthis = " << pthis << " delta " << delta << " after " << num << " iterations" << std::endl;
       
    return 0;
    
}
