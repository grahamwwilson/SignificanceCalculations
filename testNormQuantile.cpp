//
// Test 1 of MyQuantile
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
 
    CLI::App app{"Test MyQuantile with small values of p (Test 1)"};
    
    double p=2.0;
    double divisor=2.0;
    long int nmax=1080;
    
    app.add_option("--p", p, "Starting p-value");
    app.add_option("--d", divisor, "Divisor value"); 
    app.add_option("--n", nmax, "Max number of evaluations"); 
              
    CLI11_PARSE(app, argc, argv);
    
    std::cout << "p " << p << std::endl;
    std::cout << "divisor " << divisor << std::endl;
    std::cout << "nmax " << nmax << std::endl;

    double pthis = p;
    int num = 0;
    while (num < nmax){ 
        num += 1;
        pthis = pthis/divisor;
        double Z = MyQuantile(pthis);
    }
    std::cout << "Ending program with pthis = " << pthis << " after " << num << " iterations" << std::endl;
       
    return 0;
    
}
