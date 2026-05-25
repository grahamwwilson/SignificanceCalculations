// Calculate expected bias of the mean and variance of the Z-score estimate 
// using the pure Poisson case with no background uncertainty treatment (or sigma_mubhat=0)
//
// So instead of running toys can simply do appropriate (taking care of numerical issues) summing 
// over all of the discrete outcomes, where k ~ Po(mu).
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
#include "ComputeZscoreLocationAndScaleCorrection.h"   // Contains std::pair<double,double> Compute(double mu, int kmin, int kmax, int printlevel=0)

int main(int argc, char** argv){

// For pure Poisson-only case for a given Poisson-mean, mu, 
// we evaluate the mean and standard deviation that results from using the symmetrized Z-score. 
// This is "Poisson-only" equivalent to setting sigma_mub_hat = 0.

// So all we need to do is 
// 1. loop over all possible observed values k
// 2.    for each k, compute Z_k
// 3.    weight the Z_k value by the corresponding Poisson probability, Po(k; mu)
// 4. Sum over all k.

// We need to be careful about various numerical gotchas given the in principle need to 
// do infinite sums. So instead of summing from k=0, k=infinity, we sum from k=kmin, k=kmax.
// For small values of mu we of course can start at k=0.
 
    CLI::App app{"Compute Z-score expected values, Version 7"};
    
    double mumin=1.0;
    double mumax=300.0;
    double dmu = 1.0;
    bool kcontrol = false;
    int kmin=0;
    int kmax=40;
    int printlevel=0;
    int ztreatment=1;
    std::string part = "All";
    
    app.add_option("--mumin", mumin, "Poisson mean minimum in scan");
    app.add_option("--mumax", mumax, "Poisson mean maximum in scan"); 
    app.add_option("-d,--dmu",dmu, "Poisson mean increment in scan");
    app.add_option("-k,--kcontrol", kcontrol, "Boolean for using [kmin, kmax]");
    app.add_option("-f,--kmin", kmin, "Minimum count, kmin");      
    app.add_option("-l,--kmax", kmax, "Maximum count, kmax");   
    app.add_option("-s,--part", part, "Part string for file name"); 
    app.add_option("-p,--printlevel", printlevel, "Print level");
    app.add_option("-z,--ztreatment", ztreatment, "Z treatment [0 = none, 1 = switch (default)");
              
    CLI11_PARSE(app, argc, argv);
    
    std::cout << "mumin " << mumin << std::endl;
    std::cout << "mumax " << mumax << std::endl;  
    std::cout << "dmu   " << dmu << std::endl; 
    std::cout << "kcontrol " << kcontrol << std::endl;      
    std::cout << "kmin  " << kmin << std::endl;
    std::cout << "kmax  " << kmax << std::endl;
    std::cout << "part  " << part << std::endl;
    std::cout << "printlevel " << printlevel << std::endl;
    std::cout << "ztreatment " << ztreatment << std::endl;
    
    std::ofstream foutmean("ZScore-Calculator-Mean-"+part+".tpl"); 
    std::ofstream foutsd("ZScore-Calculator-SD-"+part+".tpl");           

    double mu = mumin;
    int num = 0;

// Set default value and declare result
    double zbias = 0.0;
    double zscalebias = 1.0;
    auto result = std::make_pair(zbias, zscalebias);

    while (mu < mumax){ 
        num += 1;
        if(num != 1)mu += dmu;
        if(kcontrol){
// Use specified kmin, kmax
            result = Compute(mu,kmin,kmax,printlevel);
        }
        else{
// Use kmin, kmax as hard-coded in ComputeZscoreLocationAndScaleCorrection. This is the default.
            result = ComputeZscoreLocationAndScaleCorrection(mu, ztreatment, printlevel);
        }
        foutmean << mu << " " << result.first << " " << result.second << std::endl;
        foutsd   << mu << " " << result.second << " " << result.first << std::endl;                
        std::cout << " " << std::endl;
    }
       
    foutmean.close();
    foutsd.close();   
       
    return 0;
    
}
