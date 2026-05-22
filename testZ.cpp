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

std::pair<double,double> Compute(double mu,int kmin, int kmax){

// Input quantities
// mu:   Chosen Poisson mean
// kmin: Minimum k
// kmax: Maximum k
//
// For non-negative integers k, satisfying kmin <= k <= kmax, 
// evaluate the Poisson probability pk = Poisson(k; mu)
// and compute the mean value and standard deviation for the symmetric Z-score variable 
// through summation

   std::cout << "mu " << mu << " kmin " << kmin << " kmax " << kmax << std::endl;
   
   double emu  = TMath::Exp(-mu);
   double bess = TMath::BesselI(0,2.0*mu);
   
// Historical debugging check? not used in current calculation.
   std::cout << "I0 (2 mu) " << bess << std::endl;
// Apparently Sum(k=0, infty) pk**2 = Exp(-2 mu)*I0(2mu).
   
   std::vector <std::pair<int,double>> vp;
   std::vector <std::pair<int,double>> vQ;
   std::vector <std::pair<int,double>> vP;
   
   double sumps=0.0;
   double sumpsps = 0.0;
   double sumz = 0.0;
   double sumzz = 0.0;
   
   double partialsum = 0.0;
   
   double pk;
   
// Check elements of the initialization calculation. 
// Looks like emu can underflow for big enough mu. For example mu=300, gives 5.14e-131. 
// Apparently smallest positive double is 2.225e-307.
   std::cout << " Exponential factor " << emu << std::endl;
   std::cout << " pow(mu,k)   factor " << std::pow(mu,kmin) << std::endl;
   std::cout << " Gamma(k+1)  factor " << TMath::Gamma(double(kmin+1)) << std::endl;   
   
   for (int k=kmin; k<=kmax; k++){
       if(k == kmin){
           pk = emu*std::pow(mu,k)/TMath::Gamma(double(k+1));
       }
       else{
           pk = vp.back().second*(mu/double(k)); // Use recursive definition
       }

// p_{S} calculation. This runs into floating point issues at around 2e-16 if first few terms are of order one.
// May be better to switch to a tail calculation whenever eg p(k) < 1.0e-12;
       double ps = 1.0;
       for (auto & el : vp){
           ps -= el.second;
       }
       ps -= 0.5*pk;
       
       auto pairk = std::make_pair(k, pk); vp.push_back(pairk); // Only push back new value AFTER ps calculation
       auto pairQ = std::make_pair(k, ps); vQ.push_back(pairQ);       
   }
   
// Now calculate tail probability directly using whatever range of (kmin, kmax) we decided to use
    
   for (int j=kmin; j<=kmax; j++){
       double P = 0.0;
       for (auto & el : vp){
            auto k  = el.first;
            auto pk = el.second;
            if( k==j) P += 0.5*pk;
            if( k>j) P += pk;
       }
       auto pairP = std::make_pair(j, P); vP.push_back(pairP);       
   }
   
// Now loop through all values   
   
   for (int j=kmin; j<=kmax; j++){   
      
       int idx = j-kmin;
       auto k = vp[idx].first;
       auto pk = vp[idx].second;
       auto psQ = vQ[idx].second;
       auto psP = vP[idx].second;
      
       double zQ = MyQuantile(psQ);
       double zP = MyQuantile(psP);
       
       double z;

// Choose representation with better numerical precision
       if(psQ > 0.5){
// cumulative subtraction for large probabilities
          z = zQ;
       }
       else{
// explicit tail summation for small probabilities
          z = zP;
       }
       
       double psum = psP + 0.5*pk;
       
       std::cout << "k: " << k << " p(k) " << pk << " psQ " << psQ 
                 << " psP " << psP << " ZQ " << zQ << " ZP " << zP << " Chosen z " << z 
                 << " qsum " << std::scientific << std::setprecision(16) << 1.0-psum << std::endl;
       
       sumps += pk*psQ;
       sumpsps += pk*psQ*psQ;
       if(abs(z)<10.0){
          sumz += pk*z;
          sumzz += pk*z*z;
          partialsum += pk;
       }   
   }
   
   std::cout << "partial sum = " << std::fixed << std::setprecision(16) << partialsum << std::endl;
   
   double varx = sumpsps - sumps*sumps;
   double varz = sumzz - sumz*sumz;
   
   std::cout << " For x = p_{S} " << std::endl;  
   std::cout << " E(x)   = " << std::fixed << std::setprecision(16) << sumps << std::endl;
   std::cout << " E(x^2) = " << std::fixed << std::setprecision(16) << sumpsps << std::endl;
   std::cout << " Var(x) = " << std::fixed << std::setprecision(16) << sumpsps - sumps*sumps << std::endl;
   std::cout << " SD(x)  = " << std::fixed << std::setprecision(16) << sqrt(varx) << std::endl;
   std::cout << " NSD(x) = " << std::fixed << std::setprecision(16) << sqrt(12.0*varx) << std::endl;
   
   std::cout << " ---------------------------" << std::endl;
   
   std::cout << " For x = Z_{S} " << std::endl;  
   std::cout << " E(x)   = " << std::fixed << std::setprecision(16) << sumz << std::endl;
   std::cout << " E(x^2) = " << std::fixed << std::setprecision(16) << sumzz << std::endl;
   std::cout << " Var(x) = " << std::fixed << std::setprecision(16) << varz << std::endl;
   std::cout << " SD(x)  = " << std::fixed << std::setprecision(16) << sqrt(varz) << std::endl;
   
   auto result = std::make_pair(sumz, sqrt(varz));
   return result;
      
}

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
 
    CLI::App app{"Compute Z-score expected values, Version 6"};
    
    double mumin=1.0;
    double mumax=5.0;
    double dmu = 0.1;
    int kmin=0;
    int kmax=40;
    std::string part = "P1";
    
    app.add_option("--mumin", mumin, "Poisson mean minimum in scan");
    app.add_option("--mumax", mumax, "Poisson mean maximum in scan"); 
    app.add_option("-d,--dmu",dmu, "Poisson mean increment in scan");
    app.add_option("-f,--kmin", kmin, "Minimum count, kmin");      
    app.add_option("-l,--kmax", kmax, "Maximum count, kmax");   
    app.add_option("-p,--part", part, "Part string for file name"); 
              
    CLI11_PARSE(app, argc, argv);
    
    std::cout << "mumin " << mumin << std::endl;
    std::cout << "mumax " << mumax << std::endl;  
    std::cout << "dmu   " << dmu << std::endl;        
    std::cout << "kmin  " << kmin << std::endl;
    std::cout << "kmax  " << kmax << std::endl;
    std::cout << "part  " << part << std::endl;    
    
    std::ofstream foutmean("ZScore-Calculator-Mean-"+part+".tpl"); 
    std::ofstream foutsd("ZScore-Calculator-SD-"+part+".tpl");           

    double mu = mumin;
    int num = 0;
    while (mu < mumax){ 
        num += 1;
        if(num != 1)mu += dmu;
        auto result = Compute(mu,kmin,kmax);
        foutmean << mu << " " << result.first << " " << result.second << std::endl;
        foutsd   << mu << " " << result.second << " " << result.first << std::endl;                
        std::cout << " " << std::endl;
    }
       
    foutmean.close();
    foutsd.close();   
       
    return 0;
    
}
