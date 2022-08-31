#include <iostream> 
#include <algorithm> //std::sort
#include <cmath>
#include <TRandom3.h>
#include <TMath.h>   //TMath::Prob
#include <TH1D.h>
#include <TFile.h>
#include <vector>    
#include <cstdlib>
//using namespace std;

int main(int argc, char *argv[]){

   unsigned int seed = 4359;
   if(argc==2){
      seed = atoi(argv[1]);
   }
   std::cout << "Using seed " << seed << std::endl;

   const unsigned int N = 100;
   const unsigned int NGENERATED = 125000000;   // Set so that pvalue MC fractional uncertainty is 1% for default assumptions that lead to pvalue circa 8.2e-5.
   
// Note depending on the actual p-value one may need a lot of toys, especially for very small p-values.  

   TRandom3 *rg = new TRandom3(seed);
   
   const double MUB = 3.38903;
   const double FRACERROR = 0.10;     // 10% error
   const double SIGMAB = FRACERROR*MUB;  
   const int NOBS = 13;
   int ncount = 0;
   
   std::cout << "SIGMAB " << SIGMAB << std::endl;

// Simulate a Poisson distribution with the underlying mean being a Gaussianly distributed random number.
   for (int i=0; i<NGENERATED; i++){
        double mu = MUB;
        if(FRACERROR > 1.0e-4){
            mu = rg->Gaus(MUB, SIGMAB);  // Choose Poisson mean, mu, from Gaussian with mean and rms of MUB and SIGMAB
        }
        int n = rg->Poisson(mu);         // Generate Poisson distributed random number, n, based on Poisson mean of mu  
        if(n >= NOBS)ncount++;           // Count toys with n exceeding or equal to the observed counts
   }
   
   double pvalue = double(ncount)/double(NGENERATED);
   double dp=sqrt(pvalue*(1.0-pvalue)/double(NGENERATED));  //binomial error

   std::cout << "nobs =           "  << NOBS << std::endl;
   std::cout << "muB =            "  << MUB << std::endl;
   std::cout << "Error fraction = "  << FRACERROR << std::endl;
   std::cout << "Generating toys  "  << NGENERATED << std::endl;
   std::cout << "ncount " << ncount << " p-value: " << pvalue << " +- " << dp << " probability of observing nobs or more events " << std::endl;
   
   std::cout << "Translating to a Gaussian z-score " << std::endl;
   std::cout << "z-score (ie. signed pull equivalent): " << TMath::NormQuantile(1.0-pvalue) << " sigma " << std::endl;

}
