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
   unsigned int ielement = 12345;
   if(argc==3){
      seed = atoi(argv[1]);
   }
   std::cout << "Using seed " << seed << std::endl;

   const unsigned int N = 100;
   const unsigned int NGENERATED = 10000000;

   TRandom3 *rg = new TRandom3(seed);
//   TFile *f = new TFile("out.root","RECREATE");
//   TH1D* hist = new TH1D("hist","hist",100,0.0,200.0);
//   TH1D* hpvalue = new TH1D("hpvalue","hpvalue",100,0.0,1.0);
   double mu;
   int n;
   
   const double MUB = 3.38903;
   const double FRACERROR = 0.10;     // 10% error
   const double SIGMAB = FRACERROR*MUB;  
   const int NOBS = 13;
   int ncount = 0;
   
   std::cout << "SIGMAB " << SIGMAB << std::endl;

// Simulate a Poisson distribution with the underlying mean being a Gaussianly distributed random number.
   for (int i=0; i<NGENERATED; i++){
        mu = MUB;
        if(SIGMAB/MUB > 1.0e-4){
           mu = rg->Gaus(MUB, SIGMAB);
        }
        n = rg->Poisson(mu);
        if(n >= NOBS)ncount++;        //Count toys with n exceeding or equal to the observed counts
   }
//   hist->Draw();                    // Should be distributed like chi^2_N  
//   hpvalue->Draw();                 // Should be uniform from [0,1]
//   f->Write();
   
   double pvalue = double(ncount)/double(NGENERATED);

   std::cout << "nobs =           "  << NOBS << std::endl;
   std::cout << "muB =            "  << MUB << std::endl;
   std::cout << "Error fraction = "  << FRACERROR << std::endl;
   std::cout << "Generating toys  "  << NGENERATED << std::endl;
   std::cout << "ncount " << ncount << " p-value: " << pvalue << " probability of observing nobs or more events " << std::endl;
   
   std::cout << "Translating to a Gaussian z-score " << std::endl;
   std::cout << "z-score (ie. signed pull equivalent): " << TMath::NormQuantile(1.0-pvalue) << " sigma " << std::endl;

}
