// Demonstrate chi-squared values
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
      ielement = atoi(argv[2]);
   }
   std::cout << "Using seed " << seed << std::endl;
   std::cout << "Using element " << ielement << std::endl;

   const unsigned int N = 100;
   const unsigned int NGENERATED = 100000;
   std::cout << "Chi-squared with " << N << " degrees-of-freedom" << std::endl;
   std::cout << "Number of random variates generated " << NGENERATED << std::endl;
   std::vector<double> v;

   TRandom3 *rg = new TRandom3(seed);
   TFile *f = new TFile("out.root","RECREATE");
   TH1D* hist = new TH1D("hist","hist",100,0.0,200.0);
   TH1D* hpvalue = new TH1D("hpvalue","hpvalue",100,0.0,1.0);
   double zi;

// Simulate a chi-squared distribution with N degrees of freedom
   for (int i=0; i<NGENERATED; i++){
      double chisq = 0.0;
      for (int i=0;i<N;i++){
          zi = rg->Gaus(0.0,1.0);    // standardized normal random variates
          chisq += zi*zi;            // Add zi^2 N times
      }
      v.push_back(chisq);
      double pval = TMath::Prob(chisq,N);           
      hist->Fill(chisq);            // Histogram the chi-squared value
      hpvalue->Fill(pval);          // Histogram the upper-tail integral (int(chisq, infinity))
                                    // often called the p-value to assess the 
                                    // significance of the observed chi-squared - specific to this 
                                    // value of N - the number of degrees of freedom
   }
   hist->Draw();                    // Should be distributed like chi^2_N  
   hpvalue->Draw();                 // Should be uniform from [0,1]
   f->Write();

   std::cout << "Element " << ielement << " " << v[ielement] << std::endl;
  
   double x = v[ielement];
   int ndof = N;
   double pvalue=TMath::Prob(x,ndof);
   std::cout << "p-value of " << pvalue << " for this particular variate" << std::endl;

}
