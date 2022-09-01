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
   const unsigned int NGENERATED = 10000000;   // Set so that pvalue MC fractional uncertainty is 1% for default assumptions that lead to pvalue circa 8.2e-5.
   
// Note depending on the actual p-value one may need a lot of toys, especially for very small p-values.  

   TRandom3 *rg = new TRandom3(seed);
   TFile *f = new TFile("out.root","RECREATE");
   TH1D* hist = new TH1D("hist","hist",2000,-0.5,1999.5);   
   
   const double MUB = 904.684;
//   const double FRACERROR = 0.032176;     // 10% error
   const double FRACERROR = 0.03627;
   const double SIGMAB = FRACERROR*MUB;  
   const int NOBS = 787;
   int nuppertail = 0;
   int nlowertail = 0;
   
   std::cout << "SIGMAB " << SIGMAB << std::endl;

// Simulate a Poisson distribution with the underlying mean being a Gaussianly distributed random number.
   for (int i=0; i<NGENERATED; i++){
        double mu = MUB;
        if(FRACERROR > 1.0e-4){
            mu = rg->Gaus(MUB, SIGMAB);  // Choose Poisson mean, mu, from Gaussian with mean and rms of MUB and SIGMAB
        }
        int n = rg->Poisson(mu);         // Generate Poisson distributed random number, n, based on Poisson mean of mu  
        if(n >= NOBS)nuppertail++;       // Count toys with n exceeding or equal to the observed counts
        if(n <= NOBS)nlowertail++;       // Count toys with n less than or equalt to the observed counts
        hist->Fill(double(n));
   }
   
   hist->Draw();
   f->Write();
   f->Close();
   double pvalueu = double(nuppertail)/double(NGENERATED);
   double dpu=sqrt(pvalueu*(1.0-pvalueu)/double(NGENERATED));  //binomial error
   double qvaluel = double(nlowertail)/double(NGENERATED);
   double dql=sqrt(qvaluel*(1.0-qvaluel)/double(NGENERATED));  //binomial error   

   std::cout << "nobs =           "  << NOBS << std::endl;
   std::cout << "muB =            "  << MUB << std::endl;
   std::cout << "Error fraction = "  << FRACERROR << std::endl;
   std::cout << "Generating toys  "  << NGENERATED << std::endl;
   std::cout << "nuppertail " << nuppertail << " p-value: " << pvalueu << " +- " << dpu << " probability of observing nobs or more events " << std::endl;
   std::cout << "nlowertail " << nlowertail << " q-value: " << qvaluel << " +- " << dql << " probability of observing nobs or less events " << std::endl;   
   
   std::cout << "Translating to a Gaussian z-score " << std::endl;
   if(NOBS >= MUB){
       std::cout << "z-score upper (ie. signed pull equivalent): " << TMath::NormQuantile(1.0-pvalueu) << " sigma " << std::endl;
   }    
   else{
       std::cout << "z-score lower (ie. signed pull equivalent): " << -TMath::NormQuantile(1.0-qvaluel) << " sigma " << std::endl;
   }    
}
