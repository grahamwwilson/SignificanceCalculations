#include <iostream> 
#include <algorithm> //sort
#include <cmath>
#include <TRandom3.h>
#include <TMath.h>   //TMath::Prob
#include <TH1D.h>
#include <TFile.h>
#include <vector>    
#include <cstdlib>
using namespace std;

int main(int argc, char *argv[]){

   int MILLION = 1000000;
   int NGENERATED = 1*MILLION;   // Set so that pvalue MC fractional uncertainty is better than 2.5% (ie 1600 successes) for even the worst bins.
   
// Set initial values
   double MUB = 904.684;
   int NOBS = 787;
   double FRACERROR = 0.03627;
   double pullf, pull;   
   double ratio;
   double varfit;
   int NFACTOR=1;
   
   unsigned int seed = 4359;
// Read values from input   
   if(argc==6){
       NFACTOR = atoi(argv[1]);              // Number of millions
       NGENERATED = NFACTOR*MILLION;
       NOBS = atoi(argv[4]);                 // Number of observed events
       MUB = double(atof(argv[5]));          // Expected background
       pullf = double(atof(argv[2]));        // (O-E)/sqrt(E+VF) pull 
       pull = double(atof(argv[3]));         // (O-E)/sqrt(E) pull
// Calculate the fractional error (sqrt(VF)/E) from the above two pulls
       ratio = pull*pull/(pullf*pullf);
       varfit = (ratio-1.0)*MUB;
       FRACERROR = sqrt(varfit)/MUB;         // Used to define the smearing of the expected value
       seed = seed + NOBS;                   // Change the seed between runs.
   }
   
   cout << "Generating " << NFACTOR << " million toys" << endl;
   cout << "Using seed " << seed << endl;
   cout << "MUB  = "      << MUB << endl;
   cout << "NOBS = "      << NOBS << endl;
   cout << "FRACERROR = " << FRACERROR << endl;
   double SIGMAB = FRACERROR*MUB;
   cout << "SIGMAB " << SIGMAB << endl;   
   
// Note depending on the actual p-value one may need a lot of toys, especially for very small p-values
// Plan for at least 1600 events in the tails corresponding to 2.5% errors.  

   TRandom3 *rg = new TRandom3(seed);
   TFile *f = new TFile("out.root","RECREATE");
   TH1D* hist = new TH1D("hist","hist",2000,-0.5,1999.5);   
   
   int nuppertail = 0;
   int nlowertail = 0;

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
   
//   hist->Draw();
   hist->Write();
   f->Write();
   f->Close();
   double pvalueu = double(nuppertail)/double(NGENERATED);
   double dpu=sqrt(pvalueu*(1.0-pvalueu)/double(NGENERATED));  //binomial error
   double qvaluel = double(nlowertail)/double(NGENERATED);
   double dql=sqrt(qvaluel*(1.0-qvaluel)/double(NGENERATED));  //binomial error   

   cout << "nuppertail " << nuppertail << " p-value: " << pvalueu << " +- " << dpu << " probability of observing nobs or more events " << endl;
   cout << "nlowertail " << nlowertail << " q-value: " << qvaluel << " +- " << dql << " probability of observing nobs or less events " << endl;   
   
   cout << "Translating to a Gaussian z-score " << endl;
   if(NOBS >= MUB){
       cout << "z-score upper (ie. signed pull equivalent): " << TMath::NormQuantile(1.0-pvalueu) << " sigma " << endl;
   }    
   else{
       cout << "z-score lower (ie. signed pull equivalent): " << -TMath::NormQuantile(1.0-qvaluel) << " sigma " << endl;
   }
   cout << " " << endl;
}
