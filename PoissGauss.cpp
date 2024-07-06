// You may choose to use #include "CLI11.hpp" to pick up the repo version header file 
#include "CLI11.hpp"
//#include <CLI11.hpp>   // CLI11 command line interface stuff (V2.3.1). See https://github.com/CLIUtils/CLI11
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

void ToyChucker(int NFACTOR, int NOBS, double MUB, double FRACERROR, unsigned long int seed ){

   int MILLION = 1000000;
   int NGENERATED = 1*MILLION;   // Set so that pvalue MC fractional uncertainty is better than 2.5% (ie 1600 successes) for even the worst bins.
   
   double pullf, pull;   
   double ratio;
   double varfit;
   
   NGENERATED = NFACTOR*MILLION;
   
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
   int nequal = 0;

// Simulate a Poisson distribution with the underlying mean being a Gaussianly distributed random number.
   for (int i=0; i<NGENERATED; i++){
        double mu = MUB;
        if(FRACERROR > 1.0e-4){
            mu = rg->Gaus(MUB, SIGMAB);  // Choose Poisson mean, mu, from Gaussian with mean and rms of MUB and SIGMAB
        }
        int n = rg->Poisson(mu);         // Generate Poisson distributed random number, n, based on Poisson mean of mu  
        if(n >= NOBS)nuppertail++;       // Count toys with n exceeding or equal to the observed counts
        if(n <= NOBS)nlowertail++;       // Count toys with n less than or equalt to the observed counts
        if(n==NOBS)nequal++;        
        hist->Fill(double(n));
   }
   

   double pvalueu = double(nuppertail)/double(NGENERATED);
   double dpu=sqrt(pvalueu*(1.0-pvalueu)/double(NGENERATED));  //binomial error
   double qvaluel = double(nlowertail)/double(NGENERATED);
   double dql=sqrt(qvaluel*(1.0-qvaluel)/double(NGENERATED));  //binomial error   
   
   double pequal = double(nequal)/double(NGENERATED);
   double dpequal=sqrt(pequal*(1.0-pequal)/double(NGENERATED));  //binomial error   

   cout << "nuppertail " << nuppertail << " p-value: " << pvalueu << " +- " << dpu << " probability of observing nobs or more events " << endl;
   cout << "nlowertail " << nlowertail << " q-value: " << qvaluel << " +- " << dql << " probability of observing nobs or less events " << endl;   
   std::cout << "nequal     " << nequal     << " c-value: " << pequal  << " +- " << dpequal << " probability of observing exactly nobs events " << std::endl;    
   cout << "Translating to a Gaussian z-score " << endl;
// FIXME Add uncertainties to z-scores.
   double pvalueu_UP = pvalueu - dpu;
   double pvalueu_DOWN = pvalueu + dpu;
   double qvaluel_UP = qvaluel + dql;
   double qvaluel_DOWN = qvaluel - dql;   
   
   double uperror;
   double downerror;
   double estimate;
   double averror;
   
   if(NOBS >= MUB){
       estimate = TMath::NormQuantile(1.0-pvalueu);
       uperror = TMath::NormQuantile(1.0-pvalueu_UP)- estimate;
       downerror = estimate - TMath::NormQuantile(1.0-pvalueu_DOWN);
       averror = 0.5*(uperror + downerror);
       cout << "NOBS >= MUB:  z-score upper (ie. signed pull equivalent): " << fixed << setprecision(4) << estimate << " +- "  
            << averror << " sigma   (+: " << uperror << " )  (-:  " << downerror << " ) " << endl;
   }    
   else{
       estimate = -TMath::NormQuantile(1.0-qvaluel);
       uperror = -TMath::NormQuantile(1.0-qvaluel_UP)- estimate;
       downerror = estimate - (- TMath::NormQuantile(1.0-qvaluel_DOWN));
//       cout << " Three values " << -TMath::NormQuantile(1.0-qvaluel) << " " << -TMath::NormQuantile(1.0-qvaluel_UP) << " " << -TMath::NormQuantile(1.0-qvaluel_DOWN) << endl;
       averror = 0.5*(uperror + downerror);
       cout << "NOBS <= MUB:  z-score lower (ie. signed pull equivalent): " << fixed << setprecision(4) << estimate << " +- "  
            << averror << " sigma   (+: " << uperror << " )  (-:  " << downerror << " ) " << endl;   
   }
// Alternative formulation based on measured probabilities 
   if(pvalueu < qvaluel){
       estimate = TMath::NormQuantile(1.0-pvalueu);
       uperror = TMath::NormQuantile(1.0-pvalueu_UP)- estimate;
       downerror = estimate - TMath::NormQuantile(1.0-pvalueu_DOWN);
       averror = 0.5*(uperror + downerror);
       cout << "pu < ql    :  z-score upper (ie. signed pull equivalent): " << fixed << setprecision(4) << estimate << " +- "  
            << averror << " sigma   (+: " << uperror << " )  (-:  " << downerror << " ) " << endl;   
   }
   else if(qvaluel < pvalueu){
       estimate = -TMath::NormQuantile(1.0-qvaluel);
       uperror = -TMath::NormQuantile(1.0-qvaluel_UP)- estimate;
       downerror = estimate - (- TMath::NormQuantile(1.0-qvaluel_DOWN));
//       cout << " Three values " << -TMath::NormQuantile(1.0-qvaluel) << " " << -TMath::NormQuantile(1.0-qvaluel_UP) << " " << -TMath::NormQuantile(1.0-qvaluel_DOWN) << endl;
       averror = 0.5*(uperror + downerror);
       cout << "ql < pu    :  z-score lower (ie. signed pull equivalent): " << fixed << setprecision(4) << estimate << " +- "  
            << averror << " sigma   (+: " << uperror << " )  (-:  " << downerror << " ) " << endl;      
   }
   else{
       cout << "Impossible ? " << endl;
   }
   
   if(pvalueu>0.5&&qvaluel>0.5){
       cout << "Both upper and lower tails exceed 50% so average the two? " << TMath::NormQuantile(1.0-pvalueu) << " " << -TMath::NormQuantile(1.0-qvaluel) 
                                                                 << " AV: " << 0.5*(TMath::NormQuantile(1.0-pvalueu) - TMath::NormQuantile(1.0-qvaluel) ) <<   endl;
                                                                 
// New options.
       double estimate1 = TMath::NormQuantile(1.0- (pvalueu -0.5*pequal)) ;
       double estimate2 = -TMath::NormQuantile(1.0- (qvaluel -0.5*pequal)) ;       
       cout << "Estimate 1 " << estimate1 << endl;
       cout << "Estimate 2 " << estimate2 << endl;
       cout << "Average    " << 0.5*(estimate1+estimate2) << endl;       
                                                                 
   }
   
   // Always use p-value.
   estimate = TMath::NormQuantile(1.0-pvalueu);
   uperror = TMath::NormQuantile(1.0-pvalueu_UP)- estimate;
   downerror = estimate - TMath::NormQuantile(1.0-pvalueu_DOWN);
   averror = 0.5*(uperror + downerror);
   cout << "GENERALLY USE:  z-score upper (ie. signed pull equivalent): " << fixed << setprecision(4) << estimate << " +- "  
            << averror << " sigma   (+: " << uperror << " )  (-:  " << downerror << " ) " << endl;  
   cout << " " << endl;
     
   
   cout << "Histogram summary " << fixed << setprecision(8) << hist->GetMean() << " +- " << hist->GetMeanError() << " RMS = "  << hist->GetRMS() << endl;

//   hist->Draw();
   hist->Write();
   f->Write();
   f->Close(); 
   
}

int main(int argc, char **argv) {

    CLI::App app{"Evaluate Poisson z-score"};  
    
    int nfactor = 1;
    app.add_option("-n,--nfactor", nfactor, "Million toy multiplier");    

    unsigned long int seed = 123456L;
    app.add_option("-s,--seed", seed, "Seed");
    
    int nobs = 787;
    app.add_option("-o,--nobs", nobs, "Observed event count");    
    
    double mub = 904.684;    // mu+ lifetime
    app.add_option("-b,--mub", mub, "Background mean event count");
    
    double errfrac = 0.03627;  // Use CMS flux ratio value of 1.2766 // See CMS-PAS-MUO-10-001    
    app.add_option("-e, --errfrac", errfrac, "Fractional uncertainty on background");     
      
    CLI11_PARSE(app, argc, argv);

    std::cout << "nfactor  " << nfactor << std::endl;
    std::cout << "nobs     " << nobs << std::endl;
    std::cout << "mub      " << mub << std::endl;
    std::cout << "errfrac  " << errfrac << std::endl;
    std::cout << "seed     " << seed << std::endl;  
                       
    ToyChucker(nfactor, nobs, mub, errfrac, seed );                      
       
    return 0;
    
}

