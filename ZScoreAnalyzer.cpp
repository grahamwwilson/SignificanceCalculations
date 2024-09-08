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
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include "AnalysisBin.h"
#include "statAnalysisBin.h"
using namespace std;

int NTOSAVE {};
std::vector<std::tuple<int, unsigned long int, int, int>> vtoys;  // Design with experiment-number, seed, bin-id, nobs. Can then sort by experiment number.
std::vector<std::tuple<int, unsigned long int, int, int>> gtoys;  // Design with experiment-number, seed, bin-id, nobs. Can then sort by experiment number.

#include "ToyChucker.cpp" 

std::vector<AnalysisBin> vAnaBins;
std::vector<statAnalysisBin> vstatAnaBins;

double pBi(double non, double mub, double syst){
// Based on 
// @article{Cousins:2008zz,
//    author = "Cousins, Robert D. and Linnemann, James T. and Tucker, Jordan",
//    title = "{Evaluation of three methods for calculating 
// statistical significance when incorporating a systematic uncertainty 
// into a test of the background-only hypothesis for a Poisson process}",
//    eprint = "physics/0702156",
//    archivePrefix = "arXiv",
//    doi = "10.1016/j.nima.2008.07.086",
//    journal = "Nucl. Instrum. Meth. A",
//    volume = "595",
//    number = "2",
//    pages = "480--501",
//    year = "2008"
// }

// Use Gaussian-mean background problem equivalence
// non  = Observed events
// mub  = Background mean
// syst = Fractional systematic uncertainty on background mean
   double sb = syst*mub;
   double tau = mub/(sb*sb);
   double noff = tau*mub;

   cout << "Calculating pBi for non = " << non << " mub= " << mub 
                                        << " syst = " << syst << endl;

   double p_Bi = TMath::BetaIncomplete(1.0/(1.0+tau), non, noff+1.0);

   return p_Bi;
}

double zSingleSided(double pvalue){
   double zvalue = sqrt(2.0)*TMath::ErfInverse(1.0 - 2.0*pvalue);
   return zvalue;
}


double zBi(int nobs, double mub, double syst){
   double non = double(nobs);
   double pvalue = pBi(non, mub, syst);
   double z = zSingleSided(pvalue);
   return z;
}

void ZScoreAnalyzer(bool genFixed, int afactor, int nfactor, bool PoissonOnly, int ndivide, unsigned long int baseseed, int NOBS, double MUB, double FRACERR){
    
    TFile *f = new TFile("ZScoreAnalyzer.root","RECREATE");
    TH1D* hist = new TH1D("hist","Toy; Z-Score (Algorithm 1); Bins per 0.25",49,-6.125,6.125);    
    TH1D* histz = new TH1D("histz","Toy; Z-Score (Algorithm 0MH); Bins per 0.25",49,-6.125,6.125); 
    TH1D* histzu = new TH1D("histzu","Toy; Z-Score (Algorithm 0P); Bins per 0.25",49,-6.125,6.125); 
    TH1D* histzl = new TH1D("histzl","Toy; Z-Score (Algorithm 0N); Bins per 0.25",49,-6.125,6.125);
    
    TH1D* hist_censored = new TH1D("hist_censored","Toy; Z-Score (Algorithm 1); Bins per 0.25",49,-6.125,6.125);    
    TH1D* histz_censored = new TH1D("histz_censored","Toy; Z-Score (Algorithm 0MH); Bins per 0.25",49,-6.125,6.125); 
    
    TH1D* nhist  = new TH1D("nhist","Toy; Z-Score (Algorithm 1); Bins per 0.25",40,-5.0,5.0); 
    TH1D* nhistp  = new TH1D("nhistp","Toy; Z-Score (Algorithm 1); Bins per 0.25",40,-5.0,5.0); 
    TH1D* nhistm  = new TH1D("nhistm","Toy; Z-Score (Algorithm 1); Bins per 0.25",40,-5.0,5.0);            
    TH1D* nhistz = new TH1D("nhistz","Toy; Z-Score (Algorithm 0MH); Bins per 0.25",40,-5.0,5.0);
    TH1D* nhistzu = new TH1D("nhistzu","Toy; Z-Score (ZN Upper); Bins per 0.25",40,-5.0,5.0); 
    TH1D* nhistzl = new TH1D("nhistzl","Toy; Z-Score (ZN Lower); Bins per 0.25",40,-5.0,5.0);
    TH1D* nhistzPull = new TH1D("nhistzPull","Toy; Z-Score (Raw Pull); Bins per 0.25",40,-5.0,5.0); 
    
    TH1D* mhist  = new TH1D("mhist","Toy; Z-Score (Algorithm 1); Bins per 0.025",400,-5.0,5.0);
    TH1D* mhistp  = new TH1D("mhistp","Toy; Z-Score (Algorithm 1); Bins per 0.025",400,-5.0,5.0);
    TH1D* mhistm  = new TH1D("mhistm","Toy; Z-Score (Algorithm 1); Bins per 0.025",400,-5.0,5.0);            
    TH1D* mhistz = new TH1D("mhistz","Toy; Z-Score (Algorithm 0MH); Bins per 0.025",400,-5.0,5.0);
    TH1D* mhistzu = new TH1D("mhistzu","Toy; Z-Score (ZN Upper); Bins per 0.025",400,-5.0,5.0); 
    TH1D* mhistzl = new TH1D("mhistzl","Toy; Z-Score (ZN Lower); Bins per 0.025",400,-5.0,5.0);
    TH1D* mhistzPull = new TH1D("mhistzPull","Toy; Z-Score (Raw Pull); Bins per 0.025",400,-5.0,5.0);   
    
    TH1D* phist  = new TH1D("phist","Toy; Z-Score (Algorithm 1); Bins per 0.08",125,-5.0,5.0); 
    TH1D* phistp  = new TH1D("phistp","Toy; Z-Score (Algorithm 1); Bins per 0.08",125,-5.0,5.0); 
    TH1D* phistm  = new TH1D("phistm","Toy; Z-Score (Algorithm 1); Bins per 0.08",125,-5.0,5.0);            
    TH1D* phistz = new TH1D("phistz","Toy; Z-Score (Algorithm 0MH); Bins per 0.08",125,-5.0,5.0);
    TH1D* phistz1 = new TH1D("phistz1","Toy; Z-Score (Algorithm 0MH); Bins per 0.08",125,-5.0,5.0);
    TH1D* phistz2 = new TH1D("phistz2","Toy; Z-Score (Algorithm 0MH); Bins per 0.08",125,-5.0,5.0);
    TH1D* phistz3 = new TH1D("phistz3","Toy; Z-Score (Algorithm 0MH); Bins per 0.08",125,-5.0,5.0);
    TH1D* phistz4 = new TH1D("phistz4","Toy; Z-Score (Algorithm 0MH); Bins per 0.08",125,-5.0,5.0);            
    TH1D* phistzu = new TH1D("phistzu","Toy; Z-Score (ZN Upper); Bins per 0.08",125,-5.0,5.0); 
    TH1D* phistzl = new TH1D("phistzl","Toy; Z-Score (ZN Lower); Bins per 0.08",125,-5.0,5.0);
    TH1D* phistzPull = new TH1D("phistzPull","Toy; Z-Score (Raw Pull); Bins per 0.08",125,-5.0,5.0);                    
    
    TH1D* hpval = new TH1D("hpval","Toy; p-value; Bins per 0.05",20,0.0,1.00); 
    TH1D* hqval = new TH1D("hqval","Toy; q-value; Bins per 0.05",20,0.0,1.00);
    TH1D* heval = new TH1D("heval","Toy; e-value; Bins per 0.05",20,0.0,1.00);
    TH1D* hmval = new TH1D("hmval","Toy; m-value; Bins per 0.05",20,0.0,1.00); 
    
    TH1D* hpval2 = new TH1D("hpval2","Toy; p-value; Bins per 0.01",100,0.0,1.00); 
    TH1D* hqval2 = new TH1D("hqval2","Toy; q-value; Bins per 0.01",100,0.0,1.00);
    TH1D* heval2 = new TH1D("heval2","Toy; e-value; Bins per 0.01",100,0.0,1.00);
    TH1D* hmval2 = new TH1D("hmval2","Toy; m-value; Bins per 0.01",100,0.0,1.00); 
    
    TH1D* hpemean = new TH1D("hpemean","Toy; mean Z-score; Bins per 0.01",40,-0.10,0.10);
    TH1D* hpesd   = new TH1D("hpesd","Toy; sd of Z-score; Bins per 0.005",40,0.9,1.1);    
    
// Recapitulate input parameters
    cout << "Inputs: NOBS, MUB, FRACERR " << NOBS << " " << MUB << " " << FRACERR << std::endl;

    cout << "start toys ... " << endl;
    
    int nbin = 0;
    int id = -1;
    unsigned long int seed = baseseed;

         nbin +=1;
         auto ndata = NOBS;
         auto postfitMean = MUB;
         auto postfitError = MUB*FRACERR;
        
         double satChisq;  double n = double(ndata); double y = postfitMean;
      // Use Baker and Cousins Poisson formula.  NIM 221 (1984) 437.
      // This ignores completely the postfit uncertainty.
         if(NOBS == 0){
             satChisq = 2.0*(y - n);
         }
         else{
             satChisq = 2.0*(y - n + n*log(n/y) );
         }
        
         
         double fracError = postfitError/postfitMean;
         
         cout << " " << endl;
         
         // Estimate number of standard deviations.
         double stdev = sqrt(satChisq);
         
         cout << "Bin " << nbin << " " << stdev*(ndata-postfitMean)/std::abs(ndata-postfitMean) << endl;         
        
// Standard factor is 100         
         int ntouse = nfactor;
         if (stdev > 2.5) ntouse =   1000/ndivide;         
         if (stdev > 3.0) ntouse =   2500/ndivide;
         if (stdev > 3.5) ntouse =  10000/ndivide;
         if (stdev > 4.0) ntouse =  40000/ndivide;
         if (stdev > 4.5) ntouse = 100000/ndivide;
         if (stdev > 5.0) ntouse = 200000/ndivide;         
         ntouse = ntouse*nfactor;
         
         if(genFixed)ntouse = afactor;
         
 // Throw toys if a good use of our time.
         double znew, zupper, zlower, zold, zscoreError;
         double pval,qval,eval,mval;
         double zpull;
         
         if(stdev > 9999.0){                    // never should happen
             if(ndata > postfitMean){
                znew = 5.01;
                zupper = 5.01;
                zlower = 5.01;
                zold = 5.01;
                zscoreError = -1.0;
             }
             else{
                znew = -5.01;
                zupper = -5.01;
                zlower = -5.01;
                zold = -5.01;      
                zscoreError = -1.0;         
             }
             // Fill these defaults for the two atrocious bins
             histz->Fill(znew); 
             hist->Fill(zold);
             nhistz->Fill(znew); 
             nhist->Fill(zold);             
             histzl->Fill(zlower);
             histzu->Fill(zupper);              
             // Even in this case we want to throw some toys (at least NTOSAVE) but not save the z-scores.

             int ntogen = NTOSAVE;
             int NTOUSE = 1;
             if(NTOSAVE>1000)NTOUSE=1 + (NTOSAVE/1000); 
             cout << "Defaulting to " << NTOUSE << " thousand toys for pathological bins " << endl;

             if(PoissonOnly)fracError=1.0e-6;
             std::tuple<int, double, unsigned int, double, double, double, double, double, double, double, double, double> t = 
                       ToyChucker(NTOUSE, ndata, ndata, postfitMean, fracError, seed, id );             
             
         }
         else{
          
             if(PoissonOnly)fracError=1.0e-6;
          
             std::tuple<int, double, unsigned int, double, double, double, double, double, double, double, double, double> t = 
                       ToyChucker(ntouse, ndata, ndata, postfitMean, fracError, seed, id );
         
             cout << "ZUN for bin " << nbin << " " << id << " " << std::get<6>(t)  << " +- " << std::get<7>(t) << endl;         
             cout << "ZLN for bin " << nbin << " " << id << " " << std::get<8>(t)  << " +- " << std::get<9>(t) << endl;         
             cout << "ZSN for bin " << nbin << " " << id << " " << std::get<10>(t)  << " +- " << std::get<11>(t) << endl;
             cout << "ZS  for bin " << nbin << " " << id << " " << (double(ndata)-postfitMean)/sqrt(postfitMean) << endl;
             cout << "ZBI for bin " << nbin << " " << id << " " << zBi(ndata, postfitMean, FRACERR) << endl;             
         
             pval = std::get<3>(t);
             qval = std::get<4>(t);
             eval = std::get<5>(t);
             mval = pval -0.5*eval;
         
             zpull = (double(ndata) - postfitMean)/sqrt(postfitMean);         
             znew = std::get<10>(t);
             zscoreError = std::get<11>(t);
             zupper = std::get<6>(t);
             zlower = std::get<8>(t);
             zold = zupper;
             if(ndata < postfitMean)zold = zlower;
             histz->Fill(znew); 
             hist->Fill(zold);

             nhistz->Fill(znew); 
             nhist->Fill(zold);
             nhistzu->Fill(zupper);
             nhistzl->Fill(zlower);
             nhistzPull->Fill(zpull);                          

             mhistz->Fill(znew); 
             mhist->Fill(zold);
             mhistzu->Fill(zupper);
             mhistzl->Fill(zlower);
             mhistzPull->Fill(zpull);
             
             if(postfitMean < 10.46){
                 phistz1->Fill(znew);
             }
             else if(postfitMean < 32.95){
                 phistz2->Fill(znew);        
             } 
             else if(postfitMean < 133.0){
                 phistz3->Fill(znew);        
             }
             else{
                 phistz4->Fill(znew);         
             }                 
             
             phistz->Fill(znew); 
             phist->Fill(zold);
             phistzu->Fill(zupper);
             phistzl->Fill(zlower);
             phistzPull->Fill(zpull); 
             
             if (double(ndata) >= postfitMean){
                 nhistp->Fill(zold);    
                 mhistp->Fill(zold);
                 phistp->Fill(zold);                              
             }
             else{
                 nhistm->Fill(zold);    
                 mhistm->Fill(zold);
                 phistm->Fill(zold);                         
             }             

             histzl->Fill(zlower);
             histzu->Fill(zupper);
             histz_censored->Fill(znew); 
             hist_censored->Fill(zold); 
             double EPS=1.0e-14;
//             hpval->Fill(min(pval, 1.0-EPS) );
             hpval->Fill(pval);
             hqval->Fill(qval);
             heval->Fill(eval);
             hmval->Fill(mval);
//             hpval2->Fill(min(pval, 1.0-EPS) );
             hpval2->Fill(pval);
             hqval2->Fill(qval);
             heval2->Fill(eval);
             hmval2->Fill(mval);       
             
/*
             if(pval > 0.999){
                 cout << "High p-value of " << pval << " for binid " << id << " " << ndata << " " << regionName << endl;
             }
             if(pval < 0.001){
                 cout << "Low p-value of " << pval << " for binid " << id << " " << ndata << " " << regionName << endl;
             }             
*/             
                              
         }
          
//    }
// And now sort again and print.
    
// Stick sample mean and rms in pseudo-experiment histogram                                        
    
    hpemean->Fill(phistz->GetMean());
    hpesd->Fill(phistz->GetRMS());    
    
    f->Write();
    f->Close();         

}

int main(int argc, char** argv){

    CLI::App app{"Basic Z-score analysis program"};
    
    int NOBS = 20;
    app.add_option("-n,--nobs", NOBS, "Observed data");   
    
    double MUB = 10.2;
    app.add_option("-m,--mub", MUB, "Background mean");
    
    double FRACERR = 0.035;
    app.add_option("-b,--bfrac", FRACERR, "Background mean fractional error");
    
    int nfactor = 1;
    app.add_option("-f,--nfactor", nfactor, "Thousand toy multiplier"); 
    
    bool poissonOnly = false;
    app.add_option("-p,--poissonOnly", poissonOnly, "Boolean for Poisson Only");
    
    bool genFixed = true;
    app.add_option("-g,--genFixed", genFixed, "Boolean for Fixed Number to Generate");    
    
    int afactor = 4000000;
    app.add_option("-a,--afactor", afactor, "Thousand toy fixed multiplier");    
    
    int ndivide = 1;
    app.add_option("-d,--ndivide", ndivide, "Divisor for quick test");
    
    unsigned long int baseseed = 133456L;
//    unsigned long int baseseed = 125900L;    
    app.add_option("-s,--baseseed", baseseed, "Base seed for toys");    
              
    CLI11_PARSE(app, argc, argv);
    
    std::cout << "NOBS        " << NOBS << std::endl;
    std::cout << "MUB         " << MUB << std::endl;
    std::cout << "FRACERR     " << FRACERR << std::endl;
    std::cout << "nfactor     " << nfactor << std::endl;    
    std::cout << "poissonOnly " << poissonOnly << endl;
    std::cout << "ndivide     " << ndivide << endl;
    std::cout << "baseseed    " << baseseed << endl; 
    std::cout << "afactor     " << afactor << endl;
    std::cout << "genFixed    " << genFixed << endl;
    
    baseseed += NOBS;
    
    std::cout << "Used seed = " << baseseed << std::endl;
    
    ZScoreAnalyzer(genFixed, afactor, nfactor, poissonOnly, ndivide, baseseed, NOBS, MUB, FRACERR);
       
    return 0;
    
}
