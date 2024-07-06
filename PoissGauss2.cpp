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
#include <vector>
#include <tuple>
using namespace std;

void GeneralToyChucker(int NFACTOR, int NOBSMIN, int NOBSMAX, double MUB, double FRACERROR, unsigned long int seed ){

   int THOUSAND = 1000;
   int NGENERATED = 1*THOUSAND;   // Set so that pvalue MC fractional uncertainty is better than 2.5% (ie 1600 successes) for even the worst bins.
   
   NGENERATED = NFACTOR*THOUSAND;
   
   cout << "Generating " << NFACTOR << " thousand toys" << endl;
   cout << "Using seed " << seed << endl;
   cout << "MUB  = "      << MUB << endl;
   cout << "NOBSMIN = "      << NOBSMIN << endl;
   cout << "NOBSMAX = "      << NOBSMAX << endl;   
   cout << "FRACERROR = " << FRACERROR << endl;
   double SIGMAB = FRACERROR*MUB;
   cout << "SIGMAB " << SIGMAB << endl;   
   
// Note depending on the actual p-value one may need a lot of toys, especially for very small p-values
// Plan for at least 1600 events in the tails corresponding to 2.5% errors.  

   int NBINS = NOBSMAX - NOBSMIN + 1;

   TRandom3 *rg = new TRandom3(seed);
   TFile *f = new TFile("out.root","RECREATE");
   TH1D* hist = new TH1D("hist","hist",2000,-0.5,1999.5);
   TH1D* hist0 = new TH1D("hist0","Toys; NOBS; Z-Score (Algorithm 0)",NBINS,double(NOBSMIN)-0.5,double(NOBSMAX)+0.5); 
   TH1D* hist0N = new TH1D("hist0N","Toys; NOBS; Z-Score (Algorithm 0N)",NBINS,double(NOBSMIN)-0.5,double(NOBSMAX)+0.5); 
   TH1D* hist0M = new TH1D("hist0M","Toys; NOBS; Z-Score (Algorithm 0M)",NBINS,double(NOBSMIN)-0.5,double(NOBSMAX)+0.5);  
   TH1D* hist0MH = new TH1D("hist0MH","Toys; NOBS; Z-Score (Algorithm 0MH)",NBINS,double(NOBSMIN)-0.5,double(NOBSMAX)+0.5);          
   TH1D* hist1 = new TH1D("hist1","Toys; NOBS; Z-Score (Algorithm 1)",NBINS,double(NOBSMIN)-0.5,double(NOBSMAX)+0.5); 
   TH1D* hist2 = new TH1D("hist2","Toys; NOBS; Z-Score (Algorithm 2)",NBINS,double(NOBSMIN)-0.5,double(NOBSMAX)+0.5);
   TH1D* hist3 = new TH1D("hist3","Toys; NOBS; Z-Score (Algorithm 3)",NBINS,double(NOBSMIN)-0.5,double(NOBSMAX)+0.5); 
   TH1D* hist4 = new TH1D("hist4","Toys; NOBS; Z-Score (Algorithm 4)",NBINS,double(NOBSMIN)-0.5,double(NOBSMAX)+0.5);    
   
   TH1D* histz0 = new TH1D("histz0","Toys; Z-Score (Algorithm 0); Toys per 0.25",49,-6.125,6.125); 
   TH1D* histz0N = new TH1D("histz0N","Toys; Z-Score (Algorithm 0N); Toys per 0.25",49,-6.125,6.125); 
   TH1D* histz0M = new TH1D("histz0M","Toys; Z-Score (Algorithm 0M); Toys per 0.25",49,-6.125,6.125); 
   TH1D* histz0MH = new TH1D("histz0MH","Toys; Z-Score (Algorithm 0MH); Toys per 0.25",49,-6.125,6.125);        
   TH1D* histz1 = new TH1D("histz1","Toys; Z-Score (Algorithm 1); Toys per 0.25",49,-6.125,6.125); 
   TH1D* histz2 = new TH1D("histz2","Toys; Z-Score (Algorithm 2); Toys per 0.25",49,-6.125,6.125);  
   TH1D* histz3 = new TH1D("histz3","Toys; Z-Score (Algorithm 3); Toys per 0.25",49,-6.125,6.125);
   TH1D* histz4 = new TH1D("histz4","Toys; Z-Score (Algorithm 4); Toys per 0.25",49,-6.125,6.125);   
   
   TH1D* histzz0 = new TH1D("histzz0","Toys; Z-Score (Algorithm 0); Toys per 0.01 bin",1201,-6.005,6.005); 
   TH1D* histzz0N = new TH1D("histzz0N","Toys; Z-Score (Algorithm 0N); Toys per 0.01 bin",1201,-6.005,6.005);  
   TH1D* histzz0M = new TH1D("histzz0M","Toys; Z-Score (Algorithm 0M); Toys per 0.01 bin",1201,-6.005,6.005);  
   TH1D* histzz0MH = new TH1D("histzz0MH","Toys; Z-Score (Algorithm 0MH); Toys per 0.01 bin",1201,-6.005,6.005);           
   TH1D* histzz1 = new TH1D("histzz1","Toys; Z-Score (Algorithm 1); Toys per 0.01 bin",1201,-6.005,6.005); 
   TH1D* histzz2 = new TH1D("histzz2","Toys; Z-Score (Algorithm 2); Toys per 0.01 bin",1201,-6.005,6.005);  
   TH1D* histzz3 = new TH1D("histzz3","Toys; Z-Score (Algorithm 3); Toys per 0.01 bin",1201,-6.005,6.005); 
   TH1D* histzz4 = new TH1D("histzz4","Toys; Z-Score (Algorithm 4); Toys per 0.01 bin",1201,-6.005,6.005); 
   
   TH1D* histp = new TH1D("histp","Toys; Upper-tail p-value; Toys per 0.01",100,0.0,1.0);  
   TH1D* histq = new TH1D("histq","Toys; Lower-tail q-value; Toys per 0.01",100,0.0,1.0); 
   TH1D* histh = new TH1D("histh","Toys; Upper-tail p-value - 0.5*pequal; Toys per 0.01",100,0.0,1.0); 
   
   TH1D* histcp = new TH1D("histcp","Toys; Upper-tail p-value; Toys per 0.1",10,0.0,1.0);  
   TH1D* histcq = new TH1D("histcq","Toys; Lower-tail q-value; Toys per 0.1",10,0.0,1.0); 
   TH1D* histch = new TH1D("histch","Toys; Upper-tail p-value - 0.5*pequal; Toys per 0.1",10,0.0,1.0);                              
   
   std::vector<int> vupper, vlower, vequal;
   for (int i=NOBSMIN; i<= NOBSMAX; i++){
       vupper.push_back(0);
       vlower.push_back(0);
       vequal.push_back(0);
   }
   
   std::vector<int> vgen;
   
   int nlow = 0;
   int nhigh = 0;

// Simulate a Poisson distribution with the underlying mean being a Gaussianly distributed random number.
   for (int i=0; i<NGENERATED; i++){
        double mu = MUB;
        if(FRACERROR > 1.0e-4){
            mu = rg->Gaus(MUB, SIGMAB);  // Choose Poisson mean, mu, from Gaussian with mean and rms of MUB and SIGMAB
        }
        int n = rg->Poisson(mu);         // Generate Poisson distributed random number, n, based on Poisson mean of mu  
        for (int NOBS = NOBSMIN; NOBS <=NOBSMAX; NOBS++) {
             int j = NOBS-NOBSMIN;
             if(n >= NOBS)vupper[j] +=1;       // Count toys with n exceeding or equal to the observed counts
             if(n <= NOBS)vlower[j] +=1;       // Count toys with n less than or equalt to the observed counts
             if(n==NOBS)vequal[j] +=1;
             if(n<NOBSMIN)nlow++;
             if(n>NOBSMAX)nhigh++;
        }     
        hist->Fill(double(n));
        vgen.push_back(n);
   }
   
   cout << "nlow = " << nlow << "  nhigh " << nhigh << endl; 
   
// Create summary statistics post-generation
   std::vector<std::pair<double,double>> vpvalue, vqvalue, vpequal;
   std::vector<std::tuple<int, double, int, double, double, double, double, double, double, double, double, double >> vtup;
   
   for (int NOBS = NOBSMIN; NOBS <=NOBSMAX; NOBS++) {
        int j = NOBS-NOBSMIN;
        
        double pvalue = double(vupper[j])/double(NGENERATED);
        double dp=sqrt(pvalue*(1.0-pvalue)/double(NGENERATED));  //binomial error        
        std::pair<double,double> p = std::make_pair( pvalue, dp );
        vpvalue.push_back(p);
        
        double qvalue = double(vlower[j])/double(NGENERATED);
        double dq=sqrt(qvalue*(1.0-qvalue)/double(NGENERATED));  //binomial error        
        std::pair<double,double> q = std::make_pair( qvalue, dq );
        vqvalue.push_back(q);        
        
        double pequal = double(vequal[j])/double(NGENERATED);
        double de=sqrt(pequal*(1.0-pequal)/double(NGENERATED));  //binomial error        
        std::pair<double,double> e = std::make_pair( pequal, de );
        vpequal.push_back(e);        
        
        cout << " " << endl;
        cout << "NOBS " << NOBS << endl;
        cout << "Upper-tail " << vpvalue[j].first << " +- " << vpvalue[j].second << endl;
        cout << "Lower-tail " << vqvalue[j].first << " +- " << vqvalue[j].second << endl;
        cout << "Equal      " << vpequal[j].first << " +- " << vpequal[j].second << endl; 
        
// Rather than apply the validity conditions of the estimators, let's first calculate them all 
// and defer figuring out the former.
        double zscoreu      =  TMath::NormQuantile(1.0 - vpvalue[j].first);
        double zscoreu_UP   =  TMath::NormQuantile(1.0 - (vpvalue[j].first - vpvalue[j].second));
        double zscoreu_DOWN =  TMath::NormQuantile(1.0 - (vpvalue[j].first + vpvalue[j].second));
        double zscoreu_ERR = 0.5*(zscoreu_UP - zscoreu_DOWN); 
        cout << "zscoreu: " << zscoreu << " " << zscoreu_UP << " " << zscoreu_DOWN << " " << zscoreu_ERR << endl;
               
        double zscorel      = -TMath::NormQuantile(1.0 - vqvalue[j].first);
        double zscorel_UP   = -TMath::NormQuantile(1.0 - (vqvalue[j].first + vqvalue[j].second));
        double zscorel_DOWN = -TMath::NormQuantile(1.0 - (vqvalue[j].first - vqvalue[j].second));
        double zscorel_ERR = 0.5*(zscorel_UP - zscorel_DOWN);
        cout << "zscorel: " << zscorel << " " << zscorel_UP << " " << zscorel_DOWN << " " << zscorel_ERR << endl; 
        
        double zscorem1  =  TMath::NormQuantile(1.0 - ( vpvalue[j].first -0.5*vpequal[j].first ) ); 
        double zscorem2  =  -TMath::NormQuantile(1.0 - ( vqvalue[j].first -0.5*vpequal[j].first ) );                      
                
        cout << "Z-scores " << zscoreu << " +- " << zscoreu_ERR << " " << zscorel << " +- " << zscorel_ERR << " " << 0.5*(zscoreu + zscorel) << " " << zscorem1 << " " << zscorem2 << endl;
        
        std::tuple<int, double, int, double, double, double, double, double, double, double, double, double> t = std::make_tuple(NOBS, MUB, NGENERATED, 
                                     vpvalue[j].first, vqvalue[j].first, vpequal[j].first, zscoreu, zscoreu_ERR, zscorel, zscorel_ERR, zscorem1, 0.5*(zscorel_ERR + zscoreu_ERR));
        vtup.push_back(t);
     
   }      

// Apply the various algorithms and fill histograms with the choices.
   
   for (auto & el : vtup){
        auto NOBS = std::get<0>(el);
        auto MUB  = std::get<1>(el);
        auto p = std::get<3>(el);
        auto q = std::get<4>(el);
        auto e = std::get<5>(el);
        auto zscoreu = std::get<6>(el);
        auto zscoreu_ERR = std::get<7>(el);
        auto zscorel = std::get<8>(el);
        auto zscorel_ERR = std::get<9>(el);
        auto zscorem = std::get<10>(el);
        auto zscorem_ERR = std::get<11>(el);
        int ibin = NOBS - NOBSMIN + 1;
        double z1,z2,z3,z4;
        double z0 = zscoreu;
        
        hist0->Fill(NOBS,zscoreu);
        hist0->SetBinError(ibin, zscoreu_ERR); 
        hist0N->Fill(NOBS,zscorel);
        hist0N->SetBinError(ibin, zscorel_ERR);
        hist0M->Fill(NOBS,0.5*(zscorel + zscoreu));
        hist0M->SetBinError(ibin, 0.5*(zscorel_ERR + zscoreu_ERR));
        hist0MH->Fill(NOBS,zscorem);
        hist0MH->SetBinError(ibin, zscorem_ERR);                                
        
// Algorithm 1
        if(NOBS >= MUB){
            z1 = zscoreu;
            hist1->Fill(NOBS,zscoreu);
            hist1->SetBinError(ibin, zscoreu_ERR);
        }
        else{
            z1 = zscorel;
            hist1->Fill(NOBS,zscorel);
            hist1->SetBinError(ibin, zscorel_ERR);
        }
// Algorithm 2
        if(p <= q){
            z2 = zscoreu;
            hist2->Fill(NOBS,zscoreu);
            hist2->SetBinError(ibin, zscoreu_ERR);            
        }
        else{
            z2 = zscorel;
            hist2->Fill(NOBS,zscorel);
            hist2->SetBinError(ibin, zscorel_ERR); 
        }
// Algorithm 3/4
        if(p <= q){     // Upper-tail probability is smallest
             if(p<0.5){
                 z3 = zscoreu;
//                 z4 = zscoreu;
                 hist3->Fill(NOBS, zscoreu);
                 hist3->SetBinError(ibin, zscoreu_ERR);
//                 hist4->Fill(NOBS, zscoreu);
//                 hist4->SetBinError(ibin, zscoreu_ERR);                   
             }
             else{       // Both p and q are > 0.5
                 z3 = 0.5*(zscoreu+zscorel);
                 hist3->Fill(NOBS,0.5*(zscoreu + zscorel));
                 hist3->SetBinError(ibin, 0.5*(zscorel_ERR + zscoreu_ERR)); 
             }
        }
        else{
             if(q<0.5){
                 z3 = zscorel;
                 hist3->Fill(NOBS,zscorel);
                 hist3->SetBinError(ibin, zscorel_ERR);                  
             }
             else{
                 z3 = 0.5*(zscoreu + zscorel);
                 hist3->Fill(NOBS,0.5*(zscoreu + zscorel));
                 hist3->SetBinError(ibin, 0.5*(zscorel_ERR + zscoreu_ERR));                  
             }               
        }
        cout << "NOBS = " << NOBS << " Alg 0,1,2,3 z-scores " << z0 << " " << z1 << " " << z2 << " " << z3 << endl;              
   }
   
// Make distributions based on the generated sample
   for (auto & elg : vgen){
       for (auto & el : vtup){
           auto NOBS = std::get<0>(el);
           if (NOBS == elg) {                    // Corresponds to a generated point
               auto MUB  = std::get<1>(el);
               auto p = std::get<3>(el);
               auto q = std::get<4>(el);
               auto e = std::get<5>(el);
               auto zscoreu = std::get<6>(el);
               auto zscorel = std::get<8>(el);
               auto zscorem = std::get<10>(el); 
               
               histp->Fill(p);
               histq->Fill(q);
               histh->Fill(p-0.5*e);
               histcp->Fill(p);
               histcq->Fill(q);
               histch->Fill(p-0.5*e);               
                             
               histz0->Fill(zscoreu);
               histzz0->Fill(zscoreu);
               
               histz0N->Fill(zscorel);
               histzz0N->Fill(zscorel);   
               
               histz0M->Fill(0.5*(zscoreu+zscorel));
               histzz0M->Fill(0.5*(zscoreu+zscorel)); 
               
               histz0MH->Fill(zscorem);
               histzz0MH->Fill(zscorem);                                            
   
// Algorithm 1
               if(NOBS >= MUB){
                   histz1->Fill(zscoreu);
                   histzz1->Fill(zscoreu);                   
               }
               else{
                   histz1->Fill(zscorel);
                   histzz1->Fill(zscorel);                   
               }
// Algorithm 2
               if(p <= q){
                   histz2->Fill(zscoreu);
                   histzz2->Fill(zscoreu);                   
               }
               else{
                   histz2->Fill(zscorel);
                   histzz2->Fill(zscorel);                   
               }
// Algorithm 3
               if(p <= q){     // Upper-tail probability is smallest
                   if(p<0.5){
                       histz3->Fill(zscoreu);
                       histzz3->Fill(zscoreu);                       
                   }
                   else{       // Both p and q are > 0.5
                       histz3->Fill(0.5*(zscoreu + zscorel));
                       histzz3->Fill(0.5*(zscoreu + zscorel));                       
                   }
               }
               else{
                   if(q<0.5){
                       histz3->Fill(zscorel);
                       histzz3->Fill(zscorel);                       
                   }
                   else{
                       histz3->Fill(0.5*(zscoreu + zscorel));
                       histzz3->Fill(0.5*(zscoreu + zscorel));                       
                   }               
               }
               
           }
       }        
   }

   cout << "Histogram summary " << fixed << setprecision(8) << hist->GetMean() << " +- " << hist->GetMeanError() << " RMS = "  << hist->GetRMS() << endl;

//   hist->Draw();
//   hist->Write();
   f->Write();
   f->Close(); 
   
}

int main(int argc, char **argv) {

    CLI::App app{"Evaluate Poisson Z-score"};  
    
    int nfactor = 1000;
    app.add_option("-n,--nfactor", nfactor, "Thousand toy multiplier");    

    unsigned long int seed = 123456L;
    app.add_option("-s,--seed", seed, "Seed");
    
    int nobsmin = 5;
    app.add_option("-l,--nobsmin", nobsmin, "Minimum observed event count"); 
    
    int nobsmax = 15;
    app.add_option("-u,--nobsmax", nobsmax, "Maximum observed event count");           
    
    double mub = 10.1;  
    app.add_option("-b,--mub", mub, "Background mean event count");
    
    double errfrac = 0.035;    
    app.add_option("-e, --errfrac", errfrac, "Fractional uncertainty on background");     
      
    CLI11_PARSE(app, argc, argv);

    std::cout << "nfactor  " << nfactor << std::endl;
    std::cout << "nobsmin  " << nobsmin << std::endl;
    std::cout << "nobsmax  " << nobsmax << std::endl;    
    std::cout << "mub      " << mub << std::endl;
    std::cout << "errfrac  " << errfrac << std::endl;
    std::cout << "seed     " << seed << std::endl;  
                       
    GeneralToyChucker(nfactor, nobsmin, nobsmax, mub, errfrac, seed );                      
       
    return 0;
    
}

