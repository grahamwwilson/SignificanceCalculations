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

   const unsigned int N = 2443;    // Current number of bins
   const unsigned int NGENERATED = 1000000;

   TRandom3 *rg = new TRandom3(seed);
   TFile *f = new TFile("tailtest.root","RECREATE");
   TH1D* hist = new TH1D("hist","hist",12,0.0,6.0);
   TH1D* hpvalue = new TH1D("hpvalue","hpvalue",100,0.0,1.0);
   TH1D* hx0 = new TH1D("hx0","hx0",1299,-0.5,1298.5);
   TH1D* hx1 = new TH1D("hx1","hx1",1299,-0.5,1298.5);
   TH1D* hx2 = new TH1D("hx2","hx2",1299,-0.5,1298.5);
   TH1D* hx3 = new TH1D("hx3","hx3",1299,-0.5,1298.5);
   TH1D* hx4 = new TH1D("hx4","hx4",1299,-0.5,1298.5);   
   TH1D* hx5 = new TH1D("hx5","hx5",1299,-0.5,1298.5);
   TH1D* hx6 = new TH1D("hx6","hx6",1299,-0.5,1298.5);
   TH1D* hx7 = new TH1D("hx7","hx7",1299,-0.5,1298.5);
   TH1D* hx8 = new TH1D("hx8","hx8",1299,-0.5,1298.5);      
   TH1D* hx9 = new TH1D("hx9","hx9",1299,-0.5,1298.5);
   TH1D* hmax = new TH1D("hmax","hmax",600,0.0,6.0);                             
   double zi;
   
   const int NZCUT = 10;
   double nobs[NZCUT] = {  0,   0,   0,   0,   0,   0,   0,   0,   0,   0};
   double zcut[NZCUT] = {1.0, 2.0, 2.5, 2.75, 3.0, 3.25, 3.5, 4.0, 4.5, 5.0};      

// Simulate many Gaussians
   for (int i=0; i<NGENERATED; i++){
      double chisq = 0.0;
      for (int k=0; k<NZCUT; k++){
          nobs[k] = 0;
      }
      double zmax = -1.0;
      for (int i=0;i<N;i++){
          zi = rg->Gaus(0.0,1.0);    // standardized normal random variates
          if(abs(zi) > zmax) zmax = abs(zi);
          for (int j=0; j<NZCUT; j++){
              if(abs(zi) > zcut[j]) nobs[j]++;
          }
      }
      hx0->Fill(nobs[0]);
      hx1->Fill(nobs[1]);
      hx2->Fill(nobs[2]);
      hx3->Fill(nobs[3]);
      hx4->Fill(nobs[4]);
      
      hx5->Fill(nobs[5]);
      hx6->Fill(nobs[6]);
      hx7->Fill(nobs[7]);
      hx8->Fill(nobs[8]);
      hx9->Fill(nobs[9]);
      
      hmax->Fill(zmax);                        

   }
   f->Write();

}
