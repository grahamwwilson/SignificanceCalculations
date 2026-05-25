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
#include "MyQuantile.h"
#include "ComputeZscoreLocationAndScaleCorrection.h"
using namespace std;

int NTOSAVE {};
std::vector<std::tuple<int, unsigned long int, int, int>> vtoys;  // Design with experiment-number, bin-id, nobs, seed. Can then sort by experiment number.

#include "ToyChucker.cpp" 

std::vector<AnalysisBin> vAnaBins;
std::vector<statAnalysisBin> vstatAnaBins;

bool sort_by_abszscore( const statAnalysisBin & lhs, const statAnalysisBin & rhs )
{
   return abs(lhs.zscore) > abs(rhs.zscore);  // From most to least deviating
 //  return abs(lhs.zscore) < abs(rhs.zscore);    // From least to most deviating
}


void Parser(int nfactor, bool PoissonOnly, int ndivide, bool lcorr, unsigned long int baseseed, std::string infile){

    // cout << "Hello from Parser " << endl;
    
    TFile *f = new TFile("Analysis.root","RECREATE");

// For data inspection it is important that the range of the histogram is wide enough to see far tails  
    TH1D* histz  = new TH1D("histz","Data; Z-Score (Algorithm 0MH); Bins per 0.25",48,-6.0,6.0); 
    TH1D* histzFinal  = new TH1D("histzFinal","Data; Z-Score (Latest Algorithm); Bins per 0.25",48,-6.0,6.0);
    TH1D* histzu = new TH1D("histzu","Data; Z-Score (Algorithm 0P); Bins per 0.25",48,-6.0,6.0); 
    TH1D* histzl = new TH1D("histzl","Data; Z-Score (Algorithm 0N); Bins per 0.25",48,-6.0,6.0);
    TH1D* histzo = new TH1D("histzo","Data; Z-Score (Old Algorithm); Bins per 0.25",48,-6.0,6.0); 

    TH1D* histzdu = new TH1D("histzdu","Data; Upper Z-Score (no symmetrization); Bins per 0.10",40,2.0,6.0);
    TH1D* histzdl = new TH1D("histzdl","Data; Lower Z-Score (no symmetrization); Bins per 0.10",40,-6.0,-2.0);

// Summary plots with paper binning   
    TH1D* mhistz  = new TH1D("mhistz","Data; Z-Score (Algorithm 0MH); Bins per 0.25",32,-4.0,4.0);
    TH1D* mhistzFinal  = new TH1D("mhistzFinal","Data; Z-Score (Latest Algorithm); Bins per 0.25",32,-4.0,4.0);
    TH1D* mhistzo = new TH1D("mhistzo","Data; Z-Score (Old Algorithm); Bins per 0.25",32,-4.0,4.0); 
    
// Region specific bins with paper binning
    TH1D* mhistz0L = new TH1D("mhistz0L","Data; Z-Score (Algorithm 0MH); Bins per 0.25",32,-4.0,4.0); 
    TH1D* mhistz1L = new TH1D("mhistz1L","Data; Z-Score (Algorithm 0MH); Bins per 0.25",32,-4.0,4.0); 
    TH1D* mhistz2L = new TH1D("mhistz2L","Data; Z-Score (Algorithm 0MH); Bins per 0.25",32,-4.0,4.0); 
    TH1D* mhistz3L = new TH1D("mhistz3L","Data; Z-Score (Algorithm 0MH); Bins per 0.25",32,-4.0,4.0); 
    TH1D* mhistzSV = new TH1D("mhistzSV","Data; Z-Score (Algorithm 0MH); Bins per 0.25",32,-4.0,4.0);    
    TH1D* mhistzBronze = new TH1D("mhistzBronze","Data; Z-Score (Algorithm 0MH); Bins per 0.25",32,-4.0,4.0);
    TH1D* mhistznonBronze   = new TH1D("mhistznonBronze","Data; Z-Score (Algorithm 0MH); Bins per 0.25",32,-4.0,4.0); 
    TH1D* mhistz1LnonBronze = new TH1D("mhistz1LnonBronze","Data; Z-Score (Algorithm 0MH); Bins per 0.25",32,-4.0,4.0); 
    TH1D* mhistz2LnonBronze = new TH1D("mhistz2LnonBronze","Data; Z-Score (Algorithm 0MH); Bins per 0.25",32,-4.0,4.0); 
    TH1D* mhistz3LnonBronze = new TH1D("mhistz3LnonBronze","Data; Z-Score (Algorithm 0MH); Bins per 0.25",32,-4.0,4.0);
    
    TH1D* mhistz0LSR = new TH1D("mhistz0LSR","Data; Z-Score (Algorithm 0MH); Bins per 0.25",32,-4.0,4.0);     
    TH1D* mhistz1LSR = new TH1D("mhistz1LSR","Data; Z-Score (Algorithm 0MH); Bins per 0.25",32,-4.0,4.0); 
    TH1D* mhistz2LSR = new TH1D("mhistz2LSR","Data; Z-Score (Algorithm 0MH); Bins per 0.25",32,-4.0,4.0); 
    TH1D* mhistz3LSR = new TH1D("mhistz3LSR","Data; Z-Score (Algorithm 0MH); Bins per 0.25",32,-4.0,4.0);               

    TH1D* mhistzSRD0 = new TH1D("mhistzSRD0","Data; Z-Score (Algorithm 0MH); Bins per 0.25",32,-4.0,4.0);
    TH1D* mhistzSRD1 = new TH1D("mhistzSRD1","Data; Z-Score (Algorithm 0MH); Bins per 0.25",32,-4.0,4.0);
    TH1D* mhistzNSRD1 = new TH1D("mhistzNSRD1","Data; Z-Score (Algorithm 0MH); Bins per 0.25",32,-4.0,4.0);    
    
    ifstream myfile(infile);
    
    int f1;    // nobs
    string f2; // regionName
    int f3;    // binNumber
    double f4, f5, f6, f7, f8, f9, f10; // prePull etc.

    string line;

    if( myfile.is_open() ){
 
        while( getline( myfile, line ) ) 
        {

             stringstream ss(line);

             ss >> f1 >> f2 >> f3 >> f4 >> f5 >> f6 >> f7 >> f8 >> f9 >> f10;
//             cout << f1 << " " << f2 << " " << f3 << " " << f4 << " " << f5 << " " 
//                  << f6 << " " << f7 << " " << f8 << " " << f9 << " " << f10 << endl;
                  
             struct AnalysisBin aBin;
             
             // First calculate the satChisq contribution from this bin
             int nobs = int(f7+0.01);     // ndata
             double y = f9;               // expected
             double satChisq;
             double n = double(nobs);
      // Use Baker and Cousins Poisson formula.  NIM 221 (1984) 437.
      // This ignores completely the postfit uncertainty.
             if(nobs == 0){
                 satChisq = 2.0*(y - n);
             }
             else{
                 satChisq = 2.0*(y - n + n*log(n/y) );
             }
             
             aBin = {f1, f2, nobs, int(f3+0.01), f9, f10, satChisq};
             vAnaBins.push_back(aBin);       
             
        }
        myfile.close();
    }
    else{
        cout << "Unable to open file" << endl;
    }
    
// Sort vector using saturated chisq statistic (see struct implementation)
    std::sort(vAnaBins.begin(), vAnaBins.end() );     
    
    int nbin = 0;
    
    double satChisqTot = 0.0;
// Print a header
//   cout << "Rank       ID        regionName      subBin     ndata    bMean   bErr    satChisq     satChisqTotal " << endl;    
    cout << "   Rank      ID                                           regionName   sbin    ndata       bMean        bErr       satChisq   satChisqTotal" << endl;
    for (auto & el : vAnaBins){
        nbin += 1;
        auto id = el.id;
        auto regionName = el.regionName;
        auto ndata = el.ndata;
        auto subBin = el.subBin;
        auto postfitMean = el.postfitMean;
        auto postfitError = el.postfitError;
        auto satChisq = el.satChisq;
        satChisqTot += satChisq;
        cout << setw(7) << nbin << " " << setw(7) << id << setw(53) << regionName << setw(7) << subBin 
                          << setw(9) << ndata << setw(12) << fixed << setprecision(4) << postfitMean << setw(12) <<fixed << setprecision(4) 
                          << postfitError << "   " << fixed << setprecision(5) << setw(12) << satChisq << "  " << fixed << setprecision(5) << setw(14) << satChisqTot << endl;
    }      
    cout << "Saturated chi-squared total " << satChisqTot << endl;
    
// Now start throwing toys for each bin

    cout << "start toys ... " << endl;
    
    nbin = 0;
    unsigned long int seed = baseseed;
    for (auto & el : vAnaBins){
         nbin +=1;
         auto id = el.id;
         auto regionName = el.regionName;
         auto ndata = el.ndata;
         auto subBin = el.subBin;
         auto postfitMean = el.postfitMean;
         auto postfitError = el.postfitError;
         auto satChisq = el.satChisq;
         
         double fracError = postfitError/postfitMean;
         seed += 1;
         
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
         if (stdev > 5.0) ntouse = 250000/ndivide;
         if (stdev > 6.0) ntouse = 1000000/ndivide;         
         
 // Throw toys if a good use of our time.
         double znew, zupper, zlower, zold, zscoreError;
         
         if(stdev > 20.0){
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
             histzo->Fill(zold);
             histz->Fill(znew);
             histzl->Fill(zlower);
             histzu->Fill(zupper); 

             mhistzo->Fill(zold); 
             mhistz->Fill(znew); 
           
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
         
             cout << "New raw z-score for bin " << nbin << " " << id << " " << std::get<10>(t)  << " +- " << std::get<11>(t) << endl; 

             double zbias = 0.0;
             double zscaleSD = 1.0;
             auto result = std::make_pair(zbias, zscaleSD);

             if(lcorr){
// Compute zscore bias and scale corrections. These DO NOT DEPEND on random numbers.
                 result = ComputeZscoreLocationAndScaleCorrection( postfitMean );
                 zbias = result.first;
                 zscaleSD = result.second;
             }
         
             znew = std::get<10>(t);
             zscoreError = std::get<11>(t);
             zupper = std::get<6>(t);
             zlower = std::get<8>(t);
             zold = zupper;
             if(ndata < postfitMean)zold = zlower;

// Standardization corrections for znew
             double znewCorr = (znew - zbias)/zscaleSD;
// Also scale zscoreError too.
             double zscoreErrorCorr = zscoreError/zscaleSD;
             if(lcorr){
                 std::cout << "Corrected z-score for bin " << nbin << " " << id << " " 
                           << znewCorr  << " +- " << zscoreErrorCorr << " (corrs = " << zbias << " " << zscaleSD << ")" << std::endl; 
             }
// As a sanity check also print more standard statistics too
             double zPull = ComputeNaivePull( ndata, postfitMean );
             std::cout << "Naive standardized pull (nObs - mu)/sqrt(mu) = " << zPull << std::endl;
             if ( postfitMean < 300.0 ){
                 auto poissResult = ComputePoissonOnlyStats( ndata, postfitMean );
                 std::cout << "Poisson-only standard p-value " << poissResult.first 
                           << " and corresponding Z-score " << poissResult.second << std::endl; 
             }

// Current choice.
// 1. |znew| < 1.25       Apply corrections for both bias and scale
// 2. zupper > 2.0.       Use conservative choice of just zscoreu with no bias corrections
// 3. zlower < -2.0       Use conservative choice of jsut zscorel with no bias corrections
// else                   Use znew and ony apply bias correction not the scale correction (to avoid eg a Z-value of 2.0 getting scaled up erroneously to 3.0)
//
// Note bias scale corrections can be severe espeically for small mu. For mu=0.174 (the lowest value in SUS-23-003) 
// Z is scaled up by 1.0/0.6217 = 1.61. So choice of 1.5 keeps a Z=1.5 < 2.5 ...

             double znewValue = znew;
 
             if(zupper >= 2.5){
// switch to zupper and its error
                 znew = zupper;
                 zscoreError = std::get<7>(t);
             }
             else if(zlower <= -2.5){
// switch to zlower and its error
                 znew = zlower;
                 zscoreError = std::get<9>(t);
             }
             else if(std::abs(znewValue) < 1.5){
                 znew = znewCorr;
                 zscoreError = zscoreErrorCorr;
             }
             else{
                 znew = znew - zbias;
             }


             std::cout << "Assigned Z-value of " << znew << " +- " << zscoreError << std::endl;
 
             histzo->Fill(zold);
             histz->Fill(znewValue); 
             histzl->Fill(zlower);
             histzu->Fill(zupper);
             histzFinal->Fill(znew);

             histzdu->Fill(zupper);
             histzdl->Fill(zlower);

             mhistzo->Fill(zold);
             mhistz->Fill(znewValue);
             mhistzFinal->Fill(znew);
             
             bool SR=false;
             
             if(el.is0L())mhistz0L->Fill(znewValue);
             if(el.is1L())mhistz1L->Fill(znewValue);
             if(el.is2L())mhistz2L->Fill(znewValue);
             if(el.is3L())mhistz3L->Fill(znewValue);
             if(el.isSV())mhistzSV->Fill(znewValue);             
             if(el.isBronze()){
                 mhistzBronze->Fill(znewValue);
             }
             else{
                 mhistznonBronze->Fill(znewValue);
                 if(el.is1L())mhistz1LnonBronze->Fill(znewValue);
                 if(el.is2L())mhistz2LnonBronze->Fill(znewValue);
                 if(el.is3L())mhistz3LnonBronze->Fill(znewValue);
                 if( (el.is0L() && subBin > 3) || (el.is1L() && subBin > 2) || (el.is2L() && subBin >4) || el.is3L() ){
                      SR = true;
                      mhistzSRD1->Fill(znewValue);
                 }  
                 if(subBin > 2){
                     mhistzSRD0->Fill(znewValue);
                     if(el.is0L())mhistz0LSR->Fill(znewValue);                     
                     if(el.is1L())mhistz1LSR->Fill(znewValue);
                     if(el.is2L())mhistz2LSR->Fill(znewValue);
                     if(el.is3L())mhistz3LSR->Fill(znewValue);                     
                 }  
             }
             if(!SR)mhistzNSRD1->Fill(znewValue);                           
         }

// Push bin into new struct.
         struct statAnalysisBin aBin;
         aBin = {id, regionName, ndata, subBin, postfitMean, postfitError, satChisq, znew, zscoreError };
         vstatAnaBins.push_back(aBin);                  
    }
// And now sort again and print.

// Sort vector using |zscore| (see struct implementation). Hopefully this works for the inherited one ...
    std::sort(vstatAnaBins.begin(), vstatAnaBins.end()  ); 
    std::sort(vstatAnaBins.begin(), vstatAnaBins.end(), sort_by_abszscore );    
    
    cout << "   Rank      ID                                           regionName   sbin    ndata       bMean        bErr       satChisq   satChisqTotal   zscore   zscoreError" << endl; 
    nbin = 0;
    satChisqTot = 0.0;
    for (auto & el : vstatAnaBins){
         nbin += 1;
         auto id = el.id;
         auto regionName = el.regionName;
         auto ndata = el.ndata;
         auto subBin = el.subBin;
         auto postfitMean = el.postfitMean;
         auto postfitError = el.postfitError;
         auto satChisq = el.satChisq;
         auto zscore = el.zscore;
         auto zscoreError = el.zscoreError;
 
         satChisqTot += satChisq;
         cout << setw(7) << nbin << " " << setw(7) << id << setw(53) << regionName << setw(7) << subBin 
                         << setw(9) << ndata << setw(12) << fixed << setprecision(4) << postfitMean << setw(12) << fixed << setprecision(4) 
                         << postfitError << "   " << fixed << setprecision(5) << setw(12) << satChisq << "  " 
                         << fixed << setprecision(5) << setw(14) << satChisqTot  
                         << fixed << setprecision(5) << setw(10) << zscore << fixed << setprecision(5) << setw(10) << zscoreError << endl; 
    }                                     
        
    f->Write();
    f->Close();         

}

int main(int argc, char** argv){

    CLI::App app{"Calculate z-scores using toys from fit result file"}; 
    
    int nfactor = 100;
    app.add_option("-n,--nfactor", nfactor, "Thousand toy multiplier"); 
    
    bool poissonOnly = false;
    app.add_option("-p,--poissonOnly", poissonOnly, "Boolean for Poisson Only");
    
    int ntosave = 100;
    app.add_option("-s,--ntosave", ntosave, "Number of toys to save"); 
    
    int ndivide = 1;
    app.add_option("-d,--ndivide", ndivide, "Divisor for quick test");

    bool lcorr = true;
    app.add_flag("--corrections,!--no-corrections",lcorr,"Enable/disable Zscore corrections");
    
    std::string infilename = "B135_7-7-24_all_Modified_V1.txt";
    app.add_option("-i,--infilename", infilename, "Input file from fit to parse");  
    
    unsigned long int baseseed = 123456L;
//    unsigned long int baseseed = 125900L;    
    app.add_option("-b,--baseseed", baseseed, "Base seed for toys");    
              
    CLI11_PARSE(app, argc, argv);
    
    std::cout << "nfactor  " << nfactor << std::endl;    
    std::cout << "poissonOnly " << poissonOnly << std::endl;
    std::cout << "ntosave " << ntosave << std::endl;
    std::cout << "ndivide " << ndivide << std::endl;
    std::cout << "lcorr   " << lcorr << std::endl;
    std::cout << "baseseed " << baseseed << std::endl;
    std::cout << "infilename " << infilename << std::endl;   
    
// Global value
    NTOSAVE = ntosave;

    Parser(nfactor, poissonOnly, ndivide, lcorr, baseseed, infilename);
    
// Sort global toys vector
    std::sort(vtoys.begin(), vtoys.end());

    std::ofstream fout;
    fout.open("GeneratedToys.dat");   
    for (auto & el : vtoys){
//         cout << " Toy " << std::get<0>(el) << " " << std::get<1>(el) << " " << std::get<2>(el) << " " << std::get<3>(el) << endl;
         fout << std::get<0>(el) << " " << std::get<1>(el) << " " << std::get<2>(el) << " " << std::get<3>(el) << endl;         
    } 
    fout.close();   
       
    return 0;
    
}

