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
using namespace std;

// May be best to make a vector of structs. 
//
// Can a new struct inherit from an existing struct?

struct AnalysisBin {
    int id {};
    std::string regionName {};
    int ndata {};
    int subBin {};
    double postfitMean {};
    double postfitError {};
    double satChisq {};
 // Standard sorting criterion    
    bool operator < (const AnalysisBin & aBin) const
    {
        return satChisq > aBin.satChisq;
    }
};

struct statAnalysisBin : AnalysisBin{
    double zscore {};
    double zscoreError {};
    bool operator < (const statAnalysisBin & aBin) const
    {
         return abs(zscore) > abs(aBin.zscore);
    }    
       
};

bool sort_by_abszscore( const statAnalysisBin & lhs, const statAnalysisBin & rhs )
{
   
   return abs(lhs.zscore) > abs(rhs.zscore);
}

int NTOSAVE {};
std::vector<std::tuple<int, unsigned long int, int, int>> vtoys;  // Design with experiment-number, bin-id, nobs, seed. Can then sort by experiment number.

std::vector<std::tuple<int, unsigned long int, int, int>> gtoys;  // Design with experiment-number, bin-id, nobs, seed. Can then sort by experiment number.

#include "ToyChucker.cpp" 

std::vector<AnalysisBin> vAnaBins;
std::vector<statAnalysisBin> vstatAnaBins;

void Parser(int nfactor, bool PoissonOnly, int ndivide, unsigned long int baseseed){

    // cout << "Hello from Parser " << endl;
    
    TFile *f = new TFile("Analysis.root","RECREATE");
    TH1D* hist = new TH1D("hist","Data; Z-Score (Algorithm 1); Bins per 0.25",49,-6.125,6.125);    
    TH1D* histz = new TH1D("histz","Data; Z-Score (Algorithm 0MH); Bins per 0.25",49,-6.125,6.125); 
    TH1D* histzu = new TH1D("histzu","Data; Z-Score (Algorithm 0P); Bins per 0.25",49,-6.125,6.125); 
    TH1D* histzl = new TH1D("histzl","Data; Z-Score (Algorithm 0N); Bins per 0.25",49,-6.125,6.125);
    
    TH1D* hist_censored = new TH1D("hist_censored","Data; Z-Score (Algorithm 1); Bins per 0.25",49,-6.125,6.125);    
    TH1D* histz_censored = new TH1D("histz_censored","Data; Z-Score (Algorithm 0MH); Bins per 0.25",49,-6.125,6.125); 
    
    TH1D* nhist = new TH1D("nhist","Data; Z-Score (Algorithm 1); Bins per 0.25",31,-3.875,3.875);    
    TH1D* nhistz = new TH1D("nhistz","Data; Z-Score (Algorithm 0MH); Bins per 0.25",31,-3.875,3.875);                   
    
    ifstream myfile("B120_7-2-24_all-1_Cropped.txt");
    
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
             int id = f1;                 // binnumber id
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
             
      // Now we want to find the corresponding toy element in the gtoys vector for this toy experiment
      // For now we want to be careful that we match it correctly to this data bin
             for (auto & el : gtoys){
                  auto toy_expt = std::get<0>(el);
                  auto toy_seed = std::get<1>(el);
                  auto toy_id   = std::get<2>(el);
                  auto toy_nobs = std::get<3>(el);
                  if(toy_id == id){  // we have a match
                     // Overwrite observed data and corresponding statistic
                     nobs = toy_nobs;
                     double n = nobs;
                     if(nobs == 0){
                        satChisq = 2.0*(y - n);
                     }
                     else{
                        satChisq = 2.0*(y - n + n*log(n/y) );
                     }
                     cout << "Over-wrote using toy value of " << toy_nobs << " for bin-id " << toy_id << endl;
                  }
//         cout << " Toy " << std::get<0>(el) << " " << std::get<1>(el) << " " << std::get<2>(el) << " " << std::get<3>(el) << endl;        
             } 
                  
             aBin = {f1, f2, nobs, f3, f9, f10, satChisq};
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
         
 // Throw toys if a good use of our time.
         double znew, zupper, zlower, zold, zscoreError;
         
         if(stdev > 6.0){
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
             std::tuple<int, double, int, double, double, double, double, double, double, double, double, double> t = 
                       ToyChucker(NTOUSE, ndata, ndata, postfitMean, fracError, seed, id );             
             
         }
         else{
          
             if(PoissonOnly)fracError=1.0e-6;
          
             std::tuple<int, double, int, double, double, double, double, double, double, double, double, double> t = 
                       ToyChucker(ntouse, ndata, ndata, postfitMean, fracError, seed, id );
         
             cout << "New z-score for bin " << nbin << " " << id << " " << std::get<10>(t)  << " +- " << std::get<11>(t) << endl; 
         
             znew = std::get<10>(t);
             zscoreError = std::get<11>(t);
             zupper = std::get<6>(t);
             zlower = std::get<8>(t);
             zold;
             zold = zupper;
             if(ndata < postfitMean)zold = zlower;
             histz->Fill(znew); 
             hist->Fill(zold);
             nhistz->Fill(znew); 
             nhist->Fill(zold);             
             histzl->Fill(zlower);
             histzu->Fill(zupper);
             histz_censored->Fill(znew); 
             hist_censored->Fill(zold);                  
         }

// Push bin into new struct.
         struct statAnalysisBin aBin;
         aBin = {id, regionName, ndata, subBin, postfitMean, postfitError, satChisq, znew, zscoreError };
         vstatAnaBins.push_back(aBin);                  
    }
// And now sort again and print.

// Sort vector using |zscore| (see struct implementation). Hopefully this works for the inherited one ...
    std::sort(vstatAnaBins.begin(), vstatAnaBins.end()  );
    
    
    cout << "   Rank      ID                                           regionName   sbin    ndata       bMean        bErr       satChisq   satChisqTotal   zscore   zscoreError" << endl; 
    nbin = 0;
    satChisqTot = 0.0;
    for (auto & el : vstatAnaBins){
         nbin +=1;
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
                          << setw(9) << ndata << setw(12) << fixed << setprecision(4) << postfitMean << setw(12) <<fixed << setprecision(4) 
                          << postfitError << "   " << fixed << setprecision(5) << setw(12) << satChisq << "  " 
                          << fixed << setprecision(5) << setw(14) << satChisqTot  
                          << fixed << setprecision(5) << setw(10) << zscore << fixed << setprecision(5) << setw(10) << zscoreError <<
                          endl; 
    }                                     
    
// Sort vector using |zscore| (see struct implementation). Hopefully this works for the inherited one ...
//    std::sort(vstatAnaBins.begin(), vstatAnaBins.end(), sort_by_abszscore );
    
    f->Write();
    f->Close();         

}

void ReadAToy(int toynumber, bool poissonOnly){
// Read a toy from the previously generated toys into the corresponding vector

    std::ifstream fin;
    if(poissonOnly){
        fin.open("GeneratedToys-PoissonOnly.dat");    
    }
    else{
        fin.open("GeneratedToys-Standard.dat");
    }
    
    int f1; 
    unsigned long int f2;
    int f3;
    int f4;
    
    std::string line;
    
    if( fin.is_open() ){
 
        while( getline( fin, line ) ) 
        {

             stringstream ss(line);

             ss >> f1 >> f2 >> f3 >> f4;
   // Could speed this up by quitting once we've filled the vector, but likely not time critical for now         
             if(f1==toynumber){
                std::tuple<int, unsigned long int, int, int> t = std::make_tuple(f1, f2, f3, f4);
                gtoys.push_back(t);
             }
             
        }
        fin.close();        
    }
    else{
        cout << "Unable to open file" << endl;
    }             
    cout << "Toys vector filled with " << gtoys.size() << " long experiments * bins ensemble " << endl;
    
}

int main(int argc, char** argv){

    CLI::App app{"Parse fit result file"}; 
    
    int nfactor = 100;
    app.add_option("-n,--nfactor", nfactor, "Thousand toy multiplier"); 
    
    bool poissonOnly = false;
    app.add_option("-p,--poissonOnly", poissonOnly, "Boolean for Poisson Only");
    
    int ntosave = 0;
    app.add_option("-s,--ntosave", ntosave, "Number of toys to save"); 
    
    int ndivide = 1;
    app.add_option("-d,--ndivide", ndivide, "Divisor for quick test");
    
    int toynumber = 5;
    app.add_option("-t,--toynumber", toynumber, "Toy number");    
    
//   unsigned long int baseseed = 123456L;
    unsigned long int baseseed = 125900L;    
    app.add_option("-b,--baseseed", baseseed, "Base seed for toys");    
              
    CLI11_PARSE(app, argc, argv);
    
    std::cout << "nfactor  " << nfactor << std::endl;    
    std::cout << "poissonOnly " << poissonOnly << endl;
    std::cout << "ntosave " << ntosave << endl;
    std::cout << "ndivide " << ndivide << endl;
    std::cout << "baseseed " << baseseed << endl; 
    std::cout << "toynumber " << toynumber << endl;       
    
    baseseed += (toynumber+1)*2444;
    
    std::cout << "baseseed updated to " << baseseed << endl;
    
// Global value
    NTOSAVE = ntosave;
    
// Read previously generated toys into the corresponding vector. Uses a globally declared vector.
    ReadAToy(toynumber, poissonOnly); 

    Parser(nfactor, poissonOnly, ndivide, baseseed);
    
// Can access global toys vector
// First sort it.
    std::sort(vtoys.begin(), vtoys.end());

    if ( NTOSAVE > 0 ){
        std::ofstream fout;
        fout.open("GeneratedToys.dat");   
        for (auto & el : vtoys){
//         cout << " Toy " << std::get<0>(el) << " " << std::get<1>(el) << " " << std::get<2>(el) << " " << std::get<3>(el) << endl;
             fout << std::get<0>(el) << " " << std::get<1>(el) << " " << std::get<2>(el) << " " << std::get<3>(el) << endl;         
        } 
        fout.close();   
    }
       
    return 0;
    
}

