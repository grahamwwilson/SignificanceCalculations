// You may choose to use #include "CLI11.hpp" to pick up the repo version header file 
#include "CLI11.hpp"
//#include <CLI11.hpp>   // CLI11 command line interface stuff (V2.3.1). See https://github.com/CLIUtils/CLI11
#include <iostream> 
#include <algorithm> //sort
#include <cmath>
#include <TRandom3.h>
#include <TMath.h>   //TMath::Prob
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <vector>    
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
using namespace std;

// 
// Purpose. Compare B120 and B135 and add errors to B135.
// 

struct AnalysisBin {
    int id {};
    std::string regionName {};
    int ndata {};
    int subBin {};
    double prefitMean {};
    double postfitMean {};
    double postfitError {};
    double satChisq {};
    
// Standard sorting criterion    
    bool operator < (const AnalysisBin & aBin) const
       {
        return id < aBin.id;
       }
    
 // Standard sorting criterion    
 //   bool operator < (const AnalysisBin & aBin) const
 //   {
 //       return satChisq > aBin.satChisq;
 //   }
};

//std::vector<AnalysisBin> vAnaBins;

std::vector<AnalysisBin> ParserWithFileModification(std::string infile,  std::string outfile, std::vector<AnalysisBin> vOldAnaBins ){

    cout << "Parser uses " << infile << endl;
    cout << "Writing modified input file as output file " << outfile << endl;
    
    ifstream myfile(infile);
    ofstream fout;
    fout.open(outfile);
    
    int f1;    // nobs
    string f2; // regionName
    int f3;    // binNumber
    double f4, f5, f6, f7, f8, f9, f10; // prePull etc.
    
    std::vector<AnalysisBin> vAnaBins;    

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
             double ypr = f8;             // prefit
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
             
//             aBin = {f1, f2, nobs,int(f3+0.01),f8,f9,f10,satChisq};
//             vAnaBins.push_back(aBin);
             
             // Now search through the old vector to find the corresponding bin
             
             double f10_old {};
             
             int nfound = 0;
             for (auto & el: vOldAnaBins){
                 auto id = el.id;
                 if(id == f1){ // Match!
                     auto postfitError = el.postfitError;
                     f10_old = postfitError;
                     nfound +=1;
                 }
             }          
             cout << "Found " << nfound << " bins from old file " << endl;
             // Write this modified event out. Initial implementation - use same error.
             fout << f1 << " " << f2 << " " << f3 << " " << f4 << " " << f5 << " " 
                  << f6 << " " << f7 << " " << f8 << " " << f9 << " " << f10_old << endl;
                  
// Overwrite error
             aBin = {f1, f2, nobs,int(f3+0.01),f8,f9,f10_old,satChisq};
             vAnaBins.push_back(aBin);                  
                  
        }
        myfile.close();
    }
    else{
        cout << "Unable to open file" << endl;
    }
    
    fout.close();
    
// Sort vector using id (see struct implementation)
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
    
    return vAnaBins;

}

std::vector<AnalysisBin> Parser(std::string infile){

    cout << "Parser uses " << infile << endl;
    
    ifstream myfile(infile);
    
    int f1;    // nobs
    string f2; // regionName
    int f3;    // binNumber
    double f4, f5, f6, f7, f8, f9, f10; // prePull etc.
    
    std::vector<AnalysisBin> vAnaBins;    

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
             
             aBin = {f1, f2, nobs,int(f3+0.01),f8,f9,f10,satChisq};
             vAnaBins.push_back(aBin);       
             
        }
        myfile.close();
    }
    else{
        cout << "Unable to open file" << endl;
    }
    
// Sort vector using id (see struct implementation)
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
    
    return vAnaBins;

}
int main(int argc, char** argv){

    CLI::App app{"Parse fit result files"};  
    
    CLI11_PARSE(app, argc, argv);

    std::vector<AnalysisBin> vAnaBinsB120 = Parser("B120_7-2-24_all-1_Cropped.txt");
    std::vector<AnalysisBin> vAnaBinsB135 = ParserWithFileModification("B135_7-7-24_all_Cropped.txt","B135_7-7-24_all_Modified_V2.txt", vAnaBinsB120);
    
// Histogram the bin-by-bin content of the files to help identify differences. 

    TFile *f = new TFile("Compare.root","RECREATE");
    TH1D* histbMean1 = new TH1D("histbMean1","B120; Analysis Bin; Posterior Mean",3913,-0.5,3912.5);
    TH1D* histbMean2 = new TH1D("histbMean2","B135; Analysis Bin; Posterior Mean",3913,-0.5,3912.5);
    TH1D* histbErr1 = new TH1D("histbErr1","B120; Analysis Bin; Posterior Uncertainty",3913,-0.5,3912.5);
    TH1D* histbErr2 = new TH1D("histbErr2","B135; Analysis Bin; Posterior Uncertainty ",3913,-0.5,3912.5); 
    TH1D* histndata1 = new TH1D("histndata1","B120; Analysis Bin; Observed Data",3913,-0.5,3912.5);
    TH1D* histndata2 = new TH1D("histndata2","B135; Analysis Bin; Observed Data",3913,-0.5,3912.5);    
    TH1D* histfErr1 = new TH1D("histfErr1","B120; Analysis Bin; Posterior Uncertainty / Mean",3913,-0.5,3912.5);    
    TH1D* histdiffMean = new TH1D("histdiffMean","B120/B135 Comparison; Analysis Bin; Difference in Posterior Mean",3913,-0.5,3912.5);
    TH1D* histdiffData = new TH1D("histdiffData","B120/B135 Comparison; Analysis Bin; Difference in Observed Data",3913,-0.5,3912.5); 
    TH1D* histaverageMean = new TH1D("histaverageMean","B120/B135 Comparison; Analysis Bin; Average Posterior Mean",3913,-0.5,3912.5);
    TH1D* histdifffracMean = new TH1D("histdifffracMean",
    "B120/B135 Comparison; Analysis Bin; Fractional Difference in Posterior Mean 2(B135-B120)/(B135+B120)",3913,-0.5,3912.5);   
    TH1D* histdifffracMeanDistbn = new TH1D("histdifffracMeanDistbn",
     "B120/B135 Comparison; Fractional Difference in Posterior Mean 2(B135-B120)/(B135+B120)",100,-0.5,0.5);   
    TH1D* histfErrDistbn = new TH1D("histfErrDistbn",
     "B120/B135 Comparison; Fractional Error in Posterior Mean (B120)",200,0.0,0.2); 
    TH2D* histfErrvsMean = new TH2D("histfErrvsMean","B120; 1/sqrt(Posterior Mean); fErr", 50, 0.0, 1.0, 50,0.0, 0.1);
    TH1D* histPrePostPull1 = new TH1D("histPrePostPull1","B120; (Post-fit - Data)/(Pre-fit - Data)",1000,-10.0,10.0);
    TH1D* histPrePostPull2 = new TH1D("histPrePostPull2","B135; (Post-fit - Data)/(Pre-fit - Data)",1000,-10.0,10.0);  
    
    TH1D* histPull1 = new TH1D("histPull1","B120; (Data - Post-fit)/sqrt(Post-fit)",1000,-10.0,10.0);
    TH1D* histPull2 = new TH1D("histPull2","B135; (Data - Post-fit)/sqrt(Post-fit)",1000,-10.0,10.0);
    TH1D* histLbMean1 = new TH1D("histLbMean1","B120; log10(Post-fit Mean)",35,-2.0,5.0);
    TH1D* histLbMean2 = new TH1D("histLbMean2","B135; log10(Post-fit Mean)",35,-2.0,5.0);
    TH1D* histMean1 = new TH1D("histMean1","B120; Post-fit Mean",100000,0.0,100000.0);
    TH1D* histMean2 = new TH1D("histMean2","B135; Post-fit Mean",100000,0.0,100000.0);    
    TH1D* histInvTau1 = new TH1D("histInvTau1","B120; 1/Tau",100,0.0,2.0);
    TH1D* histInvTau2 = new TH1D("histInvTau2","B135; 1/Tau",100,0.0,2.0); 
    TH1D* histFracErr1 = new TH1D("histFracErr1","B120; Post-fit Fractional Uncertainty",100,0.0,0.2);
    TH1D* histFracErr2 = new TH1D("histFracErr2","B135; Post-fit Fractional Uncertainty",100,0.0,0.2); 
    TH1D* histzPull1 = new TH1D("histzPull1","B120; Z (Pull)",64,-8.0,8.0);
    TH1D* histzPull2 = new TH1D("histzPull2","B135; Z (Pull)",64,-8.0,8.0); 
    TH1D* histzPullA1 = new TH1D("histzPullA1","B120; Z (Pull)",64,-8.0,8.0);
    TH1D* histzPullA2 = new TH1D("histzPullA2","B135; Z (Pull)",64,-8.0,8.0);
    TH1D* histzBi1 = new TH1D("histzBi1","B120; Z Bi",64,-8.0,8.0);
    TH1D* histzBi2 = new TH1D("histzBi2","B135; Z Bi",64,-8.0,8.0);                             
            
    for (auto & el : vAnaBinsB120){
        auto id = el.id;
        auto prefitMean = el.prefitMean;
        auto postfitMean = el.postfitMean;
        auto postfitError = el.postfitError;
        auto ndata = el.ndata;
        histbMean1->Fill(id,postfitMean);
        histbErr1->Fill(id,postfitError);
        histndata1->Fill(id,ndata);
        histfErr1->Fill(id, postfitError/postfitMean);
        histfErrvsMean->Fill( 1.0/sqrt(postfitMean), postfitError/postfitMean);
        double prepostPull = (postfitMean - double(ndata))/(prefitMean - double(ndata));
        if(prepostPull < -10.0) prepostPull = -9.9999;
        if(prepostPull >  10.0) prepostPull =  9.9999;        
        histPrePostPull1->Fill( prepostPull );
        histPull1->Fill( (double(ndata) - postfitMean)/sqrt(postfitMean) );
        histMean1->Fill(postfitMean);
        histLbMean1->Fill(log(postfitMean)/log(10.0));
        double tau = postfitMean/(postfitError*postfitError);
        double non = double(ndata);
        double noff = tau*postfitMean;
        double pBi = TMath::BetaIncomplete( 1.0/(1.0+tau), non, noff+1.0);
        double zBi = sqrt(2.0)*TMath::ErfInverse(1.0 - 2.0*pBi);
        histzBi1->Fill(zBi);
        histInvTau1->Fill(1.0/tau);
        histFracErr1->Fill(postfitError/postfitMean);
        double zpull = (double(ndata)-postfitMean)/sqrt(postfitMean);
        histzPull1->Fill(zpull);
        double zpulla = (double(ndata)-postfitMean)/sqrt(postfitMean + postfitError*postfitError);
        histzPullA1->Fill(zpulla);       
    } 
    for (auto & el : vAnaBinsB135){
        auto id = el.id;
        auto prefitMean = el.prefitMean;        
        auto postfitMean = el.postfitMean;
        auto postfitError = el.postfitError;
        auto ndata = el.ndata;
        histbMean2->Fill(id,postfitMean);
        histbErr2->Fill(id,postfitError);
        histndata2->Fill(id,ndata);
        double prepostPull = (postfitMean - double(ndata))/(prefitMean - double(ndata));
        if(prepostPull < -10.0) prepostPull = -9.9999;
        if(prepostPull >  10.0) prepostPull =  9.9999;         
        histPrePostPull2->Fill( prepostPull );
        histPull2->Fill( (double(ndata) - postfitMean)/sqrt(postfitMean) );
        histMean2->Fill(postfitMean);        
        histLbMean2->Fill(log(postfitMean)/log(10.0)); 
        double tau = postfitMean/(postfitError*postfitError);
        double non = double(ndata);
        double noff = tau*postfitMean;
        double pBi = TMath::BetaIncomplete( 1.0/(1.0+tau), non, noff+1.0);
        double zBi = sqrt(2.0)*TMath::ErfInverse(1.0 - 2.0*pBi);
        histzBi2->Fill(zBi);        
        histInvTau2->Fill(1.0/tau); 
        histFracErr2->Fill(postfitError/postfitMean);
        double zpull = (double(ndata)-postfitMean)/sqrt(postfitMean);
        histzPull2->Fill(zpull);
        double zpulla = (double(ndata)-postfitMean)/sqrt(postfitMean + postfitError*postfitError);
        histzPullA2->Fill(zpulla);                
        if( (1.0/tau) > 1.0 ){
            cout << "Big 1/tau " << 1.0/tau << " " << el.regionName << " " << el.id << " " << el.subBin << " " << ndata << " " << postfitMean << " " << postfitError << endl;
        }                    
    }
    histdiffMean->Add(histbMean1,-1.0);
    histdiffMean->Add(histbMean2,+1.0);
    histaverageMean->Add(histbMean1,0.5);
    histaverageMean->Add(histbMean2,0.5);        
    histdiffData->Add(histndata1,-1.0);
    histdiffData->Add(histndata2,+1.0);    
    
    histdifffracMean->Add(histdiffMean,1.0);
    histdifffracMean->Divide(histaverageMean); 
    
    for (int i=1; i<=3913;i++){
         histdifffracMean->SetBinError(i, 0.0);
         double val=histdifffracMean->GetBinContent(i);
         double ferr=histfErr1->GetBinContent(i);
         if(abs(val)>1.0e-6){
             histdifffracMeanDistbn->Fill(val);
             histfErrDistbn->Fill(ferr);
         }
    }
 // FIXME.  Check 2-d distribution of fErr vs bMean   
    
    f->Write();
    f->Close();                
    return 0;
    
}

