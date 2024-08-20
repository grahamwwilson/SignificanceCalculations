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

std::vector<statAnalysisBin> vstatAnaBins;

void ParseZScore(){
    
    ifstream myfile("B135_7-7-24_all_Modified_V1_WithZScore-NoHeader.txt");
    
    int f0;     // rank
    int f1;     // id
    string f2;  // regionName
    int f3;     // binNumber
    int f4;     // ndata
    double f5;  // bMean
    double f6;  // bErr
    double f7;  // satChisq
    double f8;  // satChisqTotal
    double f9;  // zscore
    double f10; // zscoreError

    string line;

    if( myfile.is_open() ){
 
        while( getline( myfile, line ) ) 
        {

             stringstream ss(line);

             ss >> f0 >> f1 >> f2 >> f3 >> f4 >> f5 >> f6 >> f7 >> f8 >> f9 >> f10;
                  
             struct statAnalysisBin sBin;
                          
             sBin = {f1, f2, f4, f3, f5, f6, f7, f9, f10}; // See statAnalysisBin.h and AnalysisBin.h for ordering
             vstatAnaBins.push_back(sBin);       
             
        }
        myfile.close();
    }
    else{
        cout << "Unable to open file" << endl;
    }
    
    int nbin = 0;
    
    double satChisqTot = 0.0;  
    cout << "   Rank      ID                                           regionName   sbin    ndata       bMean        bErr       satChisq    zscore   zscoreError" << endl;
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
                          << setw(9) << ndata << setw(12) << fixed << setprecision(4) << postfitMean << setw(12) <<fixed << setprecision(4) 
                          << postfitError << "   " << fixed << setprecision(5) << setw(12) << satChisq << "  " 
                          << fixed << setprecision(5) << setw(10) << zscore << fixed << setprecision(5) << setw(10) << zscoreError <<
                          endl;
    }      
    cout << "Saturated chi-squared total " << satChisqTot << endl;

}

int main(int argc, char** argv){

    ParseZScore();
       
    return 0;
    
}
