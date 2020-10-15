// Check zBi and investigate methods for combining multiple channels
// Here assume for simplicity that the individual channels 
// are independent and the background estimates are uncorrelated.
//
//      Graham W. Wilson
// 
#include <iostream>
#include <TMath.h>   // ROOT's TMath
using namespace std;

double pBi(double ns, double mub, double syst){
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

// ns  = Nsignal events
// mub = Background mean
// syst = Fractional systematic uncertainty on background mean
   double sb = syst*mub;
   double tau = mub/(sb*sb);
   double non  = mub + ns;
   double noff = tau*mub;

   cout << "Calculating pBi for ns = " << ns << " mub= " << mub 
                                       << " syst = " << syst << endl;

   double p_Bi = TMath::BetaIncomplete(1.0/(1.0+tau), non, noff+1.0);

   return p_Bi;
}

double zSingleSided(double pvalue){
   double zvalue = sqrt(2.0)*TMath::ErfInverse(1.0 - 2.0*pvalue);
   return zvalue;
}

double zDoubleSided(double pvalue){
   double zvalue = sqrt(2.0)*TMath::ErfInverse(1.0 - pvalue);
   return zvalue;
}

double pLyonsChapon(double p1, double p2){
   // Eqn 17 for n=2 in arXiv:1704.05540 
   cout << "p1 = " << p1 << " p2 = " << p2 << endl; 
   double p = p1*p2*(1.0-log(p1*p2));
   return p;
}

void Analyze(double pb1, double pb2){

// Note this is single-sided
   double zb1 = zSingleSided(pb1);
   double zb2 = zSingleSided(pb2);

   cout << " pb1 " << pb1 << " zb1 " << zb1 << endl;
   cout << " pb2 " << pb2 << " zb2 " << zb2 << endl;

// Now combine these two tests assuming corresponding chi-squared 
// statistics with the same pvalue.

// Do we need this factor of two ?
// This is just chisq1 = zb1*zb1 etc.
   double chisq1 = TMath::ChisquareQuantile(1.0-2.0*pb1, 1);
   double chisq2 = TMath::ChisquareQuantile(1.0-2.0*pb2, 1);
   cout << "chisq1, chi1 " << chisq1 << " " << sqrt(chisq1) << endl;
   cout << "chisq2, chi2 " << chisq2 << " " << sqrt(chisq2) << endl;

   double chisqtot = chisq1 + chisq2;
   cout << "chisqtot = " << chisqtot << endl;

   double pvalue=TMath::Prob(chisqtot,2);
// Reduce the p-value by a factor of 4, since there 
// are four possible quadrants for the deviations in the two measurements, 
// and only one quadrant has both measurements with upward fluctuations
   cout << "Combined p-value of " << pvalue << " " << 0.25*pvalue << endl;

   double zcombined = zSingleSided(0.25*pvalue);
   cout << "Combined z-value of " << zcombined << endl;
 
// Lyons-Chapon paper uses the product of the p-values where assuming 
// these distributions are independent and uniformly distributed, one
// computes a corresponding p-value. 
// For n=2, the expression is simply P = p1 p2 * (1 - ln(p1*p2)).
   double pLC = pLyonsChapon(pb1,pb2);
   double zLC = zSingleSided(pLC);
   cout << "Lyons-Chapon pvalue: " << pLC << " z value: " << zLC << endl; 

// This is the simple method 3 in Lyons-Chapon that is unweighted 
// statistical averaging of z values.
   double zcomb = (zb1+zb2)/sqrt(2.0);
   cout << "Stouffer z comb = " << zcomb << endl;

// An equivalent formulation is to transform the pBi into a chi-squared 
// variate (without the factor 2). ie. find the chi-squared value 
// corresponding to the same upper-tail probabilities.
   double chisq1p = TMath::ChisquareQuantile(1.0-pb1, 1);
   double chisq2p = TMath::ChisquareQuantile(1.0-pb2, 1);
 
   double chisqtotp = chisq1p+chisq2p;
   cout << "Chi2 Values " << chisq1p << " " << chisq2p << " " 
                                     << chisqtotp << endl;
 
   double pvaluep = TMath::Prob(chisqtotp,2);
   double zvaluep = zSingleSided(pvaluep);

   cout << "found pvaluep = " << pvaluep << endl;
   cout << "found zvaluep = " << zvaluep << endl;

// Fisher's method. chisqa and chisqb each have 2 degrees of freedom
   double chisqa = -2.0*log(pb1);
   double chisqb = -2.0*log(pb2);
   double chisqf = chisqa + chisqb;
   double pvaluef = TMath::Prob(chisqf,4);
   double zvaluef = zSingleSided(pvaluef);
   cout << "Fisher's method " << chisqf << " " 
                              << pvaluef << " " << zvaluef << endl;

// Simple weighted average. Here the assumption is that the 
// observed z-values, measure the same deviation, Delta, but with 
// different uncertainties, sigma_i.
// ie. z_i = Delta/sigma_i. 
// The weighted average estimate of Delta then implies that the 
// combined z is simply the square root of the quadrature sum of the 
// z-values.
   double zweighted = sqrt(zb1*zb1 + zb2*zb2);

   cout << " " << endl;
   cout << " Combined z-value SUMMARY " << endl;
   cout << " Chi-square            1: " << zcombined << endl;
   cout << " Lyons-Chapon          2: " << zLC << endl;
   cout << " Stouffer              3: " << zcomb << endl;
   cout << " Alt. Chisq.           4: " << zvaluep << endl;
   cout << " Fisher                5: " << zvaluef << endl;
   cout << " Quadrature Weighted   6: " << zweighted << endl;
   cout << " " << endl;
   cout << " " << endl;

}

int main(){

// Play with a few different specific scenarios 
// and corresponding zBinomial based p-value.
// Look at various methods for combining two p-values 
// as would happen with say a search with two bins / channels 
// in the strict counting experiment sense. 
// The combination assumes uncorrelated systematics.
   
   double pb1 = pBi(69.0, 169.2, 0.1);
   double pb2 = pBi(10.4, 58.2, 0.1);
   Analyze(pb1, pb2);

   pb1 = pBi(69.0, 169.2, 0.1);
   pb2 = pBi(69.0, 169.2, 0.1);
   Analyze(pb1, pb2);

   pb1 = pBi(44.74, 198.2, 0.1);
   pb2 = pBi(32.1, 198.2, 0.1);
   Analyze(pb1, pb2);

   pb1 = pBi(81.3079, 8930.91, 0.1);
   pb2 = pBi(62.4938, 3048.57, 0.1);
   Analyze(pb1, pb2);

}
