// Notes. Consider also to return an error flag 

double MyQuantile(double p, int printlevel = 0){

// Find standard normal equivalent Z value corresponding to upper tail probability integral p

// Convert one-sided p-value into equivalent Gaussian significance.
// Clamp pathological p=0 or p=1 cases since ROOT's NormQuantile justifiably does not return ±∞.

// NormQuantile gives an error if called with p<=0||p>=1 so we need to avoid this.
//   if ((p<=0)||(p>=1)) {
//      Error("TMath::NormQuantile", "probability outside (0, 1)");
//      return 0;
//   }
//
// Small values of p such as 1.0e-20 are easily representable as doubles that are non-zero, 
// but q = 1.0 - 1.0e-20 rounds to 1.0 with roundoff error.

    double q = 1.0 - p;
    double Zq {};
    double Zp {};

    int icase = 0;
    
    if( q<=0 ){
//
// UPGRADE? Potentially we could swap between left-tail and right-tail values in this clause too if we had a high-precision 
// value for q when it is small as an additional input parameter.
// The current implementation clamping starts when q < 1.11e-16 (64-bit epsilon such that 1 + epsilon = 1).
        icase = -1;
        Zq =  -8.3;
        if ( printlevel!=0 ) std::cout << "CASE -1: Clamped Zq (infinitesimal q) " << Zq << " for p = " << p << std::endl;
    }
    else if (q >=1){
//        std::cout << "Error for MyQuantile with arguments " << p << "( " << q << ") Setting value of +38.5 " << std::endl;
        Zq =  38.5;    
        if(p>0){
// We can still compute in these cases where q=1 numerically, the complementary Z for p in the range [4.9e-324, 1.11e-16]
            icase = 1;
            Zp = TMath::NormQuantile(p);
            Zq = -Zp;
            if ( printlevel==1 ) std::cout << "CASE 1: (small p) Zp, Zq = " << Zp << " " << Zq << " for p = " << p << std::endl;
        }
        else{
            icase = 2;
            if ( printlevel!=0 ) std::cout << "CASE 2: Clamped Zq (infinitesimal p) " << Zq << " for p = " << p << std::endl;
        }
    }
    else{
//
// Standard case (icase = 0)
// Compute the quantile (ie the cut-point) for the standard normal distribution where the cdf has a value of q. 
// ie. Integral(-infty, Z) f(z) = q, where f is the probability density function. Or rather Integral(Z, infty) = p.
//
        Zq = TMath::NormQuantile(q);
    }

    if(printlevel!=0){
       std::cout << std::scientific << std::setprecision(20) << "MyQuantile " << " case " << icase << " " << p << " " << q << " " << Zq << std::endl;
    }

    return Zq;

}
