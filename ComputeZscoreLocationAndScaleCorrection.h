
std::pair<double,double> Compute(double mu, int kmin, int kmax, int printlevel=0){

// Input quantities
// mu:   Chosen Poisson mean
// kmin: Minimum k
// kmax: Maximum k
//
// For non-negative integers k, satisfying kmin <= k <= kmax, 
// evaluate the Poisson probability pk = Poisson(k; mu)
// and compute the mean value and standard deviation for the symmetric Z-score variable 
// through summation

   double CUT = 8.29; // Potentially can do this differently for +-ve and -ve Z.
   if (printlevel!=0) std::cout << "mu " << mu << " kmin " << kmin << " kmax " << kmax << " CUT " << CUT << std::endl;
   
   double emu  = TMath::Exp(-mu);
   double bess = TMath::BesselI(0,2.0*mu);
   
// Historical debugging check? not used in current calculation.
//   std::cout << "I0 (2 mu) " << bess << std::endl;

// Apparently Sum(k=0, infty) pk**2 = Exp(-2 mu)*I0(2mu).

   std::vector <std::pair<int,double>> vp;  //vector of pairs with (k, pk) where p(k) = Prob(k; mu)
   std::vector <std::pair<int,double>> vQ;  //vector of pairs with (k, ps) where ps = 1 - [Sum (k=0, k-1) {pk} + (p(k)/2) ]. It is called Q 
                                            //as it uses a formulation that is essentially p_S = 1.0 - q_S.
   std::vector <std::pair<int,double>> vP;  //vector of pairs with (k, ps) where ps = Sum(nobs+1, infinity) {pk} + (p(nobs)/2)
   
   double sumps=0.0;
   double sumpsps = 0.0;
   double sumz = 0.0;
   double sumzz = 0.0;
   
   double partialsum = 0.0;
   double qpartialsum = 0.0;
   
   double pk;
   
// Check elements of the initialization calculation. 
// Looks like emu can underflow for big enough mu. For example mu=300, gives 5.14e-131. 
// Apparently smallest positive double is 2.225e-307.
   if(printlevel!=0){
       std::cout << " Exponential factor " << emu << std::endl;
       std::cout << " pow(mu,k)   factor " << std::pow(mu,kmin) << std::endl;
       std::cout << " Gamma(k+1)  factor " << TMath::Gamma(double(kmin+1)) << std::endl;
   }
   
   for (int k=kmin; k<=kmax; k++){
       if(k == kmin){
           pk = emu*std::pow(mu,k)/TMath::Gamma(double(k+1));
       }
       else{
           pk = vp.back().second*(mu/double(k)); // Use recursive definition
       }

// p_{S} calculation. This runs into floating point issues at around 2e-16 if first few terms are of order one.
// May be better to switch to a tail calculation whenever eg p(k) < 1.0e-12;
       double ps = 1.0;
       for (auto & el : vp){
           ps -= el.second;
       }
       ps -= 0.5*pk;
       
       auto pairk = std::make_pair(k, pk); vp.push_back(pairk); // Only push back new value AFTER ps calculation
       auto pairQ = std::make_pair(k, ps); vQ.push_back(pairQ);       
   }
   
// Now calculate tail probability directly using whatever range of (kmin, kmax) we decided to use
    
   for (int j=kmin; j<=kmax; j++){
       double P = 0.0;
       for (auto & el : vp){
            auto k  = el.first;
            auto pk = el.second;
            if( k==j) P += 0.5*pk;
            if( k>j) P += pk;
       }
       auto pairP = std::make_pair(j, P); vP.push_back(pairP);       
   }
   
// Now loop through all possible experimental outcomes  
   
   for (int j=kmin; j<=kmax; j++){   
      
       int idx = j-kmin;
       auto k   = vp[idx].first;
       auto pk  = vp[idx].second;  // Pr(k)
       auto psQ = vQ[idx].second;  // Symmetrized upper tail probability using Q based approach 
       auto psP = vP[idx].second;  // Symmetrized upper tail probability using P based approach - should be the same as psQ. 
      
       double zQ = MyQuantile(psQ);
       double zP = MyQuantile(psP);
       
// Also keep track of zl and zu
       double pU = psQ + 0.5*pk;
       double pL = psQ - 0.5*pk;
       double zu = MyQuantile(pU);
       double zl = MyQuantile(pL); 

       double z;

// Choose representation with better numerical precision
       if(psQ > 0.5){
// cumulative subtraction for large probabilities
          z = zQ;
       }
       else{
// explicit tail summation for small probabilities
          z = zP;
       }
       
       double psum = psP + 0.5*pk;
       
       if(printlevel!=0){
           std::cout << "k: " << k << " p(k) " << pk << " psQ " << psQ << " psP " << psP 
                     << " ZQ " << zQ << " ZP " << zP << " Chosen z " << z 
                     << " qsum " << std::scientific << std::setprecision(16) << 1.0-psum << std::endl;
       }

       if(printlevel!=0){
           std::cout << "k: " << k << " pL " << pL << " zL  " << zl << std::endl;
           std::cout << "k: " << k << " pU " << pU << " zU  " << zu << std::endl;
       }
       std::cout << " " << std::endl;
       
       sumps   += pk*psQ;
       sumpsps += pk*psQ*psQ;

       double zthis = z;

       if(zu >= 2.5){
// Significant excesses
// switch to zupper
           zthis = zu;
       }
       else if(zl <= -2.5){
// Significant deficits
// switch to zlower
           zthis = zl;
       }
       else if(std::abs(z) < 1.5){
// Plan here for bias and scale correction
           zthis = z;
       }
       else{
// Plan just for bias correction
           zthis = z;
       }

       if( abs(z) < CUT && abs(zu) < CUT && abs(zl) < CUT ) {
          sumz  += pk*zthis;
          sumzz += pk*zthis*zthis;
          partialsum += pk;
       }
       else{
// Also keep track of fraction that is not included even in the [kmin, kmax] range.
          qpartialsum += pk;
       }   

   }
   
   double varx = sumpsps - sumps*sumps;
   double varz = sumzz - sumz*sumz;

   if( printlevel!=0) {

       std::cout << "mu, kmin, kmax " << mu << " " << kmin << " " << kmax << std::endl;
       std::cout << "partial sum = " << std::scientific << std::setprecision(16) << partialsum << std::endl;
       std::cout << "omitted part  " << std::scientific << std::setprecision(16) << qpartialsum << std::endl;
     
       std::cout << " For x = p_{S} " << std::endl;  
       std::cout << " E(x)   = " << std::fixed << std::setprecision(16) << sumps << std::endl;
       std::cout << " E(x^2) = " << std::fixed << std::setprecision(16) << sumpsps << std::endl;
       std::cout << " Var(x) = " << std::fixed << std::setprecision(16) << sumpsps - sumps*sumps << std::endl;
       std::cout << " SD(x)  = " << std::fixed << std::setprecision(16) << sqrt(varx) << std::endl;
       std::cout << " NSD(x) = " << std::fixed << std::setprecision(16) << sqrt(12.0*varx) << std::endl;
   
       std::cout << " ---------------------------" << std::endl;
   
       std::cout << " For x = Z_{S} " << std::endl;  
       std::cout << " E(x)   = " << std::fixed << std::setprecision(16) << sumz << std::endl;
       std::cout << " E(x^2) = " << std::fixed << std::setprecision(16) << sumzz << std::endl;
       std::cout << " Var(x) = " << std::fixed << std::setprecision(16) << varz << std::endl;
       std::cout << " SD(x)  = " << std::fixed << std::setprecision(16) << sqrt(varz) << std::endl;
//       std::cout << " " << std::endl;

   }
   
// Computes the mean of Z and its standard deviation.
   auto result = std::make_pair(sumz, sqrt(varz));

   return result;
      
}

double ComputeNaivePull(int nObs, double mu){
//
// Naive pull of data observation, nObs, with respect to background-only expectation and the 
// Poisson statistical uncertainty, sqrt(mu), based on this expectation.
// 
// No accounting of post-fit uncertainty here. So assumes model is absolutely correct.
// Should only give Gaussian like results in the asymptotic Gaussian regime.
//
    double data = double(nObs);
    double zPull = (data - mu)/std::sqrt(mu);

    return zPull;

}

std::pair<double,double> ComputePoissonOnlyStats(int nObs, double mu){
//
// No attempt at symmetrization. Just sum up the Poisson probability for n>=nObs and invert.
//

    double pk;
    std::vector <std::pair<int,double>> vp;  //vector of pairs with (k, pk) where p(k) = Prob(k; mu)
    double emu  = TMath::Exp(-mu);

    double Qtot = 0.0;
    for (int k=0; k<nObs; k++){  //Here don't include the observation
        if(k == 0){
            pk = emu*std::pow(mu,k)/TMath::Gamma(double(k+1));
        }
        else{
            pk = vp.back().second*(mu/double(k)); // Use recursive definition
        }
        auto pairk = std::make_pair(k, pk); vp.push_back(pairk); // Only push back new value AFTER ps calculation
        Qtot += pk;
    }
    double Ptot = 1.0 - Qtot;
// Compute corresponding quantile
    double Z = MyQuantile(Ptot, 1);

    auto result = std::make_pair(Ptot, Z);

    return result;

}

std::pair<double,double> ComputeZscoreLocationAndScaleCorrection(double mu, int printlevel=0){

// Driver function for Zscore bias computation used for standardization of Zscore statisic

// Set default value and declare result
    double zbias = 0.0;
    double zscalebias = 1.0;
    auto result = std::make_pair(zbias, zscalebias);

    if ( mu < 0.075 ){
        std::cout << "Zscore scale corrections for very small mu are LARGE." 
                  << " Using mu=0.075 for corrections, while in this case mu = "
                  << mu << std::endl;
        int kmin = 0;
        int kmax = 20;
        result = Compute(0.075, kmin, kmax, printlevel);
    }
    else if ( mu < 1.0 ){
        int kmin = 0;
        int kmax = 20;
        result = Compute(mu, kmin, kmax, printlevel);
    }
    else if ( mu >= 1.0 && mu < 5.0 ){
        int kmin = 0;
        int kmax = 50;
        result = Compute(mu, kmin, kmax, printlevel);
    }
    else if ( mu >= 5.0 && mu < 10.0 ){
        int kmin = 0;
        int kmax = 80;
        result = Compute(mu, kmin, kmax, printlevel);
    }
    else if ( mu >= 10.0 && mu < 25.0 ){
        int kmin = 0;
        int kmax = 100;
        result = Compute(mu, kmin, kmax, printlevel);
    }
    else if ( mu >= 25.0 && mu < 100.0 ){
        int kmin = 0;
        int kmax = 250;
        result = Compute(mu, kmin, kmax, printlevel);
    }
    else if ( mu >= 100.0 && mu < 300.0 ){
        int kmin = 0;
        int kmax = 1000;
        result = Compute(mu, kmin, kmax, printlevel);
    }
    else{
        std::cout << "Zscore corrections for large mu are tiny and neglected."
                  << " Setting (bias, scale) corrections to (0.0, 1.0) for mu = "
                  << mu << std::endl;
    }

    return result;

}
