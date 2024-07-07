

std::tuple<int, double, int, double, double, double, double, double, double, double, double, double> 
   ToyChucker(int NFACTOR, int NOBSMIN, int NOBSMAX, double MUB, double FRACERROR, unsigned long int seed, int id ){

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
        vgen.push_back(n);
        if(i<NTOSAVE){
// Save toy values for later synthetic data-sets for this bin
            std::tuple<int, unsigned long int, int, int> t = std::make_tuple(i, seed, id, n);
            vtoys.push_back(t);
        }
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
        cout << "Upper-tail " << fixed << setprecision(10) << vpvalue[j].first << " +- " << fixed << setprecision(10) << vpvalue[j].second << endl;
        cout << "Lower-tail " << fixed << setprecision(10) << vqvalue[j].first << " +- " << fixed << setprecision(10) << vqvalue[j].second << endl;
        cout << "Equal      " << fixed << setprecision(10) << vpequal[j].first << " +- " << fixed << setprecision(10) << vpequal[j].second << endl; 
        
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
   
   std::tuple<int, double, int, double, double, double, double, double, double, double, double, double> t;
   
   for (auto & el : vtup){
        t = el;
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
        
/*        hist0->Fill(NOBS,zscoreu);
        hist0->SetBinError(ibin, zscoreu_ERR); 
        hist0N->Fill(NOBS,zscorel);
        hist0N->SetBinError(ibin, zscorel_ERR);
        hist0M->Fill(NOBS,0.5*(zscorel + zscoreu));
        hist0M->SetBinError(ibin, 0.5*(zscorel_ERR + zscoreu_ERR));
        hist0MH->Fill(NOBS,zscorem);
        hist0MH->SetBinError(ibin, zscorem_ERR);             */
        
// Algorithm 1
        if(NOBS >= MUB){
            z1 = zscoreu;
//            hist1->Fill(NOBS,zscoreu);
//            hist1->SetBinError(ibin, zscoreu_ERR);
        }
        else{
            z1 = zscorel;
//            hist1->Fill(NOBS,zscorel);
//            hist1->SetBinError(ibin, zscorel_ERR);
        }
// Algorithm 2
        if(p <= q){
            z2 = zscoreu;
//            hist2->Fill(NOBS,zscoreu);
//            hist2->SetBinError(ibin, zscoreu_ERR);            
        }
        else{
            z2 = zscorel;
//            hist2->Fill(NOBS,zscorel);
//            hist2->SetBinError(ibin, zscorel_ERR); 
        }
// Algorithm 3/4
        if(p <= q){     // Upper-tail probability is smallest
             if(p<0.5){
                 z3 = zscoreu;
//                 z4 = zscoreu;
//                 hist3->Fill(NOBS, zscoreu);
//                 hist3->SetBinError(ibin, zscoreu_ERR);
//                 hist4->Fill(NOBS, zscoreu);
//                 hist4->SetBinError(ibin, zscoreu_ERR);                   
             }
             else{       // Both p and q are > 0.5
                 z3 = 0.5*(zscoreu+zscorel);
//                 hist3->Fill(NOBS,0.5*(zscoreu + zscorel));
//                 hist3->SetBinError(ibin, 0.5*(zscorel_ERR + zscoreu_ERR)); 
             }
        }
        else{
             if(q<0.5){
                 z3 = zscorel;
//                 hist3->Fill(NOBS,zscorel);
//                 hist3->SetBinError(ibin, zscorel_ERR);                  
             }
             else{
                 z3 = 0.5*(zscoreu + zscorel);
//                 hist3->Fill(NOBS,0.5*(zscoreu + zscorel));
//                 hist3->SetBinError(ibin, 0.5*(zscorel_ERR + zscoreu_ERR));                  
             }               
        }
        cout << "NOBS = " << NOBS << " Alg 0,1,2,3 z-scores " << z0 << " " << z1 << " " << z2 << " " << z3 << endl;              
   }
   return t;
}


