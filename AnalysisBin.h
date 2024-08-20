struct AnalysisBin {
    int id {};                  // Justin's unique bin id
    std::string regionName {};  // regionName
    int ndata {};               // observed number of data events
    int subBin {};              // (RISR,Mperp) bin number (0-9)
    double postfitMean {};      // Post-fit background mean
    double postfitError {};     // Post-fit background uncertainty
    double satChisq {};         // Barker-Cousins saturated chi-squared (added by GWW)
    bool is0L(){
        std::string s = "0L";
        bool result = false;
        if(regionName.find(s) != std::string::npos)result = true;
        return result;
    }
    bool is1L(){
        std::string s = "1L";
        bool result = false;
        if(regionName.find(s) != std::string::npos)result = true;
        return result;
    }
    bool is2L(){
        std::string s = "2L";
        bool result = false;
        if(regionName.find(s) != std::string::npos)result = true;
        return result;
    }
    bool is3L(){
        std::string s = "3L";
        bool result = false;
        if(regionName.find(s) != std::string::npos)result = true;
        return result;
    }
    bool isBronze(){
        std::string s = "_bron";
        bool result = false;
        if(regionName.find(s) != std::string::npos)result = true;
        return result;
    }
    bool isSV(){
        std::string s = "SVeta";
        bool result = false;
        if(regionName.find(s) != std::string::npos)result = true;
        return result;
    }
    bool is0L0J(){
        std::string s = "Ch0L_0_0j";
        bool result = false;
        if(regionName.find(s) != std::string::npos)result = true;
        return result;    
    }           
        
 // Standard sorting criterion    
    bool operator < (const AnalysisBin & aBin) const
    {
        return satChisq > aBin.satChisq;
    }
};
