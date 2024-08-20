struct statAnalysisBin : AnalysisBin{
    double zscore {};
    double zscoreError {};
//    double satChisqTotal {};                   // for parsing test
    bool operator < (const statAnalysisBin & aBin) const
    {
         return abs(zscore) > abs(aBin.zscore);
    }         
};
