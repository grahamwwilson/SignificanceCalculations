struct statAnalysisBin : AnalysisBin{
    double zscore {};        // symmetrized Z-score
    double zscoreError {};   // symmetrized Z-score error
//    double satChisqTotal {};                   // for parsing test
    bool operator < (const statAnalysisBin & aBin) const
    {
         return abs(zscore) > abs(aBin.zscore); // Order by descending |Z|
//       return abs(zscore) < abs(aBin.zscore); // Order by ascending  |Z|

    }         
};
