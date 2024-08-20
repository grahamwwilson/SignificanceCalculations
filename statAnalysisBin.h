struct statAnalysisBin : AnalysisBin{
    double zscore {};
    double zscoreError {};
    bool operator < (const statAnalysisBin & aBin) const
    {
         return abs(zscore) > abs(aBin.zscore);
    }         
};
