//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMEquityDiffuseStrictSecondOrder.hpp
//
//   Description : A generator of paths using stochastic rates
//                 for Equity Assets and the 'second order method'
//                 for the SDE approximation, using a smaller number
//                 of time steps.
//
//   Date        : Aug 2006
//
//----------------------------------------------------------------------------

#ifndef EDR_SRMEQUITYDIFFUSESTRICTSECONDORDER_HPP
#define EDR_SRMEQUITYDIFFUSESTRICTSECONDORDER_HPP

DRLIB_BEGIN_NAMESPACE  

class SRMEquityDiffuseStrictSecondOrder : public SRMEquityDiffuseSecondOrder
{
    vector<int> bigStepIdxs;
    vector<double> bigStepYearFracs;
    vector<double> bigStepYearFracsSqrt;
    vector<double> cfactorBigIntegrals;

    //Used for sorting in finalize()
    struct SimDateInfo
    {
        int simDateIdx;
        int weight;     //distance to bigStep
        int upperDistance;
        int lowerDistance;

        bool operator<(const SimDateInfo &other) const { return (weight < other.weight); }
    };

    struct GapInfo
    {
        int gapID;
        double weight;
        bool operator<(const GapInfo &other) const { return (weight > other.weight); }
    };

public:

    SRMEquityDiffuseStrictSecondOrder(IQMCDiffusibleInterestRateSP domIR) : SRMEquityDiffuseSecondOrder(domIR) {}

    virtual ~SRMEquityDiffuseStrictSecondOrder() {}

    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/);
    virtual void finalize(DateTimeArrayConstSP allDates);

private:
    double getGapWeight(size_t nBegin, size_t nEnd, const vector<double> &volSize);

    size_t addTimePointsRoundDown(size_t totalGapSize, int toAllocate, const vector<int> &initialBigStepLocations, const vector<double> &volSize, vector<int> &bigStepBools);
    
    //The gap size will be weighted by volSize
    size_t addTimePointsDivideLargestGap(int toAllocate, const vector<double> &volSize, vector<int> &bigStepBools);

    size_t addTimePointsInOrder(int toAllocate, vector<int> &bigStepBools);
    
    size_t addTimePointsFillLargeGaps(double maxGapYearFrac, vector<int> &bigStepBools);

};

DRLIB_END_NAMESPACE  

#endif
