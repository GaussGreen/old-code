//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMEquityDiffuseMappingMethod.hpp
//
//   Description : A generator of paths using stochastic rates
//                 for Equity Assets and the 'mapping method' with
//                 predictor-corrector for the SDE approximation
//
//   Date        : Aug 2006
//
//
//----------------------------------------------------------------------------

#ifndef EDR_SRMEQUITYDIFFUSEMM_HPP
#define EDR_SRMEQUITYDIFFUSEMM_HPP

DRLIB_BEGIN_NAMESPACE  

class SRMEquityDiffuseMappingMethod : public SRMEquityDiffuse
{
public:

    SRMEquityDiffuseMappingMethod(IQMCDiffusibleInterestRateSP domIR) : SRMEquityDiffuse(domIR) {}
    virtual ~SRMEquityDiffuseMappingMethod() {}

    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/);
    virtual void finalize(DateTimeArrayConstSP allDates);

private:

    vector<int> bigStepIdxs;
    vector<double> bigStepYearFracs;
    vector<double> bigStepYearFracsSqrt;

    vector<double>  cumulativeLnDivYield;   //sum of lnDivYield up to given (ex-)date.
    vector<double>  cumulativeCfactorArrayEQ; //integral of cts divs, borrow
    vector<double>  cumulativeSpotVol; //integral of spot EQ vol ^ 2 to time t

    vector<double>  cumulativeSpotVol2; //integral of spot EQ vol ^ 2 to time t
    vector<double> firstMomentSpotVol2;

    /* Integrals of 't' times various quantities needed for the drift calculation: */
    vector<double> firstMomentLnDiv;    
    vector<double> firstMomentCFactor;
    vector<double> firstMomentSpotVol;
    vector<double> firstMomentLnFwdEq;  //first moment of derivative of log of forward of equity

    double computeMappingMethodDrift(int simDateIdxA, int simDateIdxB, double smileEQ, double smileDashEQ, const vector<double> &cumulativeSigmaFX, bool isFirstDateIdx) const;
    double computeCorrectedMappingMethodDrift(int simDateIdxA, int simDateIdxB, double smileEQ, double smileDashEQ, double smileEQest, double smileDashEQest, double timeSoFar, double stepYearFrac, const vector<double> &cumulativeSigmaFX, const vector<double> &firstMomentSigmaFX, const vector<double> &firstMomentRt, bool isFirstDateIdx) const;
    double computeCorrectedMappingMethodDrift2(int simDateIdxA, int simDateIdxB, double smileEQ, double smileDashEQ, double smileEQest, double smileDashEQest, double timeSoFar, double stepYearFrac, const vector<double> &cumulativeSigmaFX, const vector<double> &firstMomentSigmaFX, const vector<double> &firstMomentRt, bool isFirstDateIdx) const;

    double computeMappingMethodDrift3(int bigDateIdx, double smileEQ, double smileDashEQ) const;
    double computeCorrectedMappingMethodDrift3(int bigDateIdx, double smileEQ, double smileDashEQ, double smileEQest, double smileDashEQest, double timeSoFar, double stepYearFrac) const;

    void calculateBigStepIntegral(const vector<double> &f, bool square, const vector<int> &bigStepBools, vector<double> &bsIntegral);
    void calculateBigStepFirstMoment(const vector<double> &f, bool square, const vector<int> &bigStepBools, vector<double> &bsFirstMoment);
};

DRLIB_END_NAMESPACE  

#endif
