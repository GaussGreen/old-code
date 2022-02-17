//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : SRMCreditDiffuse.hpp
//
//   Description : Base class for credit path generation
//
//
//----------------------------------------------------------------------------

#ifndef SRMCREDITDIFFUSE_HPP
#define SRMCREDITDIFFUSE_HPP

#include "edginc/QMCCreditDiffuse.hpp"

DRLIB_BEGIN_NAMESPACE

/** base class for SRM credit path generation */
class SRMCreditDiffuse : public QMCCreditDiffuse  {
public:

    /** constructor */
    SRMCreditDiffuse(IQMCDiffusibleInterestRateSP pRatesAsset);  // the underlying rates

    /** destructor */
    virtual ~SRMCreditDiffuse();

    /** returns pointer to internal storage for MaxDiff(usion)Date / MaxCurveMaturity */
    virtual QMCHelperBoundedDiffusion* getDiffusionBound();

    /** return underlying IR asset */
    virtual IQMCDiffusibleInterestRateSP getUnderlyingIRAsset();

    virtual double getWholeTimelineLogSurvProb(size_t i);
    /** add all dates that will be needed later, extends base class definition */
    virtual void addAggregatedDates(
        DateTimeArrayConstSP _sdfDates,
        DateTimeArrayConstSP _esdfRequestedDates,
        DateTimeArrayConstSP _esdfForwardDates,
        const DateTime& maxDiffDate,
        const DateTime& maxCurveMat);


protected:

    IQMCDiffusibleInterestRateSP irAsset;
    const vector<double>*   sigmaFX;      // FX vols accross dates, passed in when created in MCPathConfigSRM
    int                     randomIndex;  // index into random numbers

    double                  crFxCorr;     // correlation(CR,FX), always one single nb

    QMCCRBoundedDiffusion diffBound;

    void        calcRemapToIRAssetIdx(const DateTimeArray& esdfForwardDates);

    vector<double>  wholeTimelineLogSurvProb;

    /** returns IR forward index */
    FwdIdx          getIRfwdIdx(FwdIdx idx) const;

private:
    vector<int>     fwdIdx2IRedfIdx;    // projection onto EDFs (expected surv. factors are subset of IR'sexpDFs)

};

/** declare smartPtr versions */
DECLARE(SRMCreditDiffuse);

DRLIB_END_NAMESPACE

#endif // SRMCREDITDIFFUSE_HPP
