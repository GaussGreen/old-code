//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : QMCCreditCIDJumps.hpp
//
//   Description : An interface of internal hooks for pure jumps asset and
//                 diffused data retrieval
//
//
//----------------------------------------------------------------------------

#ifndef QLIB_SRMCREDITCIDJUMPS_HPP
#define QLIB_SRMCREDITCIDJUMPS_HPP

#include "edginc/QMCPureJumps.hpp"
#include "edginc/IQMCDiffusibleAsset.hpp"
#include <cmath>

DRLIB_BEGIN_NAMESPACE

// the idea is -- this asset encapsulates the response of the name to
// a particular source of jumps
// the result should not be a list of jump times, but rather a full-fledged
// contribution of this source to the expected survival, etc...

class MCARLO_DLL QMCCreditCIDJumps : public IQMCDiffusibleCreditSpread
{
public:


    /** this is for efficient access to SDF values*/
    virtual double* getInternalSDFArray() { return NULL; }

    /** return underlying IR asset */
    virtual IQMCDiffusibleInterestRateSP getUnderlyingIRAsset()
    { return IQMCDiffusibleInterestRateSP(); }


    /** get the simulation start date */
    virtual DateTime getBaseDate() { return today; }

    /** generates path across all dates */
    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/);

    /** finalize the timelines, allocate necessary memory */
    virtual void finalize(DateTimeArrayConstSP allDates);

    /** returns pointer to internal storage for MaxDiff(usion)Date / MaxCurveMaturity */
    virtual QMCHelperBoundedDiffusion * getDiffusionBound()  { return &diffBound; }


    /** accessing the diffused SDF, i.e. a non-default probability by a given date */
    virtual double getSurvivalDiscFactor(SpotIdx measurementDateIdx) ;
    virtual double getRecoveryRate(SpotIdx idx) { return 0.0; } // might throw exception here -- as non-meaningful
    virtual double getExpRecoveryRate(
                FwdIdx measurementDateIdx,
                FwdIdx futureDateIdx) { return 0.0; } // might throw exception here -- as non-meaningful


//    /** accessing the semi-analytical Survival Probability, for Fast MC */
//    virtual double getConditionalSurvivalDiscFactor(SpotIdx measurementDateIdx);

    /** from the jump dates and jump sizes, compute at accessing the diffused SDF, i.e. a non-default probability by a given date */
    virtual double getWholeTimelineLogSurvProb(size_t i);


    // DanielNg: need to check the reference date
    /** returns log of deterministic surv prob */
    virtual double getOriginalLnSurvivalDiscFactor(SpotIdx measurementDateIdx);

    /** returns log of deterministic forward surv prob */
    virtual double getOriginalLnExpectedSurvivalDiscFactor(FwdIdx measurementDateIdx, FwdIdx futureDateIdx);

    /** after the path is generated, access the jump sizes*/
	virtual double getJumpSizeByIdx(size_t sourceIdx, size_t jumpIdx);

    /** Accessing the expected value ExpSDF(md, fd) where md is a
        simulated measurement date and fd is some future date after the
        measurement is made. */
    virtual double  getExpectedSurvivalDiscFactor(
                                        FwdIdx measurementDateIdx,
                                        FwdIdx futureDateIdx)
    { return std::exp(QMCCreditCIDJumps::getLnExpectedSurvivalDiscFactor(measurementDateIdx, futureDateIdx));}

    /** Accessing the natural log of the expected value ExpSDF(md, fd)
        where md is a simulated measurement date and fd is some future
        date after the measurement is made. */
    virtual double getLnExpectedSurvivalDiscFactor(
                                        FwdIdx measurementDateIdx,
                                        FwdIdx futureDateIdx);



    QMCCreditCIDJumps(IQMCPureJumpsSP _jumpGenerator, bool isFullMC) : 
        jumpGenerator(_jumpGenerator), 
        diffBound(_jumpGenerator->getDiffusionBound()),
        fullMC(isFullMC) {}
    ~QMCCreditCIDJumps() {}

    void setQMCCreditCIDJumps(
        DateTime        _today,
        const vector<double>& _jumpImpact,
        const vector<double>& _jumpMeanSize,
        const vector<double>& _jumpMeanReversionSpeed);

private:
/** given a uniform (random) variable unif, give the product of a Bernoulli(impact) times
    an Exponential with mean meanSize */
	double getJumpSizeSample(double impact, double meanSize, double unif);

/** return E(exp(-\int_t^T lambda_u du / N^{market} ) 
To obtain the conditional on all the Poisson processes survival probability,
do survProba_{t,T} = prod_{market=1}^M getConditionSurvProba(t,T,market,lambda_t^{market}) */
	double getConditionSurvProba( double t, 
								double bigT, 
								size_t market,
								double lambda_t);

    double calculateConditionalSurvivalDiscFactor(double T);
	
	DateTime today;
    IQMCPureJumpsSP  jumpGenerator;
    QMCHelperBoundedDiffusion diffBound;
    //IQMCDiffusibleInterestRateSP irAsset;

    vector<double> jumpImpact;
    vector<double> jumpMeanSize;
    vector<double> jumpMeanReversionSpeed;

    vector<vector<double> > lambda_tstar; // lambda_t[jumpSrc][fwdIdx] for all {t*}
    vector<double>  survivalProb; // sp[t] for all {t}

	vector<double> sdfRequestedDatesAsDouble; // (t)
	vector<double> esdfForwardDatesAsDouble; // (t*+T) 
	vector<double> esdfRequestedDatesAsDouble;// (t*)

    vector<size_t> offsets;
	vector<double> jumpSizes;
    vector<double> allDatesAsDouble;

    bool fullMC, STOREJUMPS;
};

/** a hookup for across-the-assets diffusion collector */
DECLARE(QMCCreditCIDJumps);


DRLIB_END_NAMESPACE
#endif //QLIB_SRMCREDITCIDJUMPS_HPP

