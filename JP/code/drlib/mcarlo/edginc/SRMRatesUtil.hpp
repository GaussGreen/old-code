//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : SRMRatesUtil.hpp
//
//   Description : base class for SRM IR util classes (calibration etc)
//
//   Author      : Mark A Robson
//
//   Date        : 14 June 2004
//
//
//----------------------------------------------------------------------------

#ifndef SRMRATESUTIL_HPP
#define SRMRATESUTIL_HPP

#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/DECLARE.hpp"

#include "edginc/QMCRatesUtil.hpp"
#include "edginc/SRMRatesFactor.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/VolProcessedBSIR.hpp"
#include "edginc/DoubleMatrix.hpp"

DRLIB_BEGIN_NAMESPACE
class FXAsset;
DECLARE(FXAsset);

/** Helper class for models to get data out of market cache */        
class MCARLO_DLL SRMRatesUtil : public QMCRatesUtil
{

public:
    
    ///// constructs and populates SRMRatesUtil object
    SRMRatesUtil(
            const DateTime&      baseDate,
            int                  numFactors,
            const string&        modelParams,
            const string&        smileParams,
            CVolProcessedSP      processedVol,
            IYieldCurveConstSP   discYC,
            IYieldCurveConstSP   diffYC,
            bool                 skipFlag,
            double               flatVolIr,
            const string&        cutoffChoice,
            double               constCutoffValue,
            const string&        corrSwapStart, // eg 1Y  (offset to today)
            const string&        corrSwapMat,   // eg 10Y (offset to start)
            const string&        corrSwapDCC,   // eg Act/365F
            const string&        corrSwapFreq); // eg 6M
           
    /** Constructs and populates SRMRatesUtil object suitable using the 'VF' style
        calibration. qLeft and qRight are hard coded to 1, as is alpha. Note
        that only single factor IR is supported */
    SRMRatesUtil(
        const DateTime&      baseDate,
        FXAssetConstSP       fx, /* optional - if specified will add FX smile
                                    dates to timeline */
        const string&        modelParamsKey,
        CVolProcessedSP      processedVol,
        IYieldCurveConstSP   discYC,
        IYieldCurveConstSP   diffYC,
        bool                 skipFlag,
        const string&        corrSwapStart, // eg 1Y  (offset to today)
        const string&        corrSwapMat,   // eg 10Y (offset to start)
        const string&        corrSwapDCC,   // eg Act/365F
        const string&        corrSwapFreq); // eg 6M

    virtual ~SRMRatesUtil(){}

    /// Accessors (inline)
    
    /** returns the swaption expiry dates */
    const DateTimeArray& getSwaptionExpiryDates() const{assertInitialized(); return swaptionExpiries; }
    /** returns the swapt start dates */
    const DateTimeArray& getSwapStartDates() const{assertInitialized(); return swapStartDates; }
    /** returns the swapt maturity dates */
    const DateTimeArray& getSwapMatDates() const{assertInitialized(); return swapMatDates; }
    /** Returns the last swaption expiry date */
    const DateTime& getLastExpiry() const{assertInitialized(); return swaptionExpiries.back(); }
    /** Returns the day count convention of the underlying swaps */
    DayCountConventionSP getSwapDCC() const{ return swapDCC; }
    /** Returns the frequency of the underlying swaps (eg 3M, 1A etc) */
    MaturityPeriodSP getSwapFrequency() const{ return swapFrequency; }
    /** returns las */
    const DateTime& getLastDate() const{ return LastDate; }
    /** returns the swaption vols */
    const DoubleArray& getSwaptionVols() const{assertInitialized(); return swaptionVols; }
    const vector<double>& getSwaptionSpotVolsAtSimDates() const {assertInitialized(); return swaptionSpotVolAtSimDates; }
    /** returns the spot vols */
    const vector<double>& getSpotVols() const{assertInitialized(); return SpotVol; }
    /** returns the spot vols as originally calibrated against the swaption
        expiries */
    const DoubleArray& getCalibratedSpotVols() const{assertInitialized(); return swaptionSpotVol; }
	/** return the start date of the swap to use when correlating against other
        IR factors */
    const DateTime& getCorrSwapStart() const{ return corrSwapStart; }
    /** wrapper around extendSpotVol to extend spot vols from this object on
        supplied dates */
    vector<double> extendSpotVol(const DateTimeArray&  newTimePoints) const;

	/** number of factors in model */
	// needs to be virtual in current set up
    virtual int numFactors() const = 0;
    // called when we know simDates
    virtual void setTimeLine(DateTimeArrayConstSP simDates) = 0;
    /** this function returns product of spot vols and forward rates
        "interest rate basis point vol". The array returned is of the same
        length as the supplied array */
    virtual DoubleArray basisPointVol(
        const DateTimeArray& dates) const = 0; /* excludes today */
	/** Returns rho[rhoIndex] where rho are model parameters */
    virtual double getRho(int rhoIndex) const = 0;  // TODO: Move to SRMRatesHJMUtil?
    /** Returns the model rho parameters - length corresponds to off diagonal
        elements of a symmetric matrix (size of which is number of factors) */
    virtual const vector<double>& getRho() const = 0;  // TODO: Move to SRMRatesHJMUtil?
    // must keep in general interface, for now
    /** vol param */
    virtual double getAlpha(int factorIdx) const = 0;  // TODO: Move to SRMRatesHJMUtil?
    virtual void getAlpha(vector<double>& alphaByFactor) const = 0; // TODO: Move to SRMRatesHJMUtil?
    virtual void getBeta(vector<double>& betaByFactor) const = 0; // TODO: Move to SRMRatesHJMUtil?
    /** Returns beta for specified factor */
    virtual double getBeta(int factorIdx) const = 0;  // TODO: Move to SRMRatesHJMUtil?
    /** Calls bFactor for the swap used for correlation purposes */
    virtual vector<double> bFactor() const = 0;
    /** Calls bFactor for the swap given by the specified index */
    virtual vector<double> bFactor(int swaptionIndex) const = 0;
    /*  Determine the value of the B coefficient (see Vladimir's memo)
        this function is Flat Forward ready */
    virtual vector<double> bFactor(const DateTime&      swapStart,
                           const DateTime&      swapMat,
                           DayCountConventionSP swapDayCC,
                           MaturityPeriodSP     swapFreq) const = 0;

    void setMomentMatchingFlag(bool mm) {momentMatching = mm;}
    bool getMomentMatchingFlag() const {return momentMatching;}

	/** vol param */
    //double getQRight() const{ return qRight; }
	virtual double getQRight() const = 0;
    /** vol param */
    //double getQLeft() const{ return qLeft; }
    virtual double getQLeft() const = 0;
    /** vol param */
    //double getFwdShift() const{ return fwdShift; }
    virtual double getFwdShift() const = 0;

    /* Exponential decay function for specified factor */
    virtual double expDecay(int factor, double t) const = 0;  // TODO:  Move to SRMRatesHJMUtil?

    /* Produces the "usual" Lower triangular matrix for orthogonalising
       corrolated factors.  
       NOTE: If Nbfac < 3 then unused matrix elements are set to zero. */
    virtual DoubleMatrix triangulation() const = 0;   // TODO:  Move to SRMRatesHJMUtil?

    /*****  From util_s::Get3A **********************************************/
    /*
     *       Calculate the decomposition of the Correlation in a 2,3-factor case
     *       given the swap start date, a tenor and a frequency
     */
    virtual DoubleMatrix get3A() const = 0; // TODO:  Move to SRMRatesHJMUtil?

	// returns the instantaneous vol by factor in the form (factor, time point)
	virtual void instFactorVol(
		vector< vector<double> >& vol,				  
		const vector<double>& DeltaTime,  // (I)  passed for convenience 
		const DateTimeArray& TPDate,      // (I)
        const vector<double>& IrFwdRate,  // (I)  pre-computed for performance
		int expiryIndex,                  // (I)  index in TPDate
		int fwdMatIndex) const = 0;       // (I)  in case, the underlying has mat > expiry

    // Calculates a 'factor variance' as used by equity/FX vol calibration 
	virtual void irAssetVariance(
		vector<double>& irVariance,        // (O)  variance over [T(i-1),  T(i)] for all i <= N
		vector<double>& irAssetCovar,      // (O)  covariance in [T(i-1),  T(i)] for all i <= N
		const vector<double>& rhoEqIr,     // (I)
		const vector<double>& DeltaTime,   // (I)  passed for convenience 
		const DateTimeArray& TPDate,       // (I)  needed for simpleFwdCurve
        const vector<double>& IrFwdRate,  // (I)  pre-computed for performance
		const vector<int>& periodIndex,    // (I)  indices in TPDate
		int fwdMatIndex) const = 0;        // (I)  indices in TPDate

    /*  Calibration routine : populates ir->SpotVol and ir->tau array */
    virtual void spotVolBM(bool skipFlag) = 0;

    CVolProcessedSP    getProcessedVol(void) { return processedVol; }

protected:

    // methods

    /** Wrapper around extendSpotVol above retrieving ExpDate and TimePt
        from this object */
    void extendSpotVol(const vector<double>& SrcSpotVol);

    //// for ease until we move everything over
    void extendSpotVol(const DoubleArray& SrcSpotVol);

    void spotVol(double flatVolIr);
    
    DateTime  LastDate; // last exposure date for volrate driver

    DateTimeArray   swaptionExpiries; // swaption expiry dates
    DateTimeArray   swapStartDates; // swap start dates
    DateTimeArray   swapMatDates; // swap maturity dates
    
	DoubleArray     swaptionVols; // 
    vector<double>  SpotVol; // LOGNORMAL vols, same size as FwdRate
    vector<double>  swaptionSpotVolAtSimDates; // swaptionExpiries.size() spot vols
    DoubleArray     swaptionSpotVol; // swaptionExpiries.size() spot vols
  
    DayCountConventionSP swapDCC;
    MaturityPeriodSP     swapFrequency;
    DateTime             corrSwapStart;
    DateTime             corrSwapMat;
    DayCountConventionSP corrSwapDCC;
    string               corrSwapFreq;
    mutable DoubleMatrix corrDecomposition; // cached internally

    string               cutoffChoice;
    double               constCutoffValue;
    CVolProcessedSP      processedVol;
    double               flatVolIr;
    bool                 skipFlag;

    bool                momentMatching;
};

DECLARE(SRMRatesUtil);

DRLIB_END_NAMESPACE
#endif
