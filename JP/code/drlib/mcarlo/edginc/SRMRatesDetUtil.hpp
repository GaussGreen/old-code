#ifndef	SRMRatesDetermUtil_HPP
#define	SRMRatesDetermUtil_HPP

#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/SRMRatesUtil.hpp"
#include "edginc/SRMRatesFactor.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/VolProcessedBSIR.hpp"
#include "edginc/DoubleMatrix.hpp"


DRLIB_BEGIN_NAMESPACE


class SRMRatesDetermUtil : public SRMRatesUtil
{
public:
    SRMRatesDetermUtil(
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
    SRMRatesDetermUtil(
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


    virtual ~SRMRatesDetermUtil(void);

    virtual int numFactors() const
    {
        return 0;
    }

    // called when we know simDates
    virtual void setTimeLine(DateTimeArrayConstSP simDates);

    /** this function returns product of spot vols and forward rates
    "interest rate basis point vol". The array returned is of the same
    length as the supplied array */
    DoubleArray basisPointVol(
        const DateTimeArray& dates) const; /* excludes today */

    /** Returns rho[rhoIndex] where rho are model parameters */
    virtual double getRho(int rhoIndex) const;
    /** Returns the model rho parameters - length corresponds to off diagonal
    elements of a symmetric matrix (size of which is number of factors) */
    virtual const vector<double>& getRho() const;
    // must keep in general interface, for now
    /** vol param */
    virtual double getAlpha(int factorIdx) const;
    virtual void getAlpha(vector<double>& alphaByFactor) const;
    virtual void getBeta(vector<double>& betaByFactor) const;
    /** Returns beta for specified factor */
    virtual double getBeta(int factorIdx) const;

    /** Calls bFactor for the swap used for correlation purposes */
    virtual vector<double> bFactor() const;
    /** Calls bFactor for the swap given by the specified index */
    virtual vector<double> bFactor(int swaptionIndex) const;
    /*  Determine the value of the B coefficient (see Vladimir's memo)
    this function is Flat Forward ready */
    virtual vector<double> bFactor(
        const DateTime&      swapStart,
        const DateTime&      swapMat,
        DayCountConventionSP swapDayCC,
        MaturityPeriodSP     swapFreq) const;

    /** vol param */
    double getQRight() const{ return qRight; }
    /** vol param */
    double getQLeft() const{ return qLeft; }
    /** vol param */
    double getFwdShift() const{ return fwdShift; }

    /* Exponential decay function for specified factor */
    virtual double expDecay(int factor, double t) const { return 1.;}

    /* Produces the "usual" Lower triangular matrix for orthogonalising
    corrolated factors.  
    NOTE: If Nbfac < 3 then unused matrix elements are set to zero. */
    virtual DoubleMatrix triangulation() const;


    virtual DoubleMatrix get3A() const;

    // returns the instantaneous vol by factor in the form (factor, time point)
    virtual void instFactorVol(
        vector< vector<double> >& vol,				  
        const vector<double>& DeltaTime,  // (I)  passed for convenience 
        const DateTimeArray& TPDate,      // (I)
        const vector<double>& IrFwdRate,  // (I)  pre-computed for performance
        int expiryIndex,                  // (I)  index in TPDate
        int fwdMatIndex) const;           // (I)  in case, the underlying has mat > expiry

    // Calculates a 'factor variance' as used by equity/FX vol calibration 
    virtual void irAssetVariance(
        vector<double>& irVariance,        // (O)  variance over [T(i-1),  T(i)] for all i <= N
        vector<double>& irAssetCovar,      // (O)  covariance in [T(i-1),  T(i)] for all i <= N
        const vector<double> &rhoEqIr,            // (I)
        const vector<double>& DeltaTime,   // (I)  passed for convenience 
        const DateTimeArray& TPDate,       // (I)  needed for simpleFwdCurve
        const vector<double>& IrFwdRate,  // (I)  pre-computed for performance
        const vector<int>& periodIndex,           // (I)  indices in TPDate
        int fwdMatIndex) const;            // (I)  indices in TPDate

    /*  Calibration routine : populates ir->SpotVol and ir->tau array */
    virtual void spotVolBM(bool skipFlag);

private:
    double      qLeft,
                qRight,
                fwdShift;

};

DECLARE(SRMRatesDetermUtil);

DRLIB_END_NAMESPACE

#endif

