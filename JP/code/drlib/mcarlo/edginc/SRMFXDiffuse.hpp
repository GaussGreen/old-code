//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMFXDiffuse.hpp
//
//   Description : A generator of paths using stochastic rates
//                 for FX Assets
//
//   Author      : Mark A Robson
//
//   Date        : 13 Aug 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_SRMFX_HPP
#define EDR_SRMFX_HPP

#include <set>

#include "edginc/DECLARE.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/SVQmcImplemented.hpp"
#include "edginc/SRMFXUtil.hpp"
#include "edginc/QMCRatesDiffuse.hpp"
#include "edginc/QMCFXBaseDiffuse.hpp"
#include "edginc/QMCHelperBoundedDiffusion.hpp"


DRLIB_BEGIN_NAMESPACE


/** A generator of paths using stochastic rates for FX Assets. The
    methodology forces us to diffuse the fx asset together with the foreign
    IR. Note that all fx assets (at least as far as the diffusion is concerned)
    must always have the same domestic IR. */
class MCARLO_DLL SRMFXDiffuse : public QMCFXBaseDiffuse
{

public:

    static const string CONSTANT_SPOT_VOL;
    static const string NO_FAILURE_ALLOWED;
    static const string USE_LAST_LEVEL;

    ~SRMFXDiffuse();

    SRMFXDiffuse(IQMCDiffusibleInterestRateSP _domesticIRAsset,
                 IQMCDiffusibleInterestRateSP _foreignIRAsset);

    void setSRMFX(
          int                   randomIndex,
          const DateTime&       today,
          SRMFXUtilSP           util,
          const vector<double>& spotFX);      // historic spot FX

    /** Prepare to simulate a new path. Returns the sigmaFX for the first
        step */
    double begin(IQMCRNGManagerSP rngMgr);

    /** Diffuses the fx to the next date. Returns the sigmaFX for the next
        date. idx will be 0 to move from the first date to the second. */
    double moveToNextDate(double forLnMONEY,
                          int    simDateIdx);

    /** [Must be] called at the end of the path */
    void end();

    /** Gives low level access to path */
    vector<double>::const_iterator getSpotPath() const;
    /** Gives low level access to path */
    vector<double>::const_iterator getExpSpotPath() const;

    /** Allows quanto equity adjustment */
    // used only within SRM
    virtual const vector<double>& getSigmaFX() const;

    /** Allows struck equity adjustment: calling this may impact performance slightly. */
    // used only within SRM
    virtual const vector<double>& requestFullSpotPath();

    //////////////////////////// IQMC interface /////////////////////////

    virtual IQMCDiffusibleInterestRateSP getDomesticIRAsset() { return domesticIRAsset; }
    virtual IQMCDiffusibleInterestRateSP getForeignIRAsset()  { return foreignIRAsset; }

    /** generate path across all dates. */
    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/);
    /** finalize the timelines, allocate necessary memory */
    virtual void finalize(DateTimeArrayConstSP allDates);
    /** getting the simulation start date */
    virtual DateTime getBaseDate() { return today; }
    /** Accessing the diffused SpotFX, i.e. a non-default probability by
        a given date */
    virtual double  getSpotFX(SpotIdx measurementDateIdx) {
        return spotFX[measurementDateIdx];
    }

    /** Accessing the ExpectedFX(md, fd) where md is a
        simulated measurement date and fd is some future date after the
        measurement is made. */
    virtual double  getExpectedFX(
                                FwdIdx measurementDateIdx,
                                FwdIdx futureDateIdx) {
        return exp(SRMFXDiffuse::getLnExpectedFX(measurementDateIdx, futureDateIdx));
    }

    /** Accessing the natural log of the ExpectedFX(md, fd)
        where md is a simulated measurement date and fd is some future
        date after the measurement is made. */
    virtual double getLnExpectedFX(
                                FwdIdx measurementDateIdx,
                                FwdIdx futureDateIdx);

    /** addAggregatedDates overrides the default implementation to consistently add dates to FX and dependent IRs */
    void   addAggregatedDates(
            const DateTimeArray&    spot,
            const DateTimeArray&    fwdDates,
            const DateTimeArray&    fwdFwdDates);
//------------------------------------------------------------------------------
    /** this is for efficient access to SpotFX values, but optional */
    virtual double* getInternalSpotFXArray() { return spotFX.empty() ? NULL : & spotFX[0]; }

    /** Returns pointer to internal storage for MaxDiff(usion)Date / MaxCurveMaturity */
    virtual QMCHelperBoundedDiffusion*    getDiffusionBound() { return &diffBound;}

    virtual double getOriginalLnFwdSpot(SpotIdx dateIdx); // where date is from getSpotDateIdx list
    virtual double getOriginalLnExpectedFwdSpot(FwdIdx dateIdx); // where date is from getForwardForwardDateIdx list

private:
    QMCRatesDiffuseSP     domesticIRAsset;
    QMCRatesDiffuseSP     foreignIRAsset;
    const vector<double>* domLnMONEY; // to domestic IR LnMONEY

    int             randomIndex;  // index into random numbers
    const double*   randoms;      // across dates
    vector<double>  sigmaFX;      // across time points (typically needed by EQ)
    double          LnQ;          // log of deflated spot fx at current time pt
    int             todayIdx;     /* number of historic sim dates */
    int             stopIdx;      // when next to save values in savedLnQ
    int             spotFXPos;    // position in spotFX
    vector<double>  spotFXComplete;       // saved values of log(spot FX) at ALL dates (may not be used)
    bool            completePathRequested; //Equity might need access to all the FX path for 'struck' ccy treatment
    int             expSpotFXPos; // position in expSpotFX
    vector<double>  SpotVol;      // t-depndnt part of the LOGNORMAL spot fx vol
    vector<double>  FXSmile_a1;   // skew
    vector<double>  FXSmile_a2;   // smile
    vector<double>  FXSmile_a3;   // ?
    vector<double>  LnFwdFx;      // log of determinstic forward fx rate
    vector<double>  yearFrac;     // between sim dates
    vector<double>  sqrtYearFrac; // sqrt of yearFrac
    int             spotFXStart;  // where to start in spotFX and spotFXIndexes

    vector<int>     spotFXIndexes;// indexes for when we save spot FX
    vector<double>  spotFX;       // saved values of log(spot FX)
    BoolArray       spotFXisIrExp;
    vector<double>  spotFXCutoff;

    vector<int>     expSpotFXIndexes;// indexes for when we save spot FX
    vector<double>  expSpotFX;       // saved values of log(spot FX)
    BoolArray       expSpotFXisIrExp;
    vector<double>  expSpotFXCutoff;

    double          origLnQ;         // log of spot at today
    DoubleArray  originalLnFwdSpot;         // [Ndf], used only by moment matching as reference values
    DoubleArray  originalLnExpFwdSpot; 


    bool            isIrExploded;    // whether or not one of the two has exploded
    double          explosionHelper;
    size_t          domDiscYCIdx;// = domesticIRAsset->getDiscYCIdx(); // get Idx of disc YC
    size_t          forDiscYCIdx;// = foreignIRAsset->getDiscYCIdx();


    bool            initialized;
    bool            datesFixed;

    DateTime        today;
    vector<FwdIdx>  fwd2DomIRFwd;  // forward date index in DomIR's fwd dates
    vector<FwdIdx>  fwd2ForIRFwd;  // forward date index in ForIR's fwd dates
    SRMFXUtilSP     fxUtil;
    QMCFXBoundedDiffusion    diffBound;
};

DECLARE(SRMFXDiffuse);


DRLIB_END_NAMESPACE
#endif

