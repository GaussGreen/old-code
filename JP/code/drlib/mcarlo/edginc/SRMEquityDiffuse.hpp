//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMEquityDiffuse.hpp
//
//   Description : A generator of paths using stochastic rates
//                 for Equity Assets
//
//   Date        : Nov 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_SRMEQUITYDIFFUSE_HPP
#define EDR_SRMEQUITYDIFFUSE_HPP

#include "edginc/DividendList.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/QMCEquityDiffuse.hpp"
#include "edginc/QMCRatesDiffuse.hpp"
#include "edginc/SRMEquityUtil.hpp"
#include "edginc/QMCHelperBoundedDiffusion.hpp"
#include "edginc/HyperTrigUtils.hpp"
#include "edginc/MCPathConfigSRM.hpp"

#include <set>

DRLIB_BEGIN_NAMESPACE

class IPastValues;

/** A generator of paths using stochastic rates for EQ Assets */
class SRMEquityDiffuse   : public QMCEquityDiffuse
{

public:
    static const string CONSTANT_SPOT_VOL;
    static const string NO_FAILURE_ALLOWED;
    static const string USE_LAST_LEVEL;

    virtual ~SRMEquityDiffuse() {}
    SRMEquityDiffuse(IQMCDiffusibleInterestRateSP domIR);

    void  setSRMEquityDiffuse(
                    int                   _randomIndex,
                    const DateTime &      _today,
                    SRMEquityUtilSP       _equityUtil,
                    const vector<double>& _pastValues,
                    int timePointsPerYear); // modifies spotEq and spotEqStart

    // used only within SRM
    static void        getDivData(AssetConstSP       asset,
                                  const DateTime&    start,
                                  const DateTime&    end,
                                  DateTimeArray&     divExDates, // time adjusted to fit timeline (->SOD)
                                  bool               withDivList, // if withDivList must supply yc too
                                  DividendListSP&    divList, // optional
                                  IYieldCurveConstSP yc);     // optional

    // used only within SRM
    static DateTimeArray getCriticalDates(AssetConstSP    asset,
                                          const DateTime& start,
                                          const DateTime& end);
///////////////////////// IQMC interface  ////////////////////////////////
        /** get all the dates important for asset */
    virtual void finalize(DateTimeArrayConstSP allDates);
    /** getting the simulation start date */
    virtual DateTime getBaseDate() { return today; }

    /** Accessing the diffused SpotFX, i.e. a non-default probability by
        a given date */
    virtual double  getSpotPrice(SpotIdx measurementDateIdx)
    {
        return spotEq[measurementDateIdx];
    }

    /// Override, so we can inform domIR about asset's dates
    virtual void   addAggregatedDates(  const DateTimeArray& spot,
                                        const DateTimeArray& forward,
                                        const DateTimeArray& forwardForward);


    /** (3) Getting the simulated ExpectedFX between [measDate, futureDate] */

    /** Accessing the ExpectedFX(md, fd) where md is a
        simulated measurement date and fd is some future date after the
        measurement is made. */
    virtual double  getExpectedPrice(FwdIdx measurementDateIdx,
                                        FwdIdx futureDateIdx);

    /** Accessing the natural log of the ExpectedFX(md, fd)
        where md is a simulated measurement date and fd is some future
        date after the measurement is made. */
    virtual double getLnExpectedPrice(FwdIdx measurementDateIdx,
                                FwdIdx futureDateIdx);

    virtual IQMCDiffusibleInterestRateSP getDomesticIRAsset() {
        return domIR;
    }
    virtual double* getInternalSpotPriceArray() { return spotEq.empty() ? NULL : & spotEq[0]; }

    /** Returns pointer to internal storage for MaxDiff(usion)Date / MaxCurveMaturity */
    virtual QMCHelperBoundedDiffusion*    getDiffusionBound() { return &diffBound;}

    virtual double getOriginalLnFwdSpot(SpotIdx dateIdx); // where date is from getSpotDateIdx list
    virtual double getOriginalLnExpectedFwdSpot(FwdIdx dateIdx); // where date is from getForwardForwardDateIdx list

protected:

    //used by strict second order and mapping method:
    //those dates which must be used for big steps: e.g. when the smile changes.
    size_t setupCriticalBigStepBools(vector<int> &bigStepBools);

    void setupBigStepIndicesAndYearFracs(const vector<int> &bigStepBools, vector<int> &bsIdxs, vector<double> &bsYearFracs, vector<double> &bsYearFracsSqrt);

    void computeLocalSmile(double logNormalEq, int simDateIdx, double *smileEQ, double *smileDashEQ);
    

    QMCRatesDiffuseSP     domIR;

    bool            smileIsFlatEverywhere; // true if 'a1' and 'a2' is zero for all time
    int             randomIndex;  // index into random numbers
    const double*   randoms;      // across dates
    vector<double>  sigmaEQ;      // local vol across time points
	double          LnE;          // log of deflated spot EQ at current time pt
    int             todayIdx;     /* number of historic sim dates */
    int             stopIdx;      // when next to save values in savedLnE
    int             spotEqPos;    // position in spotEQ
    int             expSpotEqPos; // position in expSpotEQ
    CurrencyTreatment   ccyTreatment;   //vanilla, protected, or struck.
    bool            isStruck;       //Currency treatment is 'struck' for this asset.
    const vector<double>* irLnMONEY; // from base IR
    const vector<double>* sigmaFX; // for quanto/struck adjust
    const vector<double>* spotFXFullPath; //for struck adjust
    double          corrEqFX;     // for any ccy prot adjustment
    vector<double>  SpotVol;      // t-depndnt part of the LOGNORMAL spot EQ vol
    vector<double>  EqSmile_a1;   // skew
    vector<double>  EqSmile_a2;   // smile
    vector<double>  EqSmile_a3;   // ?
    vector<int>     smileChanges; // simDate indices where smile changes
    vector<HyperLocalVolState> smileStates;    //store state of smile for each combination of smile parameters
    vector<double>  LnFwdEq;      // log of deterministic forward EQ rate
    
    DoubleArray  originalLnFwdSpot;         // [Ndf], used only by moment matching as reference values
    DoubleArray  originalLnExpFwdSpot; 

    vector<double>  yearFrac;     // between sim dates
    vector<double>  sqrtYearFrac; // sqrt of yearFrac
    int             spotEqStart;  // where to start in spotEQ and spotEQIndexes
    vector<int>     spotEqIndexes;// indexes for when we save spot EQ
    vector<double>  spotEq;       // saved values of exp(log(spot EQ))
    vector<int>     expSpotEqIndexes;// indexes for when we save spot EQ
    vector<double>  expSpotEq;       // saved values of log(spot EQ)
    double          origLnE;         // log of spot at today
    vector<int>     exDivIndexes; // when divs go ex
    int             divPos;
    vector<double>  lnDivYield;     // = log(1-Y)
    vector<double>  CfactorArrayEQ; // continuous factors (borrow, cts divs[not yet])

    //For some discretization methods (e.g. second order strict) we may ignore the timeline, 
    //and only diffuse across a subset of the dates. We may use timePointsPerYear as a minimum value.           
    int             timePointsPerYear; 
       
    DateTime        today;
    vector<FwdIdx>  fwd2DomIRFwd;  // forward date index in DomIR's fwd dates

    size_t          domDiscYCIdx;
    SRMEquityUtilSP equityUtil;
    bool            initialized;
    QMCEQBoundedDiffusion    diffBound;
};

DECLARE(SRMEquityDiffuse);

DRLIB_END_NAMESPACE
#endif

