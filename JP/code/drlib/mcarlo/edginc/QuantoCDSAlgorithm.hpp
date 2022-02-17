//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : QuantoCDSAlgorithm.hpp
//
//   Description : How to apply quanto adjustment to CDS Par Spread Curves
//
//   Author      : Mark A Robson
//
//   Date        : November 30, 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_QUANTO_CDS_ALGORITHM_HPP
#define EDR_QUANTO_CDS_ALGORITHM_HPP

#include "edginc/QuantoCDSParSpreads.hpp"
#include "edginc/IRGridPointCache.hpp"

DRLIB_BEGIN_NAMESPACE

/** How to apply quanto adjustment to CDS Par Spread Curves */
class MCARLO_DLL QuantoCDSAlgorithm: public virtual QuantoCDSParSpreads::IAlgorithm {
public:
    virtual ~QuantoCDSAlgorithm();
    
    /** Return clean spreads in the new currency 
        (ie fx base/dom currency) */
    virtual DefaultRatesConstSP currencyAdjust(
        IYieldCurveConstSP      domYC,   // set to Projection curve
        IYieldCurveConstSP      forYC,   // set to Projection curve
        ICDSParSpreadsConstSP   parSpreads,
        FXAssetConstSP          fx,
        CorrelationConstSP      corrFxCDS,     // fx and cds par spreads
        CorrelationConstSP      corrFxDom,     // fx and dom IR
        CorrelationConstSP      corrFxFor,     // fx and for IR
        CorrelationConstSP      corrCDSDom,    // cds and domestic iR
        CorrelationConstSP      corrCDSFor,    // cds and foreign IR
        CorrelationConstSP      corrDomFor) const; // dom and for IR
    
    /** Returns debug info - may be null eg if the algorithm has not been 
        used. Object must have been initialised with debug flag on */
    IObjectSP getDebugInfo() const;

    /** Constructor */
    QuantoCDSAlgorithm(
        const string&      modelParamsKey, // for domestic and foreign IR
        const string&      calibrationStyle,// for domestic and foreign IR
        const string&      calibrationMaturity,// for domestic and foreign IR
        bool               skipIRBadVols,// for domestic and foreign IR
        const string&      fxVolBootstrapMode,
        double             fxCutOffLevel,
        bool               floorIntensity, // true: floor negative vols at 0
        string             corrSwapStart, // eg 1Y  (offset to today)
        string             corrSwapMat,   // eg 10Y (offset to start)
        string             corrSwapDCC,   // eg Act/365F
        string             corrSwapFreq, // eg 6M
        bool               debug); // true: store extra info

    /** The supplied IRGridPointCache will be updated as IRVols are used.
        This replaces any previous IRGridPointCache. Note no clone of the
        supplied parameter is taken */
    void setIRGridPointCache(IRGridPointCacheSP irGridPtsCache);

    /** Factor out the caching of IR vol points */
    void cacheGridPoints(const IYieldCurveConstSP domYC,
                         const IYieldCurveConstSP forYC) const;
    //------------------------------
    // Cache-related methods
    //------------------------------
    /** Hash code function - the cache needs improved performance compared
     * to the default "hashCode" function in CObjet: only the required
     * components are hashed */
    virtual int hashCodeOpt() const;

    /** Comparison function - the cache needs improved performance compared 
     * to the default "equalTo" function in CObjet: only the required
     * components are compared */
    virtual bool equalToOpt(
        const QuantoCDSParSpreads::IAlgorithm* algorithm) const;

private:
    /** Factor out the retrieval of IR vols */
    VolProcessedBSIRSP getIRProcessedVol(const IYieldCurveConstSP irYC) const;

    class Imp;
    class Debug;
    //// fields ///
    string modelParamsKey; // for domestic and foreign IR
    string calibrationStyle;// for domestic and foreign IR
    string calibrationMaturity;// for domestic and foreign IR
    bool   skipIRBadVols;// for domestic and foreign IR
    string fxVolBootstrapMode;
    double fxCutOffLevel;
    bool   floorIntensity; // true: floor negative vols at 0
    // data used to correlate IR with FX
    string corrSwapStart; // eg 1Y  (offset to today)
    string corrSwapMat;   // eg 10Y (offset to start)
    string corrSwapDCC;   // eg Act/365F
    string corrSwapFreq; // eg 6M
    mutable IRGridPointCacheSP irGridPtsCache;
    bool   debug; // true: store extra info
    mutable smartPtr<Debug> debugOutput; // holds debug info
};

typedef smartPtr<QuantoCDSAlgorithm> QuantoCDSAlgorithmSP;
typedef smartConstPtr<QuantoCDSAlgorithm> QuantoCDSAlgorithmConstSP;
DRLIB_END_NAMESPACE

#endif




