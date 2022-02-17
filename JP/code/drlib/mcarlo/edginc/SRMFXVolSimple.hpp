//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : SRMFXVolSimple.hpp
//
//   Description : Simple version of SRMFXVol
//
//
//----------------------------------------------------------------------------

#ifndef EDR_SRMFXVOLSIMPLE_HPP
#define EDR_SRMFXVOLSIMPLE_HPP

#include "edginc/FXVolBase.hpp"
#include "edginc/VolProcessed.hpp"
#include "edginc/Theta.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/TweakableWith.hpp"
#include "edginc/PointwiseTweakableWith.hpp"
#include "edginc/CompositeVegaPointwiseTweak.hpp"
#include "edginc/CompositeVegaParallelTweak.hpp"
#include "edginc/SpotVegaParallelTweak.hpp"
#include "edginc/SpotVegaPointwiseTweak.hpp"
#include "edginc/IDynamicsParameter.hpp"
#include "edginc/CriticalDateCollector.hpp"

DRLIB_BEGIN_NAMESPACE

/** Version of SRMFXVol where the parameters A1, A2, A3 are parameterized
    in the time dimension. Each term structure uses 3 parameters; namely,
    the 1-year smile param, the short term proagation power,
    the long term propagation exponent */
    
class MCARLO_DLL SRMFXVolSimple: public FXVolBase,
             public virtual IVolProcessed,
             public virtual Theta::IShift,
             public virtual Calibrator::IAdjustable
{
public:
    static CClassConstSP const VOL_SIMPLE_TYPE;

    virtual ~SRMFXVolSimple(){}

    /** calculates the trading time between two dates */
    virtual double calcTradingTime(const DateTime &date1,
                                   const DateTime &date2) const;

    /** retrieve time measure for the vol */
    virtual TimeMetricConstSP GetTimeMetric()const;

    virtual void validatePop2Object();

    /** populate from market cache */
    void getMarket(const IModel* model, const MarketData* market);

    /** Shifts the object using given shift. */
    virtual bool sensShift(Theta* shift);

    /** I guess we can do something with the spot (=>fwd ) vols? */
    virtual CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest,
        const CAsset*      asset) const;

    /** Returns name of vol */
    virtual string getName() const;

    /** I guess we can do something with the spot (=>fwd ) vols? */
    virtual CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest,
        const CAsset*      eqAsset,
        const FXAsset*     fxAsset,
        const Correlation* eqFXCorr) const;

private:
    SRMFXVolSimple();

    static IObject* defaultConstructor();

    static void load(CClassSP& clazz);

    static void acceptCriticalDateCollector(
        const SRMFXVolSimple*      vol,
        CriticalDateCollector* collector);

    /** Turn expiries into date times */
    void cacheDates();

//    friend class SRMFXDiffuse;
    string        name;    // name of the vol
    ExpiryArraySP compVolExpiry; // need some comments here ...
    ExpiryArraySP compVolMatExpiry; // optional
    DoubleArray   compVol;
    ExpiryArraySP spotVolExpiry;
    DoubleArray   spotVol;
    // the smile expiries are still needed cause 'SRMFXVolSimple' will be
    // mapped into 'Vol' at each of those expiries
    ExpiryArraySP smileExpiry;
    // the term structure of parameter A1 is parameterized by 3 parameters,
    // namely:
    double        smileA11Y;                // 1-year smile param
    double        smileA1ShortTermPower;    // short term propagation power
    double        smileA1LongTermExpo;      // long term propagation exponent
    // and likewise for smileA1 and smileA2
    double        smileA21Y;                // 1-year smile param
    double        smileA2ShortTermPower;    // short term propagation power
    double        smileA2LongTermExpo;      // long term propagation exponent
    double        smileA31Y;                // 1-year smile param
    double        smileA3ShortTermPower;    // short term propagation power
    double        smileA3LongTermExpo;      // long term propagation exponent
    // what about time metric - need to review
    // transient fields
    DateTime      today;
    DateTimeArray compVolDate;
    DateTimeArray compVolMatDate;
    DateTimeArray spotVolDate;
    DateTimeArray smileDate;

};

typedef smartPtr<SRMFXVolSimple> SRMFXVolSimpleSP;
typedef smartConstPtr<SRMFXVolSimple> SRMFXVolSimpleConstSP;

DRLIB_END_NAMESPACE
#endif


