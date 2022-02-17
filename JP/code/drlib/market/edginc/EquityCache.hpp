//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EquityCache.hpp
//
//   Description : a simple equity that uses fwd cache etc.
//
//   Author      : Ning Shen
//
//   Date        : 24 June 2002
//
//
//----------------------------------------------------------------------------

#ifndef EQUITY_CACHE_H
#define EQUITY_CACHE_H

#include "edginc/config.hpp"
#include "edginc/EquityBase.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"
#include "edginc/RateParallel.hpp"


DRLIB_BEGIN_NAMESPACE

static const int NUM_FWD = 4; // must be the same as num of fwd in enum (= 4)
static const int NUM_TWEAK = 5; // must match the num of enum TweakType

//*** EquityCache class allows setting spot and vol and use of cahced fwd ***/
class MARKET_DLL EquityCache : public EquityBase,
                    virtual public RhoBorrowParallel::RestorableShift,
                    virtual public IRestorableWithRespectTo<RateParallel>,
                    virtual public MuParallel::Shift,
                    virtual public Theta::Shift
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
    static IObject* defaultEquityCache(){
        return new EquityCache();
    }

    EquityCache(const Equity*   equity, const CVolBase* vol_inp);

    /** set spot level */
    void setSpot(double spot);

    /** set vol level */
    void setVol(double vol);

    /** override to use cache */
    double fwdValue(const DateTime& date) const;

    /** override to use cache, returns an array of forward prices which are assumed to be in 
        ascending order */
    void fwdValue(const DateTimeArray& dateList,
                  CDoubleArray&        result) const;

    /** prc-calc fwd and cache the results  */
    void preCalcFwd(EquityBase* asset, const DateTime& valueDate, const DateTime& lastDate);

    // the IMarketObservable interface for retrieving a single sample
    virtual double pastValue(const DateTime&             sampleDate,
                             const ObservationType*      obsType,
                             const ObservationSource*    source,
                             const FixingType*           fixType,
                             const IObservationOverride* overrides,
                             const SamplingConvention*   sampleRule) const;

    // IMarketObservable - retrieve a single observation date
    // Returns false if obs is to be omitted
    virtual bool observationDate(const DateTime&           sampleDate,
                                 const ObservationSource*  source,
                                 const SamplingConvention* sampleRule,
                                 DateTime*                 obsDate) const;

    // the IMarketObservable interface for retrieving past samples events
    virtual double addPastSampleEvent(const DateTime&             sampleDate,
                                    const ObservationType*      obsType,
                                    const ObservationSource*    source,
                                    const FixingType*           fixType,
                                    const IObservationOverride* overrides,
                                    const SamplingConvention*   sampleRule,
                                    PastSamplesCollector*        collector) const;

    // the IMarketObservable interface for 
    // is the given date a holiday for the relevant source
    virtual bool isHoliday(const DateTime& sampleDate,
                           const ObservationSource*   source) const;

    // over-ride clone
    IObject* clone() const;
    
    /** Returns the name of the borrow curve - used to determine
    whether to tweak the object */
    virtual string sensName(RhoBorrowParallel* shift)const;
    
    /** Shifts the object using given shift */    
    virtual bool sensShift(RhoBorrowParallel* shift);
    
    /** Restores the object to its original form */
    virtual void sensRestore(RhoBorrowParallel* shift);

    /** Returns the name of the cash swap curve - used to determine
    whether to tweak the object */
    virtual string sensName(const RateParallel* shift)const;
    
    /** Shifts the object using given shift */    
    virtual TweakOutcome sensShift(const PropertyTweak<RateParallel>& shift);
    
    /** Restores the object to its original form */
    virtual void sensRestore(const PropertyTweak<RateParallel>& shift);

    /** Shifts the object using given shift */    
    virtual bool sensShift(Theta* shift);

    /** Returns the name of the borrow curve - used to determine
    whether to tweak the object */
    virtual string sensName(MuParallel* shift)const;
    /** Shifts the object using given shift */    
    virtual bool sensShift(MuParallel* shift);
    
    enum TweakType{BASE, RHO, RHO_BORROW, THETA, MU};

private:
    EquityCache()  : EquityBase(TYPE){};
    friend class FastQuoteEnv;

    // unregistered

    // spotRef the same as in FastQuoteEnv, kept here for easy use
    double          spotRef;  // $unregistered
    // value date integer
    int                 refDate;  // $unregistered
    // daily dates from value date to last maturity
    DateTimeArray       fwdDates; // $unregistered

    // forwards cache, [0] = fwd, [1] = delta, [2] = fwdBEX, [3] = deltaBEX
    DoubleArrayArray  fwdCache[NUM_FWD]; // $unregistered

    enum {FWD, FWD_DELTA, FWD_BEX, FWD_DELTA_BEX};

    /** set ptrs to a chosen fwd array */
    void chooseFwdCache(TweakType tweakArr);

    /** compute fwd and delta values */
    void calcFwdAndDelta(EquityBase* asset,
                     const DateTimeArray& dates,
                     const DateTimeArray& datesBEX,
                     double spotRef);

    // point to either base or rho, mutable to allow recurse to set it,
    // [0] = fwd, [1] = delta, [2] = fwdBEX, [3] = deltaBEX
    mutable DoubleArray*    forwardArray[NUM_FWD];  // $unregistered
};

typedef smartPtr<EquityCache> EquityCacheSP;

DRLIB_END_NAMESPACE

#endif
