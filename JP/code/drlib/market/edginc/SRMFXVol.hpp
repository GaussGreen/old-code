#ifndef EDR_SRMFXVOL_HPP
#define EDR_SRMFXVOL_HPP

#include "edginc/SRMFXVolBase.hpp"
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

/** Not clear whether this vol will be reused by the equity components. So
    for now just put here directly. It is market data so you could argue it
    ought to be in the market directory ... */
class MARKET_DLL SRMFXVol : public SRMFXVolBase,
             public virtual IVolProcessed,
             public virtual Theta::IShift,
             public virtual Calibrator::IAdjustable,
             public virtual TweakableWith<CompositeVegaParallelTweak>,
             public virtual PointwiseTweakableWith<CompositeVegaPointwiseTweak>,
             public virtual TweakableWith<SpotVegaParallelTweak>,
             public virtual PointwiseTweakableWith<SpotVegaPointwiseTweak>,
             public virtual IDynamicsParameter
{
//    friend class Util::Imp;
//    friend class SRMFXDiffuse;
    friend class IrConverter;
public:
    /** The type of original FX vol - in particular models (via the market
        data fetcher) will need to select a particular SRM type of fx vol */
    static CClassConstSP const VOL_TYPE;

    //static CClassConstSP const TYPE; // inside SRMFXDiffuse class for ease
    virtual ~SRMFXVol(){}

    /** ctor used by SRMFXVolSimple */
    SRMFXVol(
        string        name,
        ExpiryArraySP compVolExpiry,
        ExpiryArraySP compVolMatExpiry,
        DoubleArray   compVol,
        ExpiryArraySP spotVolExpiry,
        DoubleArray   spotVol,
        ExpiryArraySP smileExpiry,
        DoubleArray   smileA1,
        DoubleArray   smileA2,
        DoubleArray   smileA3,
        DateTime      today,
        DateTimeArray compVolDate,
        DateTimeArray compVolMatDate,
        DateTimeArray spotVolDate,
        DateTimeArray smileDate);

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

    /** returns merged benchmark dates (spot vol and comp vol) */
    DateTimeArray getVolBmDates() const;

    /** I guess we can do something with the spot (=>fwd ) vols? */
    virtual CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest,
        const CAsset*      eqAsset,
        const FXAsset*     fxAsset,
        const Correlation* eqFXCorr) const{
        throw ModelException("SRMFXVol::getProcessedVol",
                             "Currency struck processed vol not supported, but used for " +name);
    }
    /** force all times of day to be the same as today - unfortunate hack */
    void forceTimeOfDay(void);

    // Accessors

    const DateTimeArray& getCompVolDate() const { return compVolDate;}
    const DateTimeArray& getSpotVolDate() const { return spotVolDate;}
    const DoubleArray& getSpotVol() const { return spotVol; }
    const DateTimeArray& getCompVolMatDate() const { return compVolMatDate;}
    const DoubleArray& getCompVol() const { return compVol;}
    const DoubleArray& getSmileA1() const { return smileA1;}
    const DoubleArray& getSmileA2() const { return smileA2; }
    const DoubleArray& getSmileA3() const { return smileA3; }
    const DateTimeArray& getSmileDate() const { return smileDate; }

    void setToday( const DateTime& date ) {today = date;}

private:
    SRMFXVol(): SRMFXVolBase(VOL_TYPE), forcedTimeOfDay(false) {}

    static IObject* defaultConstructor(){
        return new SRMFXVol();
    }


    /** Turn BM strings into date times */
    void cacheDates();

    static ExpiryArraySP getSmileParamExpiries(const IObject* obj);

    static ExpiryArraySP getCompVolExpiries(const IObject* obj);


    static void load(CClassSP& clazz);

    static void acceptCriticalDateCollector(
        const SRMFXVol*             vol,
        CriticalDateCollector* collector);


    /****************************************
     * Composite volatility - Parallel tweak
     ****************************************/
    virtual string sensName (CompositeVegaParallelTweak* shift) const;

    virtual bool sensShift (CompositeVegaParallelTweak* shift);

    /*****************************************
     * Composite volatility - Pointwise tweak
     *****************************************/
    virtual string sensName(CompositeVegaPointwiseTweak* shift) const;

    virtual bool sensShift (CompositeVegaPointwiseTweak* shift);

    virtual ExpiryArrayConstSP sensExpiries(
        CompositeVegaPointwiseTweak* shift) const;


    /***********************************
     * Spot volatility - Parallel tweak
     ***********************************/
    virtual string sensName(SpotVegaParallelTweak* shift) const;

    virtual bool sensShift(SpotVegaParallelTweak* shift);

    /************************************
     * Spot volatility - Pointwise tweak
     ************************************/
    virtual string sensName(SpotVegaPointwiseTweak* shift) const;

    virtual bool sensShift(SpotVegaPointwiseTweak* shift);

    virtual ExpiryArrayConstSP sensExpiries(
        SpotVegaPointwiseTweak* shift) const;

    /************* exported fields *************/
    string        name;            ///< name of the vol
    ExpiryArraySP compVolExpiry;   ///< need some comments here ...
    ExpiryArraySP compVolMatExpiry;///< optional
    DoubleArray   compVol;
    ExpiryArraySP spotVolExpiry;
    DoubleArray   spotVol;
    ExpiryArraySP smileExpiry;
    DoubleArray   smileA1; // FIX - change to vector<double>
    DoubleArray   smileA2; // FIX - change to vector<double>
    DoubleArray   smileA3; // FIX - change to vector<double>
    // what about time metric - need to review

    /************* transient fields *************/
    DateTime      today;
    DateTimeArray compVolDate;
    DateTimeArray compVolMatDate;
    DateTimeArray spotVolDate;
    DateTimeArray smileDate;
    bool          forcedTimeOfDay;

    /************* end of fields *************/
};

typedef smartPtr<SRMFXVol> SRMFXVolSP;
typedef smartConstPtr<SRMFXVol> SRMFXVolConstSP;

DRLIB_END_NAMESPACE
#endif


