#ifndef EDR_SRMFXVOLSPOT_HPP
#define EDR_SRMFXVOLSPOT_HPP

#include "edginc/SRMFXVolBase.hpp"
#include "edginc/VolProcessed.hpp"
#include "edginc/Calibrator.hpp"
//#include "edginc/IDynamicsParameter.hpp"
#include "edginc/CriticalDateCollector.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/** Not clear whether this vol will be reused by the equity components. So
    for now just put here directly. It is market data so you could argue it
    ought to be in the market directory ... */
class MARKET_DLL SRMFXVolSpot : public SRMFXVolBase,
             public virtual IVolProcessed,
             public virtual Calibrator::IAdjustable//,
             //public virtual IDynamicsParameter
{
//    friend class Util::Imp;
//    friend class SRMFX;
    friend class IrConverter;
public:
    /** The type of original FX vol - in particular models (via the market
        data fetcher) will need to select a particular SRM type of fx vol */
    static CClassConstSP const VOL_TYPE;

    //static CClassConstSP const TYPE; // inside SRMFX class for ease
    virtual ~SRMFXVolSpot(){}

    /** ctor used by SRMFXVolSimple */
    SRMFXVolSpot(
        string        name,
        ExpiryArraySP spotVolExpiry,
        DoubleArray   spotVol,
        ExpiryArraySP smileExpiry,
        DoubleArray   smileA1,
        DoubleArray   smileA2,
        DoubleArray   smileA3,
        DateTime      today,
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

    /** returns merged benchmark dates (spot vol and comp vol) */
    DateTimeArray getVolBmDates() const { return spotVolDate; }

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
        const Correlation* eqFXCorr) const{
        throw ModelException("SRMFXVolSpot::getProcessedVol",
                             "Currency struck processed vol not supported, but used for " +name);
    }
    /** force all times of day to be the same as today - unfortunate hack */
    void forceTimeOfDay(void);

    // Accessors

    const DateTimeArray& getSpotVolDate() const { return spotVolDate;}
    const DoubleArray& getSpotVol() const { return spotVol; }
    const DoubleArray& getSmileA1() const { return smileA1;}
    const DoubleArray& getSmileA2() const { return smileA2; }
    const DoubleArray& getSmileA3() const { return smileA3; }
    const DateTimeArray& getSmileDate() const { return smileDate; }

    void setToday( const DateTime& date ) {today = date;}

private:
    SRMFXVolSpot(): SRMFXVolBase(VOL_TYPE) {}

    static IObject* defaultConstructor(){
        return new SRMFXVolSpot();
    }



    static ExpiryArraySP getSmileParamExpiries(const IObject* obj);

    static void load(CClassSP& clazz);

    static void acceptCriticalDateCollector(
        const SRMFXVolSpot*             vol,
        CriticalDateCollector* collector);

    /************* exported fields *************/
    string        name;            ///< name of the vol
    ExpiryArraySP spotVolExpiry;
    DoubleArray   spotVol;
    ExpiryArraySP smileExpiry;
    DoubleArray   smileA1; // FIX - change to vector<double>
    DoubleArray   smileA2; // FIX - change to vector<double>
    DoubleArray   smileA3; // FIX - change to vector<double>
    // what about time metric - need to review

    /************* transient fields *************/
    DateTime      today;
    DateTimeArray spotVolDate;
    DateTimeArray smileDate;
//    bool          forcedTimeOfDay;

    /************* end of fields *************/
};

DECLARE(SRMFXVolSpot);

DRLIB_END_NAMESPACE
#endif


