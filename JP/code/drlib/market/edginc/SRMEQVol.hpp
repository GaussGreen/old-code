#ifndef EDR_SRMEQVOL_HPP
#define EDR_SRMEQVOL_HPP

#include "edginc/VolBase.hpp"
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
class MARKET_DLL SRMEQVol: public CVolBase,
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
//    friend class SRMEquityDiffuse;
    friend class IrConverter;
public:
    static CClassConstSP const TYPE;
    virtual ~SRMEQVol(){}

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
        const CVolRequest* /*volRequest*/,
        const CAsset*      /*eqAsset*/,
        const FXAsset*     /*fxAsset*/,
        const Correlation* /*eqFXCorr*/) const{
        throw ModelException("SRMEQVol::getProcessedVol", 
                             "Currency struck processed vol not supported, but used for " +name);
    }

    // Accessors

    const DateTimeArray& getCompVolDate() const { return compVolDate;}
    const DateTimeArray& getSpotVolDate() const { return spotVolDate;}
    const DoubleArray& getSpotVol() const { return spotVol; }
    const DateTimeArray& getCompVolResetDate() const { return compVolResetDate;}
    const DoubleArray& getCompVol() const { return compVol;}
    const DoubleArray& getSmileA1() const { return smileA1;}
    const DoubleArray& getSmileA2() const { return smileA2; }
    const DoubleArray& getSmileA3() const { return smileA3; }
    const DateTimeArray& getSmileDate() const { return smileDate; }

private:
    SRMEQVol(): CVolBase(TYPE){}

    static IObject* defaultConstructor(){
        return new SRMEQVol();
    }

    // This is crude - perhaps someone has done something better elsewhere? XXX
    static ExpirySP toExpiry(string& expiryStr);

    /** Turn BM strings into date times */
    void cacheDates();

    static ExpiryArraySP getSmileParamExpiries(const IObject* obj);

    static ExpiryArraySP getCompVolExpiries(const IObject* obj);

    static ExpiryArraySP getSpotVolExpiries(const IObject* obj);

    static void load(CClassSP& clazz);

    static void acceptCriticalDateCollector(
        const SRMEQVol*             vol,
        CriticalDateCollector* collector);


    /****************************************
     * Composite volatility - Parallel tweak 
     ****************************************/
    virtual string sensName(CompositeVegaParallelTweak* shift) const;

    virtual bool sensShift(CompositeVegaParallelTweak* shift);

    /*****************************************
     * Composite volatility - Pointwise tweak 
     *****************************************/
    virtual string sensName(CompositeVegaPointwiseTweak* shift) const;

    virtual bool sensShift(CompositeVegaPointwiseTweak* shift);

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

    string        name;    // name of the vol
    StringArray   compVolBMs; // -> compVolDate
    StringArray   compVolResetBMs; // -> compVolResetDate
    DoubleArray   compVol;
    StringArray   spotVolBMs; // -> spotVolDate
    DoubleArray   spotVol;
    StringArray   smileBMs; // -> smileDate
    DoubleArray   smileA1;
    DoubleArray   smileA2;
    DoubleArray   smileA3;
    // what about time metric - need to review
    // transient fields
    DateTime      today;
    DateTimeArray compVolDate; // need some comments here ...
    ExpiryArraySP compVolExpiry;
    DateTimeArray compVolResetDate; // aka compVolMatDate ~ FX
    DateTimeArray spotVolDate;
    ExpiryArraySP spotVolExpiry;
    DateTimeArray smileDate;
    ExpiryArraySP smileExpiry;    

};

typedef smartPtr<SRMEQVol> SRMEQVolSP;
typedef smartConstPtr<SRMEQVol> SRMEQVolConstSP;

DRLIB_END_NAMESPACE
#endif // EDR_SRMEQVOL_HPP
