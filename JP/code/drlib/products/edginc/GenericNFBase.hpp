//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : GenericNFBase.hpp
//
//   Description : Useful base class for generic n factor equity products
//                 which have a reference level and past values.
//
//   Author      : Mark A Robson
//
//   Date        : 24 October 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_GENERICNFBASE_HPP
#define EDR_GENERICNFBASE_HPP

#include "edginc/GenericNFactor.hpp"
#include "edginc/RefLevel.hpp"
#include "edginc/PastValues.hpp"
#include "edginc/SensitiveStrikes.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/ISDAConvention.hpp"
#include "edginc/ObservationMap.hpp"
#include "edginc/PastSamplesEvent.hpp"
DRLIB_BEGIN_NAMESPACE


/** Useful base class for generic n factor equity products
    which have a reference level and past values. */
class PRODUCTS_DLL GenericNFBase: public GenericNFactor,
                     virtual public LastSensDate,
                     virtual public ISensitiveStrikes,
                     virtual public Theta::Shift,
                     virtual public PastSamplesEvent::IEventHandler {

public:
    static CClassConstSP const TYPE;

    void GetMarket(const IModel* model, const CMarketDataSP market);

    virtual ~GenericNFBase();

    //// roll through time (setting historic values)
    virtual bool sensShift(Theta* theta);

    /** Validate instrument having aquired market data */
    void Validate();

    //// Checks for refLevel and pastValues. Overridden where they are optional.
    virtual void validatePop2Object();

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model);

    /** returns all strikes on the vol surface to which
        this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model);

    /** when to stop tweaking - default implementation assume product can
        be priced with a MC */
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    // implementation of PastSamplesEvent::IEventHandler interface
    void getEvents(const PastSamplesEvent* samples, IModel* model,
                   const DateTime& eventDate, EventResults* events) const;

protected:
    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    // Used to make sure the PastValues object has all necessary dates in it
    // before we start pricing so existing handling can work for theta etc
    // Basically pastValues will contain the historic overrides when it comes
    // in from the client and then in the Validate method we fill in all the
    // missing dates putting zero against future ones and getting missing past
    // values from the centralised asset history
    // Note this isn't what's sampled in the Monte Carlo as we could have
    // sparse monitoring and barrier adjustment etc
    virtual const DateTimeArray samplingDates() const = 0;

    // most products have a single set of sampling dates so this defaults to true
    // some products e.g. SuperRainbow have a set per asset
    virtual bool assetsShareSamplingDates() const;

    virtual const DateTimeArrayArray* getSamplingDatesPerAsset() const;

    GenericNFBase(CClassConstSP clazz);

    GenericNFBase(CClassConstSP            clazz, 
                  const DateTime&          valueDate,
                  double                   notional,         
                  InstrumentSettlementSP   instSettle,  
                  IMultiMarketFactorsSP    assets, 
                  const YieldCurveWrapper& discount,
                  IRefLevelSP              refLevel,
                  IPastValuesSP            pastValues);

    // fields
    IRefLevelSP             refLevel;      // gives 'avg in' for each asset
    IPastValuesSP           pastValues;    // all historic (overriding) spots

    // new fields for ISDA holiday adjustment and centralised sampling
    // isdaDateAdjust optional merely for backwards compatibility and migration
    SamplingConventionSP    isdaDateAdjust;    // n factor date adjustment details
    //handles mapping between sample dates and observation dates
    ObservationMapSP        obsMap;


private:

    // this method is called after the market data is retrieved.
    // basically it fills in the missing pieces of the past values object
    // so that we can guarantee that all dates we will ever need are present
    // and all dates in the past have values (if not an Exception is thrown)
    void validatePastValues();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    GenericNFBase(); // not implemented
    GenericNFBase(const GenericNFBase& rhs);
    GenericNFBase& operator=(const GenericNFBase& rhs);
};

typedef smartConstPtr<GenericNFBase> GenericNFBaseConstSP;
typedef smartPtr<GenericNFBase> GenericNFBaseSP;

DRLIB_END_NAMESPACE
#endif

