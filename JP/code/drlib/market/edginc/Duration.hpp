//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : Duration.hpp
//
//   Description : Defines an interface for MarketObjects that can return
//                 a duration calculation and a ClientRunnable means of
//                 performing the calculation
//                 Abstract base class for concrete implementations that
//                 will set mktObject
//
//   Author      : Gordon Stephens
//
//   Date        : 12 April 2005
//
//
//----------------------------------------------------------------------------

#ifndef DURATION_HPP
#define DURATION_HPP

#include "edginc/ClientRunnable.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/Model.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/ExpiryResult.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(VectorShift)

class MARKET_DLL Duration : public CObject, public virtual ClientRunnable {
public:
    static CClassConstSP const TYPE;

    // Base interface that MarketObjects need to implement
    class MARKET_DLL IParHandler : public virtual IObject {
    public:
        static CClassConstSP const TYPE;
        static void load(CClassSP& clazz);

        //return benchmarks for par instruments
        //these will determine the durations to calculate unless specified
        virtual ExpiryArrayConstSP getParBenchmarks() const = 0;
    };

    // Subtype of IParHandler:
    // Interface for MarketObjects which support closed form calculation
    class MARKET_DLL IParHandlerWithClosedForm : public IParHandler {
    public:
        static CClassConstSP const TYPE;
        static void load(CClassSP& clazz);

        //return closed form solution
        virtual ExpiryResultArraySP getDuration(const Duration* duration) const = 0;
    };

    // Subtype of IParHandler:
    // Interface for MarketObjects which do not support closed form calculation
    class MARKET_DLL IParHandlerWithoutClosedForm : public IParHandler {
    public:
        static CClassConstSP const TYPE;
        static void load(CClassSP& clazz);

        //return a means of tweaking the MarketObject in a pointwise manner
        //assumption is that the results of this is a VectorResult...
        virtual VectorRiskPropertySensitivityConstSP getPointwiseTweaker() const = 0;

        //return a par instrument for the specified maturity
        virtual InstrumentSP getParInstrument(const ExpirySP maturity) const = 0;
    };

    // Subtype of IParHandler:
    // Interface for MarketObjects which do not support closed form calculation
    // DEPRECATED --- will go away when we port RhoPointwise to the
    // IRiskQuantityFactory sensitivity framework

    class MARKET_DLL IParHandlerWithoutClosedForm_VectorShift : public IParHandler {
    public:
        static CClassConstSP const TYPE;
        static void load(CClassSP& clazz);

        //return a means of tweaking the MarketObject in a pointwise manner
        //assumption is that the results of this is a VectorResult...
        virtual VectorShiftSP getPointwiseTweaker() const = 0;

        //return a par instrument for the specified maturity
        virtual InstrumentSP getParInstrument(const ExpirySP maturity) const = 0;
    };

    
    // Returns the discount curve (if any)
    YieldCurveWrapper getDiscount() const;
    
    // Returns the benchmark expiries (can be an empty SP)
    const ExpiryArrayConstSP getBenchmarks() const;
    
    // Returns the model
    IModelSP getModel() const;
    
    // Returns the market
    MarketDataSP getMarket() const;

    //---------------
    // ClientRunnable
    //---------------
    virtual IObjectSP run();

private:

    //----------------
    // CObject methods
    //----------------
    static void load(CClassSP& clazz);

     //----------------
    //Duration methods
    //----------------
    static IObject* defaultDuration();

    //for reflection
    Duration();

    //fields
    MarketWrapper<IParHandler> name;           //the name of the market data item for which durations are required
    MarketDataSP               market;         //the market data
    IModelSP                   model;          //optional, can be used for retrieving market data
    YieldCurveWrapper          discount;       //optional discount curve, for CDSParSpreads
    ExpiryArrayConstSP         benchmarks;     //optionally, the points for which durations are required
    string                     type;           //optionally, the explicit type of the market data with the above name
    bool                       usedClosedForm; //optionally, return a closed form result if the data supports it
};

DRLIB_END_NAMESPACE

#endif
