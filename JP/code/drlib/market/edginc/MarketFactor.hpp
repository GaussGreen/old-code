//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MarketFactor.hpp
//
//   Description : MarketFactor interface
//
//   Author      : Mark A Robson
//
//   Date        : 14 May 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_MARKETFACTOR_HPP
#define EDR_MARKETFACTOR_HPP
#include "edginc/Object.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/GetMarket.hpp"

DRLIB_BEGIN_NAMESPACE
class IMarketFactor;
typedef smartPtr<IMarketFactor> IMarketFactorSP;
typedef smartConstPtr<IMarketFactor> IMarketFactorConstSP;
typedef array<IMarketFactorSP, IMarketFactor> IMarketFactorArray;

/** A MarketFactor covers a wide range of things that you might want to use
    as an 'underlying' in an instrument. It includes equity type assets as
    well as IR assets and yield curves (the latter would be the natural
    choice for a swaption) */
class MARKET_DLL IMarketFactor: public virtual IGetMarket{
public:
    static CClassConstSP const TYPE;

    IMarketFactor();

    virtual ~IMarketFactor();

    /** Returns the name of this 'market factor' */
    virtual string getName() const = 0;
    
    /**
       For consideration: getYCName() and getProcessedVol(). Others?
     */

    /** If supplied market factor wrapper is using the market data cache, then
        retrieves, and makes if necessary, a market factor composed of the
        underlying asset together with the requested currency treatment.
        The different ccyTreatments are as per CAsset. Currently only
        CAssets can have anything other 'N' currency treatment */
    static void getMarketData(const IModel*                 model, 
                              const MarketData*             market,
                              const string&                 ccyTreatment,
                              const string&                 payOutYCName,
                              MarketWrapper<IMarketFactor>& factor);

private:
    class Action;
    static void load(CClassSP& clazz);
};

// support for wrapper class
typedef MarketWrapper<IMarketFactor> IMarketFactorWrapper;

/** specialisations of arrayObjectCast */
template <> class MARKET_DLL arrayObjectCast<IMarketFactorWrapper>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const IMarketFactorWrapper& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(IMarketFactorWrapper& value);

    /** Turns the IObjectSP into a IMarketFactor */
    static IMarketFactorWrapper fromIObject(IObjectSP& value);
};
// arrays of wrappers (note array of structures)
typedef array<IMarketFactorWrapper,
              IMarketFactorWrapper> IMarketFactorWrapperArray;

typedef smartPtr<IMarketFactorWrapperArray> IMarketFactorWrapperArraySP;

DRLIB_END_NAMESPACE
#endif
