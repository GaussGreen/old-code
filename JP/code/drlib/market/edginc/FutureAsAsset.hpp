//----------------------------------------------------------------------------
//
//   Group       : Equity Derivatives Research

//
// What a future instrument has to implement in order to be used as
// an asset (via FutureAsset)
//
//----------------------------------------------------------------------------

#ifndef FUTURE_AS_ASSET_HPP
#define FUTURE_AS_ASSET_HPP

#include "edginc/Instrument.hpp"


DRLIB_BEGIN_NAMESPACE
/** What an instrument has to implement in order to be used as an
    asset (via FutureAsset) */
class MARKET_DLL IFutureAsAsset: public virtual IInstrument{
public:
    static CClassConstSP const TYPE; // defined in FutureAsset.cpp
    ~IFutureAsAsset();     // defined in FutureAsset.cpp
    IFutureAsAsset();      // defined in FutureAsset.cpp

    /** Returns a date after which the instrument can no longer be used as
        an asset (eg after the maturity date of a bond) */
    virtual DateTime maturityDate() const = 0;

    /** Returns the yield curve used for discounting */
    virtual YieldCurveConstSP getDiscount() const = 0;
   
private:
    static void load(CClassSP& clazz);
};
typedef smartPtr<IFutureAsAsset> IFutureAsAssetSP;


DRLIB_END_NAMESPACE

#endif
