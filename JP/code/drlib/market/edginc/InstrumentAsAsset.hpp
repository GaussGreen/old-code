//----------------------------------------------------------------------------
//
//   Group       : Equity Derivatives Research
//
//   Filename    : InstrumentAsAsset.hpp
//
//   Description : What an instrument has to implement in order to be used as
//                 an asset (via InstrumentAsset)
//
//   Author      : Mark A Robson
//
//   Date        : 23 Mar 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_InstrumentAsset_HPP
#define EDR_InstrumentAsset_HPP
#include "edginc/Instrument.hpp"


DRLIB_BEGIN_NAMESPACE
/** What an instrument has to implement in order to be used as an
    asset (via InstrumentAsset) */
class MARKET_DLL IInstrumentAsAsset: public virtual IInstrument{
public:
    static CClassConstSP const TYPE; // defined in InstrumentAsset.cpp
    ~IInstrumentAsAsset(); // defined in InstrumentAsset.cpp
    IInstrumentAsAsset(); // defined in InstrumentAsset.cpp

    /** Returns a date after which the instrument can no longer be used as
        an asset (eg after the maturity date of a bond) */
    virtual DateTime maturityDate() const = 0;

    /** Returns the 'coupons' or payments that this instrument will make during
        its lifetime. This can include historic payments. */
    virtual CashFlowArraySP getCoupons() const = 0;

    /** Returns the dates on which the instrument has to be held in order to
        hold the right to the corresponding coupon as returned by 
        getCoupons() */
    virtual DateTimeArraySP getExCouponDates() const = 0;

    /** Returns the yield curve used for discounting */
    virtual YieldCurveConstSP getDiscount() const = 0;

    /** Returns the accured interest (if any) to date */
    virtual double getAccrued() const = 0;
private:
    static void load(CClassSP& clazz);
};
typedef smartPtr<IInstrumentAsAsset> IInstrumentAsAssetSP;


DRLIB_END_NAMESPACE
#endif
