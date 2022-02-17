//----------------------------------------------------------------------------
//
//   Group       : New York Credit QR
//
//   Author      : Charles Morcom
//
//   Filename    : ICreditVanillaInstrument.hpp
//
//   Description : Interface to define a simple credit instrument. This is
//                 something that can be priced just with a risky discount
//                 curve.
//
//   Date        : 27 January 2006
//
//----------------------------------------------------------------------------

#ifndef QR_ICREDITVANILLAINSTRUMENT_HPP
#define QR_ICREDITVANILLAINSTRUMENT_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Object.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/DayCountConvention.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(ICreditVanillaInstrument)
FORWARD_DECLARE(IDiscountCurveRisky)
FORWARD_DECLARE(IForwardRatePricer)
FORWARD_DECLARE(IDecretionCurve)
FORWARD_DECLARE(IBadDayAdjuster)

/**
 * An ICreditVanillaInstrument is a simple risky instrument that may be priced
 * with nothing more than a risky discount curve (IDiscountCurveRisky).
 * This is designed to make curve bootstrapping simpler, and also to make
 * common calculations for credit underlyings (such as CDS par forward spreads)
 * easier to handle.
 * In addition to the pv methods, there are also notions of accrued interest,
 * defaulted value, maturity, and notional (base notional if the instrument is
 * amortizing)
 */
class MARKET_DLL ICreditVanillaInstrument : public virtual IObject,
                                 public virtual IInstrument {

public:
    /**Accrued interest (cash amount, not percentage) for settlement on settlementDate.*/
    virtual double getAccruedInterest(
        const DateTime&                 settlementDate,
        IForwardRatePricerSP            model
        ) const = 0;

    /**Calculates PV of the instrument, given a valuation date. 
       Settlement is whatever the instrument determines it to be (e.g. T+1 CDS; T+3 bond, etc.)
       relative to the valuation date. Value should be conditional on no new defaults 
       before valuationDate.*/
    virtual double getPV(
        const DateTime&              valuationDate, 
        const IDiscountCurveRisky&   crv,
        const IDecretionCurveConstSP prepay,
        IForwardRatePricerSP         model,
        IBadDayAdjusterConstSP       bda
        ) const = 0;

    /**Calculates the PV of the instrument on a given valuationDate for unconditional settlement
       and payment on paymentDate. Used to calculate forward PVs and differentiate between
       conditional and unconditional settlement. Value should be conditional on no new defaults
       before valuationDate if valuationDate is in the future.*/
    virtual double getPV(
        const DateTime&              valuationDate, 
        const DateTime&              paymentDate, 
        const IDiscountCurveRisky&   crv,
        const IDecretionCurveConstSP prepay,
        IForwardRatePricerSP         model,
        IBadDayAdjusterConstSP       bda
        ) const = 0;

    /**Calculates the PV of an instrument on defaultDate, given that default has just occured.
       Note that the curve is needed in case there is a delay between default and recovery
       payment. This method assumes a complete default, so should be interpreted with care
       when the instrument/curve has multiple names.
       Leg valuations are discounted to valuation date */
    virtual double getPVGivenDefault(
        const DateTime&                 valuationDate, 
        const DateTime&                 defaultDate,
        const IDiscountCurveRisky&      crv,
        IForwardRatePricerSP            model
        ) const = 0;

    /**Returns true if the instrument has a finite maturity. This will be false if
       it is perpetual.*/
    virtual bool hasFiniteMaturity(
        ) const = 0;

    /**Returns the maturity of the instrument. This is the last possible payment
       date, whether contingent on survival or the last possible date that default
       can cause a payment on default (i.e. if there is a recovery delay, 
       getMaturity() on a contingent leg should be the end of the protection period. 
       If (!hasFiniteMaturity()), an exception should be thrown.*/
    virtual DateTime getMaturity(
        ) const = 0;

    /**Returns the earliest possible "interesting" date: this will be the earliest of the
       start of the first accrual period, the start of the contingent leg, the first
       cash-flow, etc.*/
    virtual DateTime getStartDate(
        ) const = 0;

    /**Return the notional of the instrument. If the instrument is amortizing, this
       will be the base notional. The idea is that getPV(...)/getNotional() should be
       the price of the instrument, in some meaningful sense.*/
    virtual double getNotional(
        ) const = 0;

    /**Set the base notional of an instrument. This is here because some derivatives
       need to be able to control the notional of their underlyings to guarantee
       that the valuations make sense. In general the instrument will be loaded with
       the correct notional, so this method will not be commonly used.*/
    virtual void setNotional(
        double                          newNotional
        ) = 0;

    virtual CashFlowArraySP getInstrumentCashFlows(
        IForwardRatePricerSP            model) const = 0;

    /** Returns the day count convention used for accruals */
    virtual DayCountConventionSP getAccrualDcc() const = 0;

    static CClassConstSP const TYPE;
    virtual ~ICreditVanillaInstrument();
private:
    static void load(CClassSP& clazz);
};

DRLIB_END_NAMESPACE
#endif
