//----------------------------------------------------------------------------
//
//   Group       : EDR
//
//   Filename    : DefaultRates.hpp
//
//   Description : Essentially a clean spread curve with explicit dates
//                 Moved here from CDSHelper
//
//   Author      : Andr�Segger
//
//   **CAUTIONS**: DefaultRates are immutable - see explanation in Expiry.hpp.
//                 Do not add any non-const methods!
//
//   Date        : 24 August 2004
//
//

#ifndef EDR_DEFAULTRATES_HPP
#define EDR_DEFAULTRATES_HPP

#include "edginc/YieldCurve.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE
/** To a credit default swap all that is needed in addition to the
    instrument definition is default rates. Default rates can either
    be implied from CDS par spreads or can come from a firm asset
    model (They can be found with other methods too, but this is all
    I've implemented). This abstract class defines the interface for an
    object from which you can compute the odds of default between two
    dates.  */
class MARKET_DLL DefaultRates {
public:
    virtual ~DefaultRates();

    /**************************************************************************
     //  CalcDefaultPV - calculates the odds you will receive a payment on the
     //  toDate from a risky counterpary conditioned on the fact that they are
     //  not in default on the fromDate. Defaults
     //  are assumed to occur overnight only (i.e. not continuously).
     *******************************************************************/
    virtual double calcDefaultPV(
        const DateTime&    fromDate,
        const DateTime&    toDate) const = 0;

    /** Interface to optimise repeated calculations of 'amounts' (eg
        default probabilities) when the calculations are close
        to each other in terms of dates */
    class MARKET_DLL IKey{
    public:
        /** Calculates the appropriate rate/factor between the two dates */
        virtual double calc(const DateTime&  loDate,
                            const DateTime&  hiDate) = 0;
        virtual ~IKey(){}
    };

    /** Returns a key used to optimise repeated calculations of
        default probabilities. The calc method for this key returns
        the natural logarithm of 1 - the default probability between
        the two dates (or equivalently the product of the forward
        default rate (continuous, Act/365F) and the negative year
        fraction (Act/365F) betweeen the two dates. Or another way of
        saying this is to say that it returns minus the integral of
        the forward default rate. The default implementation has no
        performance improvements. */
    virtual IKey* logOfDefaultPVKey() const;

    /** calculates the probability of a default due to a jump-to-zero event.
        Default implementation here returns 1 */
    virtual double calcJumpPV(
        const DateTime&      fromDate,
        const DateTime&      toDate) const;

    /** A bit undocumented - seems to combine calcDefaultPV and calcJumpPV */
    virtual double calcTotalDefaultPV(
        const DateTime&     fromDate,
        const DateTime&     toDate) const;

    /** Returns the [forward] 'clean spreads' which, loosely speaking are
        related to probability of default by
        exp(- integral of clean spread curve across time) */
    virtual CashFlowArraySP getCleanSpreadCurve() const = 0;

    /** Returns just the dates aspect of getCleanSpreadCurve */
    virtual DateTimeArraySP getDefaultDates() const = 0;

    /** Calculates an implied par spread for a spot or forward CDS */
    double cdsParSpread(
        const DateTime&                effDate,
        int                            frequency,
        double                         notional,
        double                         recovery,
        bool                           accrueFee,
        const DateTime&                accruedEffectiveDate,
        const DateTime&                maturity,
        YieldCurveConstSP              discount,
        bool                           isE2C,
        const DayCountConvention*      dcc) const; // should this be in DefaultRates

    /** Calculate spot rates
     * (implicitly from forward rates inside DefaultRates).
     * It would probably be preferable to explicitly have a
     * getFwdRates() and a getSpotRates() method (or a more
     * generic getRates(RateType type) methods). */
    CashFlowArraySP convertToSpot(
        const DateTime& valueDate,
        const DayCountConvention& dcc) const;

    /** Returns (forward) rates values */
    virtual const CDoubleArray& getRates() const = 0;

    /** Returns default dates values */
    virtual const DateTimeArray& getDates() const = 0;

    /** Returns a reference to the value date - allows default implementations
        to work */
    virtual const DateTime& getValueDate() const = 0;
private:
    class DefaultLogOfPVKey;
};

DECLARE_REF_COUNT(DefaultRates);

DRLIB_END_NAMESPACE
#endif
