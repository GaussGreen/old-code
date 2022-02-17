//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : IDiscountCurve.hpp
//
//   Author      : Mark A Robson + modifications by Charles Morcom
//
//   Date        : January 10, 2006
//
//
//----------------------------------------------------------------------------

#ifndef QR_IDISCOUNTCURVE_HPP
#define QR_IDISCOUNTCURVE_HPP

#include "edginc/Object.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"

#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"

#include "edginc/GetMarket.hpp"

#ifdef DEBUG
#define EXPOSE_IDISCOUNT_CURVE_ADDINS
#endif

DRLIB_BEGIN_NAMESPACE
FORWARD_DECLARE_WRAPPER(IDiscountCurve)

/** Interface for classes that can do 'discounting'. This definition includes
    risky discounting. YieldCurves implement this interface but the methods
    in this interface are a reduced subset of the full yield curve since
    not all methods on YieldCurve would make sense for risky discount curves
    for example.

    The main idea behing the class is to provide an interface for discount
    curves that can be used when risky discounting without forcing a risky
    discount curve to be a yield curve. A risky discount curve which allows
    you to price protection legs and compute default probabilities is
    IDiscountCurveRisky, which extends this class.

    Another idea behing this class is that instruments could declare
    their discount curve (traditionally of type YieldCurve) to be of
    type IDiscountCurve.  Then a model could select/build an
    appropriate IDiscountCurve if desired (during the call to retrieve the
    market data) which would provide a way
    to support "riskified streams" which are priced in a normal way except a
    risky discount curve is used. Of course, this only handles the case where
    the default causes a simple knock-out of the instrument with no default 
    contingent payments.

    Extra methods might yet be needed on this interface (in particular
    a getName() may well be useful(?) although the definition may be
    unclear - the danger is that it would then be used to return
    results to the client (and if the name is a composite the client will not
    understand it). For now we are leaving this out.*/

class MARKET_DLL IDiscountCurve: public virtual IObject,
                      public virtual IGetMarket {
public:
    static CClassConstSP const TYPE;

    virtual ~IDiscountCurve();

    /** @return Discounting curve's currency - typically this is the
        ISO code eg "USD" (although it is up to the client - you may
        assume however that yield curves in the same currency return the
        same value in getCcy(). Note it is NOT the name of the yield
        curve (which might be eg "GBP-LIBOR"). */
    virtual string getCcy() const = 0;

    /** Compute discount factor between two dates. This is the price
     * at time date1 of a zero-coupon (discount) bond that pays 1
     * at time date2 if it still exists then. 
     * Note that, if this bond is risky, this method
     * in IDiscountCurveRisky means the PV of such a bond which knocks
     * out on default with no recovery at all, and the price is 
     * contingent on no default before date1.
     * @param date1 payment settlement date (conditional on existence at date1)
     * @param date2 payment/zero-coupon bond maturity date
     * @return Discount factor between date1 & date2
     */
    virtual double pv(const DateTime& date1, 
                      const DateTime& date2) const = 0;
    
    /** Compute price for settlement today of a zero-coupon bond
     * maturing at date. Note settlement is to TODAY, not to
     * the curve spot date. This is because some curves
     * may have ambiguous spot-dates - for example should a combined
     * credit and rates curve have spot date T+1 or T+2?
     * @param date To get discount factor/PV for
     * @return Discount factor between today & given date
     */
    virtual double pv(const DateTime& date) const = 0;
    
    /** Calculates present value to baseDate of supplied cash flows conditional
        on continued existence (i.e. no default for a risky curve)
        Cash-flows on or before baseDate are ignored. No
        ordering of the cashflows is assumed. */
    virtual double pv(const CashFlowArray& cashFlows,
                      const DateTime&      baseDate) const = 0;

    /** Interface to optimise repeated calculations of 'amounts' (eg
        rates, discount factors etc) when the calculations are close
        to each other in terms of dates */
    class MARKET_DLL IKey{
    public:
        /** Calculates the appropriate rate/factor between the two dates */
        virtual double calc(const DateTime&  loDate,
                            const DateTime&  hiDate) = 0;
        virtual ~IKey(){}
    };

    /** Returns a key used to optimise repeated calculations of discount
        factors/forward rate. The calc method for this key returns the 
        natural logarithm of the discount factor (or equivalently the
        product of the forward rate (continuous, Act/365F) and the negative
        year fraction (Act/365F) betweeen the two dates.
        The default implementation has no performance improvements. */
    virtual IKey* logOfDiscFactorKey() const = 0;

    /** return the bootstrapped dates */
    virtual DateTimeArray zeroDates() const = 0;

private:
    static void load(CClassSP& clazz);
};

DRLIB_END_NAMESPACE
#endif

/** Other changes:
1. Change IYieldCurve to derive [virtually] from IDiscountCurve as well as
[virtually] from IMarketFactor
2. Remove methods in IDiscountCurve from IYieldCurve
3. Remove IKey from IYieldCurve and replace with a typedef:
    typedef IDiscountCurve::IKey IKey;
4. Similary change similar typedef in YieldCurve
*/
