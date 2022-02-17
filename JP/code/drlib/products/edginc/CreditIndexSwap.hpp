//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : CreditIndexSwap.hpp
//
//   Description : Credit default Index Swap: Instrument equivalent to
//                 a CredDefSwap, but where the underlier is an index 
//                 (CDX or iTraxx style), rather than a single name.
//
//   Author      : Jose Hilera
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_CREDITINDEXSWAP_HPP
#define QLIB_CREDITINDEXSWAP_HPP

#include "edginc/Instrument.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/ClosedFormCDSBasket.hpp"
#include "edginc/Theta.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"

DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart) 
 * pointers and therefore their include files are not required here - they
 * will be in the source files though)*/
FORWARD_DECLARE_WRAPPER(CreditIndex);
FORWARD_DECLARE_WRAPPER(YieldCurve);
FORWARD_DECLARE(ICreditFeeLeg);
FORWARD_DECLARE(BadDayConvention);
FORWARD_DECLARE(ICreditEventOverride);
FORWARD_DECLARE(IForwardRatePricer);


/** Credit Index Swap (CIS) is an instrument used to buy/sell protection
 * on an index (ie, in all the component names of the index).
 * The swap is priced as a basket of CDS, where the underlying name in 
 * the CDSs are the (index-basis adjusted) underlying names in the index.
 * The index basis computation and adjustment is performed within the 
 * creditIndex field of the CIS.
 * The adjCurvesCache is used to cache the price of each (index-basis 
 * adjusted) CDS */

class PRODUCTS_DLL CreditIndexSwap: public CInstrument,
                       virtual public LastSensDate,
                       virtual public ClosedFormCDSBasket::IIntoProduct,
                       virtual public Theta::Shift
{
    /** Class to handle the interaction between model and instrument.
     * Note there are more friend declarations further down in this
     * class. */
    friend class CreditIndexSwapClosedFormCDSBasket;

public:
    static CClassConstSP const TYPE;

    virtual ~CreditIndexSwap();

    /** Allow the instrument to retrieve its market data */
    virtual void GetMarket(const IModel* model, const CMarketDataSP market);

    /** Implementation of ClosedFormCDSBasket::IntoProduct interface */
    virtual ClosedFormCDSBasket::IProduct* createProduct(ClosedFormCDSBasket* model) const;

    /** Notification that (some) underlying fields have changed */
    void fieldsUpdated(const CFieldArray& fields);

    /** Price the instrument */
    virtual void price(CResults* results, Control* control, IForwardRatePricerSP model) const;

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

    /** called after market data has been retrieved */
    virtual void Validate();

    /** Required part of CInstrument */
    virtual DateTime getValueDate() const;

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

    /** When to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    /** Required part of Theta::Shift */
    virtual bool sensShift(Theta* shift);

private:
    /** For reflection */
    CreditIndexSwap();
    static void load (CClassSP& clazz);
    static IObject* defaultCreditIndexSwap();

    // Used to store the name of the adjCurvesCache field (in case it changes)
    static const string adjCurvesCacheFieldName;

    /** Deals with the CIS output requests */
    void addOutputRequests(Control*        control,
                           Results*        results,
                           IForwardRatePricerSP model) const;

    /** When do payments occur */
    DateTimeArraySP paymentDates(IForwardRatePricerSP model) const;

    /** When do riskless payments occur */
    DateTimeArraySP risklessPaymentDates() const;

    /** When do risky payments occur */
    DateTimeArraySP riskyPaymentDates() const;

    /** Returns all known cash flows.
     * Caution: Due to the use of CashFlow::merge, this is very slow */
    CashFlowArraySP knownCashFlows(IForwardRatePricerSP model) const;


    // Fields
    DateTime protectionStartDate;
    DateTime protectionEndDate;
    bool     payAccruedFee; // Whether to pay accrued fee on default or not
    string   dcc;           // Day Count Convention
    double   notional;      // Original CIS notional, ie, before any defaults
    DateTime valueDate;     // Valuation date
    YieldCurveWrapper      discount;            // Currency to discount payoff
    CreditIndexWrapper     creditIndex;         // The index to swap
    ICreditFeeLegSP        feeLeg;              // Fee payments
    ICreditEventOverrideSP creditEventOverride; // Default parameters overrides - JLH should be ConstSP
    CIntSP triggerDelay; // Delay between default and eventDeterminationDate, in days
    CIntSP defaultToSettlementDelay; // Delay between credit event and settlement
    DateTime lastTriggerDate; // Last date when a default occurred during 
                              // the protection period can be triggered

    class PriceCache;


    // Now, this is a bit of a pain:
    // PriceCache should be private (to avoid other classes using it) but
    // PriceCacheArray needs to be public to be able to register the type.
    // Also, MsVC needs to know about the class before declaring it as
    // friend, so the friend declaration for priceCache has been moved down
    // here.
public:
    DECLARE(PriceCache);

    // PriceCache is a class used to cache the price of each of 
    // the (index basis adjusted) curves in the index we are swapping.
    // It needs access to the CIS instrument in order to know how to price each 
    // of the curves as a CDS, so it needs to be a friend of CIS - in fact, the
    // only reason why this is done in a class other than CIS itself is that
    // we need to know (via the fieldsUpdated method) when any of the curves 
    // have been updated; if the array of curves was held in CIS we could
    // not know exactly which element of the array (ie which curve) had
    // changed...
    friend class PriceCache;

private:
    // Transient fields
    mutable DayCountConventionSP swapDcc; // used if payAccruedFee = true 
    mutable PriceCacheArraySP    adjCurvesCache; // cache for curves prices
};

typedef smartPtr<CreditIndexSwap>      CreditIndexSwapSP;
typedef smartConstPtr<CreditIndexSwap> CreditIndexSwapConstSP;

DRLIB_END_NAMESPACE

#endif
