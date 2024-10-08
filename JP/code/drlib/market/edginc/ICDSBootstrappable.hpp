//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : ICDSBootstrappable.hpp
//
//   Description : Interface representing 'bootstrappable' CDSParSpreads
//
//   Date        : 18 August 2005
//
//----------------------------------------------------------------------------

#ifndef ICDSBOOTSTRAPPABLE_HPP
#define ICDSBOOTSTRAPPABLE_HPP

#include "edginc/Expiry.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/ParSpreadCurve.hpp"
#include "edginc/WrapperNameCollector.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/ICDSParSpreads.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface for objects representing 'bootstrappable' CDSParSpreads:
 * allows us to write a generic 'bootstrap' method 
 * (see CDSHelper::bootstrap(ICDSBootstrappable)) for 
 * different implementations (eg. CDSParSpreads, AdjustedCDSParSpreads).
 * In particular, 'par spreads' can easily be 'adjusted' by
 * redefining getParSpreads() method. */

class MARKET_DLL ICDSBootstrappable : public virtual ICDSParSpreads {

public:
    virtual ~ICDSBootstrappable() {}
    static CClassConstSP const TYPE;

    /** Returns the corresponding discount curve */
    virtual YieldCurveConstSP getYieldCurve() const = 0;
    
    /** Returns the value date */
    virtual const DateTime& getValueDate() const = 0;
    
    /** Returns the 'effective date' = date used as reference to
     * compute cash flow dates */
    virtual const DateTime getEffDate() const = 0;
    
    /** Returns the expiries dates (absolute: eg. 27/02/2007) */
    virtual const DateTimeArray getExpiryDates() const = 0;
    
    /** Returns the expiries (relative: eg. 1M) */
    virtual ExpiryArrayConstSP getParSpreadsExpiries() const = 0;
    
    /** Returns cash flows corresponding to expiry = index */
    virtual CashFlowArraySP getCashFlowArray(int index) const = 0;

    /** Returns cash flows corresponding to expiry = index, adjusting
     * the spreads using the index basis */
    virtual CashFlowArraySP getCashFlowArray(CreditIndexBasisConstSP indexBasis,
                                             int index) const = 0;

    /** Returns the protection end date corresponding to expiry = index */
    virtual const DateTime getProtectionEndDate(int index) const = 0;

    /** Returns cash flow dates corresponding to expiry = index
     * (before bad days and holidays adjustment) */
    virtual const DateTimeArray getCashFlowDates(int index) const = 0;

    // Since bootstrapping is slow, all classes derived from ICDSBootstrappable
    // need to implement a cache. The following 3 methods deal with this cache

    /**
     * The par spreads kludged in a CashFlowArray suitable for storing in
     * Results
     *
     * Each (expiry, par spread) point on the ParSpreadCurve becomes a
     * (CashFlow::date, CashFlow::amount [sic]) pair in the CashFlowArray.  The
     * dates are generated by adding getParSpreadsExpiries() to getEffDate()
     * under getDayCountConv(), as as getExpiryDates().
     *
     * This is just a hack for supporting OutputRequest::PAR_SPREAD_CURVE.  If
     * the latter had a proper class to represent it, that might be a more
     * appropriate home for this method.  Or, it could be a function.
     */
    CashFlowArraySP asCashFlowArray() const;
};

typedef smartConstPtr<ICDSBootstrappable> ICDSBootstrappableConstSP;
typedef smartPtr<ICDSBootstrappable>      ICDSBootstrappableSP;
typedef MarketWrapper<ICDSBootstrappable> ICDSBootstrappableWrapper;

#ifndef QLIB_ICDSBOOTSTRAPPABLE_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<ICDSBootstrappable>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<ICDSBootstrappable>);
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<ICDSBootstrappable>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<ICDSBootstrappable>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<ICDSBootstrappable>);
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<ICDSBootstrappable>);
#endif

DRLIB_END_NAMESPACE

#endif
