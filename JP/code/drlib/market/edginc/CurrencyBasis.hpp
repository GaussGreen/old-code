//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CurrencyBasis.hpp
//
//   Description : Captures ccy basis curve and conventions
//
//   Author      : Andrew J Swain
//
//   Date        : 6 July 2004
//
//
//----------------------------------------------------------------------------

#ifndef CURRENCY_BASIS_HPP
#define CURRENCY_BASIS_HPP

#include "edginc/MarketObject.hpp"
#include "edginc/WrapperNameCollector.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/TweakableWith.hpp"
#include "edginc/PointwiseTweakableWith.hpp"
#include "edginc/CurrencyBasisRhoParallelTweak.hpp"
#include "edginc/CurrencyBasisRhoPointwiseTweak.hpp"
#include "edginc/CurrencyBasisSpreadLevelTweak.hpp"


using namespace std;  // string

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL CurrencyBasis: 
        public MarketObject,
        virtual public TweakableWith<CurrencyBasisRhoParallelTweak>,
        virtual public PointwiseTweakableWith<CurrencyBasisRhoPointwiseTweak>,
        virtual public TweakableWith<CurrencyBasisSpreadLevelTweak> 
{
public:
    static CClassConstSP const TYPE;

    static const string SWAP_DAYS_ACTUAL;
    static const string SWAP_DAYS_BOND;

    virtual void validatePop2Object();

    // what's my name?
    virtual string getName() const;

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

    // return basis spread to add to a given cash/swap rate and benchmark
    double basis(double        rate,
                 const Expiry* expiry,
                 const string& rateType,
                 int           mmrtDenom,
                 int           mmrtFreq,
                 const string& swapDays,
                 int           swapDenom,
                 int           swapFreq) const;  // see YieldCurve

    // return basis spread to add to a given benchmark
    double cashBasis(const Expiry& expiry) const;
    double futuresBasis(const Expiry& expiry) const;
    double swapBasis(double                    rate,
                     const Expiry&             expiry,
                     const DayCountConvention& mmrtDcc,
                     int                       mmrtFreq,
                     const DayCountConvention& swapDcc,
                     int                       swapFreq) const;

    // is this the same as another CurrencyBasis
    bool equals(const CurrencyBasis* compare) const;
    
    /** Optimized hashCode for performance : use for caching only */
    int hashCodeOpt() const;

    // Rho Parallel tweak 
    virtual string sensName(CurrencyBasisRhoParallelTweak* shift) const;
    virtual bool sensShift(CurrencyBasisRhoParallelTweak* shift);

    // Rho Pointwise tweak
    virtual string sensName(CurrencyBasisRhoPointwiseTweak* shift) const;
    virtual bool sensShift(CurrencyBasisRhoPointwiseTweak* shift);
    virtual ExpiryArrayConstSP sensExpiries(CurrencyBasisRhoPointwiseTweak* shift) const;

    // Spread Level tweak 
    virtual string sensName(CurrencyBasisSpreadLevelTweak* shift) const;
    virtual bool sensShift(CurrencyBasisSpreadLevelTweak* shift);

	const HolidayWrapper& getHolidays() const;

private:
    friend class CurrencyBasisHelper;

    CurrencyBasis();
    CurrencyBasis(const CurrencyBasis &rhs);
    CurrencyBasis& operator=(const CurrencyBasis& rhs);

    // fields
    DateTime         today;    // will be populated / validated against the market cache
    string           name;     // identifier for market cache
    ExpiryArraySP    expiries;
    DoubleArray      spreads;
    HolidayWrapper   hols;     // introduced for zc3 - optional, no default
    bool             adjusted; // introduced for zc3 - optional, default false

};

typedef smartConstPtr<CurrencyBasis> CurrencyBasisConstSP;
typedef smartPtr<CurrencyBasis> CurrencyBasisSP;
#ifndef QLIB_CURRENCYBASIS_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<CurrencyBasis>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<CurrencyBasis>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<CurrencyBasis>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<CurrencyBasis>);
#endif

// support for wrapper class
typedef MarketWrapper<CurrencyBasis> CurrencyBasisWrapper;
#ifndef QLIB_CURRENCYBASIS_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<CurrencyBasis>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<CurrencyBasis>);
#endif

DRLIB_END_NAMESPACE
#endif
