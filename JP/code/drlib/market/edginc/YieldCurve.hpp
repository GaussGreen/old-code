//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : YieldCurve.hpp
//
//   Description : Prototype yield curve interface
//
//   Author      : Andrew J Swain
//
//   Date        : 5 November 2000
//
//
//----------------------------------------------------------------------------


#ifndef YIELDCURVE_HPP
#define YIELDCURVE_HPP

#include "edginc/MarketFactor.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/IDiscountCurve.hpp"

DRLIB_BEGIN_NAMESPACE
class IVolProcessed;
class CVolRequest;
class CreditSpreadCurve;
class ICDSParSpreads;
class DefaultRates;
class Stub;
typedef smartPtr<CreditSpreadCurve> CreditSpreadCurveSP;

class IYieldCurve;
typedef smartConstPtr<IYieldCurve> IYieldCurveConstSP;
typedef smartPtr<IYieldCurve> IYieldCurveSP;

/** There are many types of yield curve - cash/swap, FRA and so on,
all with different data requirements, but fundamentally they can all
do discounting and rate interpolation
    <BR><BR>
*/
class MARKET_DLL IYieldCurve: public IMarketFactor,
                    virtual public IDiscountCurve {
public:
    static CClassConstSP const TYPE;

    virtual ~IYieldCurve();

    /** @return Yield curve's currency */
    virtual string getCcy() const = 0;
    
    /** @return Yield curve's name - used to identify sensitivities */
    virtual string getName() const = 0;

    /** @return Yield curve's spot date */
    virtual DateTime getSpotDate() const = 0;

    /** @return Yield curve's today date */
    virtual DateTime getToday() const = 0;

    /** Useful accessor methods */
    virtual ExpiryArrayConstSP getExpiries() const = 0;
    virtual StringArrayConstSP getInstruments() const = 0;

    /** Returns tradeDate + n business days where n = spot offset */
    virtual DateTime settles(const DateTime& tradeDate) const = 0;

    /** Optimized equalTo for performance: use for caching only */
    virtual bool zeroCurveEquals(const IYieldCurve* yc2) const = 0;

    /** Optimized hashCode for performance: use for caching only */
    virtual int zeroCurveHash() const = 0;

    /** Compute discount factor between two dates
     * @param lodate Lower date
     * @param hidate Upper date
     * @return Discount factor between lodate & hidate
     */
    virtual double pv(const DateTime& lodate, 
                      const DateTime& hidate) const = 0;
    
    /** Compute discount factor between value date and a date
     * @param date To get discount factor for
     * @return Discount factor between value date & given date
     */
    virtual double pv(const DateTime& date) const = 0;

    
    double discountFactor(const DateTime& date) const 
    {
        return pv(date);
    }
    
    double pv(const CashFlowArray& cashFlows) const
    {
        return pv(cashFlows, getToday());
    }

    /** Calculates present value to baseDate of supplied cash flows. Cash
        flows on or before baseDate (eg value date) are ignored. No
        ordering of the cashflows is assumed */
    virtual double pv(const CashFlowArray& cashFlows,
                      const DateTime&      baseDate) const = 0;

    /** Calculates present value to baseDate of supplied cash flows. Cash
        flows on or before baseDate (eg value date) are ignored depending
        on the optional parameter ignoreCashFlowsOnBaseDate. No
        ordering of the cashflows is assumed */
    virtual double pv(const CashFlowArray& cashFlows,
                      const DateTime& baseDate,
                      const bool ignoreCashFlowsOnBaseDate) const = 0;

    /** Interpolate zero coupon rate at a date
     * @param date Interpolate zero coupon rate here
     * @return Zero coupon rate at given date
     */
    virtual double zero(const DateTime& date) const = 0; 

    /** Returns rateMaturity->toDate(rateStart) bad day adjusted
        either using bad day convention supplied (if non null)
        together with the yield curve holidays or the yield curve's
        bad day convention and holidays */
    virtual DateTime rateMaturity(
        const DateTime&         rateStart,
        const Expiry*           rateMaturity,
        const BadDayConvention* bdc) const = 0; // optional

    /** Interpolate a forward rate between two dates 
     * see CompoundBasis for basis values
     */
    virtual double fwd(const DateTime&           lodate, 
                       const DateTime&           hidate,
                       const DayCountConvention* dcc,
                       int                       basis) const = 0;

    /** Same as fwd method above but dates for rates take into account
        spot offset, holidays etc. */
    virtual double fwd(const DateTime&           refixDate, 
                       const Expiry*             rateMaturity,
                       const BadDayConvention*   bdc, // optional
                       const DayCountConvention* dcc,
                       int                       basis) const = 0;
    
    double parSwapRate(
        const DateTime&           startDate, 
        const DateTime&           endDate,
        const MaturityPeriod&     period, 
        const DayCountConvention& dcc,
        const Stub&               stubType,
        bool                      stubAtEnd,
        const BadDayConvention&   accBadDayConv,
        const BadDayConvention&   payBadDayConv,
        const Holiday&            holidays) const ;

	/** Calculates a par swap rate for the swap starting at startDate,
     * maturing at maturityDate, with fixed day count convention 
     * fixedDayCountConv and fixed payments occuring at time 
     * intervals defined by interval. In other words, the routine calculates 
     * the fixed rate such that the present value of the fixed and 
     * floating sides are equal. The floating side is assumed to be at par.
     */
    virtual double couponRate(
        const DateTime&           startDate,   //(I) Date instrument begins at
        const DateTime&           maturityDate,//(I) Date instrument matures at
        const MaturityPeriod&     interval,    // (I) Time between payments 
        bool                      stubAtEnd,
        const DayCountConvention* dcc) const = 0;

	/** Calculates a par swap rate for the swap starting at curve base date,
     * maturing at maturity date, using the curve's fixed leg and floating 
     * leg conventions.  The floating leg is assumed to be at par, except if 
     * the estimating zero curve differs from the discounting, in which case 
     * the floating leg is explictly valued.
     */
    virtual double parSwapRate(const DateTime& maturity) const;

    /** Interpolate a futures rate between two dates 
     * see CompoundBasis for basis values
     */
    virtual double future(const DateTime&           lodate, 
                          const DateTime&           hidate,
                          const DayCountConvention* dcc,
                          int                       basis,
                          double                    irVol) const = 0;

    /** Interface to optimise repeated calculations of 'amounts' (eg
        rates, discount factors etc) when the calculations are close
        to each other in terms of dates */
//    class IKey{
//    public:
//        /** Calculates the appropriate rate/factor between the two dates */
//        virtual double calc(const DateTime&  loDate,
//                            const DateTime&  hiDate) = 0;
//        virtual ~IKey(){}
//    };

    /** Returns a key used to optimise repeated calculations of discount
        factors/forward rate. The calc method for this key returns the 
        natural logarithm of the discount factor (or equivalently the
        product of the forward rate (continuous, Act/365F) and the negative
        year fraction (Act/365F) betweeen the two dates.
        The default implementation has no performance improvements. */
//    virtual IKey* logOfDiscFactorKey() const = 0;

    /** grab the dates used in the zero curve */
    virtual DateTimeArray zeroDates() const = 0;

    /** grab the dates used in the zero curve and corresponding rates */
    virtual CashFlowArraySP getRatesAndDates() const = 0;

    /** 
     *  Return dates and rates from zero curve, specifying additional date
     *  intervals to be used in order to add extra points to the output list
     *  of dates and rates.  
     *  
     *  If requiredDates is NULL the critical dates are included.
     *  startPeriods and addIntervals must be the same size (0 is permitted).
     */
    CashFlowArraySP getRatesAndDates(
        const ExpiryArray&         startPeriods, 
        const MaturityPeriodArray& addIntervals, 
        const DateTimeArray*       requiredDates) const;

    /** drive which style of zero curve is used. Default is discounting */
    virtual void setProjectionCurve(bool useEstimatingCurve = true) const = 0;

   /** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest */
    virtual IVolProcessed* getProcessedVol(
        const CVolRequest* volRequest) const = 0;

	/** make a risky curve from a credit spread curve */
    virtual IYieldCurveSP makeRiskyCurve(
        const CreditSpreadCurve& spreadCurve,
        const DateTime*          maturityDate = NULL) const = 0;

	virtual CreditSpreadCurveSP makeCreditSpreadCurve(
		const string&        name,
		const CashFlowArray& defaultRates,
		double               recovery) const = 0;

    virtual IYieldCurveSP createForwardCurve(const DateTime& forwardDate) const = 0;

private:
    static void load(CClassSP& clazz);
};

#ifndef QLIB_YIELDCURVE_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<IYieldCurve>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<IYieldCurve>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<IYieldCurve>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<IYieldCurve>);
#endif

class MARKET_DLL YieldCurve: public MarketObject,
                  public virtual IYieldCurve {
public:
    static CClassConstSP const TYPE;
    friend class YieldCurveHelper;
    using IYieldCurve::fwd;

    // identifies money market or swap rates
    static const string MMRT_RATE;
    static const string MMRT_RATE2;
    static const string SWAP_RATE;
    static const string SWAP_RATE2;
    static const string FUTURE_RATE;
    static const string TURN_RATE;
    static const string BRAZIL_FUTURE_PRICE;

    virtual ~YieldCurve();

    ///////////////////////////////////////////////
    //// Note many methods now on IYieldCurve ////
    ///////////////////////////////////////////////

    /** @return Yield curve's name - used to identify sensitivities */
    virtual string getName() const = 0;

    /** Optimized equalTo for performance: use for caching only */
    virtual bool zeroCurveEquals(const IYieldCurve* yc2) const;

    /** Optimized hashCode for performance: use for caching only */
    virtual int zeroCurveHash() const;

    virtual double pv(const DateTime& lodate, 
                      const DateTime& hidate) const = 0;
    
    /** Compute discount factor between value date and a date
     * @param date To get discount factor for
     * @return Discount factor between value date & given date
     */
    virtual double pv(const DateTime& date) const = 0;
    
    /** Calculates present value to baseDate of supplied cash flows. Cash
        flows on or before baseDate (eg value date) are ignored. No
        ordering of the cashflows is assumed */
    virtual double pv(const CashFlowArray& cashFlows,
                      const DateTime&      baseDate) const;

    /** Calculates present value to baseDate of supplied cash flows. Cash
        flows on or before baseDate (eg value date) are ignored depending
        on the optional parameter ignoreCashFlowsOnBaseDate. No
        ordering of the cashflows is assumed */
    virtual double pv(const CashFlowArray& cashFlows,
                      const DateTime& baseDate,
                      const bool ignoreCashFlowsOnBaseDate) const;

    /** Calculates a par swap rate for the swap starting at startDate,
     * maturing at maturityDate, with fixed day count convention 
     * fixedDayCountConv and fixed payments occuring at time 
     * intervals defined by interval. In other words, the routine calculates 
     * the fixed rate such that the present value of the fixed and 
     * floating sides are equal. The floating side is assumed to be at par.
     */
    virtual double couponRate(
        const DateTime&           startDate,   //(I) Date instrument begins at
        const DateTime&           maturityDate,//(I) Date instrument matures at
        const MaturityPeriod&     interval,    // (I) Time between payments 
        bool                      stubAtEnd,
        const DayCountConvention* dcc) const;
                              
    /** Interpolate a futures rate between two dates 
     * see CompoundBasis for basis values
     */
    virtual double future(const DateTime&           lodate, 
                          const DateTime&           hidate,
                          const DayCountConvention* dcc,
                          int                       basis,
                          double                    irVol) const;

    //typedef IYieldCurve::IKey IKey;

    /** Returns a key used to optimise repeated calculations of discount
        factors/forward rate. The calc method for this key returns the 
        natural logarithm of the discount factor (or equivalently the
        product of the forward rate (continuous, Act/365F) and the negative
        year fraction (Act/365F) betweeen the two dates.
        The default implementation has no performance improvements. */
    virtual IDiscountCurve::IKey* logOfDiscFactorKey() const;

	// these to move to RiskyCurve
    /** this function calculates the discount factor based on the assumption
        that the on which the recovery is based is provided externally. 
        This allows to use different methodologies 
        (PV, face value + accrued etc.) to be included easily */
    virtual double riskyPV(const DateTime& lodate,
                           const DateTime& hidate,
                           double          cashFlow,
                           double          recoveryNotional) const = 0;

    /** this function calculates the discount factor based on the assumption
        that the on which the recovery is based is provided externally. 
        This allows to use different methodologies (PV, face value + accrued etc.)
        to be included easily  - this function will use the externally given
        recovery rather than the underlying risky curve's recovery */
    virtual double riskyPV(const DateTime& lodate,
                           const DateTime& hidate,
                           double          cashFlow,
                           double          recoveryNotional,
                           bool            useAssetRecovery,
                           double          assetRecovery) const;

    /** drive which style of zero curve is used. Default is discounting */
    // nasty default for time being
    virtual void setProjectionCurve(bool useEstimatingCurve = true) const = 0;

    /** Just does fwd(settles(start), rateMaturity(settles(start), end, bdc),
        dcc, basis) */
    virtual double fwd(const DateTime&           refixDate, 
                       const Expiry*             rateMat,
                       const BadDayConvention*   bdc, // optional
                       const DayCountConvention* dcc,
                       int                       basis) const;

    virtual double fwd(const DateTime&           refixDate, 
                       const Expiry*             rateMat,
                       const BadDayConvention*   bdc, // optional
                       const DayCountConvention* dcc,
                       const bool                isCMS) const = 0;

    virtual double fwd(const DateTime&           payDate,
                       const DateTime&           refixDate, 
                       const Expiry*             rateMat,
                       const BadDayConvention*   bdc, // optional
                       const DayCountConvention* dcc,
                       const bool                isCMS) const = 0;


protected:
    YieldCurve(CClassConstSP clazz);
private:
    static void load(CClassSP& clazz);
    
};

typedef smartConstPtr<YieldCurve> YieldCurveConstSP;
typedef smartPtr<YieldCurve> YieldCurveSP;
#ifndef QLIB_YIELDCURVE_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<YieldCurve>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<YieldCurve>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<YieldCurve>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<YieldCurve>);
#endif

// support for wrapper class
typedef MarketWrapper<YieldCurve> YieldCurveWrapper;
typedef smartPtr<YieldCurveWrapper> YieldCurveWrapperSP;
#ifndef QLIB_YIELDCURVE_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<YieldCurve>);
EXTERN_TEMPLATE(IObjectSP MARKET_DLL FieldGetInLine<YieldCurveWrapper>(
                    YieldCurveWrapper* t));
EXTERN_TEMPLATE(void MARKET_DLL FieldSetInLine<YieldCurveWrapper>(YieldCurveWrapper* t,
                                                       IObjectSP o));
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<YieldCurve>);
INSTANTIATE_TEMPLATE(IObjectSP MARKET_DLL FieldGetInLine<YieldCurveWrapper>(
                         YieldCurveWrapper* t));
INSTANTIATE_TEMPLATE(void MARKET_DLL FieldSetInLine<YieldCurveWrapper>(
                         YieldCurveWrapper* t, IObjectSP o));
#endif

// support for array of yield curves
typedef array<YieldCurveSP, YieldCurve> YieldCurveArray;
typedef smartPtr<YieldCurveArray> YieldCurveArraySP;
#ifndef QLIB_YIELDCURVE_CPP
EXTERN_TEMPLATE(class MARKET_DLL array<YieldCurveSP _COMMA_ YieldCurve>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<YieldCurveArray>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL array<YieldCurveSP _COMMA_ YieldCurve>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<YieldCurveArray>);
#endif

// arrays of wrappers (note array of structures)
typedef array<YieldCurveWrapperSP, YieldCurveWrapper> YieldCurveWrapperArray;
typedef smartPtr<YieldCurveWrapperArray> YieldCurveWrapperArraySP;
#ifndef QLIB_YIELDCURVE_CPP
EXTERN_TEMPLATE(class MARKET_DLL array<YieldCurveWrapperSP _COMMA_ YieldCurveWrapper>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<YieldCurveWrapperArray>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL array<YieldCurveWrapperSP _COMMA_ YieldCurveWrapper>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<YieldCurveWrapperArray>);
#endif

DRLIB_END_NAMESPACE

#endif
