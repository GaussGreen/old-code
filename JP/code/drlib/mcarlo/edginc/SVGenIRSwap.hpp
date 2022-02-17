//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenIRSwap.hpp
//
//   Description : A Generator of MC IR Swap State Variables
//
//   Date        : 3 June 2004
//
//
//----------------------------------------------------------------------------

#ifndef SVGenIRSwap_HPP
#define SVGenIRSwap_HPP

#include "edginc/StateVariableClient.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/SVGenExpectedDiscFactor.hpp"

DRLIB_BEGIN_NAMESPACE
class MCARLO_DLL SVGenIRSwap: public virtual IStateVariableGen,
                public virtual IStateVariableClient,
                public virtual VirtualDestructorBase
{
public:
    /** Interface for the state variable that SVGenIRSwap produces. This is
        the type that products deal with in the payoff. The payoff obtains
        it by calling the getSVDiscFactor() method below. Note support here
        for only one yield curve per state variable (unlike SVGenSpot). May
        need to reconsider at some point */
    class MCARLO_DLL IStateVar: public virtual IStateVariable{
    public:
        virtual ~IStateVar();
        /** Calculates the parYield for this swap together with a
            'spot annuity'. Uses the formula (Z0-Zn)/A, or
            sum(FwdLibor*Z)/A, depending on useGeneralFormula flag passed to
            SVGenIRSwap constructor. The parYield is what the fixed rate leg
            of the swap would be for the swap to value to par. The annuity
            is the discounted (to swap start) sum of the fixed coupons where
            a rate of 1.0 is used to calculate the coupons. */
        virtual void parYield(double& parYield, double& annuity) const = 0;
    };
    typedef smartPtr<IStateVar> IStateVarSP;
	typedef vector<IStateVarSP> StateVarArray;
	typedef refCountPtr<StateVarArray> StateVarArraySP;

    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector. Implementations typically call
        IStateVariableCollector::append */
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const;

    /** Constructor. If discountCurve is not null then the formula
        sum(FwdLibor*Z)/A is used otherwise (Z0-Zn)/A is used */
    SVGenIRSwap(
        YieldCurveConstSP couponCurve,   // for estimating coupons etc
        YieldCurveConstSP discountCurve,   /* optional, if non null then
                                              general [currency basis] formula
                                              is used */
        const DateTime&   swapStartDate, // start date of underlying swap
        const DateTime&   swapMatDate,   // maturity date of underlying swap
        string    fixedPayInterval,    // payment interval of fixed leg
        string    fixedDCC,            // day count convention of fixed leg
        string    stubType,            // type of front or back stub
        bool      stubAtEnd,           // only matters if stubType != NONE
        string    accrueBadDayConv,    //
        string    payBadDayConv,       //
        HolidaySP hols,
        bool      isCashSettled);      // T=PV wrt ytm, F=PV wrt zcurve

    /** Constructor. 
    Same as above, except allow calcDate to be before swapStartDate */
    SVGenIRSwap(
        YieldCurveConstSP couponCurve,   // for estimating coupons etc
        YieldCurveConstSP discountCurve,   /* optional, if non null then
                                           general [currency basis] formula
                                           is used */
        const DateTime&   calcDate,      // observation date of the swap yield
        const DateTime&   swapStartDate, // start date of underlying swap
        const DateTime&   swapMatDate,   // maturity date of underlying swap
        string    fixedPayPeriod,       // payment period of fixed leg, e.g. 1M, 2M, 1Q, 1Y
        string    fixedDCC,            // day count convention of fixed leg
        string    stubType,            // type of front or back stub
        bool      stubAtEnd,           // only matters if stubType != NONE
        string    accrueBadDayConv,    // 
        string    payBadDayConv,       // 
        HolidaySP hols,
        bool      isCashSettled);      // T=PV wrt ytm, F=PV wrt zcurve

    /** Constructor for precomputed swap dates.
        If useGeneralFormula is true then the formula
        sum(FwdLibor*Z)/A is used otherwise (Z0-Zn)/A is used */
    SVGenIRSwap(
        YieldCurveConstSP    couponCurve,   // for estimating coupons etc
        YieldCurveConstSP    discountCurve,   /* optional, if non null then
                                                 general [currency basis]
                                                 formula is used */
        const DateTime&      swapStartDate, // start date of swap
        const DateTimeArray& couponDates,   // coupon dates of swap
        DayCountConventionSP dcc);          // day count convention

    /** Constructor for precomputed swap dates.
    Same as above, except allow calcDate before swapStartDate */
    SVGenIRSwap(
        YieldCurveConstSP    couponCurve,   // for estimating coupons etc
        YieldCurveConstSP    discountCurve,   /* optional, if non null then
                                              general [currency basis] 
                                              formula is used */
        const DateTime&      calcDate,      // observation date of the swap yield
        const DateTime&      swapStartDate, // start date of swap
        const DateTimeArray& couponDates,   // coupon dates of swap
        DayCountConventionSP dcc);          // day count convention

    /** Create the corresponding State Variable for this State
        Variable Generator (from IStateVariableGen interface). The
        previous IStateVariableSP (may be null) should be passed in.  The
        return object may or may not be the same as oldStateVar. */
    virtual IStateVariableSP create(IStateVariableSP              oldStateVar,
                                 IStateVariableGen::IStateGen*  pathGen) const;

    /** Returns a MC IRSwap state variable which then
        provides access to simulated values etc. This is the method that
        products should call to get an SVGenIRSwap::IStateVar. */
    IStateVarSP getIRSwapSV(IStateVariableGen::IStateGen* pathGen) const;

private:
    void computeRatio(const DateTime&      swapStartDate,
                      const DateTimeArray& couponDates,
                      YieldCurveConstSP    couponCurve,
                      YieldCurveConstSP    discountCurve);
    class StateVarSimple;
    class StateVarGeneral;
    CashFlowArray  coupons;
    bool           useGeneralFormula;
    vector<double> ratio; // when useGeneralFormula is true
    SVGenExpectedDiscFactorSP mcDiscFactor; // this is our StateVar generator
    SVGenExpectedDiscFactorSP mcDFSwapStart; // df from swap start to calc date

    // Init coupon stream, create elementary sv gen, precomputation, etc
    // called by constructor
    void initialize(
        YieldCurveConstSP couponCurve,   // for estimating coupons etc
        YieldCurveConstSP discountCurve,   // optional
        const DateTime&   calcDate,      // observation date of the swap yield
        const DateTime&   swapStartDate, // start date of underlying swap
        const DateTime&   swapMatDate,   // maturity date of underlying swap
        string    fixedPayInterval,    // payment interval of fixed leg
        string    fixedDCC,            // day count convention of fixed leg
        string    stubType,            // type of front or back stub
        bool      stubAtEnd,           // only matters if stubType != NONE
        string    accrueBadDayConv,    // 
        string    payBadDayConv,       // 
        HolidaySP hols,
        bool      isCashSettled);

    // Populate coupon stream, create elementary sv gen, precomputation, etc
    // called by constructor
    void initialize(        
        YieldCurveConstSP    couponCurve,   // for estimating coupons etc
        YieldCurveConstSP    discountCurve,   // optional
        const DateTime&      calcDate,      // observation date of the swap yield
        const DateTime&      swapStartDate, // start date of swap
        const DateTimeArray& couponDates,   // coupon dates of swap
        DayCountConventionSP dcc);
};

typedef refCountPtr<SVGenIRSwap> SVGenIRSwapSP;
typedef vector<SVGenIRSwapSP> SVGenIRSwapArray;
typedef refCountPtr<SVGenIRSwapArray> SVGenIRSwapArraySP;

DRLIB_END_NAMESPACE

#endif
