//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenIRSwap.cpp
//
//   Description : A Generator of MC IR Swap State Variables
//
//   Date        : 3 June 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SVGenIRSwap.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/StubFactory.hpp"
#include "edginc/SwapTool.hpp"


DRLIB_BEGIN_NAMESPACE
SVGenIRSwap::IStateVar::~IStateVar(){}

/** Appends 'true' (ie non derived) state variable generators
    required to the supplied collector. Implementations typically call
    IStateVariableCollector::append */
void SVGenIRSwap::collectStateVars(IStateVariableCollectorSP svCollector) const{
    svCollector->append(mcDiscFactor.get());
    svCollector->append(mcDFSwapStart.get());
}

/** Constructor. If useGeneralFormula is true then the formula
    sum(FwdLibor*Z)/A is used otherwise (Z0-Zn)/A is used */
SVGenIRSwap::SVGenIRSwap(
    YieldCurveConstSP couponCurve,   // for estimating coupons etc
    YieldCurveConstSP discountCurve,   /* optional, if non null then
                                          general [currency basis] formula
                                          is used */
    const DateTime&   swapStartDate, // start date of underlying swap
    const DateTime&   swapMatDate,   // maturity date of underlying swap
    string    fixedPayInterval,      // payment interval of fixed leg
    string    fixedDCC,              // day count convention of fixed leg
    string    stubType,              // type of front or back stub
    bool      stubAtEnd,             // only matters if stubType != NONE
    string    accrueBadDayConv,      //
    string    payBadDayConv,         //
    HolidaySP hols,
    bool      isCashSettled):        // T=PV wrt ytm, F=PV wrt zcurve
    useGeneralFormula(discountCurve.get() != 0){
    try  {
        MaturityPeriod fixedPayPeriod(1, fixedPayInterval);
        // Init coupon stream, create elementary sv gen and precomputation...
        initialize(couponCurve, discountCurve, swapStartDate, swapStartDate, 
            swapMatDate, fixedPayPeriod.toString(), fixedDCC, stubType, 
            stubAtEnd, accrueBadDayConv, payBadDayConv, hols, isCashSettled);
    } 
    catch (exception& e){
        throw ModelException(e, "SVGenIRSwap::SVGenIRSwap",
            "Failed to build coupons for swap");
    }
}

SVGenIRSwap::SVGenIRSwap(
    YieldCurveConstSP couponCurve,   // for estimating coupons etc
    YieldCurveConstSP discountCurve,   /* optional, if non null then
                                        general [currency basis] formula
                                        is used */
    const DateTime&   calcDate,      // observation date of the swap yield
    const DateTime&   swapStartDate, // start date of underlying swap
    const DateTime&   swapMatDate,   // maturity date of underlying swap
    string    fixedPayPeriod,        // payment interval of fixed leg, e.g. 1S, 1Q, 2Q
    string    fixedDCC,              // day count convention of fixed leg
    string    stubType,              // type of front or back stub
    bool      stubAtEnd,             // only matters if stubType != NONE
    string    accrueBadDayConv,      //
    string    payBadDayConv,         //
    HolidaySP hols,
    bool      isCashSettled):        // T=PV wrt ytm, F=PV wrt zcurve
    useGeneralFormula(discountCurve.get() != 0){
    try  {
        // init coupon stream, elementary sv gen and precomputation...
        initialize(couponCurve, discountCurve, calcDate, swapStartDate, 
            swapMatDate, fixedPayPeriod, fixedDCC, stubType, stubAtEnd, 
            accrueBadDayConv, payBadDayConv, hols, isCashSettled);
    } 
    catch (exception& e){
        throw ModelException(e, "SVGenIRSwap::SVGenIRSwap",
            "Failed to build coupons for swap");
    }
}

void SVGenIRSwap::computeRatio(const DateTime&      swapStartDate,
                            const DateTimeArray& couponDates,
                            YieldCurveConstSP    couponCurve,
                            YieldCurveConstSP    discountCurve){
    // cache some ratios of discount factors
    if (useGeneralFormula){
        // from irdiffuse.c:SwapYield_s
        ratio = vector<double>(couponDates.size());
        for (unsigned int i = 0; i < ratio.size(); i++){
            const DateTime& start = i == 0? swapStartDate: couponDates[i-1];
            ratio[i] = discountCurve->pv(start, couponDates[i])/
                couponCurve->pv(start, couponDates[i]);
        }
    }
}

/** Constructor for precomputed swap dates.
    If a discountCurve is supplied then the formula
    sum(FwdLibor*Z)/A is used otherwise (Z0-Zn)/A is used */
SVGenIRSwap::SVGenIRSwap(
    YieldCurveConstSP    couponCurve,   // for estimating coupons etc
    YieldCurveConstSP    discountCurve,   /* optional, if non null then
                                             general [currency basis] formula
                                             is used */
    const DateTime&      swapStartDate, // start date of swap
    const DateTimeArray& couponDates,   // coupon dates of swap
    DayCountConventionSP dcc):          // day count convention
    coupons(couponDates.size()),
        useGeneralFormula(discountCurve.get() != 0){
    try {
        // Populate coupon stream, create elementary sv gen, precomputation, etc
        initialize(couponCurve, discountCurve, swapStartDate, 
            swapStartDate, couponDates, dcc);
    }
    catch (exception& e){
        throw ModelException(e, "SVGenIRSwap::SVGenIRSwap",
            "Failed to build coupons for swap");
    }
}

/** Constructor for precomputed swap dates.
If a discountCurve is supplied then the formula
sum(FwdLibor*Z)/A is used otherwise (Z0-Zn)/A is used */
SVGenIRSwap::SVGenIRSwap(
    YieldCurveConstSP    couponCurve,   // for estimating coupons etc
    YieldCurveConstSP    discountCurve,   /* optional, if non null then
                                            general [currency basis] formula
                                            is used */
    const DateTime&      calcDate,      // observation date of the swap yield
    const DateTime&      swapStartDate, // start date of swap
    const DateTimeArray& couponDates,   // coupon dates of swap
    DayCountConventionSP dcc):          // day count convention
        useGeneralFormula(discountCurve.get() != 0){
    try {
        // Populate coupon stream, create elementary sv gen, precomputation, etc
        initialize(couponCurve, discountCurve, calcDate, 
            swapStartDate, couponDates, dcc);
    }
    catch (exception& e){
        throw ModelException(e, "SVGenIRSwap::SVGenIRSwap",
            "Failed to build coupons for swap");
    }
}

/** Create the corresponding State Variable for this State
    Variable Generator (from IStateVariableGen interface). The
    previous IStateVariableSP (may be null) should be passed in.  The
    return object may or may not be the same as oldStateVar. */
IStateVariableSP SVGenIRSwap::create(IStateVariableSP              oldStateVar,
                               IStateVariableGen::IStateGen*  pathGen) const{
    return getIRSwapSV(pathGen);
}

class SVGenIRSwap::StateVarSimple: public IStateVar{
protected:
    SVExpectedDiscFactorSP discFactors;
    SVExpectedDiscFactorSP dfSwapStart;
    vector<double>            coupons; // for a rate of 1.0
public:
    virtual ~StateVarSimple(){}

    /** Returns true if this state variable is being used to 'simulate'
        the past. This is a useful method for users of state variables -
        need to see how hard it is to implement */
    bool doingPast() const{
        return discFactors->doingPast();
    }

    /** Calculates the parYield for this swap together with a
        'spot annuity'. Uses the formula (Z0-Zn)/A, not
        sum(FwdLibor*Z)/A, therefore may be inappropriate in the
        presence of ccy basis. The parYield is what the fixed rate leg
        of the swap would be for the swap to value to par. The annuity
        is the discounted (to swap start) sum of the fixed coupons where
        a rate of 1.0 is used to calculate the coupons. From
        irdiffuse.c::ParYield_s  (this is slightly simpler as the calculation
        date for the expected discount factors is [hard-coded] to the swap
        start date) */
    virtual void parYield(double& parYield, double& annuity) const{
        // get access to the expected discount factors;
        const SVPath& path = discFactors->path();
        double theAnnuity = 0.0;
        int endIdx = path.end();
        double zeroToCoupon=0;

        /*
        double firstZeroToCoupon = 1.0;
        for (int i = 0; i < endIdx; i++){
            zeroToCoupon = path[i]; // from swap start to coupon
            theAnnuity += coupons[i] * zeroToCoupon;
        }
        */
        double firstZeroToCoupon = dfSwapStart->path()[0];
        for (int i = 0; i < endIdx; i++){
            zeroToCoupon = path[i]; // from swap start to coupon
            theAnnuity += coupons[i] * zeroToCoupon;
        }

        parYield = (firstZeroToCoupon - zeroToCoupon)/theAnnuity;
        annuity = theAnnuity;
    }
    //// constructor
    StateVarSimple(SVExpectedDiscFactorSP discFactors,
                   SVExpectedDiscFactorSP dfSwapStart,
                   const CashFlowArray&      cfCoupons):
        discFactors(discFactors), 
        dfSwapStart(dfSwapStart), 
        coupons(cfCoupons.size()){
        for (int i = 0; i < cfCoupons.size(); i++){
            coupons[i] = cfCoupons[i].amount;
        }
    }
};

class SVGenIRSwap::StateVarGeneral: public StateVarSimple{
    vector<double>            ratio;   // between zeros off two curves
public:
    virtual ~StateVarGeneral(){}

    /** Calculates the parYield for this swap together with a
        'spot annuity'. Uses the formula
        sum(FwdLibor*Z)/A, which is the correct formula to use in the
        presence of ccy basis. The parYield is what the fixed rate leg
        of the swap would be for the swap to value to par. The annuity
        is the discounted (to swap start) sum of the fixed coupons where
        a rate of 1.0 is used to calculate the coupons. From
        irdiffuse.c::SwapYield_s (this is slightly simpler as the calculation
        date for the expected discount factors is [hard-coded] to the swap
        start date) */
    virtual void parYield(double& parYield, double& annuity) const{
        // get access to the expected discount factors;
        const SVPath& path = discFactors->path();
        double theAnnuity = 0.0;
        int endIdx = path.end();
        double  sumCoupon = 0.0;    /* sum(Libor*Zero) */

        /*
        double  prevZeroToCoupon = 1.0;
        for (int i = 0; i < endIdx; i++){
            double zeroToCoupon = path[i]; // from swap start to coupon
            sumCoupon  += prevZeroToCoupon * ratio[i] - zeroToCoupon;
            theAnnuity += coupons[i] * zeroToCoupon;
            prevZeroToCoupon = zeroToCoupon;
        }
        */
        double  prevZeroToCoupon = dfSwapStart->path()[0];
        for (int i = 0; i < endIdx; i++){
            double zeroToCoupon = path[i]; // from swap start to coupon
            sumCoupon  += prevZeroToCoupon * ratio[i] - zeroToCoupon;
            theAnnuity += coupons[i] * zeroToCoupon;
            prevZeroToCoupon = zeroToCoupon;
        }

        parYield = sumCoupon/theAnnuity;
        annuity = theAnnuity;
    }
    //// constructor
    StateVarGeneral(SVExpectedDiscFactorSP discFactors,
                    SVExpectedDiscFactorSP dfSwapStart,
                    const CashFlowArray&      cfCoupons,
                    const vector<double>&     ratio):
        StateVarSimple(discFactors, dfSwapStart, cfCoupons), ratio(ratio){}
};

/** Returns a MC IRSwap state variable which then
    provides access to simulated values etc. This is the method that
    products should call to get an SVGenIRSwap::IStateVar. */
SVGenIRSwap::IStateVarSP SVGenIRSwap::getIRSwapSV(
    IStateVariableGen::IStateGen* pathGen) const{
    IStateVariableSP sv = pathGen->create(mcDiscFactor.get());
    IStateVariableSP svSwapStart = pathGen->create(mcDFSwapStart.get());
    SVExpectedDiscFactorSP discFactors(& dynamic_cast<SVExpectedDiscFactor &>(*sv));
    SVExpectedDiscFactorSP dfSwapStart(& dynamic_cast<SVExpectedDiscFactor &>(*svSwapStart));

    return IStateVarSP(useGeneralFormula?
                       new StateVarGeneral(discFactors, dfSwapStart, coupons, ratio):
                       new StateVarSimple(discFactors, dfSwapStart, coupons));
}

void SVGenIRSwap::initialize(
    YieldCurveConstSP couponCurve,   // for estimating coupons etc
    YieldCurveConstSP discountCurve,   // optional
    const DateTime&   calcDate,      // observation date of the swap yield
    const DateTime&   swapStartDate, // start date of underlying swap
    const DateTime&   swapMatDate,   // maturity date of underlying swap
    string    fixedPayPeriod,       // payment period of fixed leg, e.g. 1M, 2M, 1Q, 1S
    string    fixedDCC,            // day count convention of fixed leg
    string    stubType,            // type of front or back stub
    bool      stubAtEnd,           // only matters if stubType != NONE
    string    accrueBadDayConv,    // 
    string    payBadDayConv,       // 
    HolidaySP hols,
    bool      isCashSettled)
{
    DayCountConventionSP dcc(DayCountConventionFactory::make(fixedDCC));
    StubSP               stub(StubFactory::make(stubType));
    BadDayConventionSP   accrueBDC(
        BadDayConventionFactory::make(accrueBadDayConv));
    BadDayConventionSP payBDC(BadDayConventionFactory::make(payBadDayConv));

    int count;
    string fixedPayInterval;
    MaturityPeriod fixedPayPrd(fixedPayPeriod);
    fixedPayPrd.decompose(count, fixedPayInterval);

    // build the swap coupon cash flow
    coupons = SwapTool::cashflows(swapStartDate,
        swapMatDate,
        stub.get(),
        stubAtEnd,
        accrueBDC.get(),
        payBDC.get(),
        hols.get(),
        false, // dump start date
        false,  // don't include principal
        1.0, // rate
        count, // 1
        fixedPayInterval, // e.g. Y, M, W, D //fixedPayInterval
        dcc.get());

    // create the coupon date arrays
    DateTimeArray couponDates(CashFlow::dates(coupons));

    // then build our MCDiscFactor state variable generator
    mcDiscFactor = SVGenExpectedDiscFactorSP(
        new SVGenExpectedDiscFactor(
        calcDate, // calc date
        calcDate, // pv to this date.  NOTE: we pv to calcDate as in ratesLib SwapYield_s
        //swapStartDate,
        couponCurve, 
        couponDates,
        false)); // don't compute log

    // Disc factor from swap start date to calc date
    mcDFSwapStart = SVGenExpectedDiscFactorSP(
        new SVGenExpectedDiscFactor(
        calcDate, // calc date
        calcDate, // pv to this date.  NOTE: we pv to calcDate as in ratesLib SwapYield_s
        //swapStartDate,
        couponCurve, 
        DateTimeArray(1, swapStartDate),
        false)); // don't compute log

    // precomputation...
    computeRatio(swapStartDate, couponDates, couponCurve, discountCurve);
}

void SVGenIRSwap::initialize(        
    YieldCurveConstSP    couponCurve,   // for estimating coupons etc
    YieldCurveConstSP    discountCurve,   // optional
    const DateTime&      calcDate,      // observation date of the swap yield
    const DateTime&      swapStartDate, // start date of swap
    const DateTimeArray& couponDates,   // coupon dates of swap
    DayCountConventionSP dcc)
{
    // then build our MCDiscFactor state variable generator
    mcDiscFactor = SVGenExpectedDiscFactorSP(
        new SVGenExpectedDiscFactor(
        calcDate, // calc date
        calcDate, // pv to this date.  NOTE: we pv to calcDate as in ratesLib SwapYield_s
        //swapStartDate,
        useGeneralFormula? discountCurve: couponCurve,  // Why is this inconsistent with above?
        couponDates,
        false)); // don't compute log

    // Disc factor from swap start date to calc date
    mcDFSwapStart = SVGenExpectedDiscFactorSP(
        new SVGenExpectedDiscFactor(
        calcDate, // calc date
        calcDate, // pv to this date.  NOTE: we pv to calcDate as in ratesLib SwapYield_s
        //swapStartDate,
        useGeneralFormula? discountCurve: couponCurve,  // Why is this inconsistent with above?
        DateTimeArray(1, swapStartDate),
        false)); // don't compute log

    // populate our cashflow array
    coupons = CashFlowArray(couponDates.size());
    for (int i = 0; i < couponDates.size(); i++){
        coupons[i].date = couponDates[i];
        coupons[i].amount = dcc->years(i == 0? swapStartDate: couponDates[i-1],
            couponDates[i]); // rate = 1.0
    }

    // precomputation...
    computeRatio(swapStartDate, couponDates, couponCurve, discountCurve);
}


DRLIB_END_NAMESPACE
