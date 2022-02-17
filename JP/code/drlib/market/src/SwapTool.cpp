//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SwapTool.cpp
//
//   Description : All the nasty date generation routines you need for swaps
//
//   Author      : Andrew J Swain
//
//   Date        : 7 February 2001
//
//
//----------------------------------------------------------------------------

#include <stdio.h>
#include <math.h>
#include "edginc/config.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/BadDayNone.hpp"
#include "edginc/StubFactory.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/Addin.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/StubPlacement.hpp"
#include "edginc/ZeroCurve.hpp"
#include "edginc/StubFactory.hpp"
#include "edginc/CashStream.hpp"
#include "edginc/StubSimple.hpp"
#include "edginc/NonPricingModel.hpp"

DRLIB_BEGIN_NAMESPACE


/** Calculates a par swap rate using passed zero curve for the swap starting at startDate,
 * maturing at maturityDate, with fixed day count convention 
 * fixedDayCountConv and fixed payments occuring at time 
 * intervals defined by interval. In other words, the routine calculates 
 * the fixed rate such that the present value of the fixed and 
 * floating sides are equal. The floating side is assumed to be at par.
 */
/* static */
double SwapTool::couponRate(
    const DateTime&           startDate,       // (I) Date instrument begins at
    const DateTime&           maturityDate,    // (I) Date instrument matures at
    const MaturityPeriod&     interval,        // (I) Time between payments 
    bool                      stubAtEnd,
    const DayCountConvention* dcc,
    const ZeroCurve&          zc)
{    
	static const string method = "SwapTool::couponRate";
    try {
        int    count;
        string period;

        interval.decompose(count, period);

        CashFlowArray cfl = SwapTool::cashflows(startDate,
                                                maturityDate,
                                                stubAtEnd,
                                                1.0,
                                                count,
                                                period,
                                                dcc);
        
        if (cfl.empty()) {
            throw ModelException(method, "no cashflows between " + 
                                 startDate.toString() + " and " +
                                 maturityDate.toString());
        }

        // don't want final principal repayment
        cfl.back().amount -= 1.0;

        // Get present value of 1 at startDate
        double startDatePV = zc.pv(zc.getBaseDate(), startDate);

        // Compute coupon for the given zero-coupon rates
        double couponRate = calculateCoupon(zc, cfl, startDatePV);
        return couponRate;
    }
    catch (exception&e ) {
        throw ModelException(e, method);
    }
}

/* static */
// ALIB: zr2coup.c#440
// GtoCalcCoupon
double SwapTool::calculateCoupon(const ZeroCurve& zc, const CashFlowArray& cfl, double presentValue)
{
	static const string method = "SwapTool::calculateCoupon";

    double couponsPV = zc.fv(cfl, zc.getBaseDate());

    if (!Maths::isPositive(couponsPV)) 
    {
        throw ModelException(method, "coupons with rate = 1.0 value to <= 0.0");
    }

    double lastPV = zc.pv(zc.getBaseDate(), cfl.back().date);
    return (presentValue - lastPV)/couponsPV;
}

/* static */
// ALIB: swappv.c#49
/*
 * This function is a convenience function which calls makeCFL2 followed
 * by cashFlowFV.
 */
double SwapTool::swapFixedPV(
    const ZeroCurve&          zeroCurve,
    const double              couponRate,
    const DateTime&           startDate,
    const MaturityPeriod&     couponInterval,
    const DateTime&           maturityDate,
    const DayCountConvention& dcc,
    const Stub&               stubType,
    const StubPlacement&      stubPlacement,    // stub is always short + simple
    const bool                subtractInitial,  // notional payment?
    const bool                addFinal,         // notional payment?
    const BadDayConvention&   accBadDayConv,    // accrual bad day convention
    const BadDayConvention&   payBadDayConv,    // payment bad day convention
    const Holiday&            holidays,
    const DateTime&           valueDate)
{
    static const string method ="SwapTool::swapFixedPV";

    CashFlowArray cfl = cashflows(
        couponRate,
        startDate,
        valueDate,
        maturityDate,
        stubType,
        stubPlacement.isEndStub(couponInterval, startDate, maturityDate),
        stubPlacement.isLongStub(),
        accBadDayConv,
        payBadDayConv,
        holidays,
        false,
        subtractInitial,
        addFinal,
        couponInterval,
        dcc);

    double pv = zeroCurve.pv(cfl);
    return pv;
}

/*
 * This function is a convenience function to PV swap floating leg.
 */
double SwapTool::swapFloatPV(
    const ZeroCurve&          discounting,
    const ZeroCurve&          estimating,
    const double              spread,
    const bool                isAdditive,
    const DateTime&           startDate,
    const MaturityPeriod&     interval,
    const DateTime&           maturityDate,
    const DayCountConvention& dcc,
    const Stub&               stubType,
    const StubPlacement&      stubPlacement,    // stub is always short + simple
    const bool                subtractInitial,  // notional payment?
    const bool                addFinal,         // notional payment?
    const BadDayConvention&   accBadDayConv,    // accrual bad day convention
    const BadDayConvention&   payBadDayConv,    // payment bad day convention
    const BadDayConvention&   resetBadDayConv,  // reset bad day convention
    const Holiday&            holidays,
    const bool                firstFixed,
    const double              firstFixRate,
    const DateTime&           valueDate)
{
    static const string method ="SwapTool::swapFloatPV";

        CashFlowArray cfl = cashflows(
            spread,
            estimating,
            startDate,
            valueDate,  // start date
            valueDate,  // roll date
            maturityDate,
            interval,
            dcc,
            stubType,
            stubPlacement.isEndStub(interval, startDate, maturityDate),
            stubPlacement.isLongStub(),
            subtractInitial,
            addFinal,
            accBadDayConv,
            payBadDayConv,
            resetBadDayConv,
            holidays,
            firstFixed,
            firstFixRate,
            isAdditive);

    double pv = discounting.pv(cfl);
    return pv;
}

/* static */
// ALIB: cashflow.c#809 (GtoMakeCFL2)
CashFlowArray SwapTool::cashflows(
    const double              rate,
    const DateTime&           startDate,
    const DateTime&           valueDate,
    const DateTime&           maturityDate,
    const Stub&               stubType,
    bool                      stubAtEnd,
    bool                      longStub,
    const BadDayConvention&   accrualBadDayConv,
    const BadDayConvention&   payBadDayConv,
    const Holiday&            holidays,
    bool                      keepStartDate,
    bool                      subtractInitial,
    bool                      addFinal,
    const MaturityPeriod&     period,           // as tenor, eg. 3M
    const DayCountConvention& dcc)
{
    static const string method ="SwapTool::cashflows";

    if (valueDate > startDate)
    {
        if (subtractInitial)
        {
            string msg = Format::toString(
                "Value date (%s) after start date (%s) - cannot set subtractInitial flag", 
                valueDate.toString().c_str(), 
                startDate.toString().c_str());
            throw ModelException(method, msg);
        }

        if (keepStartDate)
        {
            string msg = Format::toString( 
                "Value date (%s) after start date (%s) - cannot set keepStartDate flag", 
                valueDate.toString().c_str(), 
                startDate.toString().c_str());
            throw ModelException(method, msg);
        }
    }

    CashFlowArray cfl = cashflows(
        startDate,
        maturityDate,
        stubType,
        stubAtEnd,
        longStub,
        accrualBadDayConv,
        payBadDayConv,
        holidays,
        keepStartDate,
        subtractInitial,
        addFinal,
        rate,
        period,
        dcc);

    if (valueDate > startDate)
    {
        for (CashFlowArray::iterator iterator = cfl.begin() ; iterator!= cfl.end() ; iterator++)
        {
            if ((*iterator).date > valueDate)
            {
                cfl.erase(iterator);
            }
        }
    }

    return cfl;
}


/* static */
// ALIB: swaprate.c#120 (GtoSwapRate2)
double SwapTool::swapRate(
    const ZeroCurve&          discountCurve, 
    const DateTime&           startDate,
    const DateTime&           maturityDate,
    const MaturityPeriod&     fixedIvl,
    const DayCountConvention& fixedDayCountConv,
    const bool                valueFloating,
    const double              floatPV,
    const ZeroCurve*          estimatingCurve,
    const MaturityPeriod*     floatIvl,
    const DayCountConvention* floatDayCountConv,
    const bool                firstFixed,
    const double              firstFixRate,
    const bool                convexityDelayAdj,
    const IRVolBase*          volModelIR,
    const StubPlacement&      stubPos,
    const BadDayConvention&   accBadDayConv,    // accrual bad day convention
    const BadDayConvention&   payBadDayConv,    // payment bad day convention
    const BadDayConvention&   resetBadDayConv,  // reset bad day convention
    const Holiday&            holidays)
{
    return swapRate(
        discountCurve,
        startDate,
        startDate,  // roll date
        maturityDate,
        fixedIvl,
        fixedDayCountConv,
        valueFloating,
        floatPV,
        estimatingCurve,
        0.0,        // float spread
        floatIvl,
        floatDayCountConv,
        firstFixed,
        firstFixRate,
        convexityDelayAdj,
        volModelIR,
        NULL,       // no amortization
        stubPos,
        accBadDayConv,
        payBadDayConv,
        resetBadDayConv,
        holidays);
}


/* static */
// ALIB: swaprate.c#537 (GtoSwapRateWithRollDate)
double SwapTool::swapRate(
    const ZeroCurve&          discountCurve, 
    const DateTime&           startDate,
    const DateTime&           rollDate,
    const DateTime&           maturityDate,
    const MaturityPeriod&     fixedIvl,
    const DayCountConvention& fixedDayCountConv,
    const bool                valueFloating,
    double                    floatPV,
    const ZeroCurve*          estimatingCurve,
    const double              floatSpread,
    const MaturityPeriod*     floatIvl,
    const DayCountConvention* floatDayCountConv,
    const bool                firstFixed,
    const double              firstFixRate,
    const bool                convexityDelayAdj,
    const IRVolBase*          volModelIR,
    const CashFlowArray*      amortSched,
    const StubPlacement&      stubPos,
    const BadDayConvention&   accBadDayConv,    // accrual bad day convention
    const BadDayConvention&   payBadDayConv,    // payment bad day convention
    const BadDayConvention&   resetBadDayConv,  // reset bad day convention
    const Holiday&            holidays)
{
    static const string method ="SwapTool::swapRateWithRollDate";
    StubSimple stubType;

    double floatAdj = 0.0;
    double annuityAdj = 0.0;
    bool stubAtEnd = stubPos.isEndStub(fixedIvl, startDate, maturityDate);

    if (maturityDate.empty())
    {
        string msg = Format::toString("TBD!! support perpetual calculations [%s line %d]", __FILE__, __LINE__);
        throw ModelException(method, msg);
    }

    if (amortSched)
    {
        string msg = Format::toString("TBD!! support amortization [%s line %d]", __FILE__, __LINE__);
        throw ModelException(method, msg);
    }

    if (valueFloating && convexityDelayAdj)
    {
        string msg = Format::toString("Convexity/delay adjustments not yet supported");
        throw ModelException(method, msg);
    }

    DateTime zcFirstDate = discountCurve.firstDate();

    CashFlowArray annuityCFL = cashflows(
        startDate, 
        rollDate,
        maturityDate, 
        stubType, 
        stubAtEnd, 
        stubPos.isLongStub(), 
        accBadDayConv, 
        payBadDayConv, 
        holidays, 
        true, 
        false, 
        false, 
        1.0,    // coupon = 1 for annuity
        fixedIvl, 
        fixedDayCountConv);

    if (annuityCFL[0].date < zcFirstDate)
    {
        string msg = Format::toString("First fixed flow date (%s) < first valid zero date (%s)",
            annuityCFL[0].date.toString().c_str(), zcFirstDate.toString().c_str());
        throw ModelException(method, msg);
    }

    if (amortSched)
    {
        string msg = Format::toString("TBD!! support amortization [%s line %d]", __FILE__, __LINE__);
        throw ModelException(method, msg);
    }

    double annuityPV = discountCurve.pv(annuityCFL);

    /* Consider the case of a fixed first coupon at time t => floatPV = (1 + fix * t) * Z(t)
     *
     * What happens in the presence of amortization...?
     *
     * Floating leg value at time t:
     *   + priBalancePeriodAfterT
     *        ( resulting from payment of 1 + rt at second coupon )
     *   + priBalancePeriodBeforeT * fix * t
     *        ( fixed coupon on principal balance applying between time 0 and first coupon )
     *   + priRepaidAtTimeT
     *
     * => floatPV = (priBalancePeriodAfterT +
     *               priBalancePeriodBeforeT * fix * t +
     *               priRepaidAtTimeT) * Z(t)
     *
     * ... but: priRepaidAtTimeT = priBalancePeriodBeforeT - priBalancePeriodAfterT
     *
     * => floatPV = priBalancePeriodBeforeT * (1 + fix * t) * Z(t)
     *
     * So... the same calculation for floatPV holds providing we scale by the
     * principal balance applying between time 0 and first coupon
     *
     * This leads to our scaling rule: input floatPVs are scaled by the principal
     * balance applying between time 0 and first coupon.
     *
     * For consistency, when valuing the principal repayment schedule, we should
     * discard principal repayments before *or on* time 0.
     */
    double floatPVNoPrincipal;

    if (valueFloating)
    {
        if (!estimatingCurve || !floatIvl)
        {
            throw ModelException(method, "estimatingCurve or floating interval is NULL");
        }

        if (!floatDayCountConv)
        {
            throw ModelException(method, "Floating day count convention must be "
                "specified if valueFloating is true");
        }

        CashFlowArray floatCFL = cashflows(
            floatSpread,
            *estimatingCurve,
            startDate,
            startDate,  // value date
            rollDate,
            maturityDate,
            *floatIvl,
            *floatDayCountConv,
            stubType,
            stubAtEnd,
            stubPos.isLongStub(),
            false,
            false,
            accBadDayConv,
            payBadDayConv,
            resetBadDayConv,
            holidays,
            firstFixed,
            firstFixRate);

        if (floatCFL[0].date < zcFirstDate)
        {
            string msg = Format::toString("First floating flow date (%s) < first valid zero date (%s)",
                floatCFL[0].date.toString().c_str(), zcFirstDate.toString().c_str());
            throw ModelException(method, msg);
        }

        if (amortSched)
        {
            string msg = Format::toString("TBD!! support amortization [%s line %d]", __FILE__, __LINE__);
            throw ModelException(method, msg);
        }

        floatPVNoPrincipal = discountCurve.pv(floatCFL);
    }
    else
    {
        /*
         * User provided floatPV is the value including the principal
         * payment at maturity valued at startDate.
         *
         * floatPVNoPrincipal = floatPV@baseDate - principalPV
         */
        if (startDate < zcFirstDate)
        {
            string msg = Format::toString("Start date (%s) < first valid zero date (%s)",
                startDate.toString().c_str(), zcFirstDate.toString().c_str());
            throw ModelException(method, msg);
        }

        DateTime maturityDate = annuityCFL[annuityCFL.size() -1].date;
        double startPV = discountCurve.discountFactor(startDate);
        double maturityPV;

        if (amortSched)
        {
            string msg = Format::toString("TBD!! support amortization [%s line %d]", __FILE__, __LINE__);
            throw ModelException(method, msg);
        }
        else
        {
            maturityPV = discountCurve.discountFactor(maturityDate);
        }

        if (firstFixed)
        {
            // ensure price is 1.0
            if (!Maths::equals(1.0, floatPV))
            {
                string msg = Format::toString("Floating leg PV must equal 1.0 if "
                    "the first floating rate is fixed and valueFloating is false");
                throw ModelException(method, msg);
            }

            /*
             * We should be using valueFloating if we have different bad day
             * conventions for accrual and payment.
             */
            if (!accBadDayConv.equalTo(&payBadDayConv))
            {
                string msg = Format::toString("accBadDayConv must equal payBadDayConv "
                    "if floating rate is fixed and valueFloating is false");
                throw ModelException(method, msg);
            }

            // ensure floating dcc is specified
            if (!floatDayCountConv)
            {
                string msg = Format::toString("Floating day count convention must be specified "
                    "if the first floating rate is fixed and valueFloating is false");
                throw ModelException(method, msg);
            }

            floatPV = calculateFloatLegWithFixing(
                discountCurve,
                startDate,
                rollDate,
                maturityDate,
                *floatIvl,
                *floatDayCountConv,
                firstFixRate,
                stubAtEnd,
                accBadDayConv,
                holidays);
        }

        /*
         * For reasons discussed above, we scale floatPV by principal balance 
         * applying from startDate to first coupon date.
         */
        if (amortSched)
        {
            string msg = Format::toString("TBD!! support amortization [%s line %d]", __FILE__, __LINE__);
            throw ModelException(method, msg);
        }
        else
        {
            floatPVNoPrincipal = floatPV * startPV - maturityPV;
        }
    }

    /*
     * Adjust for unpriced fixed/floating flows (perpetual)...
     */
    floatPVNoPrincipal += floatAdj;
    annuityPV          += annuityAdj;

    /*
     * Compute coupon - this initial calculation assumes a linear accrual, and
     * will need to be adjusted in the case that we are using decompound accruals.
     */

    double couponRate = floatPVNoPrincipal / annuityPV;

    // TBD!! decompound accruals not yet supported

    return couponRate;
}


/* static */
double SwapTool::calculateFloatLegWithFixing(
    const ZeroCurve&          discountCurve,
    const DateTime&           startDate,
    const DateTime&           rollDate,
    const DateTime&           maturityDate,
    const MaturityPeriod&     floatIvl,
    const DayCountConvention& floatDcc,
    const double              firstFixRate,
    const bool                stubAtEnd,
    const BadDayConvention&   badDayConv,
    const Holiday&            holidays)
{
    DateTimeArraySP dl(SwapTool::dateArray(startDate, rollDate, maturityDate, floatIvl, stubAtEnd));
    DateTime firstCpnDate = badDayConv.adjust((*dl)[1], &holidays);
    double yearFrac = floatDcc.years(startDate, firstCpnDate);
    double fwdDiscount = discountCurve.pv(startDate, firstCpnDate);
    return (1 + firstFixRate * yearFrac) * fwdDiscount;
}


/* static */
// ALIB: cashflow.c#1227 (GtoMakeCFLFloatingRoll)
CashFlowArray SwapTool::cashflows(
    const double              spread,
    const ZeroCurve&          curve,
    const DateTime&           startDate,
    const DateTime&           valueDate,
    const DateTime&           rollDate,
    const DateTime&           maturityDate,
    const MaturityPeriod&     interval,           // as tenor, eg. 3M
    const DayCountConvention& dayCountConv,
    const Stub&               stubType,
    const bool                stubAtEnd,
    const bool                longStub,
    const bool                subtractInitial,
    const bool                addFinal,
    const BadDayConvention&   accrualBadDayConv,
    const BadDayConvention&   payBadDayConv,
    const BadDayConvention&   resetBadDayConv,
    const Holiday&            holidays,
    const bool                firstFixed,
    const double              firstFixRate,
    const bool                isAdditive)    // is spread additive?
{
    CashStream cs;
    cs.addFloatLegVanilla(
        1.0,            // notional
        spread,
        isAdditive,
        startDate,
        valueDate,
        rollDate,
        interval,
        maturityDate,
        stubType,
        stubAtEnd,
        subtractInitial,
        addFinal,
        accrualBadDayConv,
        payBadDayConv,
        resetBadDayConv,
        firstFixed,
        firstFixRate,
        dayCountConv,
        holidays);

    // Convert the cash stream into a cashflow list.
    CashFlowArraySP cfl = cs.flows(curve, NULL);
    return *cfl;
}




/* static */
// equivalent to ALIB GtoNewDateListExtended
// ALIB: datelist.c#528 (GtoNewDateListExtendedRoll)
DateTimeArray* SwapTool::dateArray(
    const DateTime&       startDate,
    const DateTime&       rollDate,
    const DateTime&       maturityDate,
    const MaturityPeriod& interval,        // as tenor, eg. 3M
    bool                  stubAtEnd,
    bool                  constrainEndDates)
{
    int count;
    string period;
    interval.decompose(count, period);

    if (startDate != rollDate)
    {
        string msg = Format::toString(
            "TBD!! handle rollDate != startDate [%s line %d]", 
            __FILE__, __LINE__);
        throw ModelException("SwapTool::dateArray", msg);
    }

    return dateArray(startDate, maturityDate, count, period, stubAtEnd, constrainEndDates);
}


/* static */
// ALIB: datelist.c#528 (GtoNewDateListExtendedRoll)
DateTimeArray* SwapTool::dateArray(
    const DateTime& startDate,
    const DateTime& rollDate,
    const DateTime& maturityDate,
    int             count,           // interval = count periods
    const string&   period,          // e.g. Y, M, W, D
    bool            stubAtEnd)
{
    static const string method = "SwapTool::dateArray";

    if (rollDate.empty() || stubAtEnd)
    {
        return dateArray(startDate, maturityDate, count, period, stubAtEnd);
    }

    int numIntervals;
    int extraDays;
    countDates(startDate, maturityDate, count, period, &numIntervals, &extraDays);

    if (extraDays == 0)
    {
        // We can use the roll date - so we count forward from this day
        if (startDate < rollDate)
        {
            string msg = Format::toString("Start date cannot be before roll date");
            throw ModelException(method, msg);
        }

        /*
         * We roll forward until we are on the startDate or if there are a
         * non-integral number of periods so that we are on the flow date
         * preceding the start date.
         */
        countDates(rollDate, startDate, count, period, &numIntervals, &extraDays);

        DateTime firstDate;
        if (extraDays == 0)
        {
            firstDate = startDate;
        }
        else
        {
            firstDate = MaturityPeriod::toDate(count * numIntervals, period, rollDate);
        }

        return dateArray(firstDate, maturityDate, count, period, true);
    }
    else
    {
        return dateArray(startDate, maturityDate, count, period, stubAtEnd);
    }
}


/* static */
CashFlowArray SwapTool::cashflows(
    const DateTime&           startDate,
    const DateTime&           maturityDate,
    const Stub&               stubType,
    bool                      stubAtEnd,
    bool                      longStub,
    const BadDayConvention&   accrualBadDayConv,
    const BadDayConvention&   payBadDayConv,
    const Holiday&            holidays,
    bool                      keepStartDate,
    bool                      subtractInitial,
    bool                      addFinal,
    double                    rate,
    const MaturityPeriod&     period,           // as tenor, eg. 3M
    const DayCountConvention& dcc)
{
    return cashflows(
        startDate, 
        startDate,  // roll date
        maturityDate,
        stubType, 
        stubAtEnd, 
        longStub, 
        accrualBadDayConv, 
        payBadDayConv, 
        holidays, 
        keepStartDate,
        subtractInitial, 
        addFinal, 
        rate, 
        period, 
        dcc);
}

/* static */
// ALIB: cashflow.c#949 (GtoMakeCFLRoll)
CashFlowArray SwapTool::cashflows(
    const DateTime&           startDate,
    const DateTime&           rollDate,
    const DateTime&           maturityDate,
    const Stub&               stubType,
    bool                      stubAtEnd,
    bool                      longStub,
    const BadDayConvention&   accrualBadDayConv,
    const BadDayConvention&   payBadDayConv,
    const Holiday&            holidays,
    bool                      keepStartDate,
    bool                      subtractInitial,
    bool                      addFinal,
    double                    rate,
    const MaturityPeriod&     period,           // as tenor, eg. 3M
    const DayCountConvention& dcc)
{
    static const string method = "SwapTool::cashflows";

    int count;
    string interval;
    period.decompose(count, interval);

    // TBD!! change SwapTool to not return array by value!!
    return cashflows(
        startDate, 
        rollDate,
        maturityDate, 
        &stubType, 
        stubAtEnd, 
        longStub, 
        &accrualBadDayConv, 
        &payBadDayConv, 
        &holidays, 
        keepStartDate, 
        subtractInitial, 
        addFinal, 
        rate, 
        count, 
        interval, 
        &dcc);
}


CashFlowArray SwapTool::cashflows(
    const DateTime&           startDate,
    const DateTime&           rollDate,
    const DateTime&           maturityDate,
    const Stub*               stub,
    bool                      stubAtEnd,
    bool                      longStub,
    const BadDayConvention*   accrualBadDayConv,
    const BadDayConvention*   payBadDayConv,
    const Holiday*            holidays,
    bool                      keepStartDate,
    bool                      subtractInitial,
    bool                      addFinal,
    double                    rate,
    int                       count,           // interval = count periods
    const string&             period,          // e.g. Y, M, W, D
    const DayCountConvention *dcc)
{
    static const string method = "SwapTool::cashflows";

    if (rollDate != startDate)
    {
        string msg = Format::toString("TBD!! support roll date [%s line %d]", __FILE__, __LINE__);
        throw ModelException(method, msg);
    }

    try {
        int             i;                  
 
        // Get dateList with adjusted coupon dates & mat date 
        DateTimeArraySP dl=DateTimeArraySP(SwapTool::dateArray(startDate,
                                                               maturityDate,
                                                               count,
                                                               period,
                                                               stubAtEnd));
        if (dl->empty()) {
            string msg = Format::toString( 
                                 "generated empty payment date list between %s and %s",
                                 startDate.toString().c_str(),
                                 maturityDate.toString().c_str());
            throw ModelException(method, msg);
        }

        bool frontStub = startDate > (*dl)[0];
        bool backStub = maturityDate < (*dl)[dl->getLength() - 1];
        if (frontStub && backStub)
        {
            throw ModelException(method, "Cannot have back & front stubs");
        }

        /* 
         * For long stubs, simply coalesce first/last two coupon periods
         * (i.e. short stub + adjacent regular period) into one...
         */
        if (longStub && (frontStub || backStub) && dl->getLength() > 2)
        {
            if (frontStub)
            {
                DateTimeArray::iterator iterator = dl->begin();
                iterator++;
                dl->erase(iterator);
            }
            else // back stub
            {
                (*dl)[dl->getLength() - 2] = (*dl)[dl->getLength() - 1];
                dl->pop_back();
            }
        }

        DateTimeArraySP accrualDates;
        DateTimeArraySP payDates;
        DateTime        accrualStart;
        DateTime        accrualMaturity;
        DateTime        payStart;
        DateTime        payMaturity;

        // get accrual dates
        if (BadDayNone::TYPE->isInstance(accrualBadDayConv)) {
            accrualDates    = dl;
            accrualStart    = startDate;
            accrualMaturity = maturityDate;
        }
        else {
            accrualDates = DateTimeArraySP(new DateTimeArray(dl->getLength()));
            for (i = 0; i < dl->getLength(); i++) {
                (*accrualDates)[i] = accrualBadDayConv->adjust((*dl)[i], holidays);
            }
            accrualStart    = accrualBadDayConv->adjust(startDate, holidays);
            accrualMaturity = accrualBadDayConv->adjust(maturityDate, holidays);
        }

        // get payment dates
        if (BadDayNone::TYPE->isInstance(payBadDayConv)) {
            payDates    = dl;
            payStart    = startDate;
            payMaturity = maturityDate;
        }
        else {
            payDates = DateTimeArraySP(new DateTimeArray(dl->getLength()));
            for (i = 0; i < dl->getLength(); i++) {
                (*payDates)[i] = payBadDayConv->adjust((*dl)[i], holidays);
            }
            payStart    = payBadDayConv->adjust(startDate, holidays);
            payMaturity = payBadDayConv->adjust(maturityDate, holidays);
        }

        CashFlowArray cfl(dl->getLength());

        cfl[0].date   = payStart;
        cfl[0].amount = 0.0;

        cfl[cfl.size()-1].date = payMaturity;

        for (i = 1; i < dl->getLength()-1; i++) {
            cfl[i].date = (*payDates)[i];
        }

        int startIdx;
        int endIdx;

        // if front stub
        if (frontStub) {
            startIdx = 2;                   // For non-stub processing 
            endIdx   = dl->getLength()-1;   // For non-stub processing 
            
            cfl[1].amount = stub->payment((*accrualDates)[0],  // Coupon start
                                          (*accrualDates)[1],  // Coupon end
                                          accrualStart,        // Accrue start
                                          (*accrualDates)[1],  // Accrue end
                                          rate,
                                          dcc);                          
        }
        else if (backStub) {
            // back stub
            startIdx = 1;                   // For non-stub processing 
            endIdx   = dl->getLength()-2;   // For non-stub processing 

            cfl[cfl.size()-1].amount = stub->payment((*accrualDates)[dl->getLength()-2],  // Coupon start
                                                     (*accrualDates)[dl->getLength()-1],  // Coupon end
                                                     (*accrualDates)[dl->getLength()-2],  // Accrue start
                                                     accrualMaturity,  // Accrue end
                                                     rate,
                                                     dcc);                               
        }
        else {
            // no stub
            startIdx = 1;
            endIdx   = dl->getLength()-1;
        }
            
        // Compute non-stub cashflows
        for (i = startIdx; i <= endIdx; i++) {
            cfl[i].amount = rate * dcc->years((*accrualDates)[i-1],
                                              (*accrualDates)[i]);
        }

        if (subtractInitial) {
            cfl[0].amount -= 1.0;
        }

        if (addFinal) {
            cfl[cfl.getLength()-1].amount += 1.0;
        }

        if (!keepStartDate && Maths::isZero(cfl[0].amount)) {
            for (i = 1; i < cfl.size(); i++) {
                cfl[i-1].date   = cfl[i].date;
                cfl[i-1].amount = cfl[i].amount;
            }
            cfl.resize(cfl.size()-1);
        }

        // sort out any multiple dates (usually if building a daily list)
        CashFlow::aggregate(cfl);
        return cfl;
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}


/** Given a frequency, turn into the number of months that separate each 
    date. Essentially 12/frequency.  */
static int freqToMonthPeriod(int freq)
{
    int prd;

    if (freq == CompoundBasis::ANNUAL){
        prd = 12;
    }
    else if (freq == CompoundBasis::SEMI_ANNUAL){
        prd = 6;
    }
    else if (freq == CompoundBasis::QUARTERLY){
        prd = 3;
    }
    else if (freq == CompoundBasis::MONTHLY){
        prd = 1;
    }
    else {
        throw ModelException("freqToMonthPeriod", "SwapTools only supports "
                             "annual, semiannual, quarterly, monthly "
                             "frequencies at present");
    }
    
    return prd;
}

// is a date on cycle ?
bool SwapTool::onCycle(
    const DateTime& valueDate,
    const DateTime& date,            // is this date on cycle ?
    int             count,           // interval = count periods
    const string&   period)          // e.g. Y, M, W, D
{
    static const string method = "SwapTool::onCycle";
    try {
        DateTime::MonthDayYear valueMDY = valueDate.toMDY();
        DateTime::MonthDayYear dateMDY  = date.toMDY();

        bool isOnCycle;

        if (valueMDY.day <= 28 && dateMDY.day <= 28)
        {
            int numIntervals;
            int extraDays;
            // We assume we can only be on cycle if date is not on or after 
            // the 29th of the month.

            // Find out if date is on cycle.
            SwapTool::countDates(valueDate,
                                 date,
                                 count, 
                                 period,
                                 &numIntervals,
                                 &extraDays);

            
            isOnCycle = (extraDays == 0);                    
        }
        else
        {
            isOnCycle = false;
        }
        return isOnCycle;
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}

// Checks if swap maturity dates are on the regular cycle
// from valueDate. Works only if floating and fixed freqs are equal.
bool SwapTool::swapDatesOnCycle(
    const DateTime&      valueDate,     
    const DateTimeArray& swapDates,      
    int                  freq) {
    static const string method = "SwapTool::swapDatesOnCycle";
    try {
        bool     onCycle = true;
        DateTime cycleDate;               
        int      swpIdx;                  
        int      coupIdx;   
             
        DateTime::MonthDayYear mdy = valueDate.toMDY();

        if (mdy.day >= 29) {
            // Exclude value dates which are on the 29th or higher,
            // since it's hard to check if we're actually "on cycle".
            onCycle = false;
        }
        else {
            cycleDate = valueDate;

            for (swpIdx=0,coupIdx=1; 
                 onCycle && swpIdx < swapDates.size(); 
                 swpIdx++)
            {
                while (cycleDate.isLess(swapDates[swpIdx]))
                {
                    int period = freqToMonthPeriod(freq);          
                    cycleDate = MaturityPeriod::toDate(period*coupIdx++, "M", valueDate);
                }
                if (!cycleDate.equals(swapDates[swpIdx])) {
                    onCycle = false;
                }
                else {
                    mdy = swapDates[swpIdx].toMDY();
                    if (mdy.day >= 29) {
                        // Exclude swap dates which are on the 29th or higher,
                        // since it's hard to check if we're actually "on cycle".
                        onCycle = false;
                    }
                }
            }
        }
        return onCycle;
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}

/*
 * Counts # intervals in a range of dates
 *
 * Caution:
 * In general countDates(dateA, dateB, count, period) produces the same results
 * as countDates(dateB, dateA, -count, period). However if dateA or dateB are 
 * the end of a month, there are differences, e.g.:
 *	startDate	endDate	result
 *	-------------------------------------
 *	29-Feb-04	28-Feb-15	11y + 0d
 *	28-Feb-15	29-Feb-04	10y + 365d
 *	
 *	28-Feb-05	29-Feb-08	3y + 1d
 *	29-Feb-08	28-Feb-05	3y + 0d
 *
 *	30-Aug-05	29-Feb-12	6y + 183d
 *	29-Feb-12	30-Aug-05	6y + 182d
 *
 *	31-Dec-2003	29-Jun-2006	4x6M + 180d
 *	29-Jun-2006	31-Dec-2003	4x6M + 181d
 *
 * Another peculiarity is that fractions of days are not included in numIntervals, but
 * ARE included in extraDays: for example, if startDate = 31-Mar-2004 (END_OF_DAY),
 * endDate = 31-Mar-2005 (START_OF_DAY), period = D and count = 1, the results are:
 * numIntervals = 364, extraDays = 1.
 *
 * This has NOT been confirmed to be correct, but is the way the code traditionally 
 * worked */
void SwapTool::countDates(
    const DateTime& fromDate,
    const DateTime& toDate,
    int             count,           // interval = count periods
    const string&   period,          // e.g. Y, M, W, D
    int*            numIntervals,     // (O) Answer (Quotient) 
    int*            extraDays)        // (O) Days left over(remainder)
{
    static const string method = "SwapTool::countDates";
    try {
        // Make sure interval has the right sign.
        if (count == 0) {
            throw ModelException(method, "zero time interval not allowed");
        }

        if (count * (toDate.getDate() - fromDate.getDate()) < 0) {
            throw ModelException(method, "cannot count from " + 
                                 fromDate.toString() + " to " +
                                 toDate.toString() + " with interval " + 
                                 Format::toString(count) + period);
        }

        switch (period[0])
        {
        case 'W':  // Weeks
            count *= 7;
        case 'D':  // Days             
            {
                int daysBetweenDates = abs(toDate.getDate() - fromDate.getDate());
                int addExtraDay = 0;

                // Check the times. If necessary adjust the number of days so that
                // we no longer need to worry about time
                if ((toDate.getTime() - fromDate.getTime()) * count < 0) {
                    // The last day is not complete, so remove it entirely
                    daysBetweenDates--;
                    // But need to add it as "extraDays", so keep it in mind
                    addExtraDay = 1;
                }

                // from now onwards only need the abs value of count 
                if (count < 0) {
                    count *= -1;
                }
                *numIntervals = daysBetweenDates / count;
                *extraDays = daysBetweenDates % count + addExtraDay;
            }
            break;
            
        default:
            // (Ugly but required for performance reasons)
            switch (period[0])
            {
            case 'A':  // Annual
            case 'Y':  // Years
                count *= 12;
                break;
            case 'S':  // Semmiannual
                count *= 6;
                break;
            case 'Q':  // Quaterly
                count *= 3;
                break;
            case 'M':  // Monthly
                break;
            default:
                throw ModelException(method, "SwapTools only supports "
                                     "Annual (Yearly), Semiannual, Quarterly, Monthly, "
                                     "Weekly and Daily frequencies at present. "
                                     "Given period of type " + period);
            }
            DateTime::MonthDayYear mdyFrom = fromDate.toMDY();

            // Check the times. If necessary adjust the dates so that
            // we no longer need to worry about time
            DateTime rolledDate;
            if ((toDate.getTime() - fromDate.getTime()) * count < 0) {
                // the last day is not complete, so remove it entirely
                rolledDate = toDate.rollDate((count > 0) ? -1 : 1);
            }
            else {
                rolledDate = toDate;
            }
            DateTime::MonthDayYear mdyTo = rolledDate.toMDY();

            int monthsBetweenDates = mdyFrom.wholeMonthsBetween(mdyTo);           
            *numIntervals = monthsBetweenDates / count;  // Integer division

            // To calculate the extra days, roll the date forwards as many 
            // months as numIntervals * count, and then calculate the 
            // days left to toDate.
            DateTime lastDate = MaturityPeriod::toDate(*numIntervals * count,
                                                       "M",
                                                       fromDate);
            *extraDays = abs(lastDate.getDate() - toDate.getDate());
        }
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}

// Returns a date list of standard euro-money swap dates.  This
// basically means dates of all coupons for a swap an integral number of
// years from now.  This differs from the coupons for the last swap, in that
// it may not be an integral number of years from the value date.
// The routine outputs a date list including all coupons for a swap 
// that's not after the given  maturity date and is after the second to
// last coupon date of the input swap.  And is an integral number of
// coupon frequencies from the value date.
DateTimeArray SwapTool::canonicalSwapDates(
    const DateTime& startDate,
    const DateTime& minDate,   // none before here
    const DateTime& endDate,
    int             frequency)
{
   static const string method = "SwapTool::canonicalSwapDates";
   try {
       int           count;           
       DateTime      swapDate; 
       DateTimeArray dates(0);
   
       int prd = freqToMonthPeriod(frequency);

       count = 1;
       do {
           swapDate = MaturityPeriod::toDate(prd * count, "M", startDate);

           count++;

           if (swapDate.isGreater(minDate)) {
               if (endDate.isGreater(swapDate)) {
                   dates.push_back(swapDate);
               }
               else {
                   dates.push_back(endDate);           
               }          
           }
       }
       while(endDate.isGreater(swapDate));

       return dates;
   }
   catch (exception &e) {
       throw ModelException(e, method);
   }
}

//// builds date list from start to end using given interval (long stub @ end)
//// Note that startDate is not included in the list of dates
void SwapTool::simpleSwapDates(
    const DateTime& startDate,
    const DateTime& endDate,
    int             frequency,
    DateTimeArray&  dates) // (O) 
{
   static const string method = "SwapTool::simpleSwapDates";
   try {
       int daysDiff = endDate.daysDiff(startDate);
       if (daysDiff <= 0){
           throw ModelException(method, "Zero or negative length interval");
       }
       dates.clear();
       // for better performance, work in MDY space
       DateTime::MonthDayYear mdy(startDate.toMDY());
       int period = freqToMonthPeriod(frequency);
       // reserve some memory, any decent guess will do
       dates.reserve(daysDiff*frequency/DateTime::DAYS_PER_YEAR+1);
       int time = startDate.getTime(); // preserve time of day
       bool notFinished; // lastDate goes out of scope before while() statement
       do {
           mdy.month += period;
           dates.push_back(mdy.toDateTime(time));
           const DateTime& lastDate = dates.back();
           if (endDate.equals(lastDate, false)){
               return; // we've finished
           }
           notFinished = endDate.isGreater(lastDate);
       } while (notFinished);
       dates.back() = endDate;
       // we are using a long stub, so remove the appropriate date 
       // Require 2 dates in order to perform this operation 
       if (dates.size() > 1) {
           dates.pop_back();              
           dates.back() = endDate;
       }
   } catch (exception& e) {
       throw ModelException(e, method);
   }
}

// Makes a cash flow list for a swap instrument.
CashFlowArray SwapTool::cashflows(
    const DateTime& valueDate,
    const DateTime& maturity,
    bool            stubAtEnd,
    double          rate,
    int             count,           // interval = count periods
    const string&   period,          // e.g. Y, M, W, D
    const DayCountConvention *dcc)
{
    static const string method = "SwapTool::cashflows";
    try {
        DateTime        prevDate; // prev date added to cash-flow list 
        int             i;                  

        // Bizarre rate==0.0 case 
        if (Maths::isZero(rate)) {
            CashFlowArray cfl(1);
            cfl[0].amount = 1.0;
            cfl[0].date   = maturity;
      
            return (cfl);
        }
    
        // Get dateList with adjusted coupon dates & mat date 
        DateTimeArraySP dl=DateTimeArraySP(SwapTool::paymentDates(valueDate,
                                                                  maturity,
                                                                  count,
                                                                  period,
                                                                  stubAtEnd));
        if (dl->empty()) {
            throw ModelException(method, 
                                 "generated empty payment date list between " +
                                 valueDate.toString() + " and " + 
                                 maturity.toString());
        }

        CashFlowArray cfl(dl->size());

        prevDate = valueDate;
        for (i = 0; i < dl->size(); i++)
        {
            DateTime cDate = (*dl)[i];       /* coupon date */

            cfl[i].amount = rate * dcc->years(prevDate, cDate);
            cfl[i].date = cDate;      /* store date */
            prevDate = cDate;
        }

        cfl.back().amount += 1.0;   /* add principal */

        return cfl;
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}

// Makes a cash flow list for a swap instrument.
CashFlowArray SwapTool::cashflows(
    const DateTime&           startDate,
    const DateTime&           maturity,
    const Stub*               stub,
    bool                      stubAtEnd,
    const BadDayConvention*   accrualBadDayConv,
    const BadDayConvention*   payBadDayConv,
    const Holiday*            holidays,
    bool                      keepStartDate,
    bool                      addPrincipal,
    double                    rate,
    int                       count,           // interval = count periods
    const string&             period,          // e.g. Y, M, W, D
    const DayCountConvention *dcc)
{
    return cashflows(
        startDate, 
        startDate,  // roll date
        maturity, 
        stub, 
        stubAtEnd, 
        false,
        accrualBadDayConv, 
        payBadDayConv, 
        holidays, 
        keepStartDate, 
        false, 
        addPrincipal, 
        rate, 
        count, 
        period, 
        dcc);
}

// build up accrue and pay dates for a swap
void SwapTool::swapDates(
    const DateTime&           startDate,
    const DateTime&           maturity,
    bool                      stubAtEnd,
    const BadDayConvention*   accrualBadDayConv,
    const BadDayConvention*   payBadDayConv,
    const Holiday*            holidays,
    int                       count,           // interval = count periods
    const string&             period,          // e.g. Y, M, W, D
    DateTimeArray&            accrualDates,    // (O) 
    DateTimeArray&            payDates)        // (O)
{
    static const string method = "SwapTool::swapDates";
    try {
        int             i;                  
 
        // start fresh
        accrualDates.resize(0);
        payDates.resize(0);

        // Get dateList with adjusted coupon dates & mat date 
        DateTimeArraySP dl=DateTimeArraySP(SwapTool::dateArray(startDate,
                                                               maturity,
                                                               count,
                                                               period,
                                                               stubAtEnd));
        if (dl->empty()) {
            throw ModelException(method, 
                                 "generated empty payment date list between " +
                                 startDate.toString() + " and " + 
                                 maturity.toString());
        }

        // get accrual dates
        accrualDates.resize(dl->size());
        for (i = 0; i < dl->size(); i++) {
            accrualDates[i] = accrualBadDayConv->adjust((*dl)[i], holidays);
        }

        // get payment dates
        payDates.resize(dl->size()-1);
        for (i = 0; i < dl->size()-1; i++) {
            payDates[i] = payBadDayConv->adjust((*dl)[i+1], holidays);
        }
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}


// create an array of swap payment dates
DateTimeArray* SwapTool::paymentDates(
    const DateTime& baseDate,        // start here
    const DateTime& maturity,        // end here
    int             count,           // interval = count periods
    const string&   period,          // e.g. Y, M, W, D
    bool            stubAtEnd)
{
    static const string method = "SwapTool::paymentDates";
    try {
        DateTimeArray* dates = SwapTool::dateArray(baseDate,
                                                   maturity,
                                                   count,
                                                   period,
                                                   stubAtEnd);

        // strip out the start date
        for (int i = 0; i < dates->size() -1; i++) {
            (*dates)[i] = (*dates)[i+1];
        }
        dates->pop_back();
        return dates;
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}

// create an array of dates separated by interval (given by count periods)
DateTimeArray* SwapTool::dateArray(
    const DateTime& baseDate,        // start here
    const DateTime& maturity,        // end here
    int             count,           // interval = count periods
    const string&   period,          // e.g. Y, M, W, D
    bool            stubAtEnd,
    bool            constrainEndDates)
{
    static const string method = "SwapTool::dateArray";
    try {
        int            numIntervals;
        int            extraDays;
        int            numDates;
        const DateTime *startDate;
        const DateTime *endDate;
        int            arrayIncrement;

        if (stubAtEnd) {
            startDate = &baseDate;
            endDate = &maturity;
            arrayIncrement = 1;
        }
        else {
            count = -count;
            startDate = &maturity;
            endDate = &baseDate;
            arrayIncrement = -1;
        }

        SwapTool::countDates(*startDate,
                             *endDate,
                             count,
                             period,
                             &numIntervals,
                             &extraDays);

        if (extraDays > 0) {
            numDates = numIntervals + 2; 
        }
        else {
            numDates = numIntervals + 1;
        }

        DateTimeArray* dates = SwapTool::dateArray(*startDate,
                                                   count,
                                                   period,
                                                   0,
                                                   arrayIncrement,
                                                   numDates - 1);
        dates->push_back(DateTime());   // allow space for extra element

        if (stubAtEnd)
        {
            (*dates)[numDates-1] = constrainEndDates ? maturity : MaturityPeriod::toDate(1, period, (*dates)[numDates - 2]);
        }
        else                                /* Stub at beginning */
        {
            for (int i = dates->size() - 1; i > 0 ; i--)
            {
                (*dates)[i] = (*dates)[i-1];
            }

            (*dates)[0] = constrainEndDates ? baseDate : MaturityPeriod::toDate(-1, period, baseDate);
        }

        return dates;
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}


// create an array of dates separated by interval (given by count periods)
DateTimeArray* SwapTool::dateArray(
    const DateTime& baseDate,        // start here
    const DateTime& endDate,         // end here
    int             count,           // interval = count periods
    const string&   period,          // e.g. Y, M, W, D
    int             startIdx,        // 0=start @ basedate, 1=start @ baseDate + interval
    bool            stubAtEnd,       // true=add dates forwards, false=add backwards
    int             numDates)        // how many dates
{
    static const string method = "SwapTool::dateArray";
    try {
        DateTimeArraySP dates(new DateTimeArray(numDates));
        int index;
        for (int i = 0; i < numDates-1; i++) {
            index = (stubAtEnd) ? i : (numDates-1-i);
            (*dates)[index] = MaturityPeriod::toDate(count * (startIdx + i), 
                                                     period, 
                                                     baseDate);
        }

        // Add the last element manually (corresponds to i = numDates-1)
        index = (stubAtEnd) ? numDates-1 : 0;
        (*dates)[index] = endDate;

        return dates.release();
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}

// create an array of dates separated by interval (given by count periods)
DateTimeArray* SwapTool::dateArray(
    const DateTime& baseDate,        // start here
    int             count,           // interval = count periods
    const string&   period,          // e.g. Y, M, W, D
    int             startIdx,        // 0=start @ basedate, 1=start @ baseDate + interval
    int             arrayIncrement,  // Usually +1 or -1
    int             numDates)        // how many dates
{
    static const string method = "SwapTool::dateArray";
    try {
        DateTimeArray* dates = new DateTimeArray(numDates);

        for (int i = 0; i < numDates; i++) {
            DateTime date = MaturityPeriod::toDate(count * (startIdx + i), period, baseDate);
            if (arrayIncrement > 0) {
                (*dates)[i] = date;
            }
            else {
                (*dates)[numDates-1 -i] = date;
            }
        }
        return dates;
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}


class CashFlowGenerator: public CObject, virtual public ClientRunnable{
    static CClassConstSP const TYPE;

    // input parameters
    DateTime             startDate;
    DateTime             maturity;
    string               stubRule;
    bool                 stubAtEnd;
    string               accrueBadDayConv;
    string               payBadDayConv;
    HolidaySP            holidays;
    bool                 keepStartDate;
    bool                 addPrincipal;
    double               rate;
    string               paymentFrequency;
    DayCountConventionSP dayCountConv;
    double               notional;
    MarketDataConstSP    marketData;

    /** for reflection */
    CashFlowGenerator():  CObject(TYPE){}

    static IObjectSP generateCashFlows(CashFlowGenerator* params) {
        static const string routine = "CashFlowGenerator::generateCashFlows";
        try {
            int     count;
            string  interval;
	    
            StubSP stub = StubSP(StubFactory::make(params->stubRule));
            BadDayConventionSP   accrualBDC(BadDayConventionFactory::make(params->accrueBadDayConv));
            BadDayConventionSP   payBDC(BadDayConventionFactory::make(params->payBadDayConv));
	    if (params->marketData.get())
            {
                // May need to fetch holidays for Business/252 DCC
                params->marketData->fetchForNonMarketObjects(params->dayCountConv, 
                                                             IModelConstSP(new NonPricingModel()),
                                                             "Business252");
            } 
            MaturityPeriod       period(params->paymentFrequency);
            period.decompose(count, interval);

            CashFlowArray cashFlows = SwapTool::cashflows(
                params->startDate,
                params->maturity,
                stub.get(),
                params->stubAtEnd,
                accrualBDC.get(),
                payBDC.get(),
                params->holidays.get(),
                params->keepStartDate,
                params->addPrincipal,
                params->rate,
                count,
                interval,
		params->dayCountConv.get());

            // scale cash flows by notional
            for (int i=0;i<cashFlows.size(); ++i) {
                cashFlows[i].amount *= params->notional;
            }

            CashFlowArraySP cashFlowsSP = CashFlowArraySP::attachToRef(&cashFlows);

            return CashFlowArraySP(cashFlowsSP.clone());
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    // EdrAction version of addin
    IObjectSP run() {
        return generateCashFlows(this);
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        clazz->setDescription("create cash flows from bond params");
        REGISTER(CashFlowGenerator, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultCashFlowGenerator);
        FIELD(startDate,          "start date");
        FIELD(maturity,           "maturity date");
        FIELD(stubRule,           "stub rule");
        FIELD(stubAtEnd,          "stub at end")
        FIELD(accrueBadDayConv,   "accrual bad day convention");
        FIELD(payBadDayConv,      "payment bad day convention");
        FIELD(holidays,                  "holidays");
        FIELD(keepStartDate,      "keep start date");
        FIELD(notional,           "notional");
        FIELD(addPrincipal,       "add principal");
        FIELD(rate,               "rate");
        FIELD(paymentFrequency,   "payment frequency");
        FIELD(dayCountConv,       "day count convention");
        FIELD(marketData, "market data");
        FIELD_MAKE_OPTIONAL(marketData);

        Addin::registerClassObjectMethod("GENERATE_CASH_FLOWS",
                                         Addin::UTILITIES,
                                         "Creates cash flows from a set of input parameters",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)generateCashFlows);

    }

    static IObject* defaultCashFlowGenerator(){
        return new CashFlowGenerator();
    }
    
};
   
CClassConstSP const CashFlowGenerator::TYPE = CClass::registerClassLoadMethod(
    "CashFlowGenerator", typeid(CashFlowGenerator), load);


class SimpleCashFlowGeneratorAddin: public CObject{
    static CClassConstSP const TYPE;

    // input parameters
    DateTime             startDate;
    DateTime             maturity;
    bool                 stubAtEnd;
    double               rate;
    MaturityPeriodSP     paymentFrequency;
    DayCountConventionSP dayCountConvention ;
    double               notional;
    bool                 includeNotional;
    MarketDataConstSP    marketData;

    /** for reflection */
    SimpleCashFlowGeneratorAddin():  CObject(TYPE){}

    static IObjectSP generateCashFlows(SimpleCashFlowGeneratorAddin* params) {
        static const string routine = "SimpleCashFlowGeneratorAddin::generateCashFlows";
        try {
            int     count;
            string  interval;

	    if (params->marketData.get())
            {
                // May need to fetch holidays for Business/252 DCC
                params->marketData->fetchForNonMarketObjects(params->dayCountConvention, 
                                                             IModelConstSP(new NonPricingModel()),
                                                             "Business252");
            }
            params->paymentFrequency->decompose(count, interval);

            CashFlowArray cashFlows = SwapTool::cashflows(
                params->startDate,
                params->maturity,
                params->stubAtEnd,
                params->rate,
                count,
                interval,
		params->dayCountConvention.get());

            // scale cash flows by notional
            for (int i=0;i<cashFlows.size()-1; ++i) {
                cashFlows[i].amount *= params->notional;
            }

            // the final payment includes the notional which possibly needs to be stripped out
            if ( !params->includeNotional && cashFlows.size() > 0 ) {
                cashFlows[cashFlows.size()-1].amount = (cashFlows[cashFlows.size()-1].amount - 1.0) * params->notional;
            } else {
                cashFlows[cashFlows.size()-1].amount *= params->notional;
            }

            CashFlowArraySP cashFlowsSP = CashFlowArraySP::attachToRef(&cashFlows);

            return CashFlowArraySP(cashFlowsSP.clone());
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }


    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(SimpleCashFlowGeneratorAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCashFlowGeneratorAddin);
        FIELD(startDate,          "start date");
        FIELD(maturity,           "maturity date");
        FIELD(stubAtEnd,          "stub at end")
        FIELD(notional,           "notional");
        FIELD(includeNotional,    "include notional payment");
        FIELD(rate,               "rate");
        FIELD(paymentFrequency,          "payment frequency");
        FIELD(dayCountConvention, "Day count convention");
        FIELD(marketData, "market data");
        FIELD_MAKE_OPTIONAL(marketData);

        Addin::registerClassObjectMethod("GENERATE_SIMPLE_CASH_FLOWS",
                                         Addin::UTILITIES,
                                         "Creates a cash flow array from a set of input parameters",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)generateCashFlows);

    }

    static IObject* defaultCashFlowGeneratorAddin(){
        return new SimpleCashFlowGeneratorAddin();
    }
    
};
   
CClassConstSP const SimpleCashFlowGeneratorAddin::TYPE = CClass::registerClassLoadMethod(
    "SimpleCashFlowGeneratorAddin", typeid(SimpleCashFlowGeneratorAddin), load);

DRLIB_END_NAMESPACE


