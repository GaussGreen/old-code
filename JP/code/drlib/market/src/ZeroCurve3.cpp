//----------------------------------------------------------------------------
//
//   Group       : QR Core Analytics Team
//
//   Filename    : ZeroCurve3.cpp
//
//   Description : Curve factory implementation for ALIB zero curve 3 method
//
//   Author      : Richard Appleton
//
//   Date        : 25th April 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ZeroCurve3.hpp"
#include "edginc/StubPlacement.hpp"
#include "edginc/SimpleZeroCurve.hpp"
#include "edginc/ZC3ZeroCurve.hpp"
#include "edginc/ZC3CurveInstrument.hpp"
#include "edginc/ZC3Interpolation.hpp"
#include "edginc/ZC3FuturesGapRule.hpp"
#include "edginc/ZC3FuturesStubRule.hpp"
#include "edginc/Atomic.hpp"

#include "edginc/Addin.hpp"
#include "edginc/Actual360.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/BadDayFollowing.hpp"
#include "edginc/BadDayNone.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/Hashtable.hpp" // for hash_string
#include <algorithm>


DRLIB_BEGIN_NAMESPACE


ZeroCurve3::ZeroCurve3()
  : IZeroCurveFactory(TYPE),
    valueFloating(false),
    mtmAdjustmentFlag(false),
    futuresStubType("None"),   // do nothing special
    futuresGapMethod("None"),
    futuresGapDays(14),
    stubType("AUTO"),
    couponInterpolationMethod("None"),
    zeroInterpolationMethod("None")
{
}


ZeroCurve3::~ZeroCurve3()
{
    //empty
}


int ZeroCurve3::hashCode() const
{
    int hCode = 0;
    hCode ^= CBool::hashCode(adjustDates);
    hCode ^= CBool::hashCode(valueFloating);
    hCode ^= CBool::hashCode(mtmAdjustmentFlag);
    hCode ^= hash_string(futuresStubType);
    hCode ^= hash_string(futuresGapMethod);
    hCode ^= futuresGapDays;
    hCode ^= hash_string(futuresMMDateBadDayConvention);
    hCode ^= hash_string(stubType);
    hCode ^= hash_string(couponInterpolationMethod);
    hCode ^= hash_string(zeroInterpolationMethod);

    if (futuresMMDate.get())
        hCode ^= futuresMMDate->hashCode();

    if (extrapDate.get())
        hCode ^= extrapDate->hashCode();

    return hCode;
}


bool ZeroCurve3::equals(const IZeroCurveFactory* zcm) const
{
    //check type
    if (!ZeroCurve3::TYPE->isInstance(zcm)) return false;

    const ZeroCurve3* zc3 = dynamic_cast<const ZeroCurve3*>(zcm);

    //check instance variables
    if (valueFloating != zc3->valueFloating) return false;
    if (adjustDates != zc3->adjustDates) return false;
    if (mtmAdjustmentFlag != zc3->mtmAdjustmentFlag) return false;
    if (futuresStubType != zc3->futuresStubType) return false;
    if (futuresGapMethod != zc3->futuresGapMethod) return false;
    if (futuresGapDays != zc3->futuresGapDays) return false;
    if (futuresMMDateBadDayConvention != zc3->futuresMMDateBadDayConvention) return false;
    if (stubType != zc3->stubType) return false;
    if (couponInterpolationMethod != zc3->couponInterpolationMethod) return false;
    if (zeroInterpolationMethod != zc3->zeroInterpolationMethod) return false;

    if (futuresMMDate.get())
    {
        if (!zc3->futuresMMDate.get() || !futuresMMDate->equalTo(zc3->futuresMMDate.get()))
        {
            return false;
        }
    }
    else if (zc3->futuresMMDate.get())
    {
        return false;
    }

    if (extrapDate.get())
    {
        if (!zc3->extrapDate.get() || !extrapDate->equalTo(zc3->extrapDate.get()))
        {
            return false;
        }
    }
    else if (zc3->extrapDate.get())
    {
        return false;
    }

    return true;
}


bool ZeroCurve3::isFloatLegValued() const
{
    return valueFloating;
}


void ZeroCurve3::validatePop2Object()
{
    static const string method = "ZeroCurve3::validatePop2Object";

    if (!ZC3CouponInterpolation::validate(couponInterpolationMethod))
    {
        throw ModelException(method, "invalid coupon interpolation type");
    }

    if (!ZC3ZeroInterpolation::validate(zeroInterpolationMethod))
    {
        throw ModelException(method, "invalid zero interpolation type");
    }
}


void ZeroCurve3::validate(
    const DateTime&                today,
    const DateTime&                spotDate,
    const ZeroCurveBenchmarkArray& benchmarks,
    const DayCountConvention*      moneyMarketDcc,
    const DayCountConvention*      swapFixedDcc,
    const DayCountConvention*      swapFloatDcc,
    const DayCountConvention*      swapBasisDcc,
    const MaturityPeriod*          swapFixedIvl,
    const MaturityPeriod*          swapFloatIvl,
    const MaturityPeriod*          swapBasisIvl,
    const BadDayConvention*        swapBadDayConvention,
    const bool                     floatRateFixed,
    const ExpiryArray*             fixDates,
    const DoubleArray*             fixRates,
    IRVolBaseConstSP               volModelIR,
    HolidayConstSP                 holidays,
    const CurrencyBasisWrapper&    ccyBasis,
    const Expiry*                  futuresMaturity
    ) const
{
    static const string method = "ZeroCurve3::validate";

    
    // verify mandatory parameters

    if (!moneyMarketDcc)
    {
        throw ModelException(method, "moneyMarketDcc not specified");
    }

    if (!swapFixedDcc)
    {
        throw ModelException(method, "swapFixedDcc not specified");
    }

    if (!swapFixedIvl)
    {
        throw ModelException(method, "swapFixedIvl not specified");
    }

    if (floatRateFixed)
    {
        if (!fixDates || !fixRates || fixDates->size() == 0 || fixRates->size() == 0)
        {
            throw ModelException(method, "No fixings provided");
        }
    }

    for (int i = 0 ; i < benchmarks.size() ; i++)
    {
        if (!benchmarks[i].get())
        {
            string msg = Format::toString("Benchmark %d is missing", i);
            throw ModelException(method, msg);
        }

        if (benchmarks[i]->isFuture() && adjustDates && !futuresMaturity)
        {
            throw ModelException(method, "Require futuresMaturity when adjusting dates for futures");
        }

        if (benchmarks[i]->isSwap() && !swapBadDayConvention)
        {
            throw ModelException(method, "Require bad day convention when using swaps");
        }

        if (benchmarks[i]->isFuture() && benchmarks[i]->getIncludeFlag() == 2)
        {
            string msg = Format::toString(
                "Invalid benchmark for futures instrument [%s to %s]",
                benchmarks[i]->getStart()->toString().c_str(),
                benchmarks[i]->getEnd() ? benchmarks[i]->getEnd()->toString().c_str() : "-");
            throw ModelException(method, msg);
        }
 
        // tests ALIB-ZC3-26 and ALIB-ZC3-111 fail, so for now disable swaps 2nd pass
        if (benchmarks[i]->isSwap() && benchmarks[i]->getIncludeFlag() == 2)
        {
            throw ModelException(method, "2-pass swap algorithm is not yet implemented in QLib");
        }
    }    
}


ZeroPairSP ZeroCurve3::bootstrap(
        const DateTime&                today,
        const DateTime&                spotDate,
        const ZeroCurveBenchmarkArray& benchmarks,
        const DayCountConvention*      moneyMarketDcc,
        const DayCountConvention*      swapFixedDcc,
        const DayCountConvention*      swapFloatDcc,
        const DayCountConvention*      swapBasisDcc,
        const MaturityPeriod*          swapFixedIvl,
        const MaturityPeriod*          swapFloatIvl,
        const MaturityPeriod*          swapBasisIvl,
        const BadDayConvention*        swapBadDayConvention,
        const bool                     floatRateFixed,
        const ExpiryArray*             fixDates,
        const DoubleArray*             fixRates,
        IRVolBaseConstSP               volModelIR,
        HolidayConstSP                 holidays,
        const CurrencyBasisWrapper&    ccyBasis,
        const Expiry*                  futuresMaturity
    ) const
{
    ZC3ZeroCurveSP szcStub = getStubCurve(spotDate, "Flat");

    // IZeroCurveFactory contract is that validate() has already been called

    DateTimeSP futuresMMDt;
    if (futuresMMDate.get())
    {
        DateTimeSP date(new DateTime(futuresMMDate->toDate(szcStub->getBaseDate())));
        
        if (futuresMMDateBadDayConvention != "")
        {
            BadDayConventionSP bdc(BadDayConventionFactory::make(futuresMMDateBadDayConvention));
            DateTime adjusted = bdc->adjust(*date, holidays.get());
            date.reset(new DateTime(adjusted));
        }

        futuresMMDt.reset(date.get());
    }

    DateTimeSP extrapDt;
    if (extrapDate.get())
    {
        extrapDt.reset(new DateTime(extrapDate->toDate(szcStub->getBaseDate())));
    }

    // ALIB: zcurveol.c#3618 (GtoZeroCurve3L)
    // ALIB: zcurveol.c#3699
    // build fixing curve (ie. SimpleZeroCurve from fixing rates & dates)
    ZeroCurveSP fixingCurve(NULL);
    if (floatRateFixed)
    {
        fixingCurve = ZeroCurveSP(
                new SimpleZeroCurve(
                    szcStub->getBaseDate(),
                    CompoundBasis::SIMPLE,
                    moneyMarketDcc,
                    ZC3ZeroInterpolation::LINEAR,
                    *fixDates,
                    *fixRates
                )
            );
    }

    HolidayConstSP hols(holidays.get());

    ZC3ZeroInterpolationArraySP interpolators = ZC3ZeroInterpolation::createArray(
        zeroInterpolationMethod, szcStub->getBaseDate(), *hols, *swapBadDayConvention);

    return bootstrap(
        today, *szcStub,
        benchmarks,
        *moneyMarketDcc, 
        futuresMMDt.get(), 
        futuresMaturity,
        *swapFixedDcc, 
        swapFloatDcc, 
        swapBasisDcc,
        *swapFixedIvl, 
        swapFloatIvl, 
        swapBasisIvl,
        *swapBadDayConvention, 
        floatRateFixed, 
        fixingCurve.get(),
        volModelIR.get(), 
        *hols, 
        ccyBasis, 
        interpolators, 
        extrapDt.get());
}


// ALIB: zcbuild3.c#1609 (GtoZeroCurve3StrucInterp)
// ALIB: zcbuild3.c#2388 (GtoZeroCurve3CBStrucInterp)
ZeroPairSP ZeroCurve3::bootstrap(
    const DateTime&                    valueDate,
    const ZC3ZeroCurve&                szcStub,
    const ZeroCurveBenchmarkArray&     benchmarks,
    const DayCountConvention&          moneyMarketDcc,
    const DateTime*                    futuresMMDt,
    const Expiry*                      futuresMaturity,
    const DayCountConvention&          swapFixedDcc,
    const DayCountConvention*          swapFloatDcc,
    const DayCountConvention*          swapBasisDcc,
    const MaturityPeriod&              swapFixedIvl,
    const MaturityPeriod*              swapFloatIvl,
    const MaturityPeriod*              swapBasisIvl,
    const BadDayConvention&            swapBadDayConvention,
    const bool                         floatRateFixed,
    const ZeroCurve*                   fixingCurve,
    const IRVolBase*                   volModelIR,
    const Holiday&                     holidays,
    const CurrencyBasisWrapper&        ccyBasis,
    const ZC3ZeroInterpolationArraySP& interpolators,
    const DateTime*                    extrapDt
    ) const
{
    static const string method = "ZeroCurve3::bootstrap";

    ZC3ZeroCurveSP outCurve;
    ZC3RateDataArray ratesDataIn;
    ZC3SwapDataArray swapsDataIn;
    ZC3TurnArray     turnsDataIn;

    try
    {
        bool hasAdjustments = false;
        for (int j = 0 ; j < benchmarks.size() ; j++)
        {
            if (benchmarks[j]->hasAdjustment())
            {
                hasAdjustments = true;
                break;
            }
        }

        // extract data from the inputs
        splitData(ratesDataIn, swapsDataIn, turnsDataIn, szcStub, benchmarks,
            moneyMarketDcc, 
            futuresMaturity,
            valueFloating, 
            swapFixedDcc, swapFloatDcc, swapBasisDcc, 
            swapFixedIvl, swapFloatIvl, swapBasisIvl,
            floatRateFixed, fixingCurve, swapBadDayConvention, 
            holidays, ccyBasis);

        ZC3ZeroCurveSP curves[2];
        curves[0] = ZC3ZeroCurveSP(const_cast<ZC3ZeroCurve*>(&szcStub));
        curves[1] = ZC3ZeroCurveSP(const_cast<ZC3ZeroCurve*>(&szcStub));
    
		ZC3RateDataArray* ratesData = &ratesDataIn;
		ZC3SwapDataArray* swapsData = &swapsDataIn;
		ZC3TurnArray*     turnsData = &turnsDataIn;

        StubPlacement stubPos(stubType);

        // successive calls - one for each segment
        DateTime lastDateUsed;
        for (int i = 0 ; i < interpolators->size() ; i++)
        {
            bool doSomething = true;
            bool lastSection = true;
            bool futuresUsed = false;

            ZC3RateDataArray ratesDataOut;
            ZC3SwapDataArray swapsDataOut;
            ZC3TurnArray     turnsDataOut;

            if (interpolators->size() > 1)
            {
                DateTime to;
                if (i < interpolators->size() - 1)
                {
                    to = (*interpolators)[i]->until();
                }

                DateTime from;
                if (i > 0)
                {
                    from = (*interpolators)[i-1]->until();
                }

                ratesData = &ratesDataOut;
                swapsData = &swapsDataOut;
                turnsData = &turnsDataOut;

                benchmarksExtractInstruments(
                    ratesDataIn, swapsDataIn, turnsDataIn, from, to, 
                    (*interpolators)[i]->isRecursive(),
                    ratesDataOut, swapsDataOut, turnsDataOut,
                    lastDateUsed, doSomething, lastSection, futuresUsed);
            }

            if (doSomething)
            {
                /* 
                 * If the instruments we are adding do not extend the curve we
                 * assume that the curve is a full discount curve rather than
                 * just a stub.
                 */
                // if !hasDisCurve else do nothing ?????
                // zcbuild3.c#2004

                ZC3CouponInterpolationConstSP couponInterpolation = 
                    ZC3CouponInterpolation::make(couponInterpolationMethod);

                const Holiday* basisHolidays = ccyBasis.get() 
                    ? ccyBasis->getHolidays().get() : NULL;

                ZC3FuturesStubRuleSP stubRule;
                ZC3FuturesGapRule gapRule(futuresGapMethod, futuresGapDays);

                if (!ccyBasis.isEmpty() && couponInterpolation.get())
                {
                    /*
				     * Note: as input instrument arrays are amended must use copies
				     * if an estimating curve is to be built as well.
                     *
                     * Also if using coupon interpolation the swaps must be done after 
                     * cash and futures, so that an estimating curve is available
                     * for coupon rate calculation.
				     */                    

                    // discounting curve cash & futures
                    stubRule = ZC3FuturesStubRule::make(futuresStubType.empty() ? "N" : futuresStubType);
                    curves[0] = ZC3ZeroCurveSP( new ZC3ZeroCurve(
                        *curves[0],
                        NULL,
                        *ratesData,
                        ZC3SwapDataArray(),    // ignore swaps at first
                        ZC3TurnArray(*turnsData),
                        hasAdjustments,
                        mtmAdjustmentFlag || futuresUsed, 
                        stubRule,
                        gapRule,
                        futuresMMDt, 
                        volModelIR,
                        stubPos, 
                        swapBadDayConvention, 
                        valueFloating, 
                        swapFloatIvl,  
                        fixingCurve,
                        couponInterpolation.get(), 
                        *(*interpolators)[i].get(), 
                        holidays, 
                        basisHolidays ? *basisHolidays : holidays,
                        lastSection ? extrapDt : NULL,
                        true));

                    // estimating curve cash & futures
                    stubRule = ZC3FuturesStubRule::make(futuresStubType.empty() ? "N" : futuresStubType);
                    curves[1] = ZC3ZeroCurveSP( new ZC3ZeroCurve(
						*curves[1],
						curves[0].get(),
                        *ratesData,
                        ZC3SwapDataArray(),    // ignore swaps at first
                        ZC3TurnArray(*turnsData),
						hasAdjustments,
						mtmAdjustmentFlag || futuresUsed, 
						stubRule,
						gapRule,
						futuresMMDt, 
						volModelIR,
						stubPos, 
						swapBadDayConvention, 
						valueFloating, 
						swapFloatIvl,  
						fixingCurve,
						couponInterpolation.get(), 
						*(*interpolators)[i].get(), 
						holidays, 
						basisHolidays ? *basisHolidays : holidays,
						lastSection ? extrapDt : NULL,
                        true));

                    /*
                     * Do coupon interpolation using current estimating curve.
                     */
                    ZC3SwapDataArraySP swapsDataInterpolated = couponInterpolation->getInterpRates(
                        swapsDataIn,
                        *curves[0],
                        curves[1].get(), 
                        swapBadDayConvention,
                        stubPos,
                        holidays,
                        basisHolidays ? *basisHolidays : holidays,
                        false,
                        volModelIR,
                        true);

                    // discounting curve swaps
                    stubRule = ZC3FuturesStubRule::make(futuresStubType.empty() ? "N" : futuresStubType);
                    curves[0] = ZC3ZeroCurveSP( new ZC3ZeroCurve(
                        *curves[0],
                        NULL,
                        ZC3RateDataArray(),
                        *swapsDataInterpolated,
                        ZC3TurnArray(),
                        hasAdjustments,
                        mtmAdjustmentFlag || futuresUsed, 
                        stubRule,
                        gapRule,
                        futuresMMDt, 
                        volModelIR,
                        stubPos, 
                        swapBadDayConvention, 
                        valueFloating, 
                        swapFloatIvl,  
                        fixingCurve,
                        NULL,   // coupon interpolation already done
                        *(*interpolators)[i].get(), 
                        holidays, 
                        basisHolidays ? *basisHolidays : holidays,
                        lastSection ? extrapDt : NULL,
                        true));

                    // estimating curve swaps
                    stubRule = ZC3FuturesStubRule::make(futuresStubType.empty() ? "N" : futuresStubType);
					curves[1] = ZC3ZeroCurveSP( new ZC3ZeroCurve(
						*curves[1],
						curves[0].get(),
                        ZC3RateDataArray(),
                        *swapsDataInterpolated,
                        ZC3TurnArray(),
						hasAdjustments,
						mtmAdjustmentFlag || futuresUsed, 
						stubRule,
						gapRule,
						futuresMMDt, 
						volModelIR,
						stubPos, 
						swapBadDayConvention, 
						valueFloating, 
						swapFloatIvl,  
						fixingCurve,
                        NULL,   // coupon interpolation already done
						*(*interpolators)[i].get(), 
						holidays, 
						basisHolidays ? *basisHolidays : holidays,
						lastSection ? extrapDt : NULL,
                        true));
                }
                else if (!ccyBasis.isEmpty())
                {
                    /*
				     * Note: as input instrument arrays are amended must use copies
				     * if an estimating curve is to be built as well.
				     */                    
                    stubRule = ZC3FuturesStubRule::make(futuresStubType.empty() ? "N" : futuresStubType);
                    curves[0] = ZC3ZeroCurveSP( new ZC3ZeroCurve(
                        *curves[0],
                        NULL,
                        ZC3RateDataArray(*ratesData),
                        ZC3SwapDataArray(*swapsData),
                        ZC3TurnArray(*turnsData),
                        hasAdjustments,
                        mtmAdjustmentFlag || futuresUsed, 
                        stubRule,
                        gapRule,
                        futuresMMDt, 
                        volModelIR,
                        stubPos, 
                        swapBadDayConvention, 
                        valueFloating, 
                        swapFloatIvl,  
                        fixingCurve,
                        couponInterpolation.get(), 
                        *(*interpolators)[i].get(), 
                        holidays, 
                        basisHolidays ? *basisHolidays : holidays,
                        lastSection ? extrapDt : NULL,
                        true));

                    stubRule = ZC3FuturesStubRule::make(futuresStubType.empty() ? "N" : futuresStubType);
					curves[1] = ZC3ZeroCurveSP( new ZC3ZeroCurve(
						*curves[1],
						curves[0].get(),
                        *ratesData,
                        *swapsData,
                        *turnsData,
						hasAdjustments,
						mtmAdjustmentFlag || futuresUsed, 
						stubRule,
						gapRule,
						futuresMMDt, 
						volModelIR,
						stubPos, 
						swapBadDayConvention, 
						valueFloating, 
						swapFloatIvl,  
						fixingCurve,
						couponInterpolation.get(), 
						*(*interpolators)[i].get(), 
						holidays, 
						basisHolidays ? *basisHolidays : holidays,
						lastSection ? extrapDt : NULL,
                        true));
                }
                else
                {
                    stubRule = ZC3FuturesStubRule::make(futuresStubType.empty() ? "N" : futuresStubType);
                    curves[0] = ZC3ZeroCurveSP( new ZC3ZeroCurve(
                        *curves[0],
                        NULL,
                        *ratesData,
                        *swapsData,
                        *turnsData,
                        hasAdjustments,
                        mtmAdjustmentFlag || futuresUsed, 
                        stubRule,
                        gapRule,
                        futuresMMDt, 
                        volModelIR,
                        stubPos, 
                        swapBadDayConvention, 
                        valueFloating, 
                        swapFloatIvl,  
                        fixingCurve,
                        couponInterpolation.get(), 
                        *(*interpolators)[i].get(), 
                        holidays, 
                        basisHolidays ? *basisHolidays : holidays,
                        lastSection ? extrapDt : NULL,
                        false));
                }
            }
        }

        // last returned curve(s) are results
        curves[0]->cleanup();

        if (!ccyBasis.isEmpty())
        {
            curves[1]->cleanup();
        }
        else
        {
            curves[1] = curves[0];
        }

        ZeroPairSP zp(new ZeroPair());
        zp->discZC = curves[0];
        zp->growZC = curves[1];
        return zp;
    }
    catch(ModelException& e)
    {
        throw ModelException(e, method);
    }
}


// ALIB: zcbuild3.c#1443
void ZeroCurve3::benchmarksExtractInstruments(
    // inputs
    const ZC3RateDataArray& ratesDataIn,
    const ZC3SwapDataArray& swapsDataIn,
    const ZC3TurnArray&     turnsDataIn,
    const DateTime&         fromDt,
    const DateTime&         to,
    bool                    recursive,
    // outputs
    ZC3RateDataArray&       ratesDataOut,
    ZC3SwapDataArray&       swapsDataOut,
    ZC3TurnArray&           turnsDataOut,
    DateTime&               lastDateUsed,
    bool&                   doSomething,
    bool&                   lastSection,
    bool&                   futuresUsed) const
{
    static const string method ="ZeroCurve3::benchmarksExtractInstruments";

    DateTime lastCurveDate = lastDateUsed;
    doSomething = false;
    lastSection = true;
    DateTime from = max(fromDt, lastCurveDate);

    // check cash / futures
    size_t i;
    for (i = 0 ; i < ratesDataIn.size() ; i++)
    {
        DateTime endDate = ratesDataIn[i]->getEndDate();
        if (endDate > from && (to.empty() || endDate <= to))
        {
            ZC3RateDataSP tmp = ratesDataIn[i];
            ratesDataOut.push_back(tmp);

            if (typeid(*ratesDataIn[i]) == typeid(ZC3FuturesData))
            {
                // if we have already added a future set stub method to zero
                futuresUsed = true;
            }

            if (endDate > lastDateUsed)
            {
                lastDateUsed = endDate;
            }
        }
        else
        {
            // isRecursiveIncludedInFlatSection() only applies to swaps 
            if (recursive && endDate > lastDateUsed)
            {
                if (endDate > lastDateUsed)
                {
                    lastDateUsed = endDate;
                }
            }
            else if (!to.empty() && endDate > to)
            {
                lastSection = false;
            }
        }
    }

    // repeat for swaps
    for (i = 0 ; i < swapsDataIn.size() ; i++)
    {
        DateTime endDate = swapsDataIn[i]->getEndDate();
        if (endDate > from && (to.empty() || endDate <= to))
        {
            swapsDataOut.push_back(swapsDataIn[i]);

            if (endDate > lastDateUsed)
            {
                lastDateUsed = endDate;
            }
        }
        else
        {
            /*
             * With recursive bootstrapping we need some way of indicating which
             * of the benchmarks are in the final flat section. We do this by
             * changing them to isRecursiveIncludedInFlatSection() = true.
             */
            if (recursive && (to.empty() || endDate > to))
            {
                if (swapsDataIn[i]->isForPass(1))
                {
                    swapsDataIn[i]->setRecursiveIncludedInFlatSection();
                }

                if (endDate > lastDateUsed)
                {
                    lastDateUsed = endDate;
                }
            }
            else if (!to.empty() && endDate > to)
            {
                lastSection = false;
            }
        }
    }

    /* 
     * The from and to dates are different for adjustments because they are tied
     * to a particular instrument.
     */
    for (i = 0 ; i < turnsDataIn.size() ; i++)
    {
        DateTime endDate = turnsDataIn[i]->getEndDate();
        if (endDate > lastCurveDate && endDate <= lastDateUsed)
        {
            turnsDataOut.push_back(turnsDataIn[i]);
        }
        else if (!to.empty() && endDate > to)
        {
            lastSection = false;
        }
    }

    /*
     * We finally check if no benchmarks are set - in which case there is 
     * nothing to do.
     */
    if (!doSomething && (ratesDataOut.size() > 0 || turnsDataOut.size() > 0))
    {
        doSomething = true;
    }

    if (!doSomething)
    {
        for (i = 0 ; i < swapsDataOut.size() ; i++)
        {
            // isRecursiveIncludedInFlatSection() only applies to swaps 
            if (!swapsDataOut[i]->isRecursiveIncludedInFlatSection())
            {
                doSomething = true;
            }
        }
    }
}


// ALIB:: zcbuild3.c#2968 (SplitData)
void ZeroCurve3::splitData(
    ZC3RateDataArray&              ratesData,
    ZC3SwapDataArray&              swapsData,
    ZC3TurnArray&                  turnsData,
    const ZC3ZeroCurve&            szcStub,
    const ZeroCurveBenchmarkArray& benchmarks,
    const DayCountConvention&      mmDcc,
    const Expiry*                  futuresMaturity,
    const bool                     valueFloating,
    const DayCountConvention&      fixedDcc,
    const DayCountConvention*      floatDcc,       // only used if value floating
    const DayCountConvention*      basisDcc,       // only used if basis rates not NULL
    const MaturityPeriod&          fixedIvl,
    const MaturityPeriod*          floatIvl,       // only used if value floating
    const MaturityPeriod*          basisIvl,       // only used if basis rates not NULL
    const bool                     floatRateFixed, // only used if value floating
    const ZeroCurve*               fixCurve,       // only used if float rate fixed
    const BadDayConvention&        badDayConv,
    const Holiday&                 holidays,
    const CurrencyBasisWrapper&    ccyBasis) const
{
    static const string method = "ZeroCurve3::splitData";

    BadDayFollowing bdf;
    const DateTime& baseDate = szcStub.getBaseDate();
    bool hasRelativeTurns = false;
    bool hasUTurns = false;

    for (int i = 0 ; i < benchmarks.size() ; i++)
    {
        DateTime startDate = benchmarks[i]->getStart() ? benchmarks[i]->getStart()->toDate(baseDate) : baseDate;

        if (benchmarks[i]->isCash())
        {
            MoneyMarketBenchmark& benchmark = static_cast<MoneyMarketBenchmark&>(*benchmarks[i]);
            DateTime endDate = benchmark.getEnd()->toDate(startDate);

            /*
             * For backwards compatibility with previous version of ZeroCurve3
             * we need to be able to specify that cash maturity dates are bad
             * day adjusted.  For ON (1D) it makes no sense to use modified
             * following, only following or none.
             */
            if (adjustDates)
            {
                if (endDate.daysDiff(baseDate) < 7 && !BadDayNone::TYPE->isInstance(&badDayConv))
                {
                    endDate = bdf.adjust(endDate, &holidays);
                }
                else
                {
                    endDate = badDayConv.adjust(endDate, &holidays);
                }
            }

            double rate = benchmark.getRate();
            double basisRate = 0.0;
            if (!ccyBasis.isEmpty() && ccyBasis.get())
            {
                basisRate = ccyBasis->cashBasis(*benchmark.getEnd());
            }

            ZC3RateData* instrument = new ZC3CashData(
                startDate,
                endDate,
                rate,
                benchmark.getDcc() ? *benchmark.getDcc() : mmDcc,
                basisRate);

            ratesData.push_back(ZC3RateDataSP(instrument));
        }
        else if (benchmarks[i]->isFuture())
        {
            FuturesBenchmark& benchmark = static_cast<FuturesBenchmark&>(*benchmarks[i]);
            DateTime endDate;

            double rate = benchmark.getRate();

            /*
             * If an end expiry has not been provided futuresMaturity is used
             * to calculate the expiry from the future's start date.
             */
            if (adjustDates || !benchmark.getEnd())
            {
                startDate = benchmark.getStart()->toDate(baseDate);
                startDate = bdf.adjust(startDate, &holidays);
                endDate = futuresMaturity->toDate(startDate);
                endDate = bdf.adjust(endDate, &holidays);
            }
            else
            {
                 endDate = benchmark.getEnd()->toDate(startDate);
            }

            if (adjustDates)
            {
                // this is for backwards compatibility with first ZC3 implementation
                rate = (FuturesBenchmark::ZERO_PERCENT - rate) / FuturesBenchmark::ZERO_PERCENT;
            }

            double basisRate = 0.0;
            if (!ccyBasis.isEmpty() && ccyBasis.get())
            {
                basisRate = ccyBasis->futuresBasis(*benchmark.getStart());

                if (adjustDates)
                {
                    // basis rate is a change in the futures PRICE
                    basisRate = -1 * basisRate / FuturesBenchmark::ZERO_PERCENT;
                }
            }

            ZC3RateData* instrument = new ZC3FuturesData(
                startDate,
                endDate,
                rate,
                benchmark.getDcc() ? *benchmark.getDcc() : mmDcc,
                benchmark.getAdjustment(),
                basisRate);
            
            ratesData.push_back(ZC3RateDataSP(instrument));
        }
        else if (benchmarks[i]->isTurn())
        {
            TurnBenchmark& benchmark = static_cast<TurnBenchmark&>(*benchmarks[i]);
            DateTime endDate = benchmark.getEnd()->toDate(startDate);

            ZC3Turn* instrument;

            if (benchmark.isAbsolute())
            {
                instrument = new ZC3AbsoluteTurn(
                    startDate, 
                    endDate, 
                    benchmark.getRate(),
                    benchmark.getDcc() ? *benchmark.getDcc() : mmDcc);
            }
            else if (benchmark.isRelative())
            {
                hasRelativeTurns = true;
                instrument = new ZC3RelativeTurn(
                    startDate, 
                    endDate, 
                    benchmark.getRate(),
                    benchmark.getDcc() ? *benchmark.getDcc() : mmDcc);
            }
            else if (benchmark.isImplied())
            {
                instrument = new ZC3ImpliedTurn(
                    startDate, 
                    endDate, 
                    benchmark.getRate(),
                    benchmark.getDcc() ? *benchmark.getDcc() : mmDcc);
            }
            else
            {
                string msg = Format::toString("Invalid instrument type: %s", 
                    benchmark.getEnd()->toString().c_str());
                throw ModelException(method, msg);
            }

            turnsData.push_back(ZC3TurnSP(instrument));

            // cf. szcrates.c#1632
            if (hasUTurns && hasRelativeTurns)
            {
                throw ModelException(method, "Cannot mix relative and U-turn adjustments");
            }
        }
        else if (benchmarks[i]->isSwap())
        {
            SwapBenchmark& benchmark = static_cast<SwapBenchmark&>(*benchmarks[i]);
            double adjustment = benchmark.getAdjustment();

            double rate = benchmark.getRate();
            double basisRate = 0.0;
            if (!ccyBasis.isEmpty() && ccyBasis.get())
            {
                basisRate = ccyBasis->swapBasis(
                    rate, 
                    *benchmark.getEnd(), 
                    mmDcc, 
                    floatIvl ? floatIvl->approxAnnualFrequency() : fixedIvl.approxAnnualFrequency(),
                    fixedDcc, 
                    fixedIvl.approxAnnualFrequency());
            }

            ZC3SwapData* instrument = new ZC3SwapData(
                baseDate,
                startDate,
                benchmark.getEnd()->toDate(startDate),
                rate,
                benchmark.hasAdjustment() ? &adjustment : NULL,
                benchmark.getPrice(),
                benchmark.getFixedIvl() ? *benchmark.getFixedIvl() : fixedIvl,
                benchmark.getFixedDcc() ? *benchmark.getFixedDcc() : fixedDcc,
                benchmark.getFloatIvl() ? benchmark.getFloatIvl() : floatIvl,
                benchmark.getFloatDcc() ? benchmark.getFloatDcc() : floatDcc,
                benchmark.getIncludeFlag(),
                valueFloating || floatRateFixed || !ccyBasis.isEmpty(),
                floatRateFixed,
                holidays,
                badDayConv,
                fixCurve,
                benchmark.getBasisIvl() ? benchmark.getBasisIvl() : basisIvl,
                benchmark.getBasisDcc() ? benchmark.getBasisDcc() : basisDcc,
                basisRate);

            swapsData.push_back(ZC3SwapDataSP(instrument));
        }
        else
        {
            string msg = Format::toString("Invalid instrument type: %s", 
                benchmarks[i]->getEnd()->toString().c_str());
            throw ModelException(method, msg);
        }
    }
}


// ALIB: zcsmooth.c#1394
ZC3ZeroCurveSP ZeroCurve3::getStubCurve(const DateTime& baseDate, const string& stubInterpType) const
{
    // TBD!! GtoSZCNewFromTCurve equivalent
    // ... for now do not permit stub curve
    ZC3ZeroCurveSP zc(new ZC3ZeroCurve(baseDate, stubInterpType));
    return zc;
}


/*
 * Reflection support.
 */

static IObject* defaultZeroCurve3()
{
    return new ZeroCurve3();
}


/** Invoked when Class is 'loaded' */
void ZeroCurve3::load(CClassSP& clazz)
{
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(ZeroCurve3, clazz);
    SUPERCLASS(IZeroCurveFactory);
    EMPTY_SHELL_METHOD(defaultZeroCurve3);

    //fields
    FIELD(valueFloating,                "value floating leg");
    FIELD(mtmAdjustmentFlag,            "MTM adjustment flag");
    FIELD(futuresStubType,              "Futures stub type");
    FIELD(futuresGapMethod,             "Futures/swaps gap method");
    FIELD(futuresGapDays,               "Futures/swaps gap days");
    FIELD       (futuresMMDate,                "MM date to keep in curve");
    FIELD(futuresMMDateBadDayConvention,"bad day convention for MM date to keep in curve");
    FIELD(couponInterpolationMethod,    "coupon interpolation style");
    FIELD(zeroInterpolationMethod,      "zero interpolation style");
    FIELD       (extrapDate,                   "Extrapolation date");
    FIELD(stubType,                     "Stub type [AUTO/FRONT/BACK]");

    FIELD_MAKE_OPTIONAL(valueFloating);
    FIELD_MAKE_OPTIONAL(mtmAdjustmentFlag);
    FIELD_MAKE_OPTIONAL(futuresStubType);
    FIELD_MAKE_OPTIONAL(futuresGapMethod);
    FIELD_MAKE_OPTIONAL(futuresGapDays);
    FIELD_MAKE_OPTIONAL(futuresMMDate);
    FIELD_MAKE_OPTIONAL(futuresMMDateBadDayConvention);
    FIELD_MAKE_OPTIONAL(couponInterpolationMethod);
    FIELD_MAKE_OPTIONAL(extrapDate);
    FIELD_MAKE_OPTIONAL(stubType);

    Addin::registerConstructor("ZERO_CURVE3",
                               Addin::MARKET,
                               "Creates factory for bootstrapping ALIB-style curves",
                               ZeroCurve3::TYPE);
}


CClassConstSP const ZeroCurve3::TYPE = 
    CClass::registerClassLoadMethod("ZeroCurve3", typeid(ZeroCurve3), load);



DRLIB_END_NAMESPACE
