//----------------------------------------------------------------------------
//
//   Group       : QR&D Core Analytics Team
//
//   Filename    : RecurseInfo.cpp
//
//   Description : Helper class for ALIB zero curve 3 bootstrapping method.
//
//   Author      : Richard Appleton
//
//   Date        : 2nd May 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ZC3RecurseInfo.hpp"
#include "edginc/ZC3ZeroCurve.hpp"
#include "edginc/ZC3CurveInstrument.hpp"
#include "edginc/ZC3FuturesGapRule.hpp"
#include "edginc/ZC3Interpolation.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SwapTool.hpp"


DRLIB_BEGIN_NAMESPACE


static const double BOOTSTRAP_ERROR_MAX = 1.0E-12;
static const int    RECURSIVE_MAX_PASSES = 40;


// ALIB: zrecurs.c#113
ZC3RecurseInfo::ZC3RecurseInfo(
        const ZC3SwapDataArray& swapsData,
        const BadDayConvention& badDayConv,
        const Holiday&          holidays)
  : cfls(0), cashStreams(0), passIndex(0), swapIndex(0)
{
    static const string method = "ZC3RecurseInfo::ZC3RecurseInfo";

    bool lastFittingDatePassed = false;
    bool firstValidFound = false;
    DateTime unadjFittingDate;

    for (size_t i = 0 ; i < swapsData.size() ; i++)
    {
        if (swapsData[i]->isForPass(1))
        {
            if (lastFittingDatePassed)
            {
                string msg = Format::toString(
                    "Swap [end date %s] benchmark of 1 after last fitting date",
                    swapsData[i]->getEndDate().toString().c_str());
                throw ModelException(method, msg);
            }

            unadjFittingDate = swapsData[i]->getEndDate();
        }

        if (swapsData[i]->isForPass(2))
        {
            string msg = "Cannot use benchmark of 2 with recursive interp type";
            throw ModelException(method, msg);
        }

        if (swapsData[i]->isRecursiveIncludedInFlatSection())
        {
            lastFittingDatePassed = true;
            swapsData[i]->setBenchmark(1);
        }
    }

    // now find common start date for swaps
    for (size_t j = 0 ; j < swapsData.size() ; j++)
    {
        if (swapsData[j]->isForPass(1))
        {
            if (firstValidFound)
            {
                if (swapsData[j]->getStartDate() != startDate)
                {
                    string msg = Format::toString(
                        "Swaps have different start dates [%s and %s]",
                        swapsData[j]->getStartDate().toString().c_str(),
                        startDate.toString().c_str());
                    throw ModelException(method, msg);
                }
            }
            else
            {
                startDate = swapsData[j]->getStartDate();
                firstValidFound = true;
            }
        }
    }

    if (!firstValidFound)
    {
        string msg = "The interp type used requires swap instruments";
        throw ModelException(method, msg);
    }

    lastSwapDate = unadjFittingDate;
    lastFittingDate = badDayConv.adjust(unadjFittingDate, &holidays);
}


bool ZC3RecurseInfo::isFirstPass() const
{
    return passIndex == 0;
}


// ALIB: szcswaps.c#992
void ZC3RecurseInfo::getCurrent(CashStream*& cs, CashFlowArray*& cfl)
{
    cs = &cashStreams[swapIndex];
    cfl = &cfls[swapIndex];
}


// ALIB: szcswaps.c#1187
void ZC3RecurseInfo::next(const CashStream& cs, const CashFlowArray& cfl)
{
    /*
     * If we are doing this recursively on the first pass we steal the
     * cfl's and cashstreams which can then be reused on later passes.
     */
    if (isFirstPass())
    {
        cashStreams.push_back(cs);
        cfls.push_back(cfl);
    }

    swapIndex++;
}


void ZC3RecurseInfo::reset()
{
    swapIndex = 0;
}


// ALIB:: zrecurs.c#275 (GtoZCRecurseEvaluateBootstrap)
bool ZC3RecurseInfo::evaluateBootstrap(
    ZC3ZeroCurve&               szc,
    const DateTimeArray&        cashDateList,
    const ZC3SwapDataArray&     swapsData,
    const BadDayConvention&     badDayConv,
    const Holiday&              holidays,
    const Holiday&              basisHolidays,
    bool                        valueFloating,
    bool                        withBasis,
    bool                        monthEndRoll,
    const MaturityPeriod*       floatIvl,
    const ZC3FuturesGapRule&    futuresGapRule,
    const ZC3ZeroInterpolation& zeroInterpType)
{
    static const string method = "ZC3RecurseInfo::evaluateBootstrap";

    DateTime lastCashDate = cashDateList[cashDateList.size() - 1];
    
    /* 
     * Find the sets of dates used for the recursive fitting. These do not
     * change every iteration, so we only need to calculate on the first pass.
     */
    if (passIndex == 0)
    {
        if (!floatIvl)
        {
            string msg = "floating leg period required for recursive interpolation type";
            throw ModelException(method, msg);
        }

        findSwapBenchmarkDates(cashDateList, lastCashDate, swapsData, valueFloating, withBasis, 
            *floatIvl, badDayConv, holidays, basisHolidays, futuresGapRule);

        findSwapCriticalDates(szc, lastCashDate, badDayConv, holidays, 
            valueFloating, withBasis, monthEndRoll, *floatIvl);
    }

    bool complete = evaluateRecursiveBootstrap(szc, zeroInterpType);
    passIndex++;

    return complete;
}


// ALIB: zcrecurs.c#568
void ZC3RecurseInfo::findSwapBenchmarkDates(
    const DateTimeArray&     cashDateList,
    const DateTime&          pLastCashDate,
    const ZC3SwapDataArray&  swapsData,
    bool                     valueFloating,
    bool                     withBasis,
    const MaturityPeriod&    floatIvl,
    const BadDayConvention&  badDayConv,
    const Holiday&           holidays,
    const Holiday&           basisHolidays,
    const ZC3FuturesGapRule& futuresGapRule)
{
    static const string method = "ZC3RecurseInfo::findSwapBenchmarkDates";

    ignoreDates.resize(0);
    benchmarkDates.resize(0);

    // populate date list from cash date list
    std::copy(cashDateList.begin(), cashDateList.end(), benchmarkDates.back_inserter());

    DateTime lastCashDate = futuresGapRule.adjust(pLastCashDate);

    for (size_t i = 0 ; i < swapsData.size() ; i++)
    {
        if (swapsData[i]->isForPass(1))
        {
            DateTime adjDate = badDayConv.adjust(swapsData[i]->getEndDate(), &holidays);

            if (adjDate > lastCashDate)
            {
                /* 
                 * If we are valuing the floating leg, the segment may actually
                 * extend to the end of the floating leg fixing. We can get 
                 * this date by counting back one floating interval from the 
                 * unadjusted maturity date, correcting for business days, 
                 * going forward one floating interval and again correcting for
                 * business days.
                 */
                if (valueFloating)
                {
                    DateTime fixStartDate = DateFwdThenAdjust(swapsData[i]->getEndDate(), floatIvl, -1, badDayConv, holidays);
                    DateTime fixEndDate = DateFwdThenAdjust(fixStartDate, floatIvl, 1, badDayConv, holidays);

                    if (fixEndDate > adjDate)
                    {
                        ignoreDates.push_back(adjDate);
                        adjDate = fixEndDate;
                    }
                }

                // The same issue can effect basis payments
                if (withBasis && !Maths::isZero(swapsData[i]->getBasisRate()))
                {
                    DateTime basisPaymentDate = badDayConv.adjust
                        (swapsData[i]->getEndDate(), &basisHolidays);

                    if (basisPaymentDate > adjDate)
                    {
                        ignoreDates.push_back(adjDate);
                        adjDate = basisPaymentDate;
                    }
                }

                benchmarkDates.push_back(adjDate);
            }
        }
    }
}


// ALIB: zcrecurs.c#736
void ZC3RecurseInfo::findSwapCriticalDates(
    const ZC3ZeroCurve&     szc,
    const DateTime&         lastCashDate,
    const BadDayConvention& badDayConv,
    const Holiday&          holidays,
    bool                    valueFloating,
    bool                    withBasis,
    bool                    monthEndRoll,
    const MaturityPeriod&   floatIvl)
{
    static const string method = "ZC3RecurseInfo::findSwapCriticalDates";

    MaturityPeriod interval(monthEndRoll ? "3F" : "3M");

    // And we add 3M segments
	DateTimeArraySP quart(SwapTool::dateArray(startDate, startDate, lastSwapDate, interval, true));

	DateTimeArray dateList;
	for (int i = 0 ; i < quart->size() ; i++)
	{
		(*quart)[i] = badDayConv.adjust((*quart)[i], &holidays);

		if ((*quart)[i] > lastCashDate && (*quart)[i] < lastFittingDate)
		{
			dateList.push_back((*quart)[i]);
		}
	}

	if (dateList.empty())
	{
		throw ModelException(method, "No swap critical dates found");
	}

    // We now do the subtraction
    DateTimeArraySP subDateList = DateTime::subtract(dateList, benchmarkDates);
	if (subDateList->empty())
	{
		throw ModelException(method, "No swap critical dates found (after removing benchmark dates)");
	}

    if ((valueFloating || withBasis) && (ignoreDates.size() > 0))
    {
        DateTimeArraySP withoutIgnoredList = 
            DateTime::subtract(*subDateList, ignoreDates);
		if (withoutIgnoredList->empty())
		{
			throw ModelException(method, "No swap critical dates found (after removing dates to ignore)");
		}

		fittingDates = *withoutIgnoredList;
    }
    else
    {
        fittingDates = *subDateList;
    }
}


// ALIB: zcrecurs.c#852
bool ZC3RecurseInfo::evaluateRecursiveBootstrap(
    ZC3ZeroCurve&               szc,
    const ZC3ZeroInterpolation& zeroInterpType)
{
    static const string method = "ZC3RecurseInfo::evaluateRecursiveBootstrap";

    if (benchmarkDates.empty() || fittingDates.empty())
    {
        throw ModelException(method, "Null inputs");
    }

	try
	{
		// We first disable the adjustments in the curve being built
		// ... only a few adjustments are expected, so copying is efficient enough
		ZC3AdjustmentArray adjustments(szc.adjustments.size());
		std::copy(szc.adjustments.begin(), szc.adjustments.end(), adjustments.begin());
		szc.adjustments.resize(0);

    
		DoubleArray benchmarkInputs(benchmarkDates.size());

		for (int i = 0 ; i < benchmarkDates.size() ; i++)
		{
			benchmarkInputs[i] = szc.discountFactor(benchmarkDates[i]);
		}

		DoubleArraySP fittedDFs = zeroInterpType.spline
			(benchmarkDates, fittingDates, benchmarkInputs);

		/* 
		 * Our criterion for fitting will be based on the normalized difference
		 * between the fitted and input discount factors.  We shall minimize based 
		 * on the sum of the squares of these differences.
		 */
		double bootstrapError = 0.0;
		for (int j = 0 ; j < fittingDates.size() ; j++)
		{
			double bootstrapDF = szc.discountFactor(fittingDates[j]);

			/* 
			 * When normalizing the difference we potentially have the choice of
			 * normalizing with bootstrapDF or fittedDF. 
			 */
			double normalizedDiff = 0.0;
			if (!Maths::isZero(bootstrapDF))
			{
				double fittedDF = (*fittedDFs)[j];
				normalizedDiff = (fittedDF - bootstrapDF) / bootstrapDF;
			}

			bootstrapError += (normalizedDiff * normalizedDiff);
		}

		bootstrapError /= fittingDates.size();

		/* 
		 * We need at least two iterations. This is because the second iteration
		 * will add the 3M critical points.
		 */
		int numPasses = passIndex + 1;
		bool complete = (bootstrapError < BOOTSTRAP_ERROR_MAX && numPasses > 1);

		if (!complete && numPasses >= RECURSIVE_MAX_PASSES)
		{
            string msg = Format::toString(
				"Curve build has not converged in %d  passes.  Bootstrap Error = %e",
				numPasses, bootstrapError);
			throw ModelException(method, msg);
		}

		/*
		 * Finally we merge the fitted discount factors with the input discount
		 * factors to make a single discount curve. Note that there should not be 
		 * any coincident days.
         *
         * FourPlusIZeroCurve is used here as we need to add pre-calculated discount
         * factors to a curve rather than bootstrap it.
		 */
		int size = benchmarkDates.size() + fittingDates.size();
		fittedCurve = FourPlusIZeroCurveSP(new FourPlusIZeroCurve(szc.getBaseDate(), size));

		int bmIndex = 0;
		int fitIndex = 0;

		while (bmIndex < benchmarkDates.size() && fitIndex < fittingDates.size())
		{
			DateTime bmDate = benchmarkDates[bmIndex];
			DateTime fitDate = fittingDates[fitIndex];

			if (bmDate < fitDate)
			{
				fittedCurve->addDiscountFactor(bmDate, benchmarkInputs[bmIndex]);
				bmIndex++;
			}
			else if (bmDate > fitDate)
			{
				fittedCurve->addDiscountFactor(fitDate, (*fittedDFs)[fitIndex]);
				fitIndex++;
			}
			else
			{
                string msg = Format::toString(
                    "Program bug. There should be no coincident days [%s] [%s line %d]",
                    bmDate.toString().c_str(), __FILE__, __LINE__);
				throw ModelException(method, msg);
			}
		}

		// Add any left over DF's
		while (bmIndex < benchmarkDates.size())
		{
			fittedCurve->addDiscountFactor
				(benchmarkDates[bmIndex], benchmarkInputs[bmIndex]);
			bmIndex++;
		}


		while (fitIndex < fittingDates.size())
		{
			fittedCurve->addDiscountFactor
				(fittingDates[fitIndex], (*fittedDFs)[fitIndex]);
			fitIndex++;
		}

		// We have finished calculating discounts without adjustments
		std::copy(adjustments.begin(), adjustments.end(), szc.adjustments.begin());

		return complete;
	}
	catch(exception& e)
	{
		throw ModelException(e, method);
	}
}


// ALIB: zcrecurs.c#383 (GtoZCRecurseAddNonBenchmarkDiscounts)
void ZC3RecurseInfo::addNonBenchmarkDiscounts(
    ZC3ZeroCurve&        discountCurve,
    const CashFlowArray& cfl) const
{
    static const string method = "ZC3RecurseInfo::addNonBenchmarkDiscounts";

	try
	{
		// We ignore all flows up to and including the last point in the curve
		DateTime lastCurveDate = discountCurve.endDate();
		DateTime lastCFLDate = cfl.back().date;

		double bootstrapPrevDiscount = discountCurve.data[discountCurve.data.size() - 1].discount;
		double fittedPrevDiscount = fittedCurve->discountFactor(lastCurveDate);
		double base = bootstrapPrevDiscount / fittedPrevDiscount;

		for (int i = 0 ; i < fittingDates.size() ; i++)
		{
			DateTime date = fittingDates[i];

			if (date <= lastCurveDate)
			{
				continue;
			}

			if (date >= lastCFLDate)
			{
				break;	// all done
			}

			double fittedDiscount = fittedCurve->discountFactor(date);
			double discount = fittedDiscount * base;
			addAdjustedDiscount(discountCurve, discount, date);
		}
	}
	catch(exception& e)
	{
		throw ModelException(e, method);
	}
}


// ALIB: szcbuild.c#500 (GtoZCRecurseAddAdjustedDiscount)
void ZC3RecurseInfo::addAdjustedDiscount(
    ZC3ZeroCurve&   curve, 
    double          discount,
    const DateTime& discDate) const
{
    static const string method= "ZC3RecurseInfo::addAdjustedDiscount";

	curve.validate(false);

	if (curve.noBootstrap)
	{
		throw ModelException(method, "Cannot bootstrap further with this curve");
	}

	ZC3CurveSegment& prevSegment = curve.data[curve.data.size() - 1];
	if (prevSegment.isSmoothed())
	{
		throw ModelException(method, "Cannot add discount factor after smooth segment");
	}

	if (curve.getBaseDate() >= discDate)
	{
		throw ModelException(method, "Cannot add discount before curve date");
	}

	if (prevSegment.date >= discDate)
	{
		throw ModelException(method, "Cannot add discount before last date in curve");
	}

	ZC3CurveSegment segment(discDate, prevSegment, discount, curve);

	// too late for this iteration, but will help us on the next iteration
	if (curve.adjAbsolute.size() > 0)
	{
		curve.flatFwdAdjAbsolute(prevSegment, segment);
	}

	curve.data.push_back(segment);
	curve.addCriticalDate(discDate);
    curve.cleanup();
}


DRLIB_END_NAMESPACE
