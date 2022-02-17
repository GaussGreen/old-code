//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : ZC3FuturesStubRule.cpp
//
//   Description : Implementation for futures stub rules.
//
//   Author      : Richard Appleton
//
//   Date        : 18th May 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ZC3FuturesStubRule.hpp"
#include "edginc/ZC3ZeroCurve.hpp"
#include "edginc/ZC3CurveInstrument.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/CompoundBasis.hpp"
#include <algorithm>



DRLIB_BEGIN_NAMESPACE

using namespace std;

static ZC3FuturesStubRuleNone     none;
static ZC3FuturesStubRuleSimple   simple;
static ZC3FuturesStubRuleInterp   interp;
static ZC3FuturesStubRuleFlatFwds flat;

// static bool tmp1 = lessThan(ZC3RateDataSP(), ZC3RateDataSP());
// static bool tmp2 = lessThan(ZC3SwapDataSP(), ZC3SwapDataSP());
// static bool tmp3 = lessThan(ZC3TurnSP(), ZC3TurnSP());


ZC3FuturesStubRuleSP ZC3FuturesStubRule::make(const string& type)
{
    if (CString::equalsIgnoreCase(type, "N")
     || CString::equalsIgnoreCase(type, "None", 4))
    {
        return ZC3FuturesStubRuleSP(new ZC3FuturesStubRuleNone());
    }
    else if (CString::equalsIgnoreCase(type, "S")
     || CString::equalsIgnoreCase(type, "Simple", 6))
    {
        return ZC3FuturesStubRuleSP(new ZC3FuturesStubRuleSimple());
    }
    else if (CString::equalsIgnoreCase(type, "I")
     || CString::equalsIgnoreCase(type, "Interp", 6))
    {
        return ZC3FuturesStubRuleSP(new ZC3FuturesStubRuleInterp());
    }
    else if (CString::equalsIgnoreCase(type, "F")
     || CString::equalsIgnoreCase(type, "Flat", 4))
    {
        return ZC3FuturesStubRuleSP(new ZC3FuturesStubRuleFlatFwds());
    }
    else
    {
        string msg = Format::toString("Unknown futures stub rule %s", type.c_str());
        throw ModelException("ZC3FuturesStubRule::make", msg);
    }
}


bool ZC3FuturesStubRule::rangesOverlap(const DateTime& date) const
{
    return !futuresStartDate.empty() && !date.empty() && (date > futuresStartDate);
}


bool ZC3FuturesStubRule::rangesTouch() const
{
    return !futuresStartDate.empty() && !lastMMDate.empty() && (lastMMDate == futuresStartDate);
}


bool ZC3FuturesStubRule::gapExists() const
{
    return !futuresStartDate.empty() && !lastMMDate.empty() && (lastMMDate < futuresStartDate);
}


// ALIB: szcrates.c#1235 (SZCRatesUsage)
ZC3RateDataArraySP ZC3FuturesStubRule::ratesUsage(
    const ZC3RateDataArray& pRatesData,
    const ZC3TurnArray&     pTurnsData,
    const DateTime*         pFuturesMMDate,
    const ZC3ZeroCurve&     pStubCurve, // base date, adjustments, adjDcc, hasUTurns
    bool                    pMtmAdjustment,
    const IRVolBase*        pVolModelIR)
{
    // initialize temporary data
    stubCurve = &pStubCurve;
    ratesData = ZC3RateDataArraySP(new ZC3RateDataArray(pRatesData));
    turnsData = &pTurnsData;
    futuresMMDate = pFuturesMMDate;
    mtmAdjustment = pMtmAdjustment;
    volModelIR = pVolModelIR;

    futuresStartDate = DateTime();
    futuresEndDate = DateTime();
    lastMMDate = DateTime();
    badMMDate = DateTime();
    deleted = ZC3RateDataArray();

    /*
     * First scan all the rates and see whether or not we have any futures and
     * also whether we have any money market rates that we may wish to exclude.
     */
    ratesInitUsage();

    /*
     * After our first scan we now have the following information:
     * (a) Whether there are any futures.
     * (b) Whether there are any MM rates (starting at base date).
     * (c) The range of relevant dates.
     *
     * Special processing requires the following conditions:
     * (a) Overlapping range of futures and MM rates.
     * (b) Non-zero futures stub type.
     *
     * An error occurs if the futures stub type is such that linear
     * interpolation on the simple rates is required, and the ranges
     * do *not* meet.
     */
    if (!futuresStartDate.empty() && !lastMMDate.empty())
    {
        process();
    }

    removeDeleted();

    return ratesData;
}


// ALIB: szcrates.c#1525 (SZCRatesInitUsage)
void ZC3FuturesStubRule::ratesInitUsage()
{
    for (size_t i = 0 ; i < ratesData->size() ; i++)
    {
        (*ratesData)[i]->ratesInitUsage(stubCurve->getBaseDate(), futuresMMDate,
            futuresStartDate, futuresEndDate, badMMDate, lastMMDate);
    }

    /*
     * Find a really bad MM date - called badMMDate.
     *
     * For example, suppose we need to preserve 6M point, and this falls on the
     * maturity date of the 1st futures contract, then the start date of the
     * 1st futures contract will (a) be the 3M point and (b) be a really bad MM
     * date, since if it is defined then we will get an error later in the code.
     *
     * This is a potentially incredibly inefficient loop to check find a really
     * bad MM date.
     *
     * We start from the start date of a futures contract. If this is also an
     * end date for another futures contract, then a really bad date would be
     * the start date of that contract, and so on. This is the basis of the
     * loop.
     *
     * However we will hardly ever enter this loop, since the initial traps are
     * pretty powerful.
     */
    if (!badMMDate.empty() && (badMMDate > futuresStartDate))
    {
        bool found;

        do
        {
            found = false;

            for (size_t j = 0 ; j < ratesData->size() ; j++)
            {
                if ((*ratesData)[j]->isReallyBadDate(badMMDate))
                {
                    found = true;
                    break;
                }
            }
        } while (found && (badMMDate > futuresStartDate));
    }
}


// ALIB: szcrates.c#1428
void ZC3FuturesStubRule::retainRate()
{
    if (futuresMMDate)
    {
        // Search for rate which defines futuresMMDate - if found retain it
        ZC3RateDataArray::iterator iterator;
        for (iterator = deleted.begin() ; iterator != deleted.end() ; iterator++)
        {
            if ((*iterator)->isCashStartingOn(stubCurve->getBaseDate())
             && (*iterator)->getEndDate() == *futuresMMDate)
            {
                deleted.erase(iterator);  // NB. iterator now invalid
                break;  // no need to search any more
            }
        }

        if (!badMMDate.empty())
        {
            /*
             * Need to remove the really bad date - this overrides any
             * stub methods.
             */
            for (iterator = ratesData->begin() ; iterator != ratesData->end() ; iterator++)
            {
                if ((*iterator)->isCashStartingOn(stubCurve->getBaseDate())
                 && (*iterator)->getEndDate() == badMMDate)
                {
                    deleted.push_back(*iterator);
                    break;  // no need to search any more
                }
             }
        }
    }
}


void ZC3FuturesStubRule::include(const ZC3RateData* rate)
{
    if (rate)
    {
        ZC3RateDataArray::iterator iterator;
        for (iterator = deleted.begin() ; iterator != deleted.end() ; iterator++)
        {
            if ((*iterator).get() == rate)
            {
                deleted.erase(iterator);    // NB. invalidates iterator
                break;
            }
        }
    }
}


// VC6 Release configuration fix for unresolved lessThan<ZC3RateData>
static inline bool lessThan_ZC3RateData(const refCountPtr< ZC3RateData >& first, const refCountPtr< ZC3RateData >& second)
{
    return lessThan(first, second);
}

void ZC3FuturesStubRule::removeDeleted()
{
    /*
     * Any instruments in the 'deleted' array are not to be used for
     * bootstrapping.  Algorithms up to now have just manipulated the 'deleted'
     * array - now these elements must be removed from the actual bootstrap
     * data input.
     */
    for (size_t i = 0 ; i < deleted.size() ; i++)
    {
        ZC3RateDataArray::iterator loc = find(ratesData->begin(), ratesData->end(), deleted[i]);
        if (loc != ratesData->end())
        {
            ratesData->erase(loc);
        }
    }

    // ensure data remains sorted
    sort(ratesData->begin(), ratesData->end(), lessThan_ZC3RateData);
}


void ZC3FuturesStubRuleNone::process()
{
    // do nothing
}


// ALIB: szcrates.c#1345
void ZC3FuturesStubRuleSimple::process()
{
    static const string method = "ZC3FuturesStubRuleSimple::process";

    if (rangesOverlap(lastMMDate))
    {
        futuresStubMarketInterp();
        retainRate();
    }
    else if (!futuresStartDate.empty() && !rangesTouch())
    {
        /*
         * Error case - we have futures and expect linear interpolation on the
         * MM rates, but this turns out not to be possible.
         */
        string msg = "Futures stub type => interpolation on MM rates. "
                     "Insufficient MM rates provided.";
        throw ModelException(method, msg);
    }
}


// ALIB: szcrates.c#1698
void ZC3FuturesStubRuleSimple::futuresStubMarketInterp()
{
    static const string method = "ZC3FuturesStubRuleSimple::futuresStubMarketInterp";

    /*
     * Algorithm:
     *
     * Loop through all rates looking for MM rates.
     * Take note of last before and first after dates.
     * Initially exclude all MM rates after.
     *
     * Interpolate between last before and first after, and include first after
     * rate.
     */
    ZC3RateData* lastBefore = NULL;
    ZC3RateData* firstAfter = NULL;

    for (size_t i = 0 ; i < ratesData->size() ; i++)
    {
        ZC3RateDataSP rate = (*ratesData)[i];

        if (rate->isCashStartingOn(stubCurve->getBaseDate()))
        {
            if (rate->getEndDate() < futuresStartDate)
            {
                if (!lastBefore || rate->getEndDate() > lastBefore->getEndDate())
                {
                    lastBefore = rate.get();
                }
            }
            else
            {
                // Exclude rate if it is before the end of the futures strip
                if (rate->getEndDate() <= futuresEndDate)
                {
                    deleted.push_back(rate);
                }

                if (!firstAfter || rate->getEndDate() < firstAfter->getEndDate())
                {
                    firstAfter = rate.get();
                }
            }
        }
    }

    /*
     * Now try to perform the linear interpolation.
     */
    if (!firstAfter)
    {
        // error (ranges do not meet) - trapped earlier
        throw ModelException(method, "Program bug. Ranges do not meet (firstAfter = NULL)");
    }
    else if (firstAfter->getEndDate() == futuresStartDate)
    {
        // exact match
        include(firstAfter);
    }
    else
    {
        include(firstAfter);

        if (!futuresMMDate || badMMDate.empty() || futuresStartDate != badMMDate)   // test moved out of retainRate as addOneRate not used
        {
            ZC3RateData* newRate = ZC3RateData::interpolate(
                lastBefore,
                *firstAfter,
                stubCurve->getBaseDate(),
                futuresStartDate);

            ratesData->push_back(ZC3RateDataSP(newRate));
        }
    }
}


// ALIB: szcrates.c#1369
void ZC3FuturesStubRuleInterp::process()
{
    if (rangesOverlap(lastMMDate))
    {
        /*
         * Regular interpolation on the curve required.  Many of the MM rates
         * are discarded.
         */
        futuresStubInterp();
        retainRate();
    }
}


// ALIB: szcrates.c#1876
void ZC3FuturesStubRuleInterp::futuresStubInterp()
{
    static const string method = "ZC3FuturesStubRuleInterp::futuresStubInterp";

    /*
     * Algorithm:
     *
     * Loop through all rates looking for MM rates.
     * Take note of first after date.
     * Initially exclude all MM rates after.
     * Go back at end and include the first after date.
     */
    ZC3RateData* firstAfter = NULL;

    for (size_t i = 0 ; i < ratesData->size() ; i++)
    {
        ZC3RateDataSP rate = (*ratesData)[i];

        if (rate->isCashStartingOn(stubCurve->getBaseDate()))
        {
            if (rate->getEndDate() >= futuresStartDate)
            {
                if (rate->getEndDate() <= futuresEndDate)
                {
                    deleted.push_back(rate);
                }

                if (!firstAfter || rate->getEndDate() < firstAfter->getEndDate())
                {
                    firstAfter = rate.get();
                }
            }
        }
    }

    /*
     * Now reset the firstAfter rate
     */
    include(firstAfter);
}


// ALIB: szcrates.c#1384
void ZC3FuturesStubRuleFlatFwds::process()
{
    /*
     * We use flat forwards using the futures MM date.
     *
     * We discard the MM rates at or after the stub (except for the futures
     * MM date).
     *
     * If the futures MM date is not provided then this will fail.
     */
    if (rangesOverlap(lastMMDate) || gapExists())
    {
        futuresStubFlatFwds();
    }
}


// ALIB: szcrates.c#1943 (includes changes made in ALIB 13.6.1 for ALIB-85)
void ZC3FuturesStubRuleFlatFwds::futuresStubFlatFwds()
{
    static const string method = "ZC3FuturesStubRuleFlatFwds::futuresStubFlatFwds";

    /*
     * We use flat forwards using the futures MM date.
     *
     * We discard the MM rates at or after the stub (except for the futures MM
     * date).
     *
     * If the futures MM date is not provided then this will fail.
     */
    const DateTime& baseDate = stubCurve->getBaseDate();
    if (!futuresMMDate || *futuresMMDate < baseDate)
    {
        throw ModelException(method, "Futures MM date is not provided, or is before base date");
    }

    ZC3RateData* found = NULL;
    ZC3RateData* expired = NULL;
    ZC3RateData* firstFuture = NULL;

    for (size_t i = 0 ; i < ratesData->size() ; i++)
    {
        ZC3RateDataSP rateData = (*ratesData)[i];

        if (rateData->isCashStartingOn(baseDate))
        {
            if (lastMMDate > futuresStartDate)
            {
                if ((rateData->getEndDate() >= futuresStartDate
                  && rateData->getEndDate() <= futuresEndDate)
                  || (rateData->getEndDate() < *futuresMMDate))
                {
                    deleted.push_back(rateData);
                }

                if (rateData->getEndDate() == *futuresMMDate)
                {
                    deleted.push_back(rateData);

                    if (found)
                    {
                        string msg = "Duplicate rates provided for " + futuresMMDate->toString();
                        throw ModelException(method, msg);
                    }
                    else
                    {
                        found = rateData.get();
                    }
                }
            }
            else    // gap exists between last MM and first future stub
            {
                if (rateData->getEndDate() == *futuresMMDate)
                {
                    if (found)
                    {
                        string msg = "Duplicate rates provided for " + futuresMMDate->toString();
                        throw ModelException(method, msg);
                    }
                    else
                    {
                        found = rateData.get();
                    }
                }
            }
        }
        else if (typeid(*rateData) == typeid(ZC3FuturesData))
        {
            if (!firstFuture || rateData->getStartDate() < firstFuture->getStartDate())
            {
                /*
                 * With contiguous futures there is the (rare) situation where
                 * the 3M date is after the first future. We handle this here.
                 * We search for a (single) future which expires before the
                 * futuresMMDate.
                 */
                if (rateData->getEndDate() < *futuresMMDate)
                {
                    if (expired)
                    {
                        string msg = "More than 1 future ends before the futuresMMDate";
                        throw ModelException(method, msg);
                    }
                    else
                    {
                        expired = rateData.get();
                    }
                }
                else
                {
                    firstFuture = rateData.get();
                }
            }
        }
    }

    if (!found)
    {
        string msg = "No rate provided for " + futuresMMDate->toString();
        throw ModelException(method, msg);
    }

    double mmRate = found->getRate();
    const DayCountConvention& mmDcc = found->getDcc();

    /*
     * Now perform the flat forwards interpolation.
     * We might have some turns in the stub curve, and some more in the list of
     * rates. So far we have only looked for those in the instruments.  So we
     * now add the adjustments from the stub curve (which are all relative
     * adjustments).
     */
    if (turnsData->size() > 0 && stubCurve->adjustments.size() > 0)
    {
        string msg = "Turns in both the list of instruments and stub curve is "
                     "currently not supported";
        throw ModelException(method, msg);
    }

    ZC3TurnArray turns = *turnsData;

    try
    {
        for (int k = 0 ; k < stubCurve->adjustments.size() ; k++)
        {
            const ZC3Adjustment& adjustment = stubCurve->adjustments[k];
            double discount = exp(adjustment.getContRate());
            double rate = RateConversion::discountToRate(discount,
                adjustment.getStartDate(),
                adjustment.getEndDate(),
                stubCurve->dcc.get(),
                CompoundBasis::SIMPLE);

            ZC3Turn* turn = new ZC3RelativeTurn(adjustment.getStartDate(),
                adjustment.getEndDate(), rate, *stubCurve->dcc);

            turns.push_back(ZC3TurnSP(turn));
        }

        double rate;
        if (rangesOverlap(*futuresMMDate))
        {
            rate = futuresStubRateFlatForwards(
                *futuresMMDate,
                mmRate,
                mmDcc,
                *firstFuture,
                turns);

            if (expired)
            {
                if (expired->getEndDate() != firstFuture->getStartDate())
                {
                    string msg = Format::toString(
                        "Future before futuresMMDate is non-contiguous [%s, not %s]",
                        expired->getEndDate().toString().c_str(),
                        firstFuture->getStartDate().toString().c_str());
                    throw ModelException(method, msg);
                }

                rate = futuresStubRateFlatForwards(
                    firstFuture->getStartDate(),
                    rate,
                    mmDcc,
                    *expired,
                    turns);

                futuresStartDate = expired->getStartDate();
            }
        }
        else
        {
            rate = futuresStubRateFlatForwards(
                *futuresMMDate,
                mmRate,
                mmDcc,
                *firstFuture,
                turns);
        }

        ZC3RateData* stub = new ZC3CashData(baseDate, futuresStartDate, rate, mmDcc, 0.0);
        ratesData->push_back(ZC3RateDataSP(stub));
    }
    catch(exception& e)
    {
        throw e;
    }
}


// ALIB: [13.6.1] futstub.c#145&275
double ZC3FuturesStubRuleFlatFwds::futuresStubRateFlatForwards(
    const DateTime&           date,
    double                    mmRate,
    const DayCountConvention& mmDcc,
    const ZC3RateData&        future,
    const ZC3TurnArray&       turns
    )
{
    static const string method = "ZC3FuturesStubRuleFlatFwds::futuresStubRateFlatForwards";

    const DateTime& valueDate = stubCurve->getBaseDate();
    double discount;

    if (rangesOverlap(*futuresMMDate))
    {
        if (date > future.getEndDate())
        {
            string msg = Format::toString(
                "FuturesMMDate (%s) > first future end date (%s)",
                date.toString().c_str(),
                future.getEndDate().toString().c_str());
            throw ModelException(method, msg);
        }

        double futDisc = RateConversion::rateToDiscount(future.getRate(),
            future.getStartDate(), future.getEndDate(), &future.getDcc(), CompoundBasis::SIMPLE);

        double stubDisc = RateConversion::rateToDiscount
            (mmRate, valueDate, date, &mmDcc, CompoundBasis::SIMPLE);

        double futStartYF;
        double cumStartDisc = turnEffectInPeriod(
            future.getStartDate(),
            date,
            turns,
            mmDcc,
            &futStartYF,
            NULL);

        double futEndYF;
        double cumEndDisc = turnEffectInPeriod(
            date,
            future.getEndDate(),
            turns,
            mmDcc,
            &futEndYF,
            NULL);

        double futContRate = -log(futDisc / (cumStartDisc * cumEndDisc));
        futContRate /= (futEndYF + futStartYF);

        double futStartDisc = exp(-futContRate * futStartYF) * cumStartDisc;
        discount = stubDisc / futStartDisc;
    }
    else    // gap exists between last MM and first future stub
    {
        if (date > future.getStartDate())
        {
            string msg = Format::toString(
                "futuresMMDate (%s) > first future start date (%s)",
                date.toString().c_str(),
                future.getStartDate().toString().c_str());
            throw ModelException(method, msg);
        }

        double mmDisc = RateConversion::rateToDiscount
            (mmRate, valueDate, date, &mmDcc, CompoundBasis::SIMPLE);

        int cashTurnDays;
        double cashTurnDisc = turnEffectInPeriod(
            valueDate,
            date,
            turns,
            mmDcc,
            NULL,
            &cashTurnDays);

        int gapTurnDays;
        double gapTurnDisc = turnEffectInPeriod(
            date,
            future.getStartDate(),
            turns,
            mmDcc,
            NULL,
            &gapTurnDays);

        double overnight = pow(mmDisc / cashTurnDisc, 1.0 / (date.daysDiff(valueDate) - cashTurnDays));
        discount = pow(overnight, (future.getStartDate().daysDiff(*futuresMMDate) - gapTurnDays)) * mmDisc * gapTurnDisc;
    }

    return RateConversion::discountToRate
        (discount, valueDate, future.getStartDate(), &mmDcc, CompoundBasis::SIMPLE);
}


// ALIB: [13.6.1] futstub.c#381
double ZC3FuturesStubRuleFlatFwds::turnEffectInPeriod(
    const DateTime&           periodStartDate,
    const DateTime&           periodEndDate,
    const ZC3TurnArray&       turns,
    const DayCountConvention& dcc,
    double*                   turnYF,
    int*                      turnDays
    ) const
{
    static const string method = "ZC3FuturesStubRuleFlatFwds::turnEffectInPeriod";

    if (periodStartDate > periodEndDate)
    {
        string msg = Format::toString(
            "start date (%s) > end date (%s)",
            periodStartDate.toString().c_str(),
            periodEndDate.toString().c_str());
        throw ModelException(method, msg);
    }

    // initialize outputs
    double discount = 1.0;

    if (turnDays)
    {
        *turnDays = 0;
    }

    if (turnYF)
    {
        *turnYF = dcc.years(periodStartDate, periodEndDate);
    }

    // iterate over turns
    for (size_t i = 0 ; i < turns.size() ; i++)
    {
        if (turns[i]->getStartDate() < periodEndDate
         && turns[i]->getEndDate() > periodStartDate)
        {
            DateTime startDate;
            DateTime endDate;

            if (turns[i]->getStartDate() > periodStartDate)
            {
                startDate = turns[i]->getStartDate();
            }
            else
            {
                startDate = periodStartDate;
            }

            if (turns[i]->getEndDate() < periodEndDate)
            {
                endDate = turns[i]->getEndDate();
            }
            else
            {
                endDate = periodEndDate;
            }

            discount *= turns[i]->turnEffectInPeriod(turnYF, turnDays);
        }
    }

    return discount;
}


DRLIB_END_NAMESPACE
