//----------------------------------------------------------------------------
//
//   Group       : QR Core Analytics Team
//
//   Filename    : LMParabCurve.cpp
//
//   Description : ALIB style local market parabolic curve
//
//   Author      : Andrew Greene
//
//   Date        : 07 April 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"

#include "edginc/LMParabCurve.hpp"

#include "edginc/Actual365F.hpp"
#include "edginc/Business252.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/Maths.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/StubSimple.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/ZC3Interpolation.hpp"

#include <algorithm>

#define LM_RELATIVE_ERROR 1e-08

#define LM_INST_RATE 0
#define LM_INST_RATE_SLOPE 1
#define LM_OVERNIGHT_RATE 2

DRLIB_BEGIN_NAMESPACE

LMParabCurve::LMParabCurve
(
    const DateTime&              baseDate,         // (I) base date for curve
    const vector<FwdRateIvl>*    frl,              // (I)
    const vector<YearEndEffect>* yel,              // (I)
    HolidayConstSP               noAccrue,         // (I) dates with no accrual
    const vector<DateTime>*      flatRateDates,    // (I) the rate will be constant in between flatRateDates,
                                                   //     last date must be <= flatParabBdryDate
    const DateTime&              flatParabBdryDate // (I) Date at which interpolations transitions
                                                   //     from flat to parabolic
)
{
    static const string method = "LMParabCurve::LMParabCurve";

    try
    {
        genFromFras(baseDate,
                    frl           != NULL ? *frl           : vector<FwdRateIvl>   (),
                    yel           != NULL ? *yel           : vector<YearEndEffect>(),
                    noAccrue,
                    flatRateDates != NULL ? *flatRateDates : vector<DateTime>     (),
                    flatParabBdryDate);
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


LMParabCurve::LMParabCurve
(
    const DateTime&           baseDate,          // (I) base date for zero curve.  Discounts from this date
    const vector<DateTime>*   mmMatDates,        // (I) vector of MM maturity dates
    const vector<double>*     mmRates,           // (I) vector of MM rates
    const vector<DateTime>*   fraStartDates,     // (I) vector of Fra start dates
    const vector<DateTime>*   fraMatDates,       // (I) vector of Fra maturity dates
    const vector<double>*     fraRates,          // (I) vector of Fra rates
    const vector<DateTime>*   yeStartDates,      // (I) vector of Year End adjustment start dates
    const vector<DateTime>*   yeEndDates,        // (I) vector of Year End adjustment end dates
    const vector<double>*     yeRates,           // (I) vector of Year End adjustment rates
    const vector<bool>*       yeIsSpread,        // (I) if TRUE the Year End adjustment is used as a
                                                 //     spread, otherwise in absolute level
    const vector<DateTime>*   relMMDates,        // (I) vector of relative MM maturity dates
    const vector<double>*     relMMSpreads,      // (I) vector of spreads for relative MM points
    const DayCountConvention& mmDayCount,        // (I) Common Day Count Convention for MM, Fras and YE
    int                       rateBasis,         // (I) Common Rate Basis for MM, Fras and YE
    const vector<DateTime>*   flatRateDates,     // (I) the rate will be constant in between flatRateDates,
                                                 //     last date must be <= flatParabBdryDate
    const DateTime&           flatParabBdryDate, // (I) Date at which interpolation transitions
                                                 //     from flat to parabolic
    HolidayConstSP            noAccrue,          // (I) No accrual dates 
    const string&             interpType         // (I) Smooth forwards or linear
)
{
    static const string method = "LMParabCurve::LMParabCurve";

    try
    {
         genFromMM(baseDate,
                   mmMatDates        != NULL ? *mmMatDates    : vector<DateTime>(),
                   mmRates           != NULL ? *mmRates       : vector<double  >(),
                   fraStartDates     != NULL ? *fraStartDates : vector<DateTime>(),
                   fraMatDates       != NULL ? *fraMatDates   : vector<DateTime>(),
                   fraRates          != NULL ? *fraRates      : vector<double  >(),
                   yeStartDates      != NULL ? *yeStartDates  : vector<DateTime>(),
                   yeEndDates        != NULL ? *yeEndDates    : vector<DateTime>(),
                   yeRates           != NULL ? *yeRates       : vector<double  >(),
                   yeIsSpread        != NULL ? *yeIsSpread    : vector<bool    >(),
                   relMMDates        != NULL ? *relMMDates    : vector<DateTime>(),
                   relMMSpreads      != NULL ? *relMMSpreads  : vector<double  >(),
                   mmDayCount,
                   rateBasis,
                   flatRateDates     != NULL ? *flatRateDates : vector<DateTime>(),
                   flatParabBdryDate,
                   noAccrue,
                   interpType        );
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


LMParabCurve::LMParabCurve
(
    const DateTime&              baseDate,          // (I) base date for zero curve.  Discounts from this date
    const vector<DateTime>*      mmMatDates,        // (I) vector of MM maturity dates
    const vector<double>*        mmRates,           // (I) vector of MM rates
    const vector<DateTime>*      fraStartDates,     // (I) vector of Fra start dates
    const vector<DateTime>*      fraMatDates,       // (I) vector of Fra maturity dates
    const vector<double>*        fraRates,          // (I) vector of Fra rates
    const vector<DateTime>*      yeStartDates,      // (I) vector of Year End adjustment start dates
    const vector<DateTime>*      yeEndDates,        // (I) vector of Year End adjustment end dates
    const vector<double>*        yeRates,           // (I) vector of Year End adjustment rates
    const vector<bool>*          yeIsSpread,        // (I) if TRUE the Year End adjustment is used as a
                                                    //     spread, otherwise in absolute level
    const vector<DateTime>*      relMMDates,        // (I) vector of relative MM maturity dates
    const vector<double>*        relMMSpreads,      // (I) vector of spreads for relative MM points
    const DayCountConvention&    mmDayCount,        // (I) common day count convention for MM, Fras and YE
    int                          rateBasis,         // (I) common rate basis for MM, Fras and YE
    const vector<DateTime>*      flatRateDates,     // (I) the rate will be constant in between flatRateDates,
    const vector<CashFlowArray>* cfl,               // (I) vector of cashflow arrays
    const vector<double>*        cflPrice,          // (I) vector of dirty prices for cashflow arrays
    const DateTime&              flatParabBdryDate, // (I) Date at which interpolation transitions
                                                    //     from flat to parabolic
    HolidayConstSP               noAccrue,          // (I) No accrual dates 
    const string&                interpType         // (I) Smooth forwards or linear
)
{
    static const string method = "LMParabCurve::LMParabCurve";

    try
    {
        genFromMM(baseDate,
                  mmMatDates        != NULL ? *mmMatDates    : vector<DateTime>(),
                  mmRates           != NULL ? *mmRates       : vector<double  >(),
                  fraStartDates     != NULL ? *fraStartDates : vector<DateTime>(),
                  fraMatDates       != NULL ? *fraMatDates   : vector<DateTime>(),
                  fraRates          != NULL ? *fraRates      : vector<double  >(),
                  yeStartDates      != NULL ? *yeStartDates  : vector<DateTime>(),
                  yeEndDates        != NULL ? *yeEndDates    : vector<DateTime>(),
                  yeRates           != NULL ? *yeRates       : vector<double  >(),
                  yeIsSpread        != NULL ? *yeIsSpread    : vector<bool    >(),
                  relMMDates        != NULL ? *relMMDates    : vector<DateTime>(),
                  relMMSpreads      != NULL ? *relMMSpreads  : vector<double  >(),
                  mmDayCount,
                  rateBasis,
                  flatRateDates     != NULL ? *flatRateDates : vector<DateTime>(),
                  flatParabBdryDate,
                  noAccrue,
                  interpType        );

        if (interpType == ZC3ZeroInterpolation::SMOOTH_FORWARDS)
        {
            addCF(cfl      != NULL ? *cfl      : vector<CashFlowArray>(),
                  cflPrice != NULL ? *cflPrice : vector<double       >(),
                  flatParabBdryDate);
        }
        else if (interpType == ZC3ZeroInterpolation::LINEAR)
        {
            linZerosAddCF(cfl      != NULL ? *cfl      : vector<CashFlowArray>(),
                          cflPrice != NULL ? *cflPrice : vector<double       >() );
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


LMParabCurveSP LMParabCurve::genFromBonds
(
    const DateTime&              baseDate,          // (I) base date for zero curve.  Discounts from this date
    const vector<DateTime>*      mmMatDates,        // (I) vector of MM maturity dates
    const vector<double>*        mmRates,           // (I) vector of MM rates
    const vector<DateTime>*      fraStartDates,     // (I) vector of Fra start dates
    const vector<DateTime>*      fraMatDates,       // (I) vector of Fra maturity dates
    const vector<double>*        fraRates,          // (I) vector of Fra rates
    const vector<DateTime>*      yeStartDates,      // (I) vector of Year End adjustment start dates
    const vector<DateTime>*      yeEndDates,        // (I) vector of Year End adjustment end dates
    const vector<double>*        yeRates,           // (I) vector of Year End adjustment rates
    const vector<bool>*          yeIsSpread,        // (I) if TRUE the Year End adjustment is used as a
                                                    //     spread, otherwise in absolute level
    const vector<DateTime>*      relMMDates,        // (I) vector of relative MM maturity dates
    const vector<double>*        relMMSpreads,      // (I) vector of spreads for relative MM points
    const DayCountConvention&    mmDayCount,        // (I) common day count convention for MM, Fras and YE
    int                          rateBasis,         // (I) common rate basis for MM, Fras and YE
    const vector<DateTime>*      flatRateDates,     // (I) the rate will be constant in between flatRateDates,
    const vector<DateTime>*      bondMat,           // (I) bond maturity dates
    const vector<double>*        bondRate,          // (I) bond rates
    int                          bondFreq,          // (I) bond number of compounding periods per year
    const DayCountConvention&    bondDayCount,      // (I) bond day count convention
    const BadDayConvention&      accBadDayConv,     // (I) bond bad day convention for interest rate accrual
    const BadDayConvention&      payBadDayConv,     // (I) bond bad day convention for payment dates
    const vector<double>*        bondPrice,         // (I) bond dirty prices
    const DateTime&              flatParabBdryDate, // (I) date at which interpolation transitions
                                                    //     from flat to parabolic
    HolidayConstSP               noAccrue,          // (I) no accrual dates 
    const string&                interpType         // (I) smooth forwards or linear
)
{
    static const string method = "LMParabCurve::genFromBonds";

    try
    {
        size_t                matSize;
        size_t                rateSize;
        size_t                idxCF;
        vector<CashFlowArray> cfl;

        matSize  = bondMat  != NULL ? bondMat ->size() : 0;
        rateSize = bondRate != NULL ? bondRate->size() : 0;

        if (rateSize != matSize)
        {
            throw ModelException("bondMat and bondRate size mismatch.");
        }

        for (idxCF = 0; idxCF < matSize; ++idxCF)
        {
            MaturityPeriod maturityPeriod(bondFreq);
            cfl.push_back(SwapTool::cashflows(baseDate,
                                              (*bondMat)[idxCF],
                                              StubSimple(),
                                              false,
                                              false,
                                              accBadDayConv,
                                              payBadDayConv,
                                              *noAccrue.get(),
                                              false,
                                              false,
                                              true,
                                              (*bondRate)[idxCF],
                                              maturityPeriod,
                                              bondDayCount             ));
        }

        return LMParabCurveSP(new LMParabCurve(baseDate,
                                               mmMatDates,
                                               mmRates,
                                               fraStartDates,
                                               fraMatDates,
                                               fraRates,
                                               yeStartDates,
                                               yeEndDates,
                                               yeRates,
                                               yeIsSpread,
                                               relMMDates,
                                               relMMSpreads,
                                               mmDayCount,
                                               rateBasis,
                                               flatRateDates,
                                               &cfl, 
                                               bondPrice, 
                                               flatParabBdryDate,
                                               noAccrue,
                                               interpType));
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


LMParabCurveSP LMParabCurve::genFromSwaps
(
    const DateTime&              baseDate,          // (I) base date for zero curve.  Discounts from this date
    const vector<DateTime>*      mmMatDates,        // (I) vector of MM maturity dates
    const vector<double>*        mmRates,           // (I) vector of MM rates
    const vector<DateTime>*      fraStartDates,     // (I) vector of Fra start dates
    const vector<DateTime>*      fraMatDates,       // (I) vector of Fra maturity dates
    const vector<double>*        fraRates,          // (I) vector of Fra rates
    const vector<DateTime>*      yeStartDates,      // (I) vector of Year End adjustment start dates
    const vector<DateTime>*      yeEndDates,        // (I) vector of Year End adjustment end dates
    const vector<double>*        yeRates,           // (I) vector of Year End adjustment rates
    const vector<bool>*          yeIsSpread,        // (I) if TRUE the Year End adjustment is used as a
                                                    //     spread, otherwise in absolute level
    const vector<DateTime>*      relMMDates,        // (I) vector of relative MM maturity dates
    const vector<double>*        relMMSpreads,      // (I) vector of spreads for relative MM points
    const DayCountConvention&    mmDayCount,        // (I) common day count convention for MM, Fras and YE
    int                          rateBasis,         // (I) common rate basis for MM, Fras and YE
    const vector<DateTime>*      flatRateDates,     // (I) the rate will be constant in between flatRateDates,
    const vector<DateTime>*      swapMat,           // (I) bond maturity dates
    const vector<double>*        swapRate,          // (I) bond rates
    int                          swapFreq,          // (I) bond number of compounding periods per year
    const DayCountConvention&    swapDayCount,      // (I) bond day count convention
    const BadDayConvention&      accBadDayConv,     // (I) bond bad day convention for interest rate accrual
    const BadDayConvention&      payBadDayConv,     // (I) bond bad day convention for payment dates
    const DateTime&              flatParabBdryDate, // (I) date at which interpolation transitions
                                                    //     from flat to parabolic
    HolidayConstSP               noAccrue,          // (I) no accrual dates
    const string&                interpType         // (I) smooth forwards or linear
)
{
    static const string method = "LMParabCurve::genFromSwaps";

    try
    {
        size_t matSize = swapMat != NULL ? swapMat ->size() : 0;
        vector<double> unitPrice(matSize, 1.);

        return genFromBonds(baseDate,
                            mmMatDates,
                            mmRates,
                            fraStartDates,
                            fraMatDates,
                            fraRates,
                            yeStartDates,
                            yeEndDates,
                            yeRates,
                            yeIsSpread,
                            relMMDates,
                            relMMSpreads,
                            mmDayCount,
                            rateBasis,
                            flatRateDates,
                            swapMat,
                            swapRate,
                            swapFreq,
                            swapDayCount,
                            accBadDayConv,
                            payBadDayConv,
                            &unitPrice,
                            flatParabBdryDate,
                            noAccrue,
                            interpType        );
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


double LMParabCurve::disc
(
    const DateTime& startDate, // (I) start date for discount factor
    const DateTime& endDate    // (I) end date for discount factor
)
const
{
    static const string method = "LMParabCurve::disc";

    try
    {
        double discStart;
        double discEnd;

        discStart = interpDisc(startDate);
        discEnd   = interpDisc(endDate);

        return discEnd/discStart;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


double LMParabCurve::mmRate
(
    const DateTime&           startDate,    // (I) start date
    const DateTime&           endDate,      // (I) end date
    const DayCountConvention& dayCountConv, // (I) day count convention
    int                       basis         // (I) rate basis
)
const
{
    static const string method = "LMParabCurve::mmRate";

    try
    {

        double df = disc(startDate,
                         endDate   );

        return RateConversion::discountToRate(df,
                                              startDate,
                                              endDate,
                                              &dayCountConv,
                                              basis         );
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


ZC3ZeroCurveSP LMParabCurve::toZeroCurve()
{
    bool   areWeekendsHolidays;
    int    daysYear;
    size_t idx;

    ZC3ZeroCurveSP zc(new ZC3ZeroCurve);
    zc->valueDate = zc->baseDate = baseDate;

    noAccrue->toALIB(areWeekendsHolidays);
    if (!areWeekendsHolidays)
    {
        zc->dcc  = DayCountConventionConstSP(new Actual365F);
        daysYear = 365;
    }
    else
    {
        zc->dcc  = DayCountConventionConstSP(new Business252(noAccrue));
        daysYear = 252;
    }

    // Set the first point at the base date with a discount factor of 1.
    zc->data.clear();
    zc->criticalDates.push_back(baseDate);
    zc->addDiscountFactor(baseDate,
                          1.,
                          *ZC3ZeroInterpolation::make("Linear"));

    if (!ratePts.empty())
    {
        zc->data[0].rate = RateConversion::discountToRate(ratePts[0].disc,
                                                          baseDate,
                                                          ratePts[0].date,
                                                          zc->dcc.get(),
                                                          zc->basis       );

        zc->data[0].averageRate(*zc, true);
    }

    for (idx = 0; idx < ratePts.size(); ++idx)
    {
        zc->criticalDates.push_back(ratePts[idx].date);

        zc->addDiscountFactor(ratePts[idx].date,
                              ratePts[idx].disc,
                              *ZC3ZeroInterpolation::make("Smooth"));

        zc->data.back().averageRate(*zc, true);

        zc->data.back().a0 = ratePts[idx].coeff[0]*daysYear;
        zc->data.back().a1 = ratePts[idx].coeff[1]*daysYear*daysYear;
        zc->data.back().a2 = ratePts[idx].coeff[2]*daysYear*daysYear*daysYear;
    }
 
    return zc;
}


void LMParabCurve::addCF
(
    const vector<CashFlowArray>& cf,               // (I) vector of cashflow arrays to be added
    const vector<double>&        cfPrice,          // (I) array of prices
    const DateTime&              flatParabBdryDate // (I) date at which interpolation transitions
                                                   //     from flat to parabolic
)
{
    static const string method = "LMParabCurve::addCF";

    try
    {
        size_t idxCF;

        bool   keepGoingFlat;
        double disc;
        double lastInstRate;

        if (cf.empty())
        {
            return;
        }
        else if (cfPrice.size() != cf.size())
        {
            throw ModelException("cf and cfPrice size mismatch.");
        }

        keepGoingFlat = !(ratePts.empty() && baseDate >= flatParabBdryDate ||
                          ratePts.back().date >= flatParabBdryDate);

        // Solve flat until flatParabBdryDate
        for (idxCF = 0; keepGoingFlat && idxCF < cf.size(); ++idxCF)
        {
            addNextCFFlat(cf[idxCF],
                          cfPrice[idxCF]);

            if (ratePts.back().date > flatParabBdryDate)
            {
                throw ModelException("flatParabBdryDate must correspond to a CF maturity.");
            }

            if (ratePts.back().date == flatParabBdryDate)
            {
                keepGoingFlat = false;
            }
        }

        // The rest is parabolic
        if (idxCF == 0 && ratePts.empty())
        {
            // For a single cashflow the curve must be flat
            if (cf.size() == 1)
            {
                addNextCFFlat(cf[0], cfPrice[0]);
                return;
            }
            // First interval must be linear for a fully parabolic curve
            else
            {
                addFirstCFLinear(cf     [0],
                                 cfPrice[0],
                                 cf     [1],
                                 cfPrice[1],
                                 lastInstRate);
            }

            ++idxCF;
        }
        else
        {
            if (idxCF == 0 && ratePts.size() > 1 && ratePts.back().date > flatParabBdryDate)
            {
                // Correct the boundary rate for MM
                // Not best solution when YE effects are in the last MM interval,
                // better to reintroduce YE effects after full curve generation
                // Remove point and reintroduce it
                CashFlowArray CFMM(1, CashFlow(ratePts.back().date, 1.));

                disc         = ratePts.back().disc;
                lastInstRate = ratePts.back().coeff[0];

                ratePts.pop_back();

                addNextCFParab(CFMM,
                               disc,
                               cf     [0],
                               cfPrice[0],
                               lastInstRate);
            }
            else
            {
                lastInstRate = instRate(ratePts.back().date,
                                        LM_INST_RATE        );
            }
        }

        for (; idxCF < cf.size() - 1; ++idxCF)
        {
            addNextCFParab(cf     [idxCF],
                           cfPrice[idxCF],
                           cf     [idxCF + 1],
                           cfPrice[idxCF + 1],
                           lastInstRate       );
        }

        // If last CF is parabolic
        if (idxCF == cf.size() - 1)
        {
            addLastCFParab(cf     [idxCF],
                           cfPrice[idxCF],
                           lastInstRate   );
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::addNextCFFlat
(
    const CashFlowArray& cashFlow, // (I) cashflow to be added
    double               cfPrice   // (I) PV of the cashflow to curve base date
)
{
    static const string method = "LMParabCurve::addNextCFFlat";

    try
    {
        DateTime prevDate;
        double   prevDisc;
        double   rate;
        double   disc;
        double   drate;

        double   sumFlows   = 0.;
        double   sumFlowsDt = 0.;
        double   flowsDur;

        double   startPV = 0.;

        int      idxCF;
        int      firstCF;

        double   price;
        double   grad;
        double   hess;
        double   error;

        if (ratePts.empty())
        {
            prevDate       = baseDate;
            prevDisc       = 1.;
            areCoeffsValid = true;
        }
        else
        {
            prevDate = ratePts.back().date;
            prevDisc = ratePts.back().disc;
        }

        // First cash flow after base date
        for (idxCF = 0; idxCF < cashFlow.size(); ++idxCF)
        {
            if (cashFlow[idxCF].date > baseDate)
            {
                break;
            }
        }

        if (idxCF == int(ratePts.size()))
        {
            throw ModelException("Cash flow expired before base date of curve.");
        }

        for (; idxCF < cashFlow.size(); ++idxCF)
        {
            if (cashFlow[idxCF].date > prevDate)
            {
                break;
            }

            disc = this->disc(baseDate,
                              ratePts[idxCF].date);

            startPV += cashFlow[idxCF].amount*disc;
        }

        if (idxCF == cashFlow.size())
        {
            throw ModelException("Cash flow expired before last date of curve.");
        }

        // Value cashflow according to initial rate guess
        firstCF = idxCF;
        vector<int> dateIvl(cashFlow.size());

        for (idxCF = firstCF; idxCF < cashFlow.size(); ++idxCF)
        {
            dateIvl[idxCF] = noAccrue->businessDaysDiff(prevDate,
                                                        cashFlow[idxCF].date);

            sumFlows   += cashFlow[idxCF].amount;
            sumFlowsDt += cashFlow[idxCF].amount*dateIvl[idxCF];
        }

        rate = 0.;
        flowsDur = sumFlowsDt/sumFlows;
        // Exact date for a single cashflow
        drate = -log((cfPrice - startPV)/sumFlows/prevDisc)/flowsDur;
        error = startPV - cfPrice;

        while (fabs(error/cfPrice)  > LM_RELATIVE_ERROR ||
               fabs(drate*flowsDur) > LM_RELATIVE_ERROR)
        {
            price = startPV;
            grad  = 0.;
            hess  = 0.;
            rate += drate;
            for (idxCF = firstCF; idxCF < cashFlow.size(); ++idxCF)
            {
                disc = prevDisc*exp(-rate*dateIvl[idxCF]);

                price += cashFlow[idxCF].amount*disc;
                grad  -= cashFlow[idxCF].amount*disc*dateIvl[idxCF];
                hess  += cashFlow[idxCF].amount*disc*dateIvl[idxCF]*dateIvl[idxCF];
            }

            error = price - cfPrice;
            drate = -grad/hess + Maths::sign(grad/hess)*sqrt(grad/hess*grad/hess - 2.*error/hess);
        }

        RatePt rp = {cashFlow.back().date,
                     disc,
                     rate,
                     0.,
                     0.                   };

        ratePts.push_back(rp);
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::addNextCFParab
(
    const CashFlowArray& currCF,      // (I) cashflows to be added
    double               currPrice,   // (I) price of the current cashflows
    const CashFlowArray& nextCF,      // (I) next cashflows to be added
    double               nextPrice,   // (I) price of the next cashflows
    double&              lastInstRate // (I/O) instantaneous forward rate at end of period
)
{
    static const string method = "LMParabCurve::addNextCFParab";

    try
    {
        double   instBdryRate;
        double   prevDisc;
        double   rate;
        double   startPV = 0.;
        double   price;
        double   grad;
        double   hess;
        double   ivlRatio;
        double   disc;
        double   error;
        double   drate;

        double   sumFlows   = 0.;
        double   sumFlowsDt = 0.;
        double   flowsDur;

        int      idxCve;
        int      idxCF;
        int      firstCF;
        int      dateIvl;
        int      currIvl;

        DateTime prevDate;

        idxCve = ratePts.size();

        // Add the next two cashflows with flat parabolic
        // and later get rid of the second one
        addNextCFFlat(currCF, currPrice);
        addNextCFFlat(nextCF, nextPrice);

        // Get the instantaneous rate at the maturity of currCF
        instBdryRate = this->instBdryRate(currCF, nextCF);

        // Variables to price currCF
        if (idxCve == 0)
        {
            prevDate = baseDate;
            prevDisc = 1.;
        }
        else
        {
            prevDate = ratePts[idxCve - 1].date;
            prevDisc = ratePts[idxCve - 1].disc;
        }

        // First cashflow after base date
        for (idxCF = 0; idxCF < currCF.size(); ++idxCF)
        {
            if (currCF[idxCF].date > baseDate)
            {
                break;
            }
        }

        for (; idxCF < currCF.size(); ++idxCF)
        {
            if (currCF[idxCF].date > prevDate)
            {
                break;
            }

            disc = this->disc(baseDate, currCF[idxCF].date);

            startPV += currCF[idxCF].amount*disc;
        }

        // Value cashflow according to initial rate guess
        firstCF = idxCF;
        vector<double> rateVar(currCF.size());
        vector<double> rateCte(currCF.size());

        currIvl = noAccrue->businessDaysDiff(prevDate, ratePts[idxCve].date);

        for (idxCF = firstCF; idxCF < currCF.size(); ++idxCF)
        {
            dateIvl = noAccrue->businessDaysDiff(prevDate, currCF[idxCF].date);

            ivlRatio = double(dateIvl)/currIvl;

            rateVar[idxCF] = -dateIvl*ivlRatio*(3. - 2.*ivlRatio);
            rateCte[idxCF] = -lastInstRate*dateIvl*(1. - ivlRatio*(2. - ivlRatio))
                             + instBdryRate*dateIvl*ivlRatio*(1. - ivlRatio);
            
            sumFlows   += currCF[idxCF].amount;
            sumFlowsDt += currCF[idxCF].amount*dateIvl;
        }
        rate = 0.;
        flowsDur = sumFlowsDt/sumFlows;
        // Exact date for a single cashflow
        drate = -log((currPrice - startPV)/sumFlows/prevDisc)/flowsDur;
        error = startPV - currPrice;

        // Newton Raphson routine (function is increasing and convex
        // so convergence is guaranteed)

        while (fabs(error/currPrice) > LM_RELATIVE_ERROR || fabs(drate*flowsDur) > LM_RELATIVE_ERROR)
        {
            price = startPV;
            grad  = 0.;
            hess  = 0.;
            rate += drate;

            for (idxCF = firstCF; idxCF < currCF.size(); ++idxCF)
            {
                disc = prevDisc*exp(rateCte[idxCF] + rateVar[idxCF]*rate);

                price += currCF[idxCF].amount*disc;
                grad  += currCF[idxCF].amount*disc*rateVar[idxCF];
                hess  += currCF[idxCF].amount*disc*rateVar[idxCF]*rateVar[idxCF];
            }
         
            error = price - currPrice;
            drate = -grad/hess + Maths::sign(grad/hess)*sqrt(grad/hess*grad/hess - 2.*error/hess);
        }

        ratePts[idxCve].disc = disc;
        ratePts[idxCve].coeff[0] = lastInstRate;
        ratePts[idxCve].coeff[1] = (6.*rate - 4.*lastInstRate - 2.*instBdryRate)/currIvl;
        ratePts[idxCve].coeff[2] = (3.*lastInstRate + 3.*instBdryRate - 6.*rate)/currIvl/currIvl;

        // Get rid of the point generated by nextCF
        ratePts.pop_back();
        lastInstRate = instBdryRate;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::addFirstCFLinear
(                                     
    const CashFlowArray& currCF,      // (I) cashflows to be added
    double               currPrice,   // (I) price of the current cashflows
    const CashFlowArray& nextCF,      // (I) next cashflow to be added
    double               nextPrice,   // (I) price of the next cashflows
    double&              lastInstRate // (O) instantaneous forward rate at end of period
)
{
    static const string method = "LMParabCurve::addFirstCFLinear";

    try
    {
        double   instBdryRate;
        double   prevDisc;
        double   rate;
        double   startPV = 0.;
        double   price;
        double   grad;
        double   hess;
        double   ivlRatio;
        double   disc;
        double   error;
        double   drate;

        double   sumFlows   = 0.;
        double   sumFlowsDt = 0.;
        double   flowsDur;

        int      idxCve;
        int      idxCF;
        int      firstCF;
        int      dateIvl;
        int      currIvl;

        DateTime prevDate;

        idxCve = ratePts.size();

        // Add the next two cashflows with flat parabolic
        // and later get rid of the second one
        addNextCFFlat(currCF, currPrice);
        addNextCFFlat(nextCF, nextPrice);

        // Get the instantaneous rate at the maturity of currCF
        instBdryRate = this->instBdryRate(currCF, nextCF);

        // Variables to price currCF
        if (idxCve == 0)
        {
            prevDate = baseDate;
            prevDisc = 1.;
        }
        else
        {
            prevDate = ratePts[idxCve - 1].date;
            prevDisc = ratePts[idxCve - 1].disc;
        }

        // First cashflow after base date
        for (idxCF = 0; idxCF < currCF.size(); ++idxCF)
        {
            if (currCF[idxCF].date > baseDate)
            {
                break;
            }
        }

        for (; idxCF < currCF.size(); ++idxCF)
        {
            if (currCF[idxCF].date > prevDate)
            {
                break;
            }

            disc = this->disc(baseDate, currCF[idxCF].date);

            startPV += currCF[idxCF].amount*disc;
        }

        // Value cashflow according to initial rate guess
        firstCF = idxCF;
        vector<double> rateVar(currCF.size());
        vector<double> rateCte(currCF.size());

        currIvl = noAccrue->businessDaysDiff(prevDate, ratePts[idxCve].date);

        for (idxCF = firstCF; idxCF < currCF.size(); ++idxCF)
        {
            dateIvl = noAccrue->businessDaysDiff(prevDate, currCF[idxCF].date);

            ivlRatio = double(dateIvl)/currIvl;

            rateVar[idxCF] = ivlRatio*(ivlRatio - 2.);
            rateCte[idxCF] = instBdryRate*dateIvl*(1. - ivlRatio);

            sumFlows   += currCF[idxCF].amount;
            sumFlowsDt += currCF[idxCF].amount*dateIvl;
        }
        rate = 0.;
        flowsDur = sumFlowsDt/sumFlows;
        // Exact date for a single cashflow
        drate = -log((currPrice - startPV)/sumFlows/prevDisc)/flowsDur;
        error = startPV - currPrice;

        // Newton Raphson routine (function is increasing and convex
        // so convergence is guaranteed)

        while (fabs(error/currPrice) > LM_RELATIVE_ERROR || fabs(drate*flowsDur) > LM_RELATIVE_ERROR)
        {
            price = startPV;
            grad  = 0.;
            hess  = 0.;
            rate += drate;

            for (idxCF = firstCF; idxCF < currCF.size(); ++idxCF)
            {
                disc = prevDisc*exp(rateCte[idxCF] + rateVar[idxCF]*rate);

                price += currCF[idxCF].amount*disc;
                grad  += currCF[idxCF].amount*disc*rateVar[idxCF];
                hess  += currCF[idxCF].amount*disc*rateVar[idxCF]*rateVar[idxCF];
            }

            error = price - currPrice;
            drate = -grad/hess + Maths::sign(grad/hess)*sqrt(grad/hess*grad/hess - 2.*error/hess);
        }

        ratePts[idxCve].disc = disc;
        ratePts[idxCve].coeff[0] = 2.*rate/currIvl - instBdryRate;
        ratePts[idxCve].coeff[1] = 2.*(instBdryRate - rate/currIvl)/currIvl;
        ratePts[idxCve].coeff[2] = 0.;

        // Get rid of the point generated by nextCF
        ratePts.pop_back();
        lastInstRate = instBdryRate;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::addLastCFParab
(
    const CashFlowArray& currCF,      // (I) cashflows to be added
    double               currPrice,   // (I) price of the current cashflows
    double               lastInstRate // (I) instantaneous forward rate at end of period
)
{
    static const string method = "LMParabCurve::addLastCFParab";

    try
    {
        double   instBdryRate;
        double   prevDisc;
        double   rate;
        double   startPV = 0.;
        double   price;
        double   grad;
        double   hess;
        double   ivlRatio;
        double   disc;
        double   error;
        double   drate;

        double   sumFlows   = 0.;
        double   sumFlowsDt = 0.;
        double   flowsDur;

        int      idxCve;
        int      idxCF;
        int      firstCF;
        int      dateIvl;
        int      currIvl;

        DateTime prevDate;

        idxCve = ratePts.size();

        // Add the next cashflows with flat parabolic
        addNextCFFlat(currCF, currPrice);

        // Get the instantaneous rate at the maturity of currCF
        instBdryRate = instBdryRateLast(currCF, lastInstRate);
        instBdryRate = 0.5*(ratePts[idxCve].coeff[0] + instBdryRate);

        // Variables to price currCF
        if (idxCve == 0)
        {
            prevDate = baseDate;
            prevDisc = 1.;
        }
        else
        {
            prevDate = ratePts[idxCve - 1].date;
            prevDisc = ratePts[idxCve - 1].disc;
        }

        // First cashflow after base date
        for (idxCF = 0; idxCF < currCF.size(); ++idxCF)
        {
            if (currCF[idxCF].date > baseDate)
            {
                break;
            }
        }

        for (; idxCF < currCF.size(); ++idxCF)
        {
            if (currCF[idxCF].date > prevDate)
            {
                break;
            }

            disc = this->disc(baseDate, currCF[idxCF].date);

            startPV += currCF[idxCF].amount*disc;
        }

        // Value cashflow according to initial rate guess
        firstCF = idxCF;
        vector<double> rateVar(currCF.size());
        vector<double> rateCte(currCF.size());

        currIvl = noAccrue->businessDaysDiff(prevDate, ratePts[idxCve].date);

        for (idxCF = firstCF; idxCF < currCF.size(); ++idxCF)
        {
            dateIvl = noAccrue->businessDaysDiff(prevDate, currCF[idxCF].date);

            ivlRatio = double(dateIvl)/currIvl;

            rateVar[idxCF] = -dateIvl*ivlRatio*(3. - 2.*ivlRatio);
            rateCte[idxCF] = -lastInstRate*dateIvl*(1. - ivlRatio*(2. - ivlRatio)) +
                             instBdryRate*dateIvl*ivlRatio*(1. - ivlRatio);
                          
            sumFlows   += currCF[idxCF].amount;
            sumFlowsDt += currCF[idxCF].amount*dateIvl;
        }
        rate = 0.;
        flowsDur = sumFlowsDt/sumFlows;
        // Exact date for a single cashflow
        drate = -log((currPrice - startPV)/sumFlows/prevDisc)/flowsDur;
        error = startPV - currPrice;

        // Newton Raphson routine (function is increasing and convex
        // so convergence is guaranteed)

        while (fabs(error/currPrice) > LM_RELATIVE_ERROR || fabs(drate*flowsDur) > LM_RELATIVE_ERROR)
        {
            price = startPV;
            grad  = 0.;
            hess  = 0.;
            rate += drate;

            for (idxCF = firstCF; idxCF < currCF.size(); ++idxCF)
            {
                disc = prevDisc*exp(rateCte[idxCF] + rateVar[idxCF]*rate);

                price += currCF[idxCF].amount*disc;
                grad  += currCF[idxCF].amount*disc*rateVar[idxCF];
                hess  += currCF[idxCF].amount*disc*rateVar[idxCF]*rateVar[idxCF];
            }

            error = price - currPrice;
            drate = -grad/hess + Maths::sign(grad/hess)*sqrt(grad/hess*grad/hess - 2.*error/hess);
        }

        ratePts[idxCve].disc = disc;
        ratePts[idxCve].coeff[0] = lastInstRate;
        ratePts[idxCve].coeff[1] = (6.*rate - 4.*lastInstRate - 2.*instBdryRate)/currIvl;
        ratePts[idxCve].coeff[2] = (3.*lastInstRate + 3.*instBdryRate - 6.*rate)/currIvl/currIvl;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


double LMParabCurve::instBdryRate
(
    const CashFlowArray& currCF, // (I) cashflows to be added
    const CashFlowArray& nextCF  // (I) next cashflow to be added
)
const
{
    static const string method = "LMParabCurve::instBdryRate";

    try
    {   
        int      idxCF;

        DateTime prevDate;
        DateTime currDate;
        double   currFlatRate;
        double   nextFlatRate;

        int      prevTime;
        int      currTime;
        int      dateIvl;

        double   currDur  = 0.;
        double   currConv = 0.;
        double   nextDur  = 0.;
        double   nextConv = 0.;

        double   aux;
        double   disc;

        double   weight;

        if (ratePts.size() < 2)
        {
            throw ModelException("Need at least two points in curve.");
        }

        prevDate = ratePts.size() == 2 ? baseDate : ratePts[ratePts.size() - 3].date;
        currDate = ratePts[ratePts.size() - 2].date;

        currFlatRate = ratePts[ratePts.size() - 2].coeff[0];
        nextFlatRate = ratePts.back().coeff[0];

        prevTime = noAccrue->businessDaysDiff(baseDate, prevDate);
        currTime = noAccrue->businessDaysDiff(baseDate, currDate);

        // Get rid of cashflows preceding prevDate
        for (idxCF = 0; idxCF < currCF.size() && currCF[idxCF].date <= prevDate; ++idxCF);
        for (; idxCF < currCF.size(); ++idxCF)
        {
            dateIvl = noAccrue->businessDaysDiff(baseDate, currCF[idxCF].date);

            disc = exp(-currFlatRate*(dateIvl - prevTime));
            aux  = currCF[idxCF].amount*(dateIvl - prevTime)*disc;

            currDur  += aux;
            currConv += aux*(dateIvl + prevTime);
        }

        // Get rid of cashflows preceding currDate
        for (idxCF = 0; idxCF < nextCF.size() && nextCF[idxCF].date <= currDate; ++idxCF);
        for (; idxCF < nextCF.size(); ++idxCF)
        {
            dateIvl = noAccrue->businessDaysDiff(baseDate, nextCF[idxCF].date);

            disc = exp(-nextFlatRate*(dateIvl - currTime));
            aux  = nextCF[idxCF].amount*(dateIvl - currTime)*disc;

            nextDur  += aux;
            nextConv += aux*(dateIvl + currTime);
        }

        weight = (nextConv/nextDur - 2.*currTime)/(nextConv/nextDur - currConv/currDur);
        return weight*currFlatRate + (1. - weight)*nextFlatRate;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


double LMParabCurve::instBdryRateLast
(
    const CashFlowArray& currCF,      // (I) cashflows to be added
    double               lastInstRate // (I) instantaneous forward rate at start of period
)
const
{
    static const string method = "LMParabCurve::instBdryRateLast";

    try
    {
        int      idxCF;

        DateTime prevDate;
        DateTime currDate;
        double   currFlatRate;

        int      prevTime;
        int      currTime;
        int      dateIvl;

        double   currDur  = 0.;
        double   currConv = 0.;

        double   aux;
        double   disc;

        double   weight;

        if (ratePts.size() < 2)
        {
            throw ModelException("Need at least two points in curve.");
        }

        prevDate = ratePts[ratePts.size() - 2].date;
        currDate = ratePts.back().date;

        currFlatRate = ratePts.back().coeff[0];

        prevTime = noAccrue->businessDaysDiff(baseDate, prevDate);
        currTime = noAccrue->businessDaysDiff(baseDate, currDate);

        // Get rid of cashflows preceding prevDate
        for (idxCF = 0; idxCF < currCF.size() && currCF[idxCF].date <= prevDate; ++idxCF);
        for (; idxCF < currCF.size(); ++idxCF)
        {
            dateIvl = noAccrue->businessDaysDiff(baseDate, currCF[idxCF].date);

            disc = exp(-currFlatRate*(dateIvl - prevTime));
            aux  = currCF[idxCF].amount*(dateIvl - prevTime)*disc;

            currDur  += aux;
            currConv += aux*(dateIvl + prevTime);
        }

        aux = 0.5*currConv/currDur;

        weight = (currTime - aux)/(prevTime - aux);

        return (1. - weight)*currFlatRate + weight*lastInstRate;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::linZerosCoeffsRatePt
(
    RatePt&         rp,        // (I/O) parabolic rate interval
    const DateTime& baseDate,  // (I) needed to compute the zero rate
    const DateTime& startDate, // (I) start date of the interval
    double          startDisc, // (I) discount factor to start of period
    double          timeDenom, // (I) time denominator for zero coupon rates
    HolidayConstSP  noAccrue   // (I) dates with no accrual
)
{
    static const string method = "LMParabCurve::linZerosCoeffsRatePt";

    try
    {
        double startZeroRate;
        double endZeroRate;
        double slopeZeroRate;

        double effRate1;
        double effRate2;
        double effRate3;

        int    daysIvl;
        int    daysToStart;

        daysToStart = noAccrue->businessDaysDiff(baseDate, startDate);

        daysIvl     = noAccrue->businessDaysDiff(startDate, rp.date);

        if (daysIvl <= 0 || daysToStart <= 0)
        {
            throw ModelException("Non-positive interval size.");
        }

        startZeroRate = pow(startDisc, -timeDenom/daysToStart) - 1.;
        endZeroRate   = pow(rp.disc, -timeDenom/(daysToStart + daysIvl)) - 1.;
        slopeZeroRate = (endZeroRate - startZeroRate)/daysIvl;

        // Match the effective rate at one third of the whole interval
        effRate1 = (daysToStart + daysIvl/3.)*
                   log(1. + startZeroRate + slopeZeroRate*daysIvl/3.)/timeDenom +
                   log(startDisc);
        effRate2 = (daysToStart + daysIvl*2./3.)*
                   log(1. + startZeroRate + slopeZeroRate*daysIvl*2./3.)/timeDenom +
                   log(startDisc);
        effRate3 = (daysToStart + daysIvl)*
                   log(1. + startZeroRate + slopeZeroRate*daysIvl)/timeDenom +
                   log(startDisc);

        rp.coeff[0] = (effRate3 - 4.5*effRate2 + 9.*effRate1)/daysIvl;
        rp.coeff[1] = -(9.*effRate3 - 36.*effRate2 + 45.*effRate1)/(daysIvl*daysIvl);
        rp.coeff[2] = 13.5*(effRate3 - 3.*effRate2 + 3.*effRate1)/(daysIvl*daysIvl*daysIvl);
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::linZerosGenFromFras
(
    const DateTime&              bd,  // (I) base date for curve
    const vector<FwdRateIvl>&    frl, // (I)
    const vector<YearEndEffect>& yel, // (I)
    HolidayConstSP               na   // (I) dates with no accrual
)
{
    static const string method = "LMParabCurve::linZerosGenFromFras";

    try
    {
        HolidaySP noAccrueYE;

        vector<FwdRateIvl> adjFraList;

        bool areWeekendsHolidays;
        double timeDenom;

        size_t idxFra;

        DateTime lastDate;
        double   lastDisc;

        ratePts.clear();

        areCoeffsValid = false;
        baseDate       = bd;
        noAccrue       = na;

        noAccrue->toALIB(areWeekendsHolidays);
        timeDenom = !areWeekendsHolidays ? 365. : 252.;

        // If curve is empty
        if (frl.empty())
        {
            return;
        }

        // Create a merged list of noAccrue with the
        // non spread year end effect dates
        yearEndAndNoAccrueCombine(yel, noAccrue, noAccrueYE);

        // Get the normalised fra list
        fraListAdjustYearEnd(yel, frl, adjFraList);

        // Check if start of first fra is after baseDate
        if (baseDate > adjFraList[0].startDate)
        {
            throw ModelException("Base date comes after first fra date.");
        }

        if (baseDate != adjFraList[0].startDate)
        {
            throw ModelException("First fra date must be baseDate.");
        }

        noAccrue = noAccrueYE;

        areCoeffsValid = true;

        lastDate = baseDate;
        lastDisc = 1.;
        idxFra = 0;

        // First period is flat
        addNextFraFlat(lastDate,
                       lastDisc,
                       adjFraList,
                       idxFra     );

        // The rest is linear zeros

        for (; idxFra < adjFraList.size(); ++idxFra)
        {
            addNextFraLinearZeros(lastDate,
                                  lastDisc,
                                  adjFraList,
                                  idxFra,
                                  timeDenom  );
        }

        addYE(yel, noAccrue);
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::addNextFraLinearZeros
(
    DateTime&                 lastDate, // (I/O) last date entered into curve
    double&                   lastDisc, // (I/O) discount up to lastDate
    const vector<FwdRateIvl>& fraList,  // (I) fras to be included
    size_t&                   idxFra,   // (I/O) index of Fra being used currently
    double                    timeDenom // (I) time denominator for zero coupon rates
)
{
    static const string method = "LMParabCurve::addNextFraLinearZeros";

    try
    {
        DateTime prevDate;
        double   prevDisc;
        RatePt   rp;

        prevDate = lastDate;
        prevDisc = lastDisc;

        nextFraDisc(fraList[idxFra],
                    rp.date,
                    rp.disc         );

        // If the Fra start is in the future, reuse the Fra
        if (rp.date == fraList[idxFra].startDate)
            --idxFra;

        linZerosCoeffsRatePt(rp,
                             baseDate,
                             prevDate,
                             prevDisc,
                             timeDenom,
                             noAccrue  );

        ratePts.push_back(rp);

        lastDate = rp.date;
        lastDisc = rp.disc;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::addNextCFLinearZeros
(
    const CashFlowArray& cashFlow, // (I) cashflow to be added
    double               cfPrice,  // (I) PV of the cashflow to curve base date        
    double               timeDenom // (I) time denominator for zero coupon rates
)
{
    static const string method = "LMParabCurve::addNextCFLinearZeros";

    try
    {
        DateTime prevDate;
        double   prevDisc;
        double   disc;
        double   startRate;
        double   rateSlope;
        double   dslope;

        double   startPV = 0.;

        int      startIvl;

        int      idxCF;
        int      firstCF;

        double   price;
        double   grad;
        double   hess;
        double   error;

        double   ivlDur   = 0.;
        double   sumFlows = 0.;

        RatePt   rp;
        double   x;

        if (ratePts.empty())
        {
            throw ModelException("Need at least one point in curve.");
        }

        prevDate = ratePts.back().date;
        prevDisc = ratePts.back().disc;

        startIvl = noAccrue->businessDaysDiff(baseDate, prevDate);

        if (startIvl <= 0)
        {
            throw ModelException ("First interval cannot be empty.");
        }

        startRate = pow(prevDisc, -timeDenom/startIvl);

        // First cash flow after base date
        for (idxCF = 0; idxCF < cashFlow.size(); ++idxCF)
        {
            if (cashFlow[idxCF].date > baseDate)
            {
                break;
            }
        }

        if (idxCF == cashFlow.size())
        {
            throw ModelException("Cash flow expired before base date of curve.");
        }

        for (; idxCF < cashFlow.size(); ++idxCF)
        {
            if (cashFlow[idxCF].date > prevDate)
            {
                break;
            }

            disc = this->disc(baseDate,
                              cashFlow[idxCF].date);

            startPV += cashFlow[idxCF].amount*disc;
        }

        if (idxCF == cashFlow.size())
        {
            throw ModelException("Cash flow expired before last date of curve.");
        }

        // Value cash flow according to initial slope guess
        firstCF = idxCF;
        vector<int   > dateIvl (cashFlow.size());
        vector<double> ctePrice(cashFlow.size());
        vector<double> cteGrad (cashFlow.size());
        vector<double> cteHess (cashFlow.size());
        vector<double> time    (cashFlow.size());

        for (; idxCF < cashFlow.size(); ++idxCF)
        {
            dateIvl [idxCF] = noAccrue->businessDaysDiff(prevDate,
                                                         cashFlow[idxCF].date);
            
            time    [idxCF] = (dateIvl[idxCF] + startIvl)/timeDenom;
            ctePrice[idxCF] = cashFlow[idxCF].amount;
            cteGrad [idxCF] = -ctePrice[idxCF]*time[idxCF]*dateIvl[idxCF];
            cteHess [idxCF] = -cteGrad[idxCF]*(1. + time[idxCF])*dateIvl[idxCF];
            ivlDur += dateIvl[idxCF]*cashFlow[idxCF].amount;
            sumFlows += cashFlow[idxCF].amount;
        }

        ivlDur/= sumFlows;
        // First guess exact for single cash flow
        rateSlope = (pow((cfPrice - startPV)/sumFlows,
                         -timeDenom/(startIvl + ivlDur)) - startRate)/ivlDur;

        price = startPV;
        grad = 0.;
        hess = 0.;

        for (idxCF = firstCF; idxCF < cashFlow.size(); ++idxCF)
        {
            disc = pow(startRate + rateSlope*dateIvl[idxCF], -time[idxCF]);
            x    = startRate + rateSlope*dateIvl[idxCF];

            price += ctePrice[idxCF]*disc;
            grad  += cteGrad[idxCF]*disc/x;
            hess  += cteHess[idxCF]*disc/(x*x);
        }

        error  = price - cfPrice;
        dslope = -grad/hess + Maths::sign(grad/hess)*sqrt(grad/hess*grad/hess - 2.*error/hess);

        // Newton Raphson routine (function is increasing and convex
        // so convergence is guaranteed)

        while (fabs(error/cfPrice) > LM_RELATIVE_ERROR || fabs(dslope*ivlDur) > LM_RELATIVE_ERROR)
        {
            price      = startPV;
            grad       = 0.;
            hess       = 0.;
            rateSlope += dslope;

            for (idxCF = firstCF; idxCF < cashFlow.size(); ++idxCF)
            {
                x    = startRate + rateSlope*dateIvl[idxCF];
                disc = pow(x, -time[idxCF]);

                price += ctePrice[idxCF]*disc;
                grad  += cteGrad[idxCF]*disc/x;
                hess  += cteHess[idxCF]*disc/(x*x);
            }
         
            error  = price - cfPrice;
            dslope = -grad/hess + Maths::sign(grad/hess)*sqrt(grad/hess*grad/hess - 2.*error/hess);
        }

        for (idxCF = firstCF; idxCF < cashFlow.size(); ++idxCF)
        {
            disc = pow(startRate + rateSlope*dateIvl[idxCF], -time[idxCF]);
            rp.date = cashFlow[idxCF].date;
            rp.disc = disc;

            linZerosCoeffsRatePt(rp,
                                 baseDate,
                                 prevDate,
                                 prevDisc,
                                 timeDenom,
                                 noAccrue  );

            ratePts.push_back(rp);
            prevDate = cashFlow[idxCF].date;
            prevDisc = disc;
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::linZerosAddCF
(
    const vector<CashFlowArray>& cf,     // (I) vector of cashflow arrays to be added
    const vector<double>&        cfPrice // (I) array of prices
)
{
    static const string method = "LMParabCurve::linZerosAddCF";

    try
    {
        size_t idxCF;

        bool   areWeekendsHolidays;
        double timeDenom;

        if (cf.empty())
        {
            return;
        }
        else if (cfPrice.size() != cf.size())
        {
            throw ModelException("cf and cfPrice size mismatch.");
        }

        noAccrue->toALIB(areWeekendsHolidays);
        timeDenom = !areWeekendsHolidays ? 365. : 252.;

        // For no previous elements starts with flat
        idxCF = 0;
        if (ratePts.empty())
        {
            addNextCFFlat(cf     [idxCF],
                          cfPrice[idxCF]);
            ++idxCF;
        }

        // The rest is linear zeros
        for (; idxCF < cf.size(); ++idxCF)
        {
            addNextCFLinearZeros(cf[idxCF],
                                 cfPrice[idxCF],
                                 timeDenom      );
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::genFromMM
(
    const DateTime&           bd,                // (I) base date for zero curve.  Discounts from this date
    const vector<DateTime>&   mmMatDates,        // (I) vector of MM maturity dates
    const vector<double>&     mmRates,           // (I) vector of MM rates
    const vector<DateTime>&   fraStartDates,     // (I) vector of Fra start dates
    const vector<DateTime>&   fraMatDates,       // (I) vector of Fra maturity dates
    const vector<double>&     fraRates,          // (I) vector of Fra rates
    const vector<DateTime>&   yeStartDates,      // (I) vector of Year End adjustment start dates
    const vector<DateTime>&   yeEndDates,        // (I) vector of Year End adjustment end dates
    const vector<double>&     yeRates,           // (I) vector of Year End adjustment rates
    const vector<bool>&       yeIsSpread,        // (I) if TRUE the Year End adjustment is used as a
                                                 //     spread, otherwise in absolute level
    const vector<DateTime>&   relMMDates,        // (I) vector of relative MM maturity dates
    const vector<double>&     relMMSpreads,      // (I) vector of spreads for relative MM points
    const DayCountConvention& mmDayCount,        // (I) Common Day Count Convention for MM, Fras and YE
    int                       rateBasis,         // (I) Common Rate Basis for MM, Fras and YE
    const vector<DateTime>&   flatRateDates,     // (I) the rate will be constant in between flatRateDates,
                                                 //     last date must be <= flatParabBdryDate
    const DateTime&           flatParabBdryDate, // (I) Date at which interpolation transitions
                                                 //     from flat to parabolic
    HolidayConstSP            na,                // (I) No accrual dates 
    const string&             interpType         // (I) Smooth forwards or linear
)
{
    static const string method = "LMParabCurve::genFromMM";

    try
    {
        double rr;
        size_t idxMM;

        ratePts.clear();

        areCoeffsValid = false;
        baseDate       = bd;
        noAccrue       = na;

        // Basic array and data validation

        validateMMData(mmMatDates,
                       mmRates,
                       fraStartDates,
                       fraMatDates,
                       fraRates,
                       yeStartDates,
                       yeEndDates,
                       yeRates,
                       yeIsSpread,
                       relMMDates,
                       relMMSpreads,
                       interpType    );

        // Create the first zero curve which only uses
        // benchmark prices and rates.  The first curve
        // has MM + fra points

        vector<DateTime>      mmStartDates(mmMatDates.size(), baseDate);
        vector<FwdRateIvl>    benchmarks;
        vector<YearEndEffect> yeList;

        appendFwdRateIvlListMM(mmStartDates,
                               mmMatDates,
                               mmRates,
                               mmDayCount,
                               rateBasis,
                               benchmarks);

        appendFwdRateIvlListMM(fraStartDates,
                               fraMatDates,
                               fraRates,
                               mmDayCount,
                               rateBasis,
                               benchmarks);

        appendYearEndEffectListMM(yeStartDates,
                                  yeEndDates,
                                  yeRates,
                                  yeIsSpread,
                                  mmDayCount,
                                  rateBasis,
                                  yeList       );

        if (interpType == ZC3ZeroInterpolation::SMOOTH_FORWARDS)
        {
            genFromFras(baseDate,
                        benchmarks,
                        yeList,
                        noAccrue,
                        flatRateDates,
                        flatParabBdryDate);
        }
        else if (interpType == ZC3ZeroInterpolation::LINEAR)
        {
            linZerosGenFromFras(baseDate,
                                benchmarks,
                                yeList,
                                noAccrue   );
        }

        // Price the relative benchmark points
        if (!relMMDates.empty())
        {
            vector<DateTime> relMMStartDates;
            vector<double>   relRates;

            for (idxMM = 0; idxMM < relMMDates.size(); ++idxMM)
            {
                rr = mmRate(baseDate,
                            relMMDates[idxMM],
                            mmDayCount,
                            rateBasis) + relMMSpreads[idxMM];

                relMMStartDates.push_back(baseDate);
                relRates       .push_back(rr      );
            }

            appendFwdRateIvlListMM(relMMStartDates,
                                   relMMDates,
                                   relRates,
                                   mmDayCount,
                                   rateBasis,
                                   benchmarks);

            if (interpType == ZC3ZeroInterpolation::SMOOTH_FORWARDS)
            {
                genFromFras(baseDate,
                            benchmarks,
                            yeList,
                            noAccrue,
                            flatRateDates,
                            flatParabBdryDate);
            }
            else if (interpType == ZC3ZeroInterpolation::LINEAR)
            {
                linZerosGenFromFras(baseDate,
                                    benchmarks,
                                    yeList,
                                    noAccrue   );
            }
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::validateMMData
(
    const vector<DateTime>& mmMatDates,    // (I) vector of MM maturity dates
    const vector<double>&   mmRates,       // (I) vector of MM rates
    const vector<DateTime>& fraStartDates, // (I) vector of Fra start dates
    const vector<DateTime>& fraMatDates,   // (I) vector of Fra maturity dates
    const vector<double>&   fraRates,      // (I) vector of Fra rates
    const vector<DateTime>& yeStartDates,  // (I) vector of Year End adjustment start dates
    const vector<DateTime>& yeEndDates,    // (I) vector of Year End adjustment end dates
    const vector<double>&   yeRates,       // (I) vector of Year End adjustment rates
    const vector<bool>&     yeIsSpread,    // (I) if TRUE the Year End adjustment is used as a
                                            //     spread, otherwise in absolute level
    const vector<DateTime>& relMMDates,    // (I) vector of relative MM maturity dates
    const vector<double>&   relMMSpreads,  // (I) vector of spreads for relative MM points
    const string&           interpType     // (I) Smooth forwards or linear
)
{
    static const string method = "LMParabCurve::validateMMData";

    try
    {
        size_t   idx;
        DateTime maxMarketDate;

        if (mmRates.size() != mmMatDates.size())
        {
            throw ModelException("Size mismatch between MM maturity dates and rates.");
        }

        for (idx = 1; idx < mmMatDates.size(); ++idx)
        {
            if (mmMatDates[idx - 1] > mmMatDates[idx])
            {
                throw ModelException("Money market maturity dates out of order.");
            }
        }

        if (!mmMatDates.empty())
        {
            maxMarketDate = mmMatDates.back();
        }

        if (fraMatDates.size() != fraStartDates.size())
        {
            throw ModelException("Size mismatch between fra start and maturity dates.");
        }

        if (fraRates.size() != fraStartDates.size())
        {
            throw ModelException("Size mismatch between fra dates and rates");
        }

        for (idx = 0; idx < fraStartDates.size(); ++idx)
        {
            if (fraStartDates[idx] > fraMatDates[idx])
            {
                throw ModelException("Fra start and maturity dates out of order.");
            }
            maxMarketDate = max(maxMarketDate, fraMatDates[idx]);
        }

        if (yeEndDates.size() != yeStartDates.size())
        {
            throw ModelException("Size mismatch between year end start and end dates");
        }

        if (yeRates.size() != yeStartDates.size())
        {
            throw ModelException("Size mismatch between year end dates and rates");
        }

        if (yeIsSpread.size() != yeStartDates.size())
        {
            throw ModelException("Size mismatch between year end dates and spread flags.");
        }

        for (idx = 0; idx < yeStartDates.size(); ++idx)
        {
            if (yeStartDates[idx] > yeEndDates[idx] ||
                idx > 0 && yeStartDates[idx] < yeEndDates[idx - 1])
            {
                throw ModelException("Year end start and maturity dates out of order.");
            }
        }

        if (relMMSpreads.size() != relMMDates.size())
        {
            throw ModelException("Size mismatch between relative MM dates and spreads.");
        }

        for (idx = 1; idx < relMMDates.size(); ++idx)
        {
            if (relMMDates[idx - 1] > relMMDates[idx])
            {
                throw ModelException("Relative MM dates out of order.");
            }
        }

        if (!relMMDates.empty() && relMMDates.back() > maxMarketDate)
        {
            throw ModelException("Relative MM date is beyond benchmark data.");
        }

        if (interpType != ZC3ZeroInterpolation::SMOOTH_FORWARDS &&
            interpType != ZC3ZeroInterpolation::LINEAR            )
        {
            throw ModelException("Wrong interpolation type " + interpType + ".");
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::appendFwdRateIvlListMM
(
    const vector<DateTime>&   startDates, // (I)
    const vector<DateTime>&   matDates,   // (I)
    const vector<double>&     rates,      // (I)
    const DayCountConvention& dayCount,   // (I)
    int                       rateBasis,  // (I)
    vector<FwdRateIvl>&       list        // (I/O) List to be added to
)
{
    static const string method = "LMParabCurve::addFwdRateIvlListMM";

    try
    {
        size_t idx;

        for (idx = 0; idx < startDates.size(); ++idx)
        {
            FwdRateIvl fri =
            {
                startDates[idx],
                matDates[idx],
                -log(RateConversion::rateToDiscount(rates[idx],
                                                    startDates[idx],
                                                    matDates[idx],
                                                    &dayCount,
                                                    rateBasis       ))
            };

            list.push_back(fri);
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::appendYearEndEffectListMM
(
    const vector<DateTime>&   startDates,   // (I)
    const vector<DateTime>&   matDates,     // (I)
    const vector<double>&     rates,        // (I)
    const vector<bool>&       isRateSpread, // (I)
    const DayCountConvention& dayCount,     // (I)
    int                       rateBasis,    // (I)
    vector<YearEndEffect>&    list          // (I/O) List to be added to
)
{
    static const string method = "LMParabCurve::appendYearEndEffectListMM";

    try
    {
        size_t idx;

        for (idx = 0; idx < startDates.size(); ++idx)
        {
            YearEndEffect yee =
            {
                {
                    startDates[idx],
                    matDates[idx],
                    -log(RateConversion::rateToDiscount(rates[idx],
                                                        startDates[idx],
                                                        matDates[idx],
                                                        &dayCount,
                                                        rateBasis       ))
                },

                isRateSpread[idx]
            };

            list.push_back(yee);
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


bool LMParabCurve::isValidFwdRateIvlList
(
    const vector<FwdRateIvl>& list // (I) List to be checked
)
{
    static const string method = "LMParabCurve::isValidFwdRateIvlList";
    size_t idx;

    try
    {
        if (list.empty())
        {
            return true;
        }
        else
        {
            for (idx = 0; idx < list.size(); ++idx)
            {
                DateTime thisStart = list[idx].startDate;
                DateTime thisMat   = list[idx].matDate;

                if (thisStart >= thisMat)
                {
                    char buffer[256];

                    sprintf(buffer, "Failed: startDate[%d] = %s\n"
                                    "      >   matDate[%d] = %s.",
                          idx, thisStart.toString().c_str(),
                          idx, thisMat.  toString().c_str());

                    throw ModelException(buffer);
                }
            }

            // The list is in order
            return true;
        }
        
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


bool LMParabCurve::isStripFwdRateIvlList
(
    const vector<FwdRateIvl>& list // (I) FwdRateIvl list to be checked
)
{
    size_t idx;

    if (list.empty())
    {
        return true;
    }
    else if (list.size() >= 2)
    {
        for (idx = 0; idx < list.size() - 2; ++idx)
        {
            if (list[idx].matDate > list[idx + 1].startDate)
            {
                return false;
            }
        }
    }

    return true;
}


bool LMParabCurve::isOrderedYearEndEffectList
(
    const vector<YearEndEffect>& list // (I) YearEndEffect list to be checked
)
{
    static const string method = "LMParabCurve::isOrderedYearEndEffectList";

    try
    {
        size_t idx;

        if (list.empty())
        {
            return true;
        }
        else
        {
            DateTime lastStart = list.back().ivl.startDate;
            DateTime lastMat   = list.back().ivl.matDate;
            char buffer[256];

            // Test the final ivl first.
            if (lastStart >= lastMat)
            {
                sprintf(buffer, "   startDate[%d] = %s\n"
                                ">=   matDate[%d] = %s.",
                          list.size() - 1, lastStart.toString().c_str(),
                          list.size() - 1, lastMat  .toString().c_str());
                throw ModelException(buffer);
            }

            // Check other intervals.  Note that if numItems is 1
            // this loop is avoided.
            for (idx = 0; idx < list.size() - 1; ++idx)
            {
                DateTime thisStart = list[idx].ivl.startDate;
                DateTime thisMat   = list[idx].ivl.matDate;
                DateTime nextStart = list[idx + 1].ivl.startDate;

                if (thisMat > nextStart)
                {
                    sprintf(buffer, "    matDate[%d] = %s\n"
                                    "> startDate[%d] = %s.",
                              idx,     thisMat  .toString().c_str(),
                              idx + 1, nextStart.toString().c_str());
                    throw ModelException(buffer);
                }
                if (thisStart >= thisMat)
                {
                    sprintf(buffer, "   startDate[%d] = %s\n"
                                    ">=   matDate[%d] = %s.",
                              idx, thisStart.toString().c_str(),
                              idx, thisMat  .toString().c_str());
                    throw ModelException(buffer);
                }
            }

            // The list is in order
            return true;
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


bool LMParabCurve::sortFwdRateIvlList
(
    vector<FwdRateIvl>& listFras
)
{
    static const string method = "LMParabCurve::sortFwdRateIvlList";

    try
    {
        size_t idx;

        // Check the imput data for validity
        if (listFras.empty())
        {
            return true;
        }

        // Sort the list
        sort(listFras.begin(), listFras.end(), compareFra);

        // Verify the list is monotonic otherwise fail
        for (idx = 0; idx < listFras.size() - 1; ++idx)
        {
            if (compareFra(listFras[idx+1], listFras[idx]))
            {
                char buffer[256];
                sprintf(buffer,
                        "list of Fras out of order at idx=%d and %d.",
                        idx, idx + 1);
                throw ModelException(buffer);
            }
        }

        return true;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


bool LMParabCurve::compareFra
(
    const LMParabCurve::FwdRateIvl& fra1, /*(I) fra1 */
    const LMParabCurve::FwdRateIvl& fra2  /*(I) fra2 */
)
{
    return fra1.matDate < fra2.matDate ||
           fra1.matDate == fra2.matDate && fra1.startDate < fra2.startDate;
}


bool LMParabCurve::sortAndNormalizeFRAList
(
    const vector<FwdRateIvl>& list,
    vector<FwdRateIvl>&       newList
)
{
    static const string method = "LMParabCurve::sortAndNormalizeFRAList";

    try
    {
        newList.clear();
        size_t currentIdx;

        // Initial validation of input list
        if (list.empty())
        {
            return false;
        }

        if (!isValidFwdRateIvlList(list))
        {
            goto failed;
        }

        // Make a copy of the original list and sort it
        // before normalising each cluster
        newList = list;
        if (!sortFwdRateIvlList(newList))
        {
            goto failed;
        }

        // Normalize each cluster starting from the maxDate.
        // Then sort the list again due to changes caused by the
        // cluster.
        currentIdx = newList.size() - 1;

        while (currentIdx > 0)
        {
           normalize1FRADate(currentIdx, newList);
           if (!sortFwdRateIvlList(newList))
           {
               goto failed;
           }
           --currentIdx;
        }

        return true;

failed:
        throw ModelException("Failed.");
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


bool LMParabCurve::normalize1FRADate
(
    size_t indexFRA,         // (I) index of FRA Date to check
    vector<FwdRateIvl>& list // (I/O) FRA dates to normalize and sort
)
{
    size_t idxLow;
    size_t idx;

    if (list.empty())
        return false;

    // Search for 1st FRA in cluster whose matDate is
    // current matDate.
    idxLow = indexFRA;
    while (idxLow >= 0 &&
           list[idxLow].matDate == list[indexFRA].matDate)
    {
        --idxLow;
    }
    ++idxLow;

    if (idxLow == indexFRA)
    {
        // This means a cluster of size 1. No normalization.
        return true;
    }
    else // idxLow < indexFRA
    {
        // Redefine the current FRA by modifying its matDate
        // and rate.
        for (idx = idxLow; idx < indexFRA; ++idx)
        {
            list[idx].matDate = list[idx + 1].startDate;
            list[idx].rate    = list[idx].rate - list[idx + 1].rate;
        }

        return true;
    }
}


void LMParabCurve::yearEndAndNoAccrueCombine
(
    const vector<YearEndEffect>& listYE, // (I)
    HolidayConstSP noAccrue,             // (I)   current no accrue
    HolidaySP&     noAccrueYE            // (I/O) new no accrue
)
{
    static const string method = "LMParabCurve::yearEndAndNoAccrueCombine";

    try
    {
        vector<YearEndEffect>::const_iterator i;
        DateTimeArray::const_iterator j;

        DateTime             currentDate;
        DateTimeArrayConstSP noAccrueList;
        DateTimeArray        dlYE;

        bool useWeekends;

        // Convert the year end dates to a date vector.
        // For an interval (start, mat) add the dates
        // (start + 1 , mat) to the datelist.
        for (i = listYE.begin(); i != listYE.end(); ++i)
        {
            if (i->isRateSpread)
                continue; // Only non Spread YE should be included in the no Accrue File

            // Current date to be added to date vector
            currentDate = i->ivl.startDate.rollDate(1);
            while (currentDate <= i->ivl.matDate)
            {
                dlYE.push_back(currentDate);
                currentDate.rollDate(1);
            }
        }

        // Merge the datelist associated with the no accrue
        // and the datelist associated with YE list.
        noAccrueList = noAccrue->toALIB(useWeekends);
        for (j = noAccrueList->begin(); j != noAccrueList->end(); ++j)
        {
            dlYE.push_back(*j);
        }

        noAccrueYE = HolidaySP(new Holiday("No Accrue", dlYE, useWeekends));
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::fraListAdjustYearEnd
(
    const vector<YearEndEffect>& yel,   
    const vector<FwdRateIvl>&    frl,
    vector<FwdRateIvl>&          newList
)
{
    static const string method = "LMParabCurve::fraListAdjustYearEnd";

    try
    {
        size_t idxFRA; // At the beginning of the loop points to
                       // the current FRA to be adjusted
        size_t idxYE;  // At the beginning of the loop points to
                       // the current YE possibly inserted into the
                       // current FRA

        // End points of current FRA and YE intervals
        DateTime yStart, yMat, fStart, fMat;

        char buffer[256];

        // Sort and normalize the FRA list.  If already so, do nothing
        if (!sortAndNormalizeFRAList(frl, newList))
        {
            goto failed;
        }

        // empty YE list means there are no adjustments
        if (yel.empty())
        {
            return;
        }

        if (!isOrderedYearEndEffectList(yel))
        {
            goto failed;
        }

        // If the FRA List is stripped, then the algorithm is much faster
        // O(numFras + numYE).  This is the first algorithm.
      
        // However, if the FRA List is NOT stripped, then the algorithm
        // cannot be faster than O(numFras*numYE).  This is implemented
        // as a nested for loop.  This is the second algorithm.

        if (isStripFwdRateIvlList(newList))
        {
            idxFRA = 0;
            idxYE  = 0;
            while (idxFRA < newList.size() &&
                   idxYE  < yel.size())
            {
                // Maintain concurrent indices of YE and FRA list separately and
                // scan each list in ascending order.
                // Test whether either current YE or FRA is past its respective
                // End of List, if we are done.
                // Test whether the current YE falls inside the current FRA, 
                // if there is partial overlap return FAILURE
                // if current YE comes AFTER current FRA, then
                // the current FRA does not need anymore adjustment so update
                // the current FRA.
                // if current YE comes BEFORE current FRA, then update current
                // YE with no adjustment.
                // if current YE falls inside current FRA, then do adjustment
                // and update current YE.

                yStart = yel[idxYE].ivl.startDate;
                yMat   = yel[idxYE].ivl.matDate;
                fStart = newList[idxFRA].startDate;
                fMat   = newList[idxFRA].matDate;

                // YE precedes FRA

                if (yMat <= fStart)
                {
                    ++idxYE;
                    continue;
                }

                // YE comes after FRA

                else if ( yStart >= fMat)
                {
                    ++idxFRA;
                    continue;
                }

                // YE contains FRA

                else if (yStart <= fStart && yMat >= fMat)
                {
                    sprintf(buffer, "YE %d contains FRA %d.",
                                    idxYE, idxFRA);
                    throw ModelException(buffer);
                }

                // Normal case: YE is inside FRA

                else if (yStart >= fStart && yMat <= fMat)
                {
                    newList[idxFRA].rate -=
                        yel[idxYE].ivl.rate;
                    ++idxYE;
                    continue;
                }

                // YE intersects FRA ( YE start inside FRA)

                else if (yStart > fStart &&
                         yStart < fMat   &&
                         yMat   > fMat     )
                {
                    sprintf(buffer, "Yearend effect idx %d partially\n"
                                    "intersects FRA idx %d.",
                                    idxYE, idxFRA);
                    throw ModelException(buffer);
                }

                // YE intersects FRA (YE mat inside FRA)

                else if (yMat   > fStart &&
                         yMat   < fMat   &&
                         yStart < fStart   )
                {
                    sprintf(buffer, "Yearend effect idx %d partially\n"
                                    "intersects FRA idx %d.",
                                    idxYE, idxFRA);
                    throw ModelException(buffer);
                }
                // ??? Any remaining bad cases ?

                else
                {
                    throw ModelException("Invalid FRA or YE interval data.");
                }
            }
        }
        else
        {
            // Second algorithm

            for (idxFRA = 0; idxFRA < newList.size(); ++idxFRA)
            {
                for (idxYE = 0; idxYE < yel.size(); ++idxYE)
                {
                    yStart = yel[idxYE].ivl.startDate;
                    yMat   = yel[idxYE].ivl.matDate;
                    fStart = newList[idxFRA].startDate;
                    fMat   = newList[idxFRA].matDate;

                    // YE precedes FRA or
                    // YE comes after FRA
                    if (yMat <= fStart || yStart >= fMat)
                    {
                        continue;
                        // Nothing
                    }
                    // YE contains FRA

                    else if (yStart <= fStart && yMat >= fMat)
                    {
                        sprintf(buffer, "YE %d  contains FRA %d.",
                                        idxYE, idxFRA);
                        throw ModelException(buffer);
                    }

                    // YE intersects FRA (YE start inside FRA)

                    else if (yStart > fStart &&
                             yStart < fMat   &&
                             yMat   > fMat     )
                    {
                        sprintf(buffer, "Yearend effect idx %d partially\n"
                                        "intersects FRA idx %d.",
                                        idxYE, idxFRA);
                        throw ModelException(buffer);
                    }

                    // Normal case: YE is inside FRA

                    else if (yStart >= fStart && yMat <= fMat)
                    {
                        newList[idxFRA].rate -=
                            yel[idxYE].ivl.rate;
                    }

                    // YE intersects FRA (YE mat inside FRA)

                    else if (yMat   > fStart &&
                             yMat   < fMat   &&
                             yStart < fStart)
                    {
                        sprintf(buffer, "Yearend effect Num %d partially\n"
                                        "intersects FRA Num %d.",
                                        idxYE + 1, idxFRA + 1);
                        throw ModelException(buffer);
                    }
                    // ???? Any remaining bad cases?
                    else
                    {
                        throw ModelException("Invalid FRA or YE interval data.");
                    }
                }
            }
        }

        return;
        
failed:
        throw ModelException("Failed.");
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


double LMParabCurve::interpDisc
(
    DateTime desiredDate // (I) where interpolation is done
)
const
{
    static const string method = "LMParabCurve::interpDisc";

    try
    {
         double prevDisc;
         int daysIvl;
         double coeff[3];

         // Forward PV = exp(-INTEGRAL(a + bs + cs^2, s=0, s= T))
         // where T = days from lower benchmark to desired date          
         //
         // Interpolated PV = benchmark PV * fwd PV
         interpParam(desiredDate, prevDisc, daysIvl, coeff);

         return prevDisc*exp(-(daysIvl*(coeff[0] + 
                               daysIvl*(coeff[1]/2.0 +
                               daysIvl* coeff[2]/3.0   ))));
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


double LMParabCurve::instRate
(
    DateTime desiredDate, // (I) where interpolation is done
    int rateType          // (I) LM_INST_RATE, LM_INST_RATE_SLOPE or LM_OVERNIGHT_RATE
)
const
{
    static const string method = "LMParabCurve::instRate";

    try
    {
        double  prevDisc;
        int     daysIvl;  // from desired date to lower benchmark date
        double  coeff[3]; // parabolic coefficient of interpolating interval

        interpParam(desiredDate,
                    prevDisc,
                    daysIvl,
                    coeff       );

        switch (rateType)
        {
            case LM_OVERNIGHT_RATE:
                return coeff[0] + coeff[1]/2. + coeff[2]/3. +
                       daysIvl*(coeff[1] + coeff[2]*(1. + daysIvl));

            case LM_INST_RATE:
                return coeff[0] + daysIvl*(coeff[1] + coeff[2]*daysIvl);

            case LM_INST_RATE_SLOPE:
                return coeff[1] + 2.*coeff[2]*daysIvl;

            default:
                throw ModelException("Wrong rate type.");
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::interpParam
(
    DateTime desiredDate, // (I) where interpolation is done.
    double& prevDisc,     // (O) Discount at start of appropriate interval
    int&    daysIvl,      // (O) days from start of appropriate interval
    double  coeff[3]      // (O) curve coefficients of appropriate interval
)
const
{
    static const string method = "LMParabCurve::interpParam";

    try
    {
        int idxLow;
        int idxHigh;
        DateTime dateLow;
        int totalIvl;

        // Check for input error

        if (ratePts.empty())
        {
            throw ModelException("Curve cannot be empty.");
        }

        if (noAccrue.get() == NULL)
        {
            throw ModelException("No accrual dates not specified.");
        }

        // Look for valid upper and lower bracketing dates
        // and valid desired date
        indexFind(desiredDate, idxLow, idxHigh);
        
        // Compute the interpolated disc

        // Check if desired date is curve base date
        if (idxLow == -1 && idxHigh == -1)
        {
            prevDisc = 1.;
            daysIvl  = 0;
            coeff[0] = 0.;
            coeff[1] = 0.;
            coeff[2] = 0.;

            return;
        }

        // Check if desired date us between base date and first
        // benchmark date
        else if (idxHigh == 0)
        {
            prevDisc = 1.;
            dateLow  = baseDate;
        }

        // Check if desired date is a benchmark date
        else
        {
            prevDisc = ratePts[idxHigh - 1].disc;
            dateLow  = ratePts[idxHigh - 1].date;
        }

        // Compute number of days from lower bracketing
        // benchmark date to desired date
        daysIvl = noAccrue->businessDaysDiff(dateLow, desiredDate);

        if (!areCoeffsValid)
        {
            // Use as default flat forwards
            totalIvl = noAccrue->businessDaysDiff(dateLow, ratePts[idxHigh].date);

            coeff[0] = log(prevDisc/ratePts[idxHigh].disc)/totalIvl;
            coeff[1] = 0.;
            coeff[2] = 0.;
        }
        else
        {
            coeff[0] = ratePts[idxHigh].coeff[0];
            coeff[1] = ratePts[idxHigh].coeff[1];
            coeff[2] = ratePts[idxHigh].coeff[2];
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::indexFind
(
    DateTime date, // (I) Date to intepolate at
    int& low,      // (O) lower bracketing point
    int& high      // (O) upper bracketing point
)
const
{
    static const string method = "LMParabCurve::indexFind";

    try
    {
        if (ratePts.empty())
        {
            throw ModelException("Curve is empty.");
        }

        if (date < baseDate || date > ratePts.back().date)
        {
            char buffer[256];
            sprintf(buffer, "Date %s is out of range [%s, %s].",
                date.toString().c_str(), baseDate.toString().c_str(),
                ratePts.back().date.toString().c_str());
            throw ModelException(buffer);
        }

        // Check the Extreme Endpoints first

        if (date == baseDate)
        {
            low  = -1;
            high = low;
            return;
        }

        if (date == ratePts.back().date)
        {
            low  = ratePts.size() - 1;
            high = low;
            return;
        }

        // Check if the date falls between baseDate and first
        // bechmark date

        if (date > baseDate &&
            date < ratePts.front().date)
        {
            low  = -1;
            high = 0;
            return;
        }

        // Do the binary search

        RatePt rp = {date};
        vector<RatePt>::const_iterator i;
        i = lower_bound(ratePts.begin(), ratePts.end(), rp, compareRatePt);
        if (i == ratePts.end())
        {
            throw ModelException("Failed.");
        }

        low  = &*i - &ratePts[0] - 1;
        high = low + 1;
        for (++i; i != ratePts.end() && i->date == date; ++i, ++high);

        // Check if date falls on a benchmark date

        if (date == ratePts[low].date)
        {
            high = low;
        }
        else if (date == ratePts[high].date)
        {
            low  = high;
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


bool LMParabCurve::compareRatePt
(
    const RatePt& rp1,
    const RatePt& rp2
)
{
    return rp1.date < rp2.date;
}


void LMParabCurve::calcParabolicCoeffsRatePt
(
    RatePt&         rp,            // (I/O) parabolic rate interval
    const DateTime& startDate,     // (I) start date of the interval
    double          startDisc,     // (I) discount factor to start of period
    double          startInstRate, // (I) instantaneous forward rate at start of period
    const DateTime& nextDate,      // (I) date for following interval
    double          nextDisc,      // (I) discount factor to the end of next period
    double&         endInstRate,   // (O) instantaneous forward rate at end of period
    HolidayConstSP  noAccrue       // (I) dates with no accrual
)
{
    static const string method = "LMParabCurve::calcParabolicCoeffsRatePt";

    try
    {
        int    daysIvl;
        int    daysIvlNext;
        double flatRate;
        double flatRateNext;
        double dTime;
        double mrFac;
        double difFac;

        if (nextDate < rp.date)
        {
            throw ModelException("nextDate must not be before end of period.");
        }

        flatRate = getFlatFwdFromDiscountHoliday(startDate,
                                                 startDisc,
                                                 rp.date,
                                                 rp.disc,
                                                 noAccrue  );

        daysIvl = noAccrue->businessDaysDiff(startDate, rp.date);

        if (nextDate == rp.date)
        {
            // Hardcode MeanReversion = 5 % for last interval;
            // This is just to establish some good shape on final interval,
            // and thus is fine to hardcode.
            dTime = (rp.date.daysDiff(startDate))*(rp.date.daysDiff(startDate))/daysIvl/365.;
            mrFac = 0.05*dTime;
            difFac = Maths::isZero(mrFac) ? 
                        1. :
                        mrFac*(1. - exp(-mrFac))/(exp(-mrFac) + mrFac - 1.) - 1.;
            endInstRate = flatRate - difFac*(startInstRate - flatRate);
        }
        else
        {
            flatRateNext = getFlatFwdFromDiscountHoliday(rp.date,
                                                         rp.disc,
                                                         nextDate,
                                                         nextDisc,
                                                         noAccrue );

            daysIvlNext = noAccrue->businessDaysDiff(rp.date, nextDate);

            if (daysIvl <=0 || daysIvlNext <= 0)
            {
                throw ModelException("Non-positive interval size.");
            }

            endInstRate = (flatRate*daysIvlNext + flatRateNext*daysIvl)/
                          (daysIvl + daysIvlNext);
        }

        rp.coeff[0] = startInstRate;
        rp.coeff[1] = (6.*flatRate - 4.*startInstRate - 2.*endInstRate)/daysIvl;
        rp.coeff[2] = (3.*startInstRate + 3.*endInstRate - 6.*flatRate)/daysIvl/daysIvl;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::calcLinearCoeffsRatePt
(
    RatePt&         rp,            // (I/O) parabolic rate interval
    const DateTime& startDate,     // (I) start date of the interval
    double          startDisc,     // (I) discount factor to start of period
    double          startInstRate, // (I) instantaneous forward rate at start of period
    double&         endInstRate,   // (O) instantaneous forwad rate at end of period
    HolidayConstSP  noAccrue       // (I) dates with no accrual
)
{
    static const string method = "LMParabCurve::calcLinearCoeffsRatePt";

    try
    {
        long daysIvl;
        double flatRate;

        flatRate = getFlatFwdFromDiscountHoliday(startDate,
                                                 startDisc,
                                                 rp.date,
                                                 rp.disc,
                                                 noAccrue  );

        daysIvl = noAccrue->businessDaysDiff(startDate, rp.date);

        if (daysIvl <= 0)
        {
            throw ModelException("Non-positive interval size.");
        }
        
        endInstRate = 2*flatRate - startInstRate;

        rp.coeff[0] = startInstRate;
        rp.coeff[1] = (endInstRate - startInstRate)/daysIvl;
        rp.coeff[2] = 0.;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::calcFlatCoeffsRatePt
(
    RatePt&         rp,        // (I/O) parabolic rate interval
    const DateTime& startDate, // (I) start date of the interval
    double          startDisc, // (I) discoiunt factor to start of period
    HolidayConstSP  noAccrue   // (I) dates with no accrual
)
{
    static const string method = "LMParabCurve::calcFlatCoeffsRatePt";

    try
    {
        double flatRate = getFlatFwdFromDiscountHoliday(startDate,
                                                        startDisc,
                                                        rp.date,
                                                        rp.disc,
                                                        noAccrue);

        rp.coeff[0] = flatRate;
        rp.coeff[1] = 0.;
        rp.coeff[2] = 0.;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::nextFraDisc
(
    const FwdRateIvl& fra,  // (I) fra rate interval 
    DateTime&         date, // (O) date for next discount point
    double&           disc  // (O) next disc factor
)
const
{
    static const string method = "LMParabCurve::nextFraDisc";

    try
    {
        DateTime prevDate;
        double   prevDisc;
        int      prevIvl;
        int      midIvl;
        int      nextIvl;

        DateTime lastDate;
        double   lastDisc;
        double   stubDisc;

        double   leftRate;
        double   midRate;
        double   rightRate;

        lastDate = !ratePts.empty() ? ratePts.back().date : baseDate;
        lastDisc = !ratePts.empty() ? ratePts.back().disc : 1.;

        if (fra.matDate <= lastDate)
        {
            throw ModelException("Fra maturity expected to come after curve end.");
        }

        if (fra.startDate == lastDate)
        {
            date = fra.matDate;
            disc = lastDisc*exp(-fra.rate);
            return;
        }

        if (fra.startDate < lastDate)
        {
            stubDisc = this->disc(fra.startDate, lastDate);

            date = fra.matDate;
            disc = lastDisc*exp(-fra.rate)/stubDisc;
            return;
        }

        // For a gap in the curve must find an inverse weighted average
        // of the last interval in the curve and the fra
        if (lastDate == baseDate)
        {
            throw ModelException("First fra rate must start at BaseDate.");
        }

        if (ratePts.size() == 1)
        {
            prevDate = baseDate;
            prevDisc = 1.;
        }
        else 
        {
            prevDate = ratePts[ratePts.size() - 2].date;
            prevDisc = ratePts[ratePts.size() - 2].disc;
        }

        prevIvl = noAccrue->businessDaysDiff(prevDate,      lastDate     );
        midIvl  = noAccrue->businessDaysDiff(lastDate,      fra.startDate);
        nextIvl = noAccrue->businessDaysDiff(fra.startDate, fra.matDate  );

        leftRate = getFlatFwdFromDiscountHoliday(prevDate,
                                                 prevDisc,
                                                 lastDate,
                                                 lastDisc,
                                                 noAccrue);

        rightRate = fra.rate/nextIvl;
        midRate   = (leftRate*nextIvl + rightRate*prevIvl)/(nextIvl + prevIvl);

        date = fra.startDate;
        disc = lastDisc*exp(-midRate*midIvl);
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::addYE
(
    const vector<YearEndEffect>& yel,
    HolidayConstSP               noAccrue
)
{
    static const string method = "LMParabCurve::addYE";

    try
    {
        size_t     idxCve;
        size_t     idxYE;
        int        daysIvl;

        DateTime   currDate;

        double     cumYeDisc;

        RatePt     newPt;
        RatePt     currPt;

        if (yel.empty())
            return;

        LMParabCurve curve = *this;
        ratePts.clear();

        cumYeDisc = 1.;

        // Find the first relevant YE interval
        for (idxYE = 0; idxYE < yel.size(); ++idxYE)
        {
            if (yel[idxYE].ivl.matDate > curve.baseDate)
                break;
        }

        currDate = curve.baseDate;

        for (idxCve = 0; idxCve < curve.ratePts.size(); ++idxCve)
        {
            currPt = curve.ratePts[idxCve];

            while (idxYE < yel.size() && yel[idxYE].ivl.matDate <= currPt.date)
            {
                if (yel[idxYE].ivl.startDate == currDate)
                {
                     newPt.date = yel[idxYE].ivl.matDate;

                     cumYeDisc*= exp(-yel[idxYE].ivl.rate);

                     daysIvl = noAccrue->businessDaysDiff(currDate,
                                                          yel[idxYE].ivl.matDate);
                     if (daysIvl == 0)
                     {
                         throw ModelException(
                             "Year End Interval corresponds to period of no Accrual.");
                     }

                     newPt.disc = cumYeDisc*disc(baseDate, yel[idxYE].ivl.matDate);

                     if (!yel[idxYE].isRateSpread)
                     {
                         newPt.coeff[0] = yel[idxYE].ivl.rate/daysIvl;
                         newPt.coeff[1] = 0.;
                         newPt.coeff[2] = 0.;
                     }
                     else
                     {
                         newPt.coeff[0] = currPt.coeff[0] + yel[idxYE].ivl.rate/daysIvl;
                         newPt.coeff[1] = currPt.coeff[1];
                         newPt.coeff[2] = currPt.coeff[2];

                         currPt.coeff[0]+= currPt.coeff[1]*daysIvl +
                                           currPt.coeff[2]*daysIvl*daysIvl;
                         currPt.coeff[1]+= 2.*currPt.coeff[2]*daysIvl;
                     }

                     currDate = yel[idxYE].ivl.matDate;

                     ratePts.push_back(newPt);
                     ++idxYE;
                     continue;
                }
                else // Start of interval is in the future
                {
                    newPt.date = yel[idxYE].ivl.startDate;

                    daysIvl = noAccrue->businessDaysDiff(currDate,
                                                         yel[idxYE].ivl.startDate);

                    newPt.disc = cumYeDisc*disc(baseDate, yel[idxYE].ivl.startDate);

                    newPt.coeff[0] = currPt.coeff[0];
                    newPt.coeff[1] = currPt.coeff[1];
                    newPt.coeff[2] = currPt.coeff[2];

                    currPt.coeff[0]+= currPt.coeff[1]*daysIvl +
                                      currPt.coeff[2]*daysIvl*daysIvl;
                    currPt.coeff[1]+= 2.*currPt.coeff[2]*daysIvl;

                    currDate = yel[idxYE].ivl.startDate;

                    ratePts.push_back(newPt);
                    ++idxYE;
                    continue;
                }
            }

            // All the relevant YE intervals have now been introduced
            if (currDate == currPt.date)
                continue;

            newPt = currPt;
            newPt.disc*= cumYeDisc;
            currDate = currPt.date;

            ratePts.push_back(newPt);
        }

        areCoeffsValid = true;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::genFromFras
(
    const DateTime&              bd,               // (I) base date for curve
    const vector<FwdRateIvl>&    frl,              // (I)
    const vector<YearEndEffect>& yel,              // (I)
    HolidayConstSP               na,               // (I) dates with no accrual
    const vector<DateTime>&      flatDates,        // (I) the rate will be constant in between flatRateDates,
                                                   //     last date must be <= flatParabBdryDate
    const DateTime&              flatParabBdryDate // (I) Date at which interpolations transitions
                                                   //     from flat to parabolic
)
{
    static const string method = "LMParabCurve::genFromFras";

    try
    {
        HolidaySP noAccrueYE;

        vector<FwdRateIvl> adjFraList;

        size_t numFlatDates;
        size_t idxFlatDate;

        size_t idxFra;

        DateTime lastDate;
        double   lastInstRate;
        double   lastDisc;

        ratePts.clear();

        areCoeffsValid = false;
        baseDate       = bd;
        noAccrue       = na;

        // If curve is empty
        if (frl.empty())
        {
            return;
        }

        // Create a merged list of noAccrue with the
        // non spread year end effect dates
        yearEndAndNoAccrueCombine(yel, noAccrue, noAccrueYE);

        // Get the normalised fra list
        fraListAdjustYearEnd(yel, frl, adjFraList);

        // Check if start of first fra is after baseDate
        if (baseDate > adjFraList[0].startDate)
        {
            throw ModelException("Base date comes after first fra date.");
        }

        if (baseDate != adjFraList[0].startDate)
        {
            throw ModelException("First fra date must be baseDate.");
        }

        // In case of only one rate create a flat curve
        if (adjFraList.size() == 1)
        {
            noAccrue = noAccrueYE;

            RatePt rp;
            rp.date = adjFraList[0].matDate;
            rp.disc = exp(-adjFraList[0].rate);

            rp.coeff[0] = getFlatFwdFromDiscountHoliday(baseDate,
                                                        1.,
                                                        rp.date,
                                                        rp.disc,
                                                        noAccrueYE);
            rp.coeff[1] = 0.;
            rp.coeff[2] = 0.;
            
            ratePts.push_back(rp);
            areCoeffsValid = true;

            addYE(yel, noAccrue);

            return;
        }

        numFlatDates = flatDates.size();

        if (numFlatDates > 0 && flatDates.back() > flatParabBdryDate)
        {
            throw ModelException(
                "The Flat Parabolic Boundary Date must be greater than all flat dates.");
        }

        // Disregard all flat dates after last fra date
        for (idxFlatDate = numFlatDates - 1; idxFlatDate > 0; --idxFlatDate)
        {
            if (flatDates[idxFlatDate] < adjFraList.back().matDate)
                break;
        }
        numFlatDates = idxFlatDate + 1;

        // Find first relevant flat date
        for (idxFlatDate = 0; idxFlatDate < numFlatDates; ++idxFlatDate)
        {
            if (flatDates[idxFlatDate] > baseDate)
                break;
        }

        noAccrue = noAccrueYE;

        areCoeffsValid = true;

        lastDate = baseDate;
        lastDisc = 1.;
        idxFra = 0;

        // Solve the flat dates first
        for (; idxFlatDate < numFlatDates; ++idxFlatDate)
        {
            addNextFlatPeriod(flatDates,
                              idxFlatDate,
                              numFlatDates,
                              lastDate,
                              lastDisc,
                              adjFraList,
                              idxFra       );
        }

        // Solve Flat until flatParabBdryDate
        for (; idxFra < adjFraList.size(); ++idxFra)
        {
            if (lastDate >= flatParabBdryDate)
                break;

            addNextFraFlat(lastDate,
                           lastDisc,
                           adjFraList,
                           idxFra     );
        }

        // The rest is parabolic
        if (idxFra == 0)
        {
            // First interval must be linear for a fully parabolic curve
            addFirstFraLinear(adjFraList,
                              lastInstRate);
            
            lastDate = ratePts.front().date;
            lastDisc = ratePts.front().disc;

            idxFra = 1;
        }
        else if (idxFra < adjFraList.size())
        {
            lastInstRate = ratePts.back().coeff[0];
            addBdryFraFlatParab((ratePts.size() > 1) ? 
                                 ratePts[ratePts.size() - 2].date :
                                 baseDate,
                                 lastDate,
                                 lastDisc,
                                 lastInstRate,
                                 adjFraList,
                                 idxFra       );
            ++idxFra;
        }
        
        for(; idxFra < adjFraList.size(); ++idxFra)
        {
            addNextFraParab(lastDate,
                            lastDisc,
                            lastInstRate,
                            adjFraList,
                            idxFra       );
        }

        addYE(yel, noAccrue);
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::addNextFlatPeriod
(
    const vector<DateTime>&   flatDates,    // (I) list of flat dates still to cover
    size_t&                   idxFlatDate,  // (I/O) current flat date index
    size_t                    numFlatDates, // (I) number of flat dates
    DateTime&                 lastFlatDate, // (I/O) last flat date covered
    double&                   lastDisc,     // (I/O) discount up to lastFlatDate
    const vector<FwdRateIvl>& fraList,      // (I) fras to be included
    size_t&                   idxFra        // (I/O) index of Fra being used currently
)
{
    static const string method = "LMParabCurve::addNextFlatPeriod";

    try
    {
        DateTime prevFlatDate;
        DateTime currFlatDate;
        DateTime currFraDate;
        DateTime nextFraDate;
        DateTime nextFraStart;

        double prevDisc;
        double currDisc;
        double currFraDisc;
        double nextFraDisc;
        double currRate;
        double nextRate;

        int currFraIvl;
        int currFlatIvl;
        int nextFraIvl;
        int overlapIvl;
        int auxMax;
        int auxMin;

        prevFlatDate = lastFlatDate;
        currFlatDate = prevFlatDate;
        prevDisc     = lastDisc;

        this->nextFraDisc(fraList[idxFra], currFraDate, currFraDisc);

        // if the Fra start is in the future, reuse the Fra
        if (currFraDate == fraList[idxFra].startDate)
            --idxFra;

        if (idxFra + 1 < fraList.size())
        {
            nextFraDate = fraList[idxFra + 1].startDate > currFraDate ?
                          fraList[idxFra + 1].startDate :
                          fraList[idxFra + 1].matDate;
        }
        else
        {
            nextFraDate = currFraDate;
        }

        // remove unnecessary flat dates
        for ( ; idxFlatDate < numFlatDates; ++idxFlatDate)
        {
            if (flatDates[idxFlatDate] == nextFraDate)
            {
                if (currFlatDate == prevFlatDate)
                {
                    throw ModelException(
                        "Flat Date Interval spans more than one Fra Interval.");
                }
                else
                {
                    break;
                }
            }

            currFlatDate = flatDates[idxFlatDate];
            if (currFlatDate >= currFraDate)
                break;
        }

        if (currFlatDate < currFraDate)
        {
            // if it is the last flat date treat it as if it were at the currFraDate
            if (idxFlatDate == numFlatDates)
            {
                currFlatDate = currFraDate;
            }
            // if there is a gap then remove it
            else if (currFlatDate < currFraDate && nextFraDate == fraList[idxFra + 1].startDate)
            {
                currFlatDate = nextFraDate;
            }
        }

        currFraIvl = noAccrue->businessDaysDiff(prevFlatDate, currFraDate);

        if (currFraIvl == 0)
        {
            throw ModelException("No Accrue Days on Interval.");
        }

        currRate = -log(currFraDisc/prevDisc)/currFraIvl;

        currFlatIvl = noAccrue->businessDaysDiff(prevFlatDate, currFlatDate);

        if (currFlatDate >= currFraDate)
        {
            currDisc = prevDisc*exp(-currFlatIvl*currRate);

            RatePt rp = {currFlatDate, currDisc};

            calcFlatCoeffsRatePt(rp,
                                 prevFlatDate,
                                 prevDisc,
                                 noAccrue);

            ratePts.push_back(rp);

            ++idxFra;
        }
        else
        {
            // The remaining case is where currFlatDate is before currFraDisc and
            // the following FlatDate[idxFlatDate+1] is nextFraDate
            ++idxFra;

            nextFraStart = max(fraList[idxFra].startDate, currFlatDate);

            nextFraIvl = noAccrue->businessDaysDiff(nextFraStart, fraList[idxFra].matDate);

            if (nextFraIvl == 0)
            {
                throw ModelException("No Accrue Days on Interval.");
            }

            overlapIvl = noAccrue->businessDaysDiff(nextFraStart, currFraDate);

            // get currRate and nextRate as total areas
            currRate*= currFraIvl;

            if (nextFraStart > currFlatDate)
            {
                nextRate = fraList[idxFra].rate;
            }
            else
            {
                this->nextFraDisc(fraList[idxFra], nextFraDate, nextFraDisc);

                nextRate = -log(nextFraDisc/prevDisc);
            }

            // Now get the rate (not area) for the current flat period
            auxMax = max(currFlatIvl + overlapIvl - currFraIvl, 0);
            auxMin = nextFraIvl - auxMax;

            currDisc = exp(-currFlatIvl*
                           (currRate*auxMin - nextRate*(currFraIvl - currFlatIvl))/
                           (currFlatIvl*auxMin - (currFraIvl - currFlatIvl)*auxMax));

            RatePt rp1 = {currFlatDate, currDisc};

            calcFlatCoeffsRatePt(rp1,
                                 prevFlatDate,
                                 prevDisc,
                                 noAccrue     );

            ratePts.push_back(rp1);

            currRate = (currRate*auxMax - nextRate*currFlatIvl)/
                       ((currFraIvl - currFlatIvl)*auxMax - currFlatIvl*auxMin);

            prevFlatDate = currFlatDate;
            prevDisc     = currDisc;

            currFlatDate = nextFraDate;
            currFlatIvl  = nextFraIvl - currFlatIvl - overlapIvl + currFraIvl;
            currDisc     = exp (-currRate*currFlatIvl);

            RatePt rp2 = {currFlatDate, currDisc};

            calcFlatCoeffsRatePt(rp2,
                                 prevFlatDate,
                                 prevDisc,
                                 noAccrue     );

            ratePts.push_back(rp2);

            ++idxFra;
        }

        lastFlatDate = currFlatDate;
        lastDisc     = currDisc;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::addNextFraFlat
(
    DateTime&                 lastDate, // (I/O) lastdate enetered into curve
    double&                   lastDisc, // (I/O) discount up to lastDate
    const vector<FwdRateIvl>& fraList,  // (I) fras to be included
    size_t&                   idxFra    // (I/O) index of Fra being used currently
)
{
    static const string method = "LMParabCurve::addNextFraFlat";

    try
    {
        DateTime prevDate;
        double   prevDisc;
        RatePt   rp;

        prevDate = lastDate;
        prevDisc = lastDisc;

        nextFraDisc(fraList[idxFra],
                    rp.date,
                    rp.disc         );

        // If the Fra start is in the future, reuse the Fra
        if (rp.date == fraList[idxFra].startDate)
            --idxFra;

        calcFlatCoeffsRatePt(rp,
                             prevDate,
                             prevDisc,
                             noAccrue);

        ratePts.push_back(rp);

        lastDate = rp.date;
        lastDisc = rp.disc;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::addFirstFraLinear
(
    const vector<FwdRateIvl>& fraList,      // (I) fras to be included
    double&                   currInstRate  // (O) instantaneous rate at the end of first period
)
{
    static const string method = "LMParabCurve::addFirstFraLinear";

    try
    {
        DateTime prevDate;
        DateTime currDate;
        DateTime nextDate;

        double   prevDisc;
        double   currRate;
        double   currDisc;
        double   nextRate;
        double   nextDisc;

        double   prevInstRate;

        int      currIvl;
        int      nextIvl;

        RatePt   rp;

        prevDate = baseDate;
        prevDisc = 1.;

        currDate = fraList.front().matDate;

        currIvl = noAccrue->businessDaysDiff(prevDate, currDate);

        if (currIvl == 0)
        {
            throw ModelException("No Accrue Days on Interval.");
        }

        currRate = fraList.front().rate/currIvl;
        currDisc = exp(-fraList.front().rate);

        rp.date = currDate;
        rp.disc = currDisc;

        calcFlatCoeffsRatePt(rp,
                             prevDate,
                             prevDisc,
                             noAccrue);

        ratePts.push_back(rp);

        nextFraDisc(fraList[1],
                    nextDate,
                    nextDisc   );

        nextIvl = noAccrue->businessDaysDiff(currDate, nextDate);

        nextRate = -log(nextDisc/currDisc)/nextIvl;

        // Need to adjust the first interval to be linear
        prevInstRate = 2*currRate -
                       (currRate*nextIvl + nextRate*currIvl)/
                       (currIvl + nextIvl);

        calcLinearCoeffsRatePt(ratePts.front(),
                               prevDate,
                               prevDisc,
                               prevInstRate,
                               currInstRate,
                               noAccrue        );
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::addNextFraParab
    (
        DateTime&                 lastDate,     // (I/O) last date entered into curve
        double&                   lastDisc,     // (I/O) discount up to lastDate
        double&                   lastInstRate, // (I/O) instantaneous forward rate at lastDate
        const vector<FwdRateIvl>& fraList,      // (I) fras to be included
        size_t&                   idxFra        // (I/O) index of Fra being used currently
    )
{
    static const string method = "LMParabCurve::addNextFraParab";

    try
    {
         DateTime prevDate;
         DateTime currDate;
         DateTime nextDate;

         double prevDisc;
         double currDisc;
         double nextDisc;

         double prevInstRate;

         int    currIvl;

         RatePt rp;

         prevDate = lastDate;
         prevDisc = lastDisc;

         prevInstRate = lastInstRate;

         nextFraDisc(fraList[idxFra],
                     currDate,
                     currDisc        );

         // If there was a gap the fra should be reused for another interval
         if (currDate == fraList[idxFra].startDate)
             --idxFra;

         currIvl = noAccrue->businessDaysDiff(prevDate, currDate);

         if (currIvl == 0)
         {
             throw ModelException("No Accrue Days on Interval.");
         }

         // Add this point as a flat forwards point
         rp.date = currDate;
         rp.disc = currDisc;

         calcFlatCoeffsRatePt(rp,
                              prevDate,
                              prevDisc,
                              noAccrue );

         ratePts.push_back(rp);

         if (size_t(idxFra + 1) < fraList.size())
         {
             nextFraDisc(fraList[idxFra + 1],
                         nextDate,
                         nextDisc            );
         }
         else
         {
             nextDate = currDate;
             nextDisc = currDisc;
         }

         calcParabolicCoeffsRatePt(rp,
                                   prevDate,
                                   prevDisc,
                                   prevInstRate,
                                   nextDate,
                                   nextDisc,
                                   lastInstRate,
                                   noAccrue     );

         ratePts.back() = rp;

         lastDate = currDate;
         lastDisc = currDisc;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void LMParabCurve::addBdryFraFlatParab
(
    const DateTime&           prevStartDate, // (I) date at which the flat period started
    DateTime&                 lastDate,      // (I/O) last date entered into curve
    double&                   lastDisc,      // (I/O) discount up to lastDate 
    double&                   lastInstRate,  // (I/O) instantaneous forward rate at lastDate
    const vector<FwdRateIvl>& fraList,       // (I) fras to be included
    size_t&                   idxFra         // (I/O) index of Fra being used currently
)
{
    static const string method = "LMParabCurve::addBdryFraFlatParab";

    try
    {
        DateTime prevDate;
        DateTime currDate;

        double   prevDisc;
        double   currDisc;

        double   prevFlatRate;
        double   currFlatRate;

        int      prevDaysIvl;
        int      currDaysIvl;

        prevDate = lastDate;
        prevDisc = lastDisc;

        prevFlatRate = lastInstRate;

        nextFraDisc(fraList[idxFra],
                    currDate,
                    currDisc        );

        currFlatRate = getFlatFwdFromDiscountHoliday(prevDate,
                                                     prevDisc,
                                                     currDate,
                                                     currDisc,
                                                     noAccrue);

        prevDaysIvl = noAccrue->businessDaysDiff(prevStartDate, prevDate);
        currDaysIvl = noAccrue->businessDaysDiff(prevDate     , currDate);

        lastInstRate = (prevFlatRate*currDaysIvl + currFlatRate*prevDaysIvl)/
                       (currDaysIvl + prevDaysIvl);

        addNextFraParab(lastDate,
                        lastDisc,
                        lastInstRate,
                        fraList,
                        idxFra       );
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


double LMParabCurve::getFlatFwdFromDiscountHoliday
(
    const DateTime& leftDate,  // (I)
    double          leftRate,  // (I)
    const DateTime& rightDate, // (I)
    double          rightRate, // (I)
    HolidayConstSP  holiday    // (I) used for time between 1st, 2nd dates
)
{
    static const string method = "LMParabCurve::getFlatFwdFromDiscountHoliday";

    try
    {
        int daysDiff = holiday->businessDaysDiff(leftDate, rightDate);

        if (daysDiff == 0)
        {
            char buffer[256];
            sprintf(buffer, "No business days between %s and %s.",
                leftDate .toString().c_str(),
                rightDate.toString().c_str());
            throw ModelException(buffer);
        }

        if (Maths::isZero(leftRate) ||
            Maths::isZero(rightRate/leftRate) ||
            rightRate / leftRate < 0.)
        {
            throw ModelException(
                "Discount Factors are negative or out of range.");
        }

        return -log(rightRate/leftRate)/daysDiff;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

DRLIB_END_NAMESPACE
