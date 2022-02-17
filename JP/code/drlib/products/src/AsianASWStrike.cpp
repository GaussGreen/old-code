//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//   Filename    : AssetSwapStrike.cpp
//
//   Description : Base class for asset swap strike objects 
//
//   Author      : André Segger
//
//   Date        : 02 June 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/AsianASWStrike.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/PayStream.hpp"


DRLIB_BEGIN_NAMESPACE

/** calculate the strike of the asset swap */
double AsianASWStrike::getStrike(const DateTime&    baseDate,
                                 ConvBondConstSP    cvb,
                                 YieldCurveConstSP  discountCurve) const
{
    static const string method = "AsianASWStrike::getStrike";

    try {
        int i;
        DayCountConventionSP dcc       = DayCountConventionSP(DayCountConventionFactory::make(swapDCC));
        BadDayConventionSP   bdc       = BadDayConventionSP(BadDayConventionFactory::make(rateBadDayConv));
        MaturityPeriodSP     matPeriod(new MaturityPeriod(swapInterval));

        DoubleArray notionals(refixDates->size());
        DoubleArray spreads(refixDates->size());
        DoubleArray weights(refixDates->size());
        BoolArray   compound(refixDates->size());
        for (i=0 ; i < refixDates->size() ; ++i) {
            notionals[i] = swapNotional;
            spreads[i]   = spread;
            weights[i]   = 1.0;
            compound[i]  = false;
        }

        FloatRateSP floatRate (new FloatRate(
                    discountCurve.get(),
                    hol.get(),
                    dcc.get(),
                    bdc.get(),
                    matPeriod.get(),
                    matPeriod.get(),
                    false));
        /*

        // create a payment stream
        PayStreamSP payStream(new PayStream(
                  dcc.get(),
                  true,         // compound flat
                  floatRate.get(),
                  notionals,
                  *(refixDates.get()),
                  *(accrualDates.get()),
                  *(payDates.get()),
                  *(fixings.get()),
                  spreads,
                  weights,
                  compound));


        DateTime swapMaturity = (*(payDates)[payDates->size()-1]);
        // see if we're passed swap maturity
        if (baseDate > swapMaturity) {
            strike = 0.;
        } else {
            double exerPct;
            // see if it's exercisable today
            try { // need a try block as an error is thrown if it's not exercisable
                exerPct = exerSched->interpolate(baseDate);
            } catch (exception& ) {
                // we should still calculate the strike as if it were exercisable today for info purposes
                // use the exer level on the first day
                exerPct = exerSched->interpolate(exerSched->firstDate()); 
            }    

            double floatLeg = getFloatLegPV(baseDate);
            double fixedLeg = cvb->bond->couponsPV(baseDate, swapMaturity, cvb->discount.getSP());
            
            fixedLeg += swapBackEndFee*cvb->discount->pv(baseDate, swapMaturity);

            *strike = exerPct*swapNotional + floatLeg + fixedLeg;

            // calculate break of funds for Credit Suisse style swaps
            if ( ocbType == "C" ) {
                double currentRate         = 0.0;
                double currentDayCountFrac = 0.0;
                int count;
                string interval;
                swapInterval->decompose(count, interval);

                CashFlowArray spreadFlows = SwapTool::cashflows(
                    baseDate.rollDate(-1),
                    swapMaturity,
                    !shortFrontStub,
                    -spread,
                    count,           // interval = count periods
                    interval,          // e.g. Y, M, W, D
                    swapDCC.get());

                if ( baseDate < swapLastFixDate ) {
                    throw ModelException(method, "last swap fixing date is in the future");
                } else if ( spreadFlows.size() > 0 && spreadFlows[0].date == baseDate) {
                    currentRate         = 0.0;
                    currentDayCountFrac = 0.0;
                } else if ( spreadFlows.size() > 0 ){
                    currentRate         = cvb->discount->fwd(baseDate,
                                                             spreadFlows[0].date,
                                                             swapDCC.get(),
                                                             CompoundBasis::SIMPLE);
                    currentDayCountFrac = swapDCC->years(baseDate, spreadFlows[0].date);
                } else {
                    currentRate         = 0.0;
                    currentDayCountFrac = 0.0;
                }

                double breakOfFunds;
                if (swapLastFixRate > currentRate) {
                    breakOfFunds = (swapLastFixRate - currentRate + swapBreakOfFundsRate) * 
                                   currentDayCountFrac                                    *
                                   swapNotional;
                }

                *strike += breakOfFunds;
            }
        }
            */
    } catch (exception& e) {
        throw ModelException(e, method);
    }
    return 0.0;
}

/** returns true if the asset swap is exercisable on the baseDate */
bool    AsianASWStrike::isExercisable(const DateTime& baseDate)
{
    return true;
}

void    AsianASWStrike::validatePop2Object()
{
    static const string method = "AsianASWStrike::validatePop2Object";

    int numFixings = refixDates->size();
    if ( accrualDates->size() != numFixings ) {
        throw ModelException(method, "Must have same number of fixings(" + Format::toString(numFixings) +
                        ") as number of accrual dates (" + Format::toString(accrualDates->size()) + ")");
    }
    if ( payDates->size() != numFixings ) {
        throw ModelException(method, "Must have same number of fixings(" + Format::toString(numFixings) +
                        ") as number of accrual dates (" + Format::toString(payDates->size()) + ")");
    }
    if ( fixings->size() != numFixings ) {
        throw ModelException(method, "Must have same number of fixings(" + Format::toString(numFixings) +
                        ") as number of accrual dates (" + Format::toString(fixings->size()) + ")");
    }
}


class AsianASWStrikeHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(AsianASWStrike, clazz);
        SUPERCLASS(AssetSwapStrike);

        FIELD(valueDate,            "current value date");
        FIELD_MAKE_OPTIONAL(valueDate);

        FIELD(refixDates,                  "refix dates for Libor leg");
        FIELD(accrualDates,                "accrual dates for Libor leg");
        FIELD(payDates,                    "payment dates for Libor leg");
        FIELD(fixings,                     "Libor fixings");

        FIELD(spread,               "Libor spread");
        FIELD(swapNotional,         "Notional of Swap");
        FIELD(exerSched,                   "exercise Schedule");
        FIELD(swapDCC,              "Day count convention");
        FIELD(swapBackEndFee,       "Back end fee");

        // FIELD(payStream,                   "internal field");
        // FIELD_MAKE_TRANSIENT(payStream);
        // FIELD_MAKE_TWEAKABLE(payStream);
    }
};

CClassConstSP const AsianASWStrike::TYPE = CClass::registerClassLoadMethod(
    "AsianASWStrike", typeid(AsianASWStrike), AsianASWStrikeHelper::load);
bool  AsianASWStrikeLoad() {
    return (AsianASWStrike::TYPE != 0);
   }


DRLIB_END_NAMESPACE
