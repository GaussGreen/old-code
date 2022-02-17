//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : InstrumentUtil.hpp
//
//   Description : Utility functions for instruments
//
//   Author      : André Segger
//
//   Date        : 04 Jun 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Algorithm.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/StruckEquity.hpp"

DRLIB_BEGIN_NAMESPACE

/** outputs all forwards at maturity */
void InstrumentUtil::recordFwdAtMat(
                           Control*         control,
                           Results*         results,
                           const DateTime&  matDate,
                           const DateTime&  valueDate,
                           const CAsset*    asset)
{
    OutputRequest* request = NULL;
    if ( matDate.isGreaterOrEqual(valueDate) )
    {
        if ( control->requestsOutput(OutputRequest::FWD_AT_MAT, request))
        {
            /* do fwd @ mat by stock */
            asset->recordFwdAtMat(request, results, matDate);
        }
    }
}

/** outputs indicative volatility */
void InstrumentUtil::recordIndicativeVol(
                         Control*     control,
                         Results*     results,
                         double       indVol)
{
    OutputRequest* request = NULL;

    if (control->requestsOutput(OutputRequest::IND_VOL, request)) {
        results->storeRequestResult(request, indVol);
    }
}

/** outputs delay Price*/
void InstrumentUtil::delayPriceHelper(
    Control*                    control,
    Results*                    results,
    const double&               fairValue,
    const DateTime&             valueDate,
    const YieldCurve*           discount,
    const Asset*                asset,
    const InstrumentSettlement* premiumSettlement)
{
    if (premiumSettlement) {
        OutputRequest* req=control->requestsOutput(OutputRequest::DELAY_PRICE);
        if (req) {
            try {
                if (premiumSettlement) {
                    DateTime paymentDate = 
                        premiumSettlement->settles(valueDate, asset);
                    double delayPrice = InstrumentUtil::calculateDelayPrice(
                        fairValue,
                        valueDate,
                        paymentDate,
                        discount);
                    results->storeRequestResult(req, delayPrice);
                } 
            }
            catch (exception& e) {
                results->
                    storeRequestResult(req, IObjectSP(new Untweakable(e)));
            }
        }
   }  
}

void InstrumentUtil::recordDelayPrice(Control*     control,
                                      Results*     results,
                                      double       delayPrice)
{
    OutputRequest* request = NULL;

    if ((request = control->requestsOutput(OutputRequest::DELAY_PRICE))) {
        results->storeRequestResult(request, delayPrice);
    }
}

void InstrumentUtil::recordPricePCTpayoff(Control*     control,
                                          Results*     results,
                                          double       pricePCTpayoff)
{
    OutputRequest* request = NULL;

    if ((request = control->requestsOutput(OutputRequest::PRICE_PCT_PAYOFF))) {
        results->storeRequestResult(request, pricePCTpayoff);
    }
}

void InstrumentUtil::recordMaxPayoff(Control*     control,
                                     Results*     results,
                                     double       maxPayoff)
{
    OutputRequest* request = NULL;

    if ((request = control->requestsOutput(OutputRequest::MAX_PAYOFF))) {
        results->storeRequestResult(request, maxPayoff);
    }
}


double InstrumentUtil::calculateDelayPrice(const double&     fairValue,
                                           const DateTime&   valueDate,
                                           const DateTime&   payDate,
                                           const YieldCurve* discount)
{
    double discFactor, delayPrice;
    if (payDate > valueDate)
    {
        discFactor = discount->pv(valueDate, payDate);
    }
    else
    {
        discFactor = 1.0;
    }

    if ( !Maths::isZero(discFactor) ) {
        delayPrice = fairValue/discFactor;
    } else {
        delayPrice = 0.0;
    }

    return delayPrice;
}

double InstrumentUtil::scalePremium(bool          oneContract,
                                    bool          fwdStarting,
                                    const double& notional,
                                    const double& fwdAtStart,
                                    const double& initialSpot)
{
    double scalingFactor = 1.0;
    if (!oneContract)
    {
        /* handle notional */
        if (fwdStarting)
        {
            if ( Maths::isZero(fwdAtStart) )
            {
                throw ModelException("InstrumentUtil::scalePremium", 
                                     "Forward at start is 0.0. Infinite premium.");
            }
            scalingFactor = notional/fwdAtStart;
        }
        else
        {
            if ( Maths::isZero(initialSpot) )
            {
                throw ModelException("InstrumentUtil::scalePremium",
                                     "initial Spot is 0.0. Infinite number of contracts.");
            }
            /* handle position */
            scalingFactor = notional/initialSpot;
        }
    }
    return scalingFactor;
}

/** Sorts the supplied performances and then takes the weight sum of them
    using supplied weights. weights.size() must equal perf.size(). The
    sorted performances are returned in perf.  */
double InstrumentUtil::rainbowPerformance(const DoubleArray& weights,
                                          DoubleArray&       perf){
    Algorithm::shellSort(perf);
    // now form the rainbow basket
    return weightedPerformance(weights, perf);
}
    

/** Takes the weighted sum of the supplied performances.
    weights.size() must equal perf.size().  */
double InstrumentUtil::weightedPerformance(const DoubleArray& weights,
                                           const DoubleArray& perf){
    double bskPerf = 0.0;
    for (int iAsset = 0; iAsset < weights.size(); iAsset++) {
        bskPerf += weights[iAsset] * perf[iAsset];
    }
    return bskPerf;
}

DRLIB_END_NAMESPACE
