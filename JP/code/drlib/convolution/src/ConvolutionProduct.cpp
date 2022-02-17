//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Description : View onto instrument as required by ConvolutionEngine
//
//   Date        : 18th Nov 2005
//
//   Author      : Mark Robson
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/ConvolutionEngine.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/ConvolutionModel.hpp"
#include "edginc/CounterPartyCredit.hpp"
#include "edginc/FixedTrancheLossCalculator.hpp"
#include "edginc/ExpiryResult.hpp"
#include "edginc/Control.hpp"
#include "edginc/Results.hpp"
#include "edginc/EffectiveCurve.hpp"
#include "edginc/IEffectiveCurveLossGen.hpp"
#include "edginc/IModelConfigMapper.hpp"
#include "edginc/ConvolutionProduct.hpp"
#include "edginc/IForwardRatePricer.hpp"

DRLIB_BEGIN_NAMESPACE

ConvolutionProduct::Output::~Output(){}
ConvolutionProduct::Output::Output(){}

const string ConvolutionProduct::CC_CONSERVATIVE("CONSERVATIVE");
const string ConvolutionProduct::CC_AGGRESSIVE("AGGRESSIVE");

ConvolutionProduct::~ConvolutionProduct(){}

//// Constructor just takes references
ConvolutionProduct::ConvolutionProduct(
        ConvolutionEngineConstSP convolutionEngine,
        YieldCurveConstSP        discount) :
    convolutionEngine(convolutionEngine), discount(discount)
{}
    

/** Returns today aka value date. Note not yield curve spot date */
DateTime ConvolutionProduct::getToday() const{
    return convolutionEngine->getValueDate();
}

/** Returns the [riskless] curve for discounting (could remove this if
    we had better methods on ICDSParSpreads - only used for calculating
    durations) */
YieldCurveConstSP ConvolutionProduct::getDiscount() const{
    return discount;
}

/** Kapital backward compatibility mode - to be removed eventually.
    Returns whether we pv to cfCutOffDate() */
bool ConvolutionProduct::pvToSpot() const{
    return convolutionEngine->pvToSpot();
}
    
/** Kapital backward compatibility mode - to be removed eventually.
    Ignore cashflows as well as protection etc before this date. */
DateTime ConvolutionProduct::cfCutOffDate() const {
    const DateTime& cfCutOffDate = convolutionEngine->cfCutOffDate();
    if (cfCutOffDate.empty()) {
        return pvToSpot()? discount->getSpotDate(): getToday();
    } 
    else if (getToday() > cfCutOffDate) {
        throw ModelException("ConvolutionProduct::cfCutOffDate",
                             "Tranche CF value date "+
                             cfCutOffDate.toString()+" < today "+
                             getToday().toString());
    }
    return cfCutOffDate;
}


/** Returns highStrike-lowStrike with highStrike and lowStrike obtained from
    trancheStrikes */
double ConvolutionProduct::trancheNotional() const{
    double lowStrike;
    double highStrike;
    trancheStrikes(lowStrike, highStrike);
    return (highStrike-lowStrike);
}

/** Returns the tranche notional currently at risk, either due to
    a) short losses, in which case the outstanding notional will increase
    b) long losses, in which case the outstanding notional will decrease.
    Risky notional is equivalent to OutstandingNotional in case there
    is no short. */
double ConvolutionProduct::trancheRiskyNotional(const bool recoverNotional) const
{
    double lowStrike;
    double highStrike;
    trancheStrikes(lowStrike, highStrike);
    double longLoss = portfolioLongLoss();
    double shortNotional = portfolioShortNotional();
    double minLoss = longLoss - shortNotional;

    if (recoverNotional)
    {
        double topLoss = 
            Maths::shiftedCollar(-portfolioLongNotional()+ portfolioLongRecoveredNotional(),
                                 -highStrike,
                                 -lowStrike);
        return highStrike - Maths::max(lowStrike,0.0) - topLoss
                          - Maths::creditCollar(minLoss, lowStrike, highStrike);
    }
    else
    {
        return highStrike - Maths::max(lowStrike,0.0) - Maths::creditCollar(minLoss, lowStrike, highStrike);
    }
}

// jlhp this is now in CreditTrancheLossConfig, can it be removed here? Note it is called from calculateEffectiveCurve
/** Returns tranche outstanding notional as seen today */
double ConvolutionProduct::trancheOutstandingNotional(const bool recoverNotional) const
{
    double lowStrike;
    double highStrike;
    trancheStrikes(lowStrike, highStrike);
    double bottomLoss = Maths::creditCollar(portfolioLoss(), lowStrike, highStrike);
        if (recoverNotional) {
        // need to consider recovered notional separately from losses to the
        // bottom of the tranche. Normally would write 
        // collar(recoveredNotional, 1-high strike, 1-low strike) 
        // where 1 = portfolio notional. But collar is invariant under 
        // a transformation where a constant is added to each parameter.
        // We can use shiftedCollar rather than creditCollar since  
        // "1-high strike" etc will be positive provided we enfore 
        // high strike <= portfolio notional
        double topLoss = 
            Maths::shiftedCollar(-portfolioLongNotional()+ portfolioLongRecoveredNotional(),
                                 -highStrike,
                                 -lowStrike);
        return highStrike - Maths::max(lowStrike,0.0) - bottomLoss - topLoss;
    } 
    else {
        return highStrike - Maths::max(lowStrike,0.0) - bottomLoss;
    }
}

/** Utility method: calculates counter party survival probabilities along
    timeline. Resizes counterPartyProb as required. If there is no counterparty
    information (ie getCounterParty() returns null) then does nothing */
void ConvolutionProduct::computeCounterPartySurvivalProb(
    const DateTimeArray& timeline,               // (I)
    DoubleArray&         counterPartyProb) const // (O) [time index]
{
    try {
        CounterPartyCreditConstSP counterParty(getCounterParty());
        if (counterParty.get()) {
            counterPartyProb.resize(timeline.size()); // reserve space
            // get hold of relevant market data
            SingleCreditAssetConstSP asset(counterParty->getAsset());
            DefaultRatesSP defaultRates(asset->
                                        getParSpreadCurve()->defaultRates());
            const DateTime& today = getToday(); /* should have method on
                                                   DefaultRates to avoid this */
            for (int t = 0; t < timeline.size(); t++) {
                counterPartyProb[t]=defaultRates->calcDefaultPV(today,
                                                                timeline[t]);
            }
        }
    } 
    catch (exception& e) {
        throw ModelException(e, "ConvolutionProduct::"
                             "computeCounterPartySurvivalProb");
    }
}

/** Computes the effective curve along the specified timeline
    storing the result in effectiveCurve. If computeCondCurve is
    true then an effective curve conditional on the counterparty
    not defaulting should be returned too (in this case
    ConvolutionProduct::getCounterParty() will not return null).
    The timeline supplied here is the same as that generated by 
    generateTimeline, which will have been called first. 
    The DoubleArraySP will initially be null. */
void ConvolutionProduct::calculateEffectiveCurve(
    IFixedTrancheLossCalculatorConstSP  lossCalculator,              /* (I) */
    IFixedTrancheLossCalculatorConstSP  recoveredNotionalCalculator, /* (I, may be null) */
    const DateTimeArray&                timeline,                    /* (I) */
    DoubleArraySP&                      ctgEffectiveCurve,           /* (O) */
    DoubleArraySP&                      feeEffectiveCurve,           /* (O) */
    bool                                computeCondCurve,            /* (I) */
    DoubleArraySP&                      ctgEffectiveCurveCond,       /* (O) */ 
    DoubleArraySP&                      feeEffectiveCurveCond)       /* (O) */ const
{
    static const string method("ConvolutionProduct::calculateEffectiveCurve");
    try {
        int numTimePoints = timeline.size(); // for ease
        ctgEffectiveCurve.reset(new DoubleArray(numTimePoints));
        feeEffectiveCurve.reset(new DoubleArray(numTimePoints));

        //initialise as empty arrays
        ctgEffectiveCurveCond.reset(new DoubleArray(numTimePoints));
        feeEffectiveCurveCond.reset(new DoubleArray(numTimePoints));

        // get information on the tranche
        bool recoverNotional = !!recoveredNotionalCalculator.get();

        double ctgLegOutstandingNotional = 
            this->trancheOutstandingNotional(false);
        double feeLegOutstandingNotional = 
            this->trancheOutstandingNotional(recoverNotional);

        /* Below, the effective curve is expressed as percentage of
         * the initial outstanding notional. In the presence of short names
         * the outstanding notional can be 0 and the future expected notional
         * can be >0, in which case the clean spread is not defined.
         * For the moment, we will fail in that case.
        /* Treat t=today separately */
        if (    (ctgLegOutstandingNotional <= 0 && this->portfolioShortNotional() > 0)
             || (feeLegOutstandingNotional <= 0 && this->portfolioShortNotional() > 0))
        {
            // fail miserably
            throw ModelException("Outstanding notional is 0 and short names are present."
                                 "This is not currently handled.", method);
        }

        (*ctgEffectiveCurve)[0] = 1.0; 
        (*feeEffectiveCurve)[0] = 1.0;

        // exactly the same for cond curves
        if (computeCondCurve) {
            (*ctgEffectiveCurveCond)[0] = 1.0;
            (*feeEffectiveCurveCond)[0] = 1.0;
        }

        double lowStrike;
        double highStrike;
        trancheStrikes(lowStrike, highStrike);
        double tranchePastLoss = Maths::creditCollar(portfolioLoss(), 
                                                     lowStrike, 
                                                     highStrike);
        double topLoss = recoverNotional ?
            (Maths::shiftedCollar(-portfolioLongNotional()+portfolioLongRecoveredNotional(),
                                  -highStrike, 
                                  -lowStrike)) :
            0;
        
        /* then loop through points  of the timeline */
        for (int t = 1; t < numTimePoints; ++t) {
            double expectedLoss     = 0.0;
            double expectedLossCond = 0.0;

            if ((!Maths::isZero(ctgLegOutstandingNotional)) || 
                (!Maths::isZero(feeLegOutstandingNotional)))
            {
                //get expected loss
                lossCalculator->loss(t, expectedLoss, expectedLossCond);            
            }

            // deal with the contingent leg curve
            if (!Maths::isZero(ctgLegOutstandingNotional)) {
                double expectedOutstandingNtl =
                    highStrike - Maths::max(lowStrike,0.0) - expectedLoss;
                if (expectedOutstandingNtl < 0.0)
                {
                    throw ModelException(
                        method,
                        "Timepoint " + Format::toString(t) + ": "
                        "Expected loss (" +
                        Format::toString(expectedLoss) +
                        ") is greater than notional (" +
                        Format::toString(ctgLegOutstandingNotional + tranchePastLoss) +
                        ")");
                }
                
                (*ctgEffectiveCurve)[t] = expectedOutstandingNtl / ctgLegOutstandingNotional;

                if (computeCondCurve){
                    double expectedOutstandingNtlCond =
                        highStrike - Maths::max(lowStrike,0.0) - expectedLossCond;
                    if (expectedOutstandingNtlCond < 0.0)
                    {
                        throw ModelException(
                            method,
                            "Timepoint " + Format::toString(t) + ": "
                            "Conditional expected loss (" +
                            Format::toString(expectedLoss) +
                            ") is greater than notional (" +
                            Format::toString(ctgLegOutstandingNotional + tranchePastLoss) +
                            ")");
                    }
                    
                    (*ctgEffectiveCurveCond)[t] = 
                        expectedOutstandingNtlCond / ctgLegOutstandingNotional;
                }

            }
            else { /* outstanding notional is 0 */ 
                (*ctgEffectiveCurve)[t] = 1.0;
                if (computeCondCurve) {
                    (*ctgEffectiveCurveCond)[t] = 1.0;
                }               
            }
    
            // deal with the fee leg curve
            if (!Maths::isZero(feeLegOutstandingNotional)) {
                double expectedRecNtl = 0.0;
                double expectedRecNtlCond = 0.0;
                //get expected recovered notional if needed
                if (recoverNotional) {
                    recoveredNotionalCalculator->loss(t, expectedRecNtl,
                                                      expectedRecNtlCond);
                }

                (*feeEffectiveCurve)[t] = 
                    (highStrike - Maths::max(lowStrike,0.0) - expectedLoss - expectedRecNtl)
                        / feeLegOutstandingNotional;

                if (computeCondCurve) {
                    (*feeEffectiveCurveCond)[t] = 
                        (highStrike - Maths::max(lowStrike,0.0) - expectedLossCond - expectedRecNtlCond)
                        / feeLegOutstandingNotional;
                }
            }
            else { /* risky notional for the fee leg is 0 */ 
                (*feeEffectiveCurve)[t] = 1.0;
                if (computeCondCurve) {
                    (*feeEffectiveCurveCond)[t] = 1.0;
                }
            }
        }
    } 
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** Invoked by model to compute the price. Note not virtual (although
    could be) as don't expect products to override */
void ConvolutionProduct::price(const IConvolutionModel* convolutionModel,
                               Control*                 control, 
                               Results*                 results,
                               IForwardRatePricerSP     model) const
{
    static const string method("ConvolutionProduct::price");
    try {
        ICreditLossConfigConstSP lossConfig = getLossConfig();
        ICreditLossGenSP lossGen = 
            convolutionEngine->lossGenerator(lossConfig,
                                             IModelConfigMapperConstSP());
        IEffectiveCurveLossGenSP effCurveLossGen =
            DYNAMIC_POINTER_CAST<IEffectiveCurveLossGen>(lossGen);
        if (!effCurveLossGen) {
            throw ModelException(method,
                                 "Internal error: The loss generator is not of "
                                 "type IEffectiveCurveLossGen");
        }

        // Build timeline using the loss generator
        const DateTime& maturity = lastObservationDate();
        DateTimeArraySP timeline = 
            effCurveLossGen->generateTimeline(maturity);

        // are we doing main pricing and has CPTY_CREDIT_CHARGE been requested?
        bool isPricing = control->isPricing();
        OutputRequest* cccRequest = 
            control->requestsOutput(OutputRequest::CPTY_CREDIT_CHARGE);

        // CCC must be calculated if we have a counterparty
        bool doCCC = !!getCounterParty();

        //does the instrument dictate that we
        //should recover notional as well ?
        bool doNotionalRecovery = recoverNotional();
        //and if so, can the model assist in optimising out
        //cases where it really makes no difference
        if (doNotionalRecovery) {
            doNotionalRecovery = effCurveLossGen->modelRecoveredNotional();
        }

        // build effective curve
        // with recovered notional, which is only applied to the fee leg
        // we get a pair of effective curves - one for the contingent
        // leg and one for the fee leg
        DoubleArraySP ctgEffectiveCurve;
        DoubleArraySP ctgEffectiveCurveCond;
        DoubleArraySP feeEffectiveCurve;
        DoubleArraySP feeEffectiveCurveCond;
        IFixedTrancheLossCalculatorConstSP conditionalLossCalc;
        IFixedTrancheLossCalculatorConstSP conditionalRecNtnlCalc;
        EffectiveCurveSP ctgEffCurve;
        EffectiveCurveSP feeEffCurve;
        EffectiveCurveSP ctgCondEffCurve;
        EffectiveCurveSP feeCondEffCurve;

        if (numNames() > 0){        
            // get hold of the loss calculator(s)
            IFixedTrancheLossCalculatorConstSP lossCalculator;
            IFixedTrancheLossCalculatorConstSP recoveredNotionalCalculator;

            try {
                effCurveLossGen->createLossCalculators(
                    *timeline, getCounterParty(), maturity,
                    control, results, doNotionalRecovery,
                    getDiscount(),
                    lossCalculator, recoveredNotionalCalculator,
                    conditionalLossCalc, conditionalRecNtnlCalc);
            } 
            catch (exception& e){
                throw ModelException(e, method,
                                     "Failed to build loss calculators");
            }

            // to do: put in check for trancheOutstandingNotional != 0.0
            calculateEffectiveCurve(lossCalculator,
                                    recoveredNotionalCalculator,
                                    *timeline,
                                    ctgEffectiveCurve,
                                    feeEffectiveCurve,
                                    doCCC,
                                    ctgEffectiveCurveCond,
                                    feeEffectiveCurveCond);
            storeCalculatorResults(results, control, lossCalculator);
        } 
        else {
            //initialise to array of correct size
            //with with 0 rates
            int numTimes = timeline->size();
            ctgEffectiveCurve.reset(new DoubleArray(numTimes));
            feeEffectiveCurve.reset(new DoubleArray(numTimes));
            ctgEffectiveCurveCond.reset(new DoubleArray(numTimes));
            feeEffectiveCurveCond.reset(new DoubleArray(numTimes));
        }

        //for convenience
        string lossInterpolation = convolutionModel->getLossInterpolation();

        //package effective curves into a useful structure
        //...contingent leg
        ctgEffCurve = EffectiveCurveSP(
            new EffectiveCurve(getToday(), discount, *timeline,
                               *ctgEffectiveCurve, lossInterpolation));
        //...fee leg
        feeEffCurve = EffectiveCurveSP(
            new EffectiveCurve(getToday(), discount, *timeline,
                               *feeEffectiveCurve, lossInterpolation));

        //get scaling factor for Kapital backwards compatability
        double forwardFactor = calculateForwardFactor(discount);

        //now call the payoff (first without CCC)
        OutputSP output = payoff(ctgEffCurve, feeEffCurve, forwardFactor, model);

        OutputSP outputWithCCC;
        double   ccc = 0.0;
        if (doCCC){
            OutputSP outputWithoutCCC = output;
            SingleCreditAssetConstSP counterParty(getCounterParty()->getAsset());

            //package conditional effective curves
            //...contingent leg conditional on counterparty
            ctgCondEffCurve = EffectiveCurveSP(
                new EffectiveCurve(getToday(), discount, counterParty, *timeline,
                                *ctgEffectiveCurveCond, lossInterpolation));

            //...fee leg conditional on counterparty
            feeCondEffCurve = EffectiveCurveSP(
                new EffectiveCurve(getToday(), discount, counterParty, *timeline,
                                *feeEffectiveCurveCond, lossInterpolation));

            if (conditionalLossCalc.get()) {
                // we have to recalculate the fair value without CCC
                // to do: put in check for trancheOutstandingNotional != 0.0
                DoubleArraySP ctgEffectiveCurveWithoutCCC;
                DoubleArraySP feeEffectiveCurveWithoutCCC;
                ctgEffectiveCurveCond.reset();
                feeEffectiveCurveCond.reset();
                // rebuild curve: to do: switch off doCCC in call to 
                // calculateEffectiveCurve above in this case
                calculateEffectiveCurve(conditionalLossCalc,
                                        conditionalRecNtnlCalc,
                                        *timeline, 
                                        ctgEffectiveCurveWithoutCCC,
                                        feeEffectiveCurveWithoutCCC,
                                        doCCC,
                                        ctgEffectiveCurveCond,
                                        feeEffectiveCurveCond);

                //package effective curves
                //...contingent leg
                EffectiveCurveSP ctgEffCurveWithoutCCC = EffectiveCurveSP(
                    new EffectiveCurve(getToday(), discount, *timeline,
                                    *ctgEffectiveCurveWithoutCCC, lossInterpolation));
                //...fee leg
                EffectiveCurveSP feeEffCurveWithoutCCC = EffectiveCurveSP(
                    new EffectiveCurve(getToday(), discount, *timeline,
                                    *feeEffectiveCurveWithoutCCC, lossInterpolation));

                //...contingent leg conditional on counterparty
                ctgCondEffCurve.reset(
                    new EffectiveCurve(getToday(), discount, counterParty, *timeline,
                                    *ctgEffectiveCurveCond, lossInterpolation));
                //...fee leg conditional on counterparty
                feeCondEffCurve.reset(
                    new EffectiveCurve(getToday(), discount, counterParty, *timeline,
                                    *feeEffectiveCurveCond, lossInterpolation));

                // and then fair value
                outputWithoutCCC = payoff(ctgEffCurveWithoutCCC,
                                          feeEffCurveWithoutCCC,
                                          forwardFactor,
                                          model);
            }
            // and then handling counterparty (and with risky discount curve).
            outputWithCCC = payoff(ctgCondEffCurve,
                                   feeCondEffCurve,
                                   forwardFactor,
                                   model);

            //// calculate counterparty credit charge (CCC)
            EffectiveCurveSP cptyCurve(
                new EffectiveCurve(
                    getToday(),
                    discount,
                    counterParty,
                    EffectiveCurve::FLAT_FORWARD));

            ccc = computeCreditCharge(cptyCurve,
                                      convolutionEngine->getCreditChargeViewType(),
                                      outputWithoutCCC, outputWithCCC);
        }
        // store our [primary] results
        storePrimaryResults(ctgEffCurve, feeEffCurve,
                            results, control, output,
                            cccRequest, doCCC, ccc);
        if (isPricing){ // store any extra results
            storeExtraResults(results, 
                              control, 
                              output, 
                              outputWithCCC, 
                              model, 
                              convolutionModel);
        }
    } 
    catch (exception& e){
        const string& modelName = convolutionModel->getClass()->getName();
        throw ModelException(e, method, "Failed whilst using " + modelName);
    }
}

/** calculates forward factor to price things as of IR spot date */
double ConvolutionProduct::calculateForwardFactor(
    YieldCurveConstSP discountCurve) const{
    return (pvToSpot()? 1.0/discountCurve->pv(cfCutOffDate()): 1.0);
}

/** ask calculator to stores its results */
void ConvolutionProduct::storeCalculatorResults(
    Results* results, Control* control, 
    IFixedTrancheLossCalculatorConstSP calculator) const
{
    calculator->storeResults(results, control);
}

/** Store results that ConvolutionProduct is capable of storing */
void ConvolutionProduct::storePrimaryResults(
    const EffectiveCurveSP ctgEffectiveCurve,
    const EffectiveCurveSP feeEffectiveCurve,
    Results*               results,
    Control*               control,
    OutputSP               output,
    OutputRequest*         cccRequest,
    bool                   doCCC,      
    double                 ccc) const
{
    double price = output->price();
    //price now takes CCC into account
    results->storePrice(price-ccc, discount->getCcy());
    if (control->isPricing()){
        // CPTY_CREDIT_CHARGE
        if (cccRequest) {
            if (!doCCC) {
                results->storeNotApplicable(cccRequest);
            } else {
                results->storeRequestResult(cccRequest, ccc);
            }
        }
        // FV_MINUS_CCC : redundant with price taking CCC into account
        // but remains for backwards compatability
        OutputRequest* request = 
            control->requestsOutput(OutputRequest::FV_MINUS_CCC);
        if (request) {
            if (!doCCC) {
                results->storeNotApplicable(request);
            } else {
                results->storeRequestResult(request, price - ccc);
            }
        }
        // VALUE_WITHOUT_CCC : price disregarding CCC
        request = control->requestsOutput(OutputRequest::VALUE_WITHOUT_CCC);
        if (request) {
            results->storeRequestResult(request, price);
        }

        request = control->requestsOutput(OutputRequest::IND_CDS_PAR_SPREAD);
        OutputRequest* request2 =
            control->requestsOutput(OutputRequest::PAR_SPREAD_CURVE);
        if (request || request2){
            storeIndicativeSpreads(request, request2, results);
        }

        request = control->requestsOutput(
            OutputRequest::TRANCHE_EXPECTED_LOSS_CURVE);
        if (request){
            storeExpectedLossCurve(ctgEffectiveCurve, request, results);
        }

        request = control->requestsOutput(
            OutputRequest::TRANCHE_EXPECTED_LOSS_CURVE_WITH_RECOVERED_NOTIONAL);
        if (request){
            storeExpectedLossCurve(feeEffectiveCurve, request, results);
        }
    }
}

/*
 * Compute the 'Effective Curve'
 * i.e. 'zero curve' with the following default conventions :
 * - Day count : Act/365F
 * - Basis : Annual
 * */
void ConvolutionProduct::storeExpectedLossCurve(
    const EffectiveCurveSP effectiveCurve,
    OutputRequest*         request,
    Results*               results) const
{
    CashFlowArraySP cashFlows = effectiveCurve->asZeroRates();
    results->storeRequestResult(request, cashFlows);
}

//// Little helper method to implement IND_CDS_PAR_SPREAD and/or
//// PAR_SPREAD_CURVE request
void ConvolutionProduct::storeIndicativeSpreads(OutputRequest* indRequest,
                                                OutputRequest* parRequest,
                                                Results*       results) const
{
    const DateTime& today = getToday(); // for ease
    const DateTime& lastPayDate = this->lastPayDate();
    int numNames = this->numNames();
    for (int i = 0; i< numNames; i++) {
        SingleCreditAssetConstSP myAsset(nameAsset(i));
        ICDSParSpreadsConstSP parSpreadCurve(myAsset->getParSpreadCurve());
        // create name for result
        OutputNameConstSP nm(new OutputName(myAsset->getName()));
        if (indRequest){
            // call dubious ICreditCurve method
            IObjectSP indSpread(
                parSpreadCurve->getCurrentSpreadOrUntweakable(today, 
                                                              lastPayDate));
            // and store
            results->storeRequestResult(indRequest, indSpread, nm);
        }
        if (parRequest){
            ExpiryResultArraySP curve(
                ExpiryResult::createExpiryResultArray(
                    *parSpreadCurve->getParSpreadsExpiries(),
                    *parSpreadCurve->getParSpreads()));
            results->storeRequestResult(parRequest, curve, nm);
        }
    }
}

/** Create risky discount curve -- but in old style structure */
/*
FlatFwdZeroCurveSP ConvolutionProduct::createCptyCleanSpreadCurve() const{
    SingleCreditAssetConstSP counterParty(getCounterParty()->getAsset());
    DefaultRatesSP defaultRates(counterParty->getParSpreadCurve()->
                                defaultRates());
    // get the raw data
    CashFlowArraySP cleanFwdSpreadCurve(defaultRates->getCleanSpreadCurve());
    return FlatFwdZeroCurveSP(new FlatFwdZeroCurve(getToday(), 
                                                   *cleanFwdSpreadCurve));
}
*/
DRLIB_END_NAMESPACE
