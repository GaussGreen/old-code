#include "edginc/config.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/TwoSidedDeriv.hpp"
#include "edginc/AssetSpotGreek.hpp"
#include "edginc/Delta.hpp"
#include "edginc/CrossGamma.hpp"
#include "edginc/MCPricing.hpp"  // Could remain in .cpp directory?
#include "edginc/MCProduct.hpp"
#include "edginc/MCProductIQuicks.hpp"

DRLIB_BEGIN_NAMESPACE

//////////////////////////////////////////////////////////////////////
// MCPricing class
//////////////////////////////////////////////////////////////////////
MCPricing::~MCPricing(){}

/** Creates IMCPrices object for 'past' simulation run */
IMCPrices* MCPricing::createPastPrices(
    const MonteCarlo*          mcarlo,  // the model
    IMCProduct*                 product)
{
    return product->createOrigPrices(1, 1, 0 /* don't save path*/);
    // no need to reset() since this is done on construction anyway
}

/** whether we should ask the path generator to save data such
    that random access to paths can be provided */
bool MCPricing::needRandomAccessToPaths(
    MonteCarlo*                mcarlo,  // the model
    IMCProduct*                 product,
    Control*                   control)
{
    return false;
}

/** Creates IMCPrices object for future simulation run.
    May have request to cache product when will need to
    use cached values for greeks, or may be simple case. */
IMCPricesSP MCPricing::createFuturePrices(
    MonteCarlo*                mcarlo,  // the model
    IMCProduct*                 product,
    Control*                   control,
    const MCPathGeneratorSP&   pathGen) // path generator
{
    IMCPricesSP p;
    if (!control->isPricing() &&
        (mcarlo->cachedQckGreeksOfType( IMCProduct::CACHE_PRODUCT_BIT ))){
        // greek call with product caching requested so reuse IMCPrices class
        p = mcarlo->getOrigPrices();
        p->reset();
    } else {
        // otherwise a new one each time is fine
        p = IMCPricesSP(
            product->createOrigPrices(mcarlo->getNbIters(),
                                      mcarlo->getNbSubSamples(),
                                      mcarlo->getCachedQckGreeksModes()));
    }
    return p;
}

/** Run the simulation */
void MCPricing::runSim(
    MonteCarlo*                mcarlo,  // the model
    IMCProduct*                 product,
    Control*                   control,
    MCPathGenerator* pathGen, // path generator
    IMCPrices&         prices){
    int nbIter = mcarlo->getNbIters();

    // loop requested number of times
    for (int idx = 0; idx < nbIter; idx++){
        // inform path generator that we're on the next path
        pathGen->generatePath(idx);
        // value this path
        product->payoff(pathGen, prices);
    }
}

/** Store the results of the simulation in the results. This
    must include the price */
void MCPricing::postResults(
    MonteCarlo*                mcarlo,  // the model
    IMCProduct*                 product,
    Control*                   control,
    IMCPrices&                  _prices,
    const string&              discountISOCode,
    Results*                   results)
{
    IMCPricesSimple&     prices= dynamic_cast<IMCPricesSimple&>(_prices);
    double thePrice;
    double stdErr;
    // get price and standard error
    prices.getResult(&thePrice, &stdErr);
    // then pv
    double pv = product->pvFromPaymentDate();
    thePrice *= pv;
    stdErr *= pv;

    results->storePrice(thePrice, discountISOCode);
    // are we pricing and has std err been requested
    if (control->isPricing()) {
        OutputRequest* request =
            control->requestsOutput(OutputRequest::VALUE_STE);
        if (request){ // yes so store it under <INSTRUMENT_PACKET, VALUE_STE>
            results->storeRequestResult(request, stdErr);
        }
    }
}

//////////////////////////////////////////////////////////////////////
// end of MCPricing class
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Pricers based on stateless payoffs
//////////////////////////////////////////////////////////////////////
void MCSlicePricing::runSim(
    MonteCarlo* mcarlo,  // the model
    IMCProduct* product,
    Control* control,
    MCPathGenerator* pathGen, // path generator
    IMCPrices& prices)
{
    static const string method("MCSlicePricing");

    // Check that the product supports the stateless payoff interface
    IMCStatelessProductClient* instrument =
        dynamic_cast<IMCStatelessProductClient*>(product);
    if ( !instrument )
        throw ModelException(method, "Product does not support stateless payoff interface");

    // Check that the path generator supports the new API
    MCStatelessSliceGen* gen = dynamic_cast<MCStatelessSliceGen*>(pathGen);
    if (!gen)
        throw ModelException(method, "PathConfig does not support stateless slice based simulation");

    // Finalize the instrument and record dates on which to call instrument
    vector<int> keyDateIdx = instrument->finalize( simDates );

    int nbIter = mcarlo->getNbIters();

    IHistoricalContextSP context = instrument->getInitialHistoricalContext();
    vector<IHistoricalContextSP> hContext( nbIter );
    if ( context.get() ) {
        for ( size_t i = 0; i < hContext.size(); ++i ) {
            hContext[i] = IHistoricalContextSP( instrument->createHistoricalContext() );
            context->deepCopyTo( hContext[i].get() );
        }
    }

    int instIdx = 0;
    for (int t = 0; t < simDates.size(); t++){
        gen->generateSlice(t);
        if ( keyDateIdx[instIdx] == t ) {
            for (int p = 0; p < nbIter; ++p) {
                instrument->statelessPayOff( instIdx, hContext[p], prices );
                gen->advance();
                instIdx += 1;
            }
        }
    }
}

DateTimeArray& MCSlicePricing::getTimeLine()
{
    return simDates;
}

void MCPathPricing::runSim(
    MonteCarlo* mcarlo,  // the model
    IMCProduct* product,
    Control* control,
    MCPathGenerator* pathGen, // path generator
    IMCPrices& prices)
{
    static const string method("MCPathPricing");

    // Check that the product supports the stateless payoff interface
    IMCStatelessProductClient* instrument =
        dynamic_cast<IMCStatelessProductClient*>(product);
    if ( !instrument )
        throw ModelException(method, "Product does not support stateless payoff interface");

    // Check that the path generator supports the new API
    MCStatelessPathGen* gen = dynamic_cast<MCStatelessPathGen*>(pathGen);
    if (!gen)
        throw ModelException(method, "PathConfig does not support stateless path based simulation");

    // Finalize the instrument and record dates on which to call the payoff
    vector<int> keyDateIdx = instrument->finalize( simDates );

    IHistoricalContextSP snapshot = instrument->getInitialHistoricalContext();
    IHistoricalContextSP context = instrument->createHistoricalContext();

    bool pathDependent = snapshot.get() != 0;
    if (pathDependent) {
        ASSERT(context.get() != 0);
        snapshot->deepCopyTo(context.get());
    }

    int nbIter = mcarlo->getNbIters();
    for (int p = 0; p < nbIter; ++p ) {
        gen->generatePath(p);
        for (size_t t = 0; t < keyDateIdx.size(); ++t) {
            instrument->statelessPayOff(t, context, prices);
            gen->advance();
        }
        if (pathDependent)
            snapshot->deepCopyTo(context.get());
    }
}

DateTimeArray& MCPathPricing::getTimeLine()
{
    return simDates;
}

void MCPastPricing::runSim(
    MonteCarlo* mcarlo,  // the model
    IMCProduct* product,
    Control* control,
    MCPathGenerator* pathGen, // path generator
    IMCPrices& prices)
{
    static const string method("MCPastPricing");

    // Check that the product supports the stateless payoff interface
    IMCStatelessProductClient* instrument =
        dynamic_cast<IMCStatelessProductClient*>(product);
    if ( !instrument )
        throw ModelException(method, "Product does not support stateless payoff interface");

    // Check that the path generator supports the new API
    IMCStatelessGen* gen = dynamic_cast<IMCStatelessGen*>(pathGen);
    if (!gen)
        throw ModelException(method, "PathConfig does not support stateless simulation");

    DateTimeArray pastDates = instrument->getPastDates();
    IHistoricalContextSP context = instrument->getInitialHistoricalContext();
    // Context could be null if there is no historical context
    // ASSERT(context.get() != 0);

    for (int i = 0; i < pastDates.size(); ++i) {
        instrument->statelessPayOff(i, context, prices);
        gen->advance();
    }
}



//////////////////////////////////////////////////////////////////////
// Specialised QGPricing class for quick greeks
//////////////////////////////////////////////////////////////////////
/** whether we should ask the path generator to save data such
    that random access to paths can be provided */
bool QGPricing::needRandomAccessToPaths(
    MonteCarlo*                mcarlo,  // the model
    IMCProduct*                 product,
    Control*                   control)
{
    // true on first pricing run
    return (control->isPricing());
}

/** Creates IMCPrices object for future simulation run. */
IMCPricesSP QGPricing::createFuturePrices(
    MonteCarlo*                mcarlo,  // the model
    IMCProduct*                 product,
    Control*                   control,
    const MCPathGeneratorSP&   pathGen){ // path generator
    IMCPricesSP p;
    if (control->isPricing()){
        // first pricing call
        p = IMCPricesSP(
            product->createOrigPrices(mcarlo->getNbIters(),
                                      mcarlo->getNbSubSamples(),
                                      mcarlo->getCachedQckGreeksModes()));
    } else {
        IMCQuickGreeks& qckGreeks =
            dynamic_cast<IMCQuickGreeks&>(*product);
        // have to switch between normal quick greeks and quick
        // delta/gamma. Could have more specialised classes going
        // forward
        SensitivitySP sens(control->getCurrentSensitivity());
        p = mcarlo->getOrigPrices();
        const ITwoSidedDeriv* twoSidedDeriv =
            dynamic_cast<const ITwoSidedDeriv*>(sens.get());
        if (twoSidedDeriv != 0){
            // only supported for IAssetSpotGreek at the moment
            if (dynamic_cast<const IAssetSpotGreek*>(sens.get())){
                qckGreeks.setPricesForTwoSidedGreek(
                    p.get(), pathGen.get(),
                    twoSidedDeriv->getComponentShifts());
            } else {
                // must not use quick greeks here
                p = IMCPricesSP(
                    product->createOrigPrices(mcarlo->getNbIters(),
                                              mcarlo->getNbSubSamples(),
                                              0)); // don't save
            }
        } else {
            // one sided case
            qckGreeks.setPricesForGreek(p.get(), pathGen.get(), sens.get());
        }
    }
    p->reset();
    return p;
}

/** Run the simulation */
void QGPricing::runSim(
    MonteCarlo*                mcarlo,  // the model
    IMCProduct*                 product,
    Control*                   control,
    MCPathGenerator* pathGen, // path generator
    IMCPrices&         _prices)
{
    MCPricesGeneral & prices = dynamic_cast<MCPricesGeneral&> (_prices);
#ifdef _DEBUG
    int numCalced = 0;
#endif
    bool isPricing = control->isPricing();
    int nbIter = mcarlo->getNbIters();
    // loop requested number of times
    for (int idx = 0; idx < nbIter; idx++){
        bool doCalc = isPricing || prices.repriceForGreek(idx);
        if (doCalc){
            // inform path generator that we're on the next path
            pathGen->generatePath(idx);
            // value this path
            product->payoff(pathGen, prices);
#ifdef _DEBUG
            numCalced++;
#endif
        }
    }
}

//////////////////////////////////////////////////////////////////////
// Specialised QXGammaPricing class for quick X gamma
//////////////////////////////////////////////////////////////////////
QXGamma::QXGamma(const IMCPricesSP& qxgPrices): qxgPrices(qxgPrices){}

/** Creates IMCPrices object for future simulation run. */
IMCPricesSP QXGamma::createFuturePrices(
    MonteCarlo*                mcarlo,  // the model
    IMCProduct*                 product,
    Control*                   control,
    const MCPathGeneratorSP&   pathGen) // path generator
{
    // we reuse the same IMCPrices object each time (for each individual
    // Cross Gamma calculation)
    qxgPrices->reset();
    return qxgPrices;
}

//////////////////////////////////////////////////////////////////////
// Specialised SubGamma class - for SubGamma (see quick X Gamma)
//////////////////////////////////////////////////////////////////////
/** Class is for calculating what we call SubGamma. SubGamma for asset
    A is the same as Gamma for asset A but it is only done on a subset
    of paths for which the cross gamma wrt assets A and B is non
    zero. Hence for each B there is a different SubGamma for A */

SubGamma::SubGamma(
    const MCPricingSP&              origPricing,
    bool                            doRealDelta,
    const MCPricesGeneralSP&        origPrices,
    int                             numAssets):
        origPricing(origPricing), doRealDelta(doRealDelta),
        gammaPrices(numAssets),
        gammaVals(numAssets, numAssets), assetIdx(0),
        firstPricingCall(true), completed(numAssets)
{
    for (int i = 0; i < numAssets; i++){
        gammaPrices[i] = MCPricesGeneralSP(dynamic_cast<MCPricesGeneral *>(origPrices->clone()));
    }
}

/** Creates IMCPrices object for future simulation run. */
IMCPricesSP SubGamma::createFuturePrices(
    MonteCarlo*                 mcarlo,  // the model
    IMCProduct*                 product,
    Control*                    control,
    const MCPathGeneratorSP&    pathGen){// path generator
    // we must use whatever would be used for Delta
    return origPricing->createFuturePrices(mcarlo, product,
                                           control, pathGen);
}

/** Run the simulation */
void SubGamma::runSim(
    MonteCarlo*                 mcarlo,  // the model
    IMCProduct*                 product,
    Control*                    control,
    MCPathGenerator*            pathGen, // path generator
    IMCPrices&                  _pricesRef)
{
    MCPricesGeneral& pricesRef = dynamic_cast<MCPricesGeneral&>(_pricesRef);
    int nbIter = mcarlo->getNbIters();
    // loop requested number of times
    for (int pathIdx = 0; pathIdx < nbIter; pathIdx++){
        double lastPrice;
        bool doneCalc = false;
        if (doRealDelta && pricesRef.repriceForGreek(pathIdx)){
            pathGen->generatePath(pathIdx);
            product->payoff(pathGen, pricesRef);
            doneCalc = true;
            lastPrice = pricesRef.lastPrice();
        }
        /* This block of code is in a nutshell what the sub gamma
           code is about. It allows us to calculate n-1 sub gammas
           'at the same time' */
        for (unsigned int i = 0; i < gammaPrices.size(); i++){
            if ((int)i != assetIdx){
                MCPricesGeneral& prices = *gammaPrices[i];
                if (prices.repriceForGreek(pathIdx)){
                    if (!doneCalc){
                        // need to price this path
                        pathGen->generatePath(pathIdx);
                        product->payoff(pathGen, prices);
#ifdef _DEBUG
                        debugCount++;
#endif
                        lastPrice = prices.lastPrice();
                        doneCalc = true;
                    } else {
                        prices.add(lastPrice);
                    }
                }
            }
        }
    }
}

/** Store the results of the simulation in the results. This
    must include the price */
void SubGamma::postResults(
    MonteCarlo*                 mcarlo,  // the model
    IMCProduct*                 product,
    Control*                    control,
    IMCPrices&                  prices,
    const string&               discountISOCode,
    Results*                    results)
{
    static const string method("SubGamma::store");
    if (completed[assetIdx]){
        throw ModelException(method, "Internal error 1");
    }
    SensitivitySP currentSens(control->getCurrentSensitivity());
    double divisor;
    if (!firstPricingCall){
        const Delta* delta =
            dynamic_cast<const Delta*>(currentSens.get());
        if (!delta){
            throw ModelException(method, "Internal error 2");
        }
        divisor = delta->divisor();
    }
    double pv = product->pvFromPaymentDate();
    for (unsigned int i = 0; i < gammaPrices.size(); i++){
        if ((int)i != assetIdx){
            double result, resultStdErrDoNotUse; // no pv done on std err
            gammaPrices[i]->getResult(&result, &resultStdErrDoNotUse);
            result *= pv;
            if (firstPricingCall){
                // just store
                gammaVals[assetIdx][i] = result;
                // and reset our IMCPrices object for next pricing call
                gammaPrices[i]->reset();
            } else {
                // calculate 'gamma'
                double shiftUpPrice = gammaVals[assetIdx][i];
                // note: 'origPrice' = 0.0 as we set IMCPrices object to
                // report change in price
                gammaVals[assetIdx][i] = (shiftUpPrice + result)/
                    (divisor * divisor);
            }
            // you could argue that we should be calling
            // recordExtraOutput for each IMCPrices but what does this
            // mean? We only have one Results object
        }
    }
    if (firstPricingCall){
        firstPricingCall = false;
    } else {
        completed[assetIdx] = true;
    }
    // handle main prices object
    if (!doRealDelta){
        // need to store a number of some sort
        results->storePrice(0.0, discountISOCode);
    } else {
        // invoke parent's method to do default for delta prices object
        MCPricing::postResults(mcarlo, product, control, prices,
                             discountISOCode, results);
    }
}

/** Prepare for calculating another set of subGammas wrt specified
    asset */
void SubGamma::resetForAsset(
    int                         assetIdx,
    IMCQuickXGamma*             qckXGamma,
    const MCPathGeneratorSP&    futurePathGen,
    const ScalarShiftArray&     theShifts)
{
    this->assetIdx = assetIdx;
#ifdef _DEBUG
    debugCount = 0;
#endif
    OutputNameArrayConstSP names1 = theShifts[0]->overrideNames();
    OutputNameArrayConstSP names2 = theShifts[1]->overrideNames();
    // theShifts passed to Quick X Gamma-says what is being shifted
    theShifts[0]->setMarketDataName((*names1)[assetIdx]);
    for (unsigned int i = 0; i < gammaPrices.size(); i++){
        if ((int)i != assetIdx){
            gammaPrices[i]->reset();
            theShifts[1]->setMarketDataName((*names2)[i]);
            qckXGamma->setPricesForXGamma(gammaPrices[i].get(),
                                          futurePathGen.get(),
                                          theShifts);
        }
    }
    firstPricingCall = true;
    completed[assetIdx] = false;
}
/** Calcuates all the sub gammas that we need */
void SubGamma::calculateSubGamma(
    MonteCarlo*                         mcarlo,
    IMCQuickXGamma*            qckXGamma,
    const MCPathGeneratorSP&  futurePathGen,
    Control*                            control,
    CInstrument*                        instrument,
    CResults*                           results,
    const ScalarShiftArray&             theShifts){
    // have to calculate the subGammas one set at a time
    // namesForTweak will hold what we want to be calculated
    OutputNameArraySP namesForTweak(new OutputNameArray(1));
    // clone the shift - as we want to alter the names inside it
    ScalarShiftSP theShift(theShifts[0].clone());
    OutputNameArrayConstSP names1 = theShifts[0]->overrideNames();
    for (int i = 0; i < names1->size(); i++){
        // this part needs fixing if not doing regular fx cross gamma
        // switch to next asset
        resetForAsset(i, qckXGamma, futurePathGen, theShifts);
        // then calculate delta/gamma + n-1 SubGammas
        (*namesForTweak)[0] = (*names1)[i];
        theShift->storeOverrideNames(namesForTweak);
        control->calculateSens(theShift, mcarlo, instrument, results);
    }
}

//////////////////////////////////////////////////////////////////////
// LRPricing class (do price + greeks using LR method)
//////////////////////////////////////////////////////////////////////
LRPricing::LRPricing(const MCPricingSP&         mainPricing):
    mainPricing(mainPricing) {}

/** whether we should ask the path generator to save data such
    that random access to paths can be provided */
bool LRPricing::needRandomAccessToPaths(
    MonteCarlo*                mcarlo,  // the model
    IMCProduct*                 product,
    Control*                   control){
    return mainPricing->needRandomAccessToPaths(mcarlo, product, control);
}

/** Creates IMCPrices object for future simulation run. */
IMCPricesSP LRPricing::createFuturePrices(
    MonteCarlo*                mcarlo,  // the model
    IMCProduct*                 product,
    Control*                   control,
    const MCPathGeneratorSP&   pathGen){ // path generator
    static const string method("LRPricing::createFuturePrices");
    // only do LR greeks on main pricing run
    if (control->isPricing()){
        // get hold of our LRGenerator
        auto_ptr<LRGenerator> lrGenerator(
            mcarlo->getPathConfig()->createLRGenerator(
                pathGen,
                mcarlo->getNbIters(),
                mcarlo->getNbSubSamples()));
        if (!lrGenerator.get()){
            throw ModelException(method,
                                 "Path Generator of type '"+
                                 mcarlo->getPathConfig()->getClass()->getName()+
                                 "' does not support calculating greeks "
                                 "using Likelihood Ratio methodology");
        }
        // then scan through greeks getting GreekCalculator for each one
        SensitivityArrayConstSP sens(control->getSens());
        for (int i = 0; i < sens->size(); i++){
            const Sensitivity* theSens = (*sens)[i].get();
            LRGenerator::GreekCalculator* calc =
                lrGenerator->greekCalculator(theSens,
                                             false/*no stderr for greeks*/);
            if (calc){
                greekCalcs.push_back(LRGenerator::GreekCalculatorSP(calc));
                if (theSens->getClass() == CrossGamma::TYPE){
                    // switch off quick x gamma
                    mcarlo->disableCachedQckGreeksMode(
                        IMCProduct::QUICK_X_GAMMA_BIT );
                }
            }
        }
    }
    return mainPricing->createFuturePrices(mcarlo, product,
                                           control, pathGen);
}

/** Run the simulation */
void LRPricing::runSim(
    MonteCarlo*                mcarlo,  // the model
    IMCProduct*                 product,
    Control*                   control,
    MCPathGenerator*           pathGen, // path generator
    IMCPrices&                  _prices)
{
    IMCPricesSimple&         prices = dynamic_cast<IMCPricesSimple&>(_prices);
    if (!control->isPricing()){
        mainPricing->runSim(mcarlo, product, control, pathGen, prices);
        return;
    }
    int nbIter = mcarlo->getNbIters();
    unsigned int numGreeks = greekCalcs.size();
    // loop requested number of times
    for (int idx = 0; idx < nbIter; idx++){
        // inform path generator that we're on the next path
        pathGen->generatePath(idx);
        // value this path
        product->payoff(pathGen, prices);
        // retrieve the price for this path
        double thePrice = prices.lastPrice();
        // now do greeks
        for (unsigned int i = 0; i < numGreeks; i++){
            greekCalcs[i]->processPath(idx, thePrice);
        }
    }
}
/** Store the results of the simulation in the results. This
    must include the price */
    // Simply calls mainPricing and then posts results for all greekCalcs
void LRPricing::postResults(
    MonteCarlo*                mcarlo,  // the model
    IMCProduct*                 product,
    Control*                   control,
    IMCPrices&         prices,
    const string&              discountISOCode,
    Results*                   results)
{
    // handle price and std err
    mainPricing->postResults(mcarlo, product, control, prices,
                             discountISOCode, results);
    // then store LR greeks
    if (control->isPricing()){
        double pv = product->pvFromPaymentDate();
        for (unsigned int i = 0; i < greekCalcs.size(); i++){
            greekCalcs[i]->storeResult(results, pv);
        }
    }
    // tidy up
    greekCalcs.resize(0);
}

//////////////////////////////////////////////////////////////////////
// end of Specialised LRPricing class
//////////////////////////////////////////////////////////////////////

DRLIB_END_NAMESPACE
