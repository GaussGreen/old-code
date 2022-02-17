//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MonteCarlo.cpp
//
//   Description : Monte Carlo Algorithm
//
//   Date        : May 2001
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/MonteCarlo.hpp"
//#include "edginc/MCPathConfig.hpp"
//#include "edginc/Format.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/IInstrumentCollection.hpp"
//#include "edginc/Results.hpp"
//#include "edginc/SimSeries.hpp"
//#include "edginc/InstrumentSettlement.hpp"
//#include "edginc/PastPathGenerator.hpp"
#include "edginc/Addin.hpp"
//#include "edginc/TwoSidedDeriv.hpp"
//#include "edginc/RhoParallel.hpp"
//#include "edginc/RhoPointwise.hpp"
#include "edginc/CommandLineParams.hpp"
//#include "edginc/OutputRequestUtil.hpp"
//#include "edginc/HandlePaymentEvents.hpp"
#include "edginc/SpotCrossDerivative.hpp"
//#include "edginc/Reprice.hpp"
//#include "edginc/Maths.hpp"
//#include "edginc/AssetSpotGreek.hpp"
#include "edginc/Delta.hpp"
//#include "edginc/Untweakable.hpp"
#include "edginc/FXCrossGamma.hpp"
//#include "edginc/CrossGamma.hpp"
//#include "edginc/CliquetVolRequest.hpp"
//#include "edginc/Timer.hpp"
//#include "edginc/SampleList.hpp"
#include "edginc/MCProductEngineClient.hpp"
#include "edginc/MCProductClient.hpp"
#include "edginc/MCPricing.hpp"
#include "edginc/MCProductIQuicks.hpp"
#include "edginc/DataDictionary.hpp"
#include "edginc/ToolkitDebug.hpp"
#include "edginc/MDFUtil.hpp"
#include "edginc/MCPathConfigCCM.hpp"

#include "edginc/VolProcessedBS.hpp"
#include "edginc/fieldRiskProperty.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/TweakFunction.hpp"
#include "edginc/FieldPath.hpp"
#include "edginc/IPDFBoundaryProb.hpp"
DRLIB_BEGIN_NAMESPACE


const int MonteCarlo::MAX_STORAGE = (50 * 1024 * 1024); // 50 Mb
const string MonteCarlo::CACHE_PATH_DEFAULT = "DEFAULT";
const string MonteCarlo::CACHE_PROD_EXTENSION = "_WITH_PROD";
const string MonteCarlo::SIM_MODE_DEFAULT = "PATHOPT";
const string MonteCarlo::QUICK_GREEK_DEFAULT = "DEFAULT";
const string MonteCarlo::LR_GREEK_DEFAULT = "DEFAULT";
const string MonteCarlo::DEBUG_GREEK_STD_ERR; // defined in 'load'
const string MonteCarlo::QUICK_GREEK_AND_X_GAMMA = "GREEKS_AND_X_GAMMA";
const string MonteCarlo::QUICK_GREEK_NO_X_GAMMA = "GREEKS_NO_X_GAMMA";
const string MonteCarlo::QUICK_X_GAMMA_ONLY = "X_GAMMA_ONLY";

MonteCarlo::~MonteCarlo(){} 

#define xDEBUG_INST_XGAMMA
/* Useful in conjuction with ScalarShift.cpp to spot issues with quick
   X Gamma. Use -loop command line option to find problem paths */
#ifdef DEBUG_INST_XGAMMA
static int debugNumIter = 0;
#endif

static bool doLRGreeks(const string& lrGreeks){
    bool useLR = false;
    if (lrGreeks == MonteCarlo::LR_GREEK_DEFAULT){
        // default is off
    } else {
        try{
            useLR = CBool::fromString(lrGreeks);
        } catch (exception& e){
            throw ModelException(e, "MonteCarlo::doLRGreeks",
                                 "Value '"+lrGreeks+
                                 "' for 'lrGreeks' not recognised");
        }
    }
    return useLR;
}


bool MonteCarlo::stateVarUsed() const{
    bool overridingUseStateVars = 
        CommandLineParams::hasParameter(CommandLineParams::UseStateVars);

    return overridingUseStateVars || useStateVars;
}

bool MonteCarlo::inStatelessMode() const{
    return CString::equalsIgnoreCase(simMode, "path")
            || CString::equalsIgnoreCase(simMode, "slice");
}

const string* MonteCarlo::getFamilyName() const{
    // delegate via pathConfig
    IModelFamily* modelFamily = dynamic_cast<IModelFamily*>(pathConfig.get());
    if (modelFamily) {
        return modelFamily->getFamilyName();
    }
    return 0;
}

void MonteCarlo::initializePathConfigPostPricing()
{
    pathConfigPostPricing = IMCPathConfigSP(pathConfig.clone());
}

bool MonteCarlo::cachedQckGreeksOfType( int type ) const
{
	// !! avoids compiler warning converting into to bool
    return !!(cachedQckGreeks & type); 
}

IMCPricesSP MonteCarlo::getOrigPrices() const
{
    return origPrices;
}

// bypass is false by default: 
// if true, or if there is no sampling, we ignore the number of paths in 
//     the sampling
int MonteCarlo::getNbIters(bool bypass) const
{
	MCPathConfigCCM *pc; 
    if (!(pc = dynamic_cast<MCPathConfigCCM *>(pathConfig.get())))
		return nbIter;

	if (bypass)
		return nbIter;

	if (!(pc->getSampling().get()))
		return nbIter;

	return pc->getSampling()->getNumSamples();
}

int MonteCarlo::getNbSubSamples() const
{
    return nbSubSamples;
}

int MonteCarlo::getCachedQckGreeksModes() const
{
    return cachedQckGreeks;
}

MCPricingSP MonteCarlo::getPricing() const
{
    return pricing;
}

bool MonteCarlo::getCachedCachePaths() const
{
    return cachedCachePaths;
}

IMCPathConfigSP MonteCarlo::getPathConfig() const
{
    return pathConfig;
}

void MonteCarlo::setOrigPrices( IMCPricesSP prices )
{
    origPrices = prices;
}

void MonteCarlo::enableCachedQckGreeksMode(int mode)
{
    cachedQckGreeks |= mode;
}

void MonteCarlo::disableCachedQckGreeksMode(int mode)
{
    cachedQckGreeks &= ~mode;
}

IMCProduct* MonteCarlo::createProduct(const Instrument*  instrument) const
{
    static const string method("MonteCarlo::createProduct");
    if (!IMCIntoProduct::TYPE->isInstance(instrument)){
        throw ModelException(method, "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support MonteCarlo");
    }
    try{
        const IMCIntoProduct& intoProd = 
            dynamic_cast<const IMCIntoProduct&>(*instrument);
        
        // Put product in auto_ptr for temporary memory management:
        // see exceptions below
        auto_ptr<IMCProduct> prod(intoProd.createProduct(this));
        
        // Check consistency with state variables flag
        MCProductClient* prodClient = dynamic_cast<MCProductClient*>(prod.get());
        if( (prodClient  && stateVarUsed()) ||
            (!prodClient && !stateVarUsed()) ) {
            // Flags are consistent
        } else {
            // Flags are inconsitent
            string msg(stateVarUsed() ? " not " : " ");
            throw ModelException(
                "Inconsistent state variables information: "
                "useStateVars flag is " + 
                Format::toString(stateVarUsed()) + " but product is" + 
                msg + "using state variables");
        }

        // Release auto_ptr otherwise product will be deleted
        return prod.release();
    } catch (exception& e){
        throw ModelException(e, method,  
                             "Failed to createProduct for instrument of type "+
                             instrument->getClass()->getName());
    }
}


/** Called on initial pricing run */
void MonteCarlo::createPricingObject(Control*   control,
                                     IMCProduct* product)
{
    if (inStatelessMode()) {
        if (cachedQckGreeks & IMCProduct::QUICK_GREEKS_BIT)
            throw ModelException("quick greeks not supported in stateless mode!");
        if (doLRGreeks(lrGreeks))
            throw ModelException("LR greeks not supported in stateless mode!");
    }
    // detect sim modes
    if (CString::equalsIgnoreCase(simMode, "slice")) {
        pricing = MCPricingSP(new MCSlicePricing());
        return;
    } else if (CString::equalsIgnoreCase(simMode, "path")) {
        pricing = MCPricingSP(new MCPathPricing());
        return;
    } else if (CString::equalsIgnoreCase(simMode, "pathOpt")) {
        /* currently either standard or quick greeks (quick x gamma uses others
           but at a later stage) 
           Product-level cache is independent of this, and should work ok
           whichever MCPricing class is built.
        */
        if (cachedQckGreeks & IMCProduct::QUICK_GREEKS_BIT){
            pricing = MCPricingSP(new QGPricing());
        } else {
            pricing = MCPricingSP(new MCPricing());
        }
        // then wrap this pricing object in LR pricing object
        if (doLRGreeks(lrGreeks)){
            pricing = MCPricingSP(new LRPricing(pricing));
        }
    } else throw ModelException("unknown simMode " + simMode);
}

static double getSumSqr(double value, double stdErr, int nbSubSamples){
    double sumSqr = nbSubSamples * 
        (stdErr * stdErr * (nbSubSamples - 1) + value * value);
    return sumSqr;
}
enum AggregateStdErrMode{
    FIRST_TIME = 0, // first time in
    SUB_BLOCK,      // doing sub block (but not first time)
    CALCULATE       // calculate final value
};

/** Tots up sumSqrSoFar for each names in packetName, sumSqrSoFar will
    be initialised if firstTime is true else it must be of the same
    length as names */
static void aggregateStdErrForGreek(
    const string&             packetName,
    AggregateStdErrMode       mode,
    int                       numSubSamples,
    Results*                  results,
    DoubleArray&              sumSqrSoFar){
    static const string method("aggregateStdErrForGreek");
    static const string errorMsg(
        "Internal error. Mismatch between number of [LR] greeks"
        " calculated in sub-blocks");
    // look to see if corresponding packet exists
    if (!results->packetExists(MonteCarlo::DEBUG_GREEK_STD_ERR+packetName)){
        if (!sumSqrSoFar.empty()){
            throw ModelException(method, errorMsg);
        }
        return;
    }
    OutputNameArraySP namesSP(results->packetContents(packetName));
    OutputNameArray& names = *namesSP; // for ease
    if (mode == FIRST_TIME){
        sumSqrSoFar = DoubleArray(names.size());
    } else if (sumSqrSoFar.size() != names.size()){
        throw ModelException(method, errorMsg);
    }
    for (int i = 0; i < names.size(); i++){
        // get the value for this greek
        IObjectConstSP value(results->retrieveGreek(packetName, names[i]));
        // get the std err for this greek
        IObjectConstSP stdErr(results->retrieveGreek(
                                  MonteCarlo::DEBUG_GREEK_STD_ERR+packetName,
                                  names[i]));
        if (CDouble::TYPE->isInstance(value) && 
            CDouble::TYPE->isInstance(stdErr)){
            double val = dynamic_cast<const CDouble&>(*value).doubleValue();
            if (mode == CALCULATE){
                double stdErr = numSubSamples > 1?
                    (sumSqrSoFar[i]/numSubSamples - val*val)/(numSubSamples-1): 
                    0.0;
                stdErr = stdErr <= 0.0? 0.0: sqrt(stdErr);
                results->storeScalarGreek(
                    stdErr, MonteCarlo::DEBUG_GREEK_STD_ERR+packetName,
                    names[i]);
            } else {
                const CDouble& stdE = dynamic_cast<const CDouble&>(*stdErr);
                sumSqrSoFar[i] += getSumSqr(val, stdE.doubleValue(),
                                            numSubSamples);
            }
        }
    }
}

/** Simple loop around aggregateStdErrForGreek */
static void aggregateStdErrForGreek(
    const SensitivityArrayConstSP& sens,
    AggregateStdErrMode            mode,
    int                            numSubSamples,
    Results*                       results,
    DoubleArrayArray&              sumSqrSoFar){ // of size sens->size()+1
    bool doGamma = false;
    for (int iSens = 0; iSens < sens->size(); iSens++){
        const string& packetName = (*sens)[iSens]->getPacketName();
        aggregateStdErrForGreek(packetName, mode, numSubSamples,
                                results, sumSqrSoFar[iSens]);
        if (packetName == Delta::NAME){
            doGamma = true; 
        }
    }
    if (doGamma){
        // we should have made gamma an explicit request ...
        aggregateStdErrForGreek(Delta::SECOND_ORDER_NAME, mode, 
                                numSubSamples, results, sumSqrSoFar.back());
    }
}

static void aggregateStdErrForGreek_multi(
    const SensitivityArrayConstSP& sens,
    AggregateStdErrMode            mode,
    int                            numSubSamples,
    CResultsArray*                  resultss,
    vector<DoubleArrayArray>&      sumSqrSoFar){

    for (int i = 0; i < resultss->size(); ++i)
        aggregateStdErrForGreek(sens, mode, numSubSamples,
                                (*resultss)[i].get(), sumSqrSoFar[i]);
}

/** Returns the maximum storage space to use for paths given the users
    input. A return value of zero means do not cache anything */
static int getMaxStorage(const string& cachePaths,
                         bool&         cacheProd){
    int maxStorage;
    try{
        cacheProd = false; // unless ...
        if (CString::equalsIgnoreCase(MonteCarlo::CACHE_PATH_DEFAULT, 
                                      cachePaths)){
            maxStorage = MonteCarlo::MAX_STORAGE;
        } else {
            // Recognise a trailing "_WITH_PROD"
            string::size_type iWithProd = cachePaths.find(MonteCarlo::CACHE_PROD_EXTENSION);
            string cp = cachePaths;
            if (iWithProd != string::npos) {
                cacheProd = true;
                cp.erase(iWithProd); // remove the trailing "_WITH_PROD"
            }
            // see if it's a number
            char* endPos;
            double dVal = strtod(cp.c_str(), &endPos);
            if (endPos == cp.c_str() || *endPos != '\0'){
                // not recognised
                maxStorage = CBool::fromString(cp)? 
                    MonteCarlo::MAX_STORAGE: 0;
            } else {
                maxStorage = int(Maths::max(dVal, 0.0) * 1024 * 1024);
            }
        }
    } catch (exception& e){
        throw ModelException(e, "getMaxStorage", "Failed to interpret '"+
                             cachePaths+"' for cachePaths field in "
                             "MonteCarlo");
    }
    return maxStorage;
}

static int productStoragePerPath(IMCProduct* prod, IMCPathConfig* conf,
                                 bool cacheProd) {
    // FIXME or is this once-only?
    int storage = conf->storagePerPath(prod);

    // KKK - something here for product caching
    // Need prices class just to know memory demands when caching
    // Better would be a static method on the specific IMCPrices class, 
    // but this should do.
    IMCPricesSP pricesAP(
        prod->createOrigPrices(1, 1, 
                               cacheProd?IMCProduct::CACHE_PRODUCT_BIT:0));
    storage += pricesAP->storagePerPath(prod);

    return storage;
}

/** override of main control - calculates price and sensitivities in
    blocks to limit memory usage. */
CResultsArraySP MonteCarlo::RunMulti(IInstrumentCollectionSP instruments,
                                     CControl* control){
    static const string method("MonteCarlo::RunMulti");
    try {
        control->startTiming(); // start clock
        control->reset();
        // do we cache paths?
        bool cacheProd;
        int maxStorage = getMaxStorage(cachePaths, cacheProd); // validate setting first
        // and then override if it makes no sense
        SensitivityArrayConstSP sens(control->getSens());
        if (sens->empty()){
            maxStorage = 0; // do calc in one loop
        }    
#ifdef DEBUG_INST_XGAMMA
        maxStorage = 0; // do calc in one loop
#endif
        cachedCachePaths = maxStorage != 0;
        // validate first
        instruments->Validate();
        if (inStatelessMode())
            cachedCachePaths = false; // turn off path caching for statelessPayoff

        // For overnight intraweek testing allow fewer iterations via
        // command line flag.
        bool doingFewIter = 
            CommandLineParams::hasParameter(CommandLineParams::FewIter);
        int nbIter = doingFewIter? 20: this->nbIter; // hide field
        int nbSubSamples = doingFewIter? 10: this->nbSubSamples; // hide field
        // then find out our memory limit - need product though
        int storagePerPath = 0; // if not caching paths
        if (cachedCachePaths){
            // FIXME ugly hack, works because instruments isa VanillaGridInstrumentCollection
            if (IMCIntoProduct::TYPE->isInstance(instruments)) {
                IMCProductSP prod(dynamic_cast<IMCIntoProduct *>(
                    instruments.get())->createProduct(this));
                storagePerPath += productStoragePerPath(
                    prod.get(), pathConfig.get(), cacheProd);
            }
            else {
                for (int i = 0; i < instruments->size(); ++i) {
                    IMCProductSP prod(createProduct((*instruments)[i].get()));
                    storagePerPath += productStoragePerPath(
                        prod.get(), pathConfig.get(), cacheProd);
                }
            }
        }
        int numPerSubSample = nbIter/nbSubSamples; // must keep this the same
        int storagePerSubSample = numPerSubSample * storagePerPath;
        if (doingFewIter) {
            maxStorage = storagePerPath * 9; /* doing 20 iter - should split
                                                into 4, 8, and 8 */
        }
        int maxNumSubSamples = storagePerSubSample == 0?
            nbSubSamples: maxStorage/storagePerSubSample;
        if (maxNumSubSamples == 0){
            cachedCachePaths = false; // switch off
            maxNumSubSamples = nbSubSamples; // do in one block
        }
        // clone model
        smartPtr<MonteCarlo> modelCopy(copy(this));// then must set nbSubSamples.
        // nbIter set from nbSubSamples & numPerSubSample
        int numBlocks = nbSubSamples/maxNumSubSamples;
        if ((nbSubSamples % maxNumSubSamples) != 0){
            // have stub - put at front
            modelCopy->nbSubSamples = nbSubSamples - numBlocks*maxNumSubSamples;
            numBlocks++;
        } else {
            modelCopy->nbSubSamples = maxNumSubSamples;
        }
#ifdef DEBUG_INST_XGAMMA
        // force through single pass with debugNumIter of paths
        debugNumIter++;
        numPerSubSample = debugNumIter;
        modelCopy->nbSubSamples = 1;
#endif

        // then get ready to price the blocks
        CResultsArraySP resultss = instruments->emptyResults();

        CControlSP controlToUse(copy(control));
        // for performance, switch off COMPUTE_INDEX and COMPUTE_ESTIMATE
        // in subblocks
        OutputRequest* request1 = 
            controlToUse->requestsOutput(OutputRequest::COMPUTE_INDEX);
        if (request1){
            controlToUse->removeRequest(OutputRequest::COMPUTE_INDEX);
        }
        OutputRequest* request2 = 
            controlToUse->requestsOutput(OutputRequest::COMPUTE_ESTIMATE);
        if (request2){
            controlToUse->removeRequest(OutputRequest::COMPUTE_ESTIMATE);
        }
        bool removePriceTime = false;
        OutputRequest* requestPriceTime = 
            controlToUse->requestsOutput(OutputRequest::PRICE_TIME);
        if (!requestPriceTime){
            // need to ensure this gets calculated (and aggregated) so that
            // COMPUTE_ESTIMATE etc can be done
            controlToUse->addRequest(
                OutputRequestSP(new OutputRequest(OutputRequest::PRICE_TIME)));
            removePriceTime = true;
        }

        OutputRequest* stdErrRequest = 
            controlToUse->requestsOutput(OutputRequest::VALUE_STE);
        vector<bool> doStdErr(instruments->size(), numBlocks > 1 && stdErrRequest);
        double sumSqrSoFar = 0.0; // for working out true std err
        // same but for greeks
        bool useLR = doLRGreeks(lrGreeks);
        vector<DoubleArrayArray> greeksSumSqrSoFar(
            instruments->size(), DoubleArrayArray(useLR? sens->size()+1: 0));

        // then loop over our blocks
        for (int b = 0; b < numBlocks; b++){
            try{
                CResultsArraySP tmpResults = b == 0 ?
                    resultss : instruments->emptyResults();
                int ourNbSubSamples = modelCopy->nbSubSamples;
                // set nbIter for this block
                modelCopy->nbIter = ourNbSubSamples * numPerSubSample;
                // do the actual calculation
                controlToUse->calculateMulti(modelCopy.get(), instruments, tmpResults);
                // then do work needed to allow us to aggregate std err correctly
                // note: doing this before we scale resultss

                for (int i = 0; i < resultss->size(); ++i)
                    if (doStdErr[i]){ // if std err is still ok (and we're doing it)
                        IObjectConstSP stdErrObj = 
                            (*tmpResults)[i]->retrieveRequestResult(
                                OutputRequest::VALUE_STE);
                        if (CDouble::TYPE->isInstance(stdErrObj)){
                            const CDouble& dbl = 
                                dynamic_cast<const CDouble&>(*stdErrObj);
                            sumSqrSoFar += getSumSqr(
                                (*tmpResults)[i]->retrievePrice(), 
                                dbl.doubleValue(), ourNbSubSamples);
                        } else {
                            doStdErr[i] = false;
                        } 
                    }

                // sort out standard error for [LR] greeks
                if (numBlocks > 1 && useLR){
                    aggregateStdErrForGreek_multi(sens,
                                                  b == 0? FIRST_TIME: SUB_BLOCK,
                                                  ourNbSubSamples,
                                                  tmpResults.get(),
                                                  greeksSumSqrSoFar);
                }
                // must calculate weighted mean of resultss
                double scaling = (double)ourNbSubSamples/(double)nbSubSamples;
                if (b != 0){
                    // merge resultss in
                    for (int i = 0; i < resultss->size(); ++i)
                        (*resultss)[i]->add((*tmpResults)[i].get(), controlToUse,
                                            scaling, true);
                } else if (numBlocks > 1){ // performance improvement
                    for (int i = 0; i < resultss->size(); ++i)
                        (*resultss)[i]->scale(controlToUse, scaling, true);
                    // set nbIter and nbSubSamples correctly for future runs
                    modelCopy->nbSubSamples = maxNumSubSamples;
                }
                if (b < numBlocks - 1){
                    // switch path config to the state as if just done pricing
                    modelCopy->pathConfig = modelCopy->pathConfigPostPricing;
                    // tell path config to create new generators using current
                    // state of random number generator
                    modelCopy->pathConfig->saveRandomGeneratorState();
                }
                modelCopy->pathConfigPostPricing.reset(); // tidy up
            } catch (exception& e){
                throw ModelException(e, method);
            }
        }
        // then work out our std err
        for (int i = 0; i < resultss->size(); ++i) {
            if (doStdErr[i]) {
                // standard error is aggregated so any error will be propagated
                double price = (*resultss)[i]->retrievePrice();
                double stdErr = nbSubSamples > 1?
                    (sumSqrSoFar/nbSubSamples - price * price)/(nbSubSamples - 1) :
                    0.0;
                stdErr = stdErr <= 0.0? 0.0: sqrt(stdErr);
                (*resultss)[i]->storeRequestResult(stdErrRequest, stdErr);
            }
        }
        if (numBlocks > 1 && useLR){
            // calculate the std error for greeks 
            aggregateStdErrForGreek_multi(sens,
                                          CALCULATE,
                                          nbSubSamples,
                                          resultss.get(),
                                          greeksSumSqrSoFar);
        }

        /* will do COMPUTE_INDEX/COMPUTE_ESTIMATE
           and sort out debug values for price/compute time */  

        control->endTiming(this, instruments, resultss);

        // finally tidy up
        if (removePriceTime){
            OutputRequest priceTimeRequest(OutputRequest::PRICE_TIME);
            for (int i = 0; i < resultss->size(); ++i)
                (*resultss)[i]->removeGreek(priceTimeRequest.getPacketName(),
                                            OutputNameSP(new OutputName(priceTimeRequest.
                                                                        getRequestName())));
        }

        return resultss;
    } 
    catch (exception& e){
        throw ModelException(e, method);
    }
}


void MonteCarlo::Price(CInstrument*  instrument,
                       CControl*     control, 
                       CResults*     results)
{
    static const string method = "MonteCarlo::Price";
    try{
        auto_ptr<IMCProduct> prodAP(createProduct(instrument));
        prodAP->validate();
        // we call ourselves recursively
        bool isPricing = control->isPricing();
        IMCQuickXGamma* qckXGamma = 0;
        if (isPricing){
            // determine quick greek settings before pricing
            IMCQuickGreeks* qckGreeks = 
                prodAP->quickGreeksSupported();
            qckXGamma = prodAP->quickXGammaSupported();
            // This assigns a value to cachedQckGreeks...
            quickGreeksChoice(control, qckGreeks, qckXGamma);
            // instantiate MCPricing object
            createPricingObject(control, prodAP.get());
        }
        MCPathGeneratorSP futurePathGen(prodAP->price(this,
                                                      control,
                                                      results));
        if (isPricing && futurePathGen.get()){
            // if on pricing run (and there is a future) do quick x gamma ...
            // do we do internal tweaks? Has it been requested and do we
            // support it?
            if (cachedQckGreeks & IMCProduct::QUICK_X_GAMMA_BIT){
                // Then, is x gamma requested?
                // to do - move into Control ?
                SensitivityArrayConstSP sens = control->getSens();
                CResultsArraySP resultss(new CResultsArray(1, CResultsSP::attachToRef(results)));
                for (int i = 0; i < sens->size(); i++){
                    const Sensitivity* s = (*sens)[i].get();
                    if (ISpotCrossDerivative::TYPE->isInstance(s)){
                        control->recordPriceTime(resultss);//doesn't overwrite
                        // must work on clone of model
                        smartPtr<MonteCarlo> tweakModel(copy(this));
                        // do it
                        tweakModel->priceXGamma(qckXGamma, futurePathGen, 
                                                instrument,
                                                control, results, s);
                    }
                }
            }
        }
    } catch (exception& e){
        throw ModelException(e, method
#ifdef DEBUG_INST_XGAMMA
                             ,"Failed to nbIter: "+
                             Format::toString(debugNumIter)
#endif
            );
    }
}

MarketDataFetcherSP MonteCarlo::createMDF() const {
    MarketDataFetcherSP mcMDF = pathConfig->marketDataFetcher();
    MDFUtil::setUseCurrencyBasis(*mcMDF, useCcyBasis);
    return mcMDF;
}

void MonteCarlo::updateMarketBoundaryProb(CInstrumentSP instr, TweakGroupSP tweakGroup, CClassConstSP clazz)
{
    IScalarRiskPropertyConstSP prop = fieldRiskProperty::scalar(
        clazz,
        IFieldTweak::bulk(
        FieldPath::SP("pdfBoundaryProb"), FieldPathSP(),
        CDouble::SP(1.),
        IFieldTweak::IOperator::numeric(
        TweakFunction::setter(),
        InfiniteRange::InfiniteRange(),
        true, false)),
        true);

    double newValue = pdfBoundaryFactor * CVolProcessedBS::LEGACY_PDF_BOUNDARY_PROB;
    prop->axisFor(OutputNameConstSP())->hypothesis(newValue)->applyTo(tweakGroup);
}

/** Invoked after instrument and model have their market data. Passes on the
    call to the MCPathConfig object (which by default is empty)*/
void MonteCarlo::getMarket(const MarketData*  market,
                           IInstrumentCollectionSP instrument) {
    pathConfig->getMarket(this, market, instrument);

    if (pdfBoundaryFactor != 1)
    {
        QLIB_VERIFY(pdfBoundaryFactor >= 1, "pdfBoundaryFactor must be at least 1");
        double inv = 1.0/CVolProcessedBS::LEGACY_PDF_BOUNDARY_PROB;
        QLIB_VERIFY(Maths::isPositive(inv - pdfBoundaryFactor), "pdfBoundaryFactor must be less than " + Format::toString(inv));
        for (int i = 0;i<instrument->size();i++)
        {
            CInstrumentSP instr = (*instrument)[i];
            TweakGroupSP tweakGroup = TweakGroupSP(new TweakGroup(instr, IModelSP::attachToRef(this)));
            updateMarketBoundaryProb(instr, tweakGroup, IPDFBoundaryProb::TYPE);
        }
    }    
}


/** when to stop tweaking */
DateTime MonteCarlo::endDate(const CInstrument*  instrument,
                             const Sensitivity*  sensitivity) const{
    auto_ptr<IMCProduct> prodAP(createProduct(instrument));
    return prodAP->endDate(sensitivity);
}

/** indicates whether the MonteCarlo infrastructure can support 
    Vega Matrix for supplied instrument */
bool MonteCarlo::vegaMatrixSupported(CInstrument* instrument) const {
    // start by checking if pathConfig is happy, e.g. not Implied
    bool supported = pathConfig->vegaMatrixSupported();
    /* Then only do it if the associated product also supports the
       IMCProductLN interface */
    if (supported &&
        IMCIntoProduct::TYPE->isInstance(instrument)){
        try{
            IMCProductSP prodAP(createProduct(instrument));
            // this is probably not the most efficient way of doing this
            if (dynamic_cast<const IMCProductLN*>(prodAP.get())){
                supported = true;
            } else {
                supported = false;
            }
        } catch (exception& e){
            throw ModelException(e, "MonteCarlo::vegaMatrixSupported", 
                                 "Failed to createProduct for "
                                 "instrument of type "+
                                 instrument->getClass()->getName());
        }
    }
    return supported;
}

IModel::WantsRiskMapping MonteCarlo::wantsRiskMapping() const {
    return pathConfig->wantsRiskMapping();
}

/** returns all strikes on the vol surface to which 
    the supplied instrument is sensitive. Use vegaMatrixSupported to see
    if this routine should work */

DoubleArraySP MonteCarlo::getSensitiveStrikes(
    CInstrument*             instrument,
    const OutputNameConstSP& outputName) const {
    static const string method("MonteCarlo::getSensitiveStrikes");
    if (!IMCIntoProduct::TYPE->isInstance(instrument)){
        throw ModelException(method, "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support MonteCarlo::IntoProduct");
    }
    try{
        DoubleArraySP sensStrikes(new DoubleArray()); // what's returned
        // 1. Create product
        auto_ptr<IMCProduct> prodAP(createProduct(instrument));
        prodAP->validate();
        // 2. Create past path generator and update product
        MCPathGeneratorSP pastPathGen(pathConfig->pastPathGenerator(prodAP.get()));
        prodAP->pathGenUpdated(pastPathGen.get());
        // 3. see if there is anything left to simulate
        if (prodAP->hasFuture()){ // no sim => no vega matrix sensitivity
            // 4. do past if needed

            // Now if we're being invoked via say COMPUTE_INDEX, we're
            // not actually pricing and this will toot, so make sure
            // we don't die at least.
            // Did throw an exception, but NT-opt cores if there's an exception
            // somewhere between comments (4) and (5)
            if (pastPathGen->hasPast() && pricing.get()) {
                 // see if our product implements the IMCQuickGreeks interface
                IMCPricesSP priceForPast(
                    pricing->createPastPrices(this, prodAP.get()));
                pastPathGen->generatePath(0); // may well do nothing
                prodAP->payoff(pastPathGen.get(), *priceForPast);
            }

            // 5. cast to IMCProductLN
            const IMCProductLN* prodLN = 
                dynamic_cast<const IMCProductLN*>(prodAP.get());
            if (!prodLN){
                throw ModelException(method, "Product does not implement "
                                     "IMCProductLN interface");
            }
            // 6. Set up ourSensitiveStrikeDescriptor
            SensitiveStrikeDescriptor sensStrikeDesc;
            sensStrikeDesc.forwardOnly = false;
            // 7. loop across assets getting vol interps
            for (int iAsset=0; iAsset < prodAP->getNumAssets(); iAsset++) {
                // for each asset have an array of vol requests
                CVolRequestLNArray volRequests(
                    prodLN->getVolInterp(pastPathGen.get(), iAsset));
                // 8. loop across each path
                for (int iPath = 0; iPath < volRequests.size(); iPath++){
                    prodAP->getMultiFactors()->
                        assetGetSensitiveStrikes(iAsset,
                                                 volRequests[iPath].get(),
                                                 outputName,
                                                 sensStrikeDesc,
                                                 sensStrikes);
                }
            }
        }
        return sensStrikes;
    } catch (exception& e){
        throw ModelException(e, method, "Failed to calculate sensitive "
                             "strikes for instrument of type "+
                             instrument->getClass()->getName());
    }
}

/** Populate a Results object ready for quick X gamma calculation */
ResultsSP MonteCarlo::createXGammaResults(
    const string&                 ccyName,
    const OutputNameSP&           name1,
    const OutputNameSP&           name2,
    double                        subGamma12,
    double                        subGamma21,
    const ScalarShiftSP&          theShift){
    ResultsSP localResults(new Results());
    // we get product to report change in price for path
    localResults->storePrice(0.0, ccyName);
    // we also need to store some delta's as otherwise they will get
    // calculated
    localResults->storeScalarGreek(0.0, // value is irrelevant
                                   Delta::NAME, name1);
    localResults->storeScalarGreek(0.0, // value is irrelevant
                                   Delta::NAME, name2);
    // then store our pair of SubGammas
    // sub gamma for asset i when doing cross gamma wrt i and j
    localResults->storeScalarGreek(subGamma12,
                                   Delta::SECOND_ORDER_NAME, name1);
    // delta shift size
    string shiftSizeID = theShift->getSensOutputName()+
        Results::SHIFT_SIZE_POSTFIX;
    localResults->storeScalarGreek(theShift->getShiftSize(),
                                   shiftSizeID, name1);
    // sub gamma for asset j when doing cross gamma wrt i and j
    localResults->storeScalarGreek(subGamma21,
                                   Delta::SECOND_ORDER_NAME, name2);
    // delta shift size
    localResults->storeScalarGreek(theShift->getShiftSize(),
                                   shiftSizeID, name2);
    return localResults;
}

/** Manual override for calculation of cross gammas. Assumes that
    quick x gamma has been requested, and that the product supports
    it, the instrument has been priced and that the createOrigPrices()
    was called for that run. 'this' is modified - so should be cloned 
    before this method is invoked */
void MonteCarlo::priceXGamma(
    IMCQuickXGamma*            qckXGamma,
    const MCPathGeneratorSP&  futurePathGen,
    CInstrument*                        instrument,
    Control*                            control, 
    CResults*                           results,
    const Sensitivity*                  xGammaSens){
    static const string method("MonteCarlo::priceXGamma");
    try{
        // Clone xGamma as we'll be altering it
        SensitivitySP tweak(copy(xGammaSens));
        // Get hold of what we need to be shifting (and also what we're going
        // to be shifting (Needed for createXGammaPrices() method). 1st do cast
        const ISpotCrossDerivative& crossDeriv = 
            dynamic_cast<const ISpotCrossDerivative&>(*xGammaSens);
        ScalarShiftArray theShifts(crossDeriv.definingShifts(instrument,
                                                             this,
                                                             control));
        if (theShifts.size() != 2){
            throw ModelException(method, "Num shifts != 2");
        }
        // are the shifts shifting the same thing (eq spot vs fx spot)
        bool sameShiftType = theShifts[0]->shiftInterface() == 
            theShifts[1]->shiftInterface();
        if (!sameShiftType){
            throw ModelException(method, "Method does not support "
                                 "different shift types yet");
        }
        /* Need to set up two loops. In inner loop need to set
           relevant pair of names to calculate in SensControl and
           then call normal calculate methods. */
        OutputNameArrayConstSP names1 = theShifts[0]->overrideNames();
        OutputNameArrayConstSP names2 = theShifts[1]->overrideNames();

        // see if delta has been requested
        SensitivitySP deltaRequested(control->
                                     sensitivityRequested(theShifts[0]));
        /* Delta & FXCrossGamma will want the deltas. FXCrossGamma will
           delete them if Delta is not requested */
        bool doRealDelta = 
            !(!deltaRequested &&
              !control->sensitivityRequested(FXCrossGamma::TYPE));
        MCPricesGeneralSP genPrices = DYNAMIC_POINTER_CAST<MCPricesGeneral>(origPrices);
        /* start by calculating what we call SubGamma for each
           stock (see SubGamma for a description) */
        refCountPtr<SubGamma> subGamma(new SubGamma(pricing,
                                                    doRealDelta, 
                                                    genPrices,
                                                    names1->size()));
        // where we will do the sub gamma calculation
        ResultsSP subGammaResults;
        if (doRealDelta){
            subGammaResults = ResultsSP::attachToRef(results);
        } else {
            // create bogus results object for SubGamma calc
            subGammaResults = ResultsSP(new Results());
            subGammaResults->storePrice(0.0, "");
        }
        // then calculate all the sub gammas - using our pricing object
        pricing = subGamma;
        subGamma->calculateSubGamma(this, qckXGamma, futurePathGen, control, 
                                    instrument, subGammaResults.get(),
                                    theShifts);

        // now set the MCPricing object for Cross Gamma
        pricing = MCPricingSP(new QXGamma(origPrices));
        // for cross gamma use same set of names (use in pairs though)
        OutputNameArraySP namesForTweak(new OutputNameArray(2));
        const string& ccyName = results->getCcyName();
        // for cross gamma use same set of names (use in pairs though)
        for (int i = 0; i < names1->size(); i++){
            // namesForTweak used by CrossGamma (only calc one at a time)
            (*namesForTweak)[0] = (*names1)[i];
            // theShifts passed to Quick X Gamma-says what is being shifted
            theShifts[0]->setMarketDataName((*names1)[i]);
            const string& packet = tweak->getSensOutputName();
            for (int j = sameShiftType? i+1: 0; j < names2->size(); j++){
                IObjectSP xGamma; // the result
                OutputNameSP outName(new OutputName((*names1)[i].get(), 
                                                    (*names2)[j].get()));
                if (!subGamma->completed[i] || !subGamma->completed[j]){
                    xGamma = IObjectSP(new Untweakable(
                        string("Gamma calculation failed")));
                } else {
                    // to do: for fx cross gamma, need to skip when we
                    // are tweaking fx spot and equity spot of the
                    // same asset (no longer cross gamma)
                    (*namesForTweak)[1] = (*names2)[j];
                    theShifts[1]->setMarketDataName((*names2)[j]);
                    qckXGamma->setPricesForXGamma(origPrices.get(),
                                                  futurePathGen.get(),
                                                  theShifts);
                    // store pair of names in sens control
                    tweak->storeOverrideNames(namesForTweak);
                    // create a results object (must not use 
                    // original:wrong gammas)
                    ResultsSP localResults(createXGammaResults(
                        ccyName, (*names1)[i], (*names2)[j],
                        subGamma->gammaVals[i][j],
                        subGamma->gammaVals[j][i],
                        theShifts[0]));
                    // then calculate
                    control->calculateSens(tweak, this, instrument,
                                           localResults.get());
                    // copy over our results - stored in standard way
                    // important to manipulate objects here in case calc 
                    // failed
                    xGamma = IObjectSP(localResults->
                                       retrieveGreek(packet, outName).clone());
                }
                results->storeGreek(xGamma, packet, outName);
                outName = OutputNameSP(new OutputName((*names2)[j].get(), 
                                                      (*names1)[i].get()));
                results->storeGreek(IObjectSP(xGamma.clone()), packet,
                                    outName);
            }
        }
    } catch (exception& e){
        // think we should fail here as something catastrophic has happened - 
        // don't want greeks to be calculated in normal way
        throw ModelException(e, method);
    }
}

MonteCarlo::MonteCarlo(CClassConstSP clazz): 
    CModel(clazz), quickGreeks(QUICK_GREEK_DEFAULT), 
    lrGreeks(LR_GREEK_DEFAULT), useStateVars(false), useCcyBasis(false),
    simMode(SIM_MODE_DEFAULT), cachePaths(CACHE_PATH_DEFAULT), cachedQckGreeks(0), 
    cachedCachePaths(false), pdfBoundaryFactor(1) {}

// for reflection
MonteCarlo::MonteCarlo():
    CModel(TYPE),     
    nbIter(0), nbSubSamples(0), quickGreeks(QUICK_GREEK_DEFAULT), 
    lrGreeks(LR_GREEK_DEFAULT), useStateVars(false), useCcyBasis(false),
    simMode(SIM_MODE_DEFAULT), cachePaths(CACHE_PATH_DEFAULT), cachedQckGreeks(0),
    cachedCachePaths(false), pdfBoundaryFactor(1){}

void MonteCarlo::validatePop2Object() {
    static const string method = "MonteCarlo::validatePop2Object";
    try {
        if (nbIter <= 1) {
            throw ModelException(method, "Number of iterations ("+
                                 Format::toString(nbIter) + ") should "
                                 "be more than 1");
        }
        if (nbSubSamples < 1) {
            throw ModelException(method, "Number of sub-samples ("+
                                 Format::toString(nbSubSamples) + ") should "
                                 "be at least 1");
        }
        if ((nbIter % (2*nbSubSamples)) != 0) {
            throw ModelException(method, "Twice number of sub-samples ("+
                                 Format::toString(nbSubSamples) + ") should "
                                 "divide the number of iterations (" + 
                                 Format::toString(nbIter) + ")");            
        }          
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

class MonteCarloHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(MonteCarlo, clazz);
        SUPERCLASS(CModel);
        EMPTY_SHELL_METHOD(defaultMonteCarlo);
        FIELD(nbIter,            "NbIter");
        FIELD(nbSubSamples,      "NbSubSamples");
        FIELD(pathConfig,               "Path generator configuration");
        FIELD(quickGreeks,       "Use quick greeks algorithm: yes, no"
                     " or default");
        FIELD_MAKE_OPTIONAL(quickGreeks);
        FIELD(lrGreeks,       "Use Likelihood Ratio algorithm for "
                     "greeks where possible: yes, no or default");
        FIELD_MAKE_OPTIONAL(lrGreeks);
        FIELD(simMode, "simulation mode: slice, path, or pathOpt");
        FIELD_MAKE_OPTIONAL(simMode);
        FIELD(cachePaths,       "Cache paths between tweaks: yes, no,"
                     " default or limit in Mb");
        FIELD_MAKE_OPTIONAL(cachePaths);
        FIELD_NO_DESC(cachedQckGreeks);
        FIELD_MAKE_TRANSIENT(cachedQckGreeks);
        FIELD_NO_DESC(cachedCachePaths);
        FIELD_MAKE_TRANSIENT(cachedCachePaths);
        FIELD(pathConfigPostPricing, "Path config after pricing");
        FIELD_MAKE_TRANSIENT(pathConfigPostPricing);
        FIELD(useStateVars, "If true, use state variables; don't otherwise");
        FIELD_MAKE_OPTIONAL(useStateVars);
        FIELD(useCcyBasis, "use currency basis (default = false)");
        FIELD_MAKE_OPTIONAL(useCcyBasis);

        FIELD(pdfBoundaryFactor, "used to set pdf lower limit in VolProcessedBS.");
        FIELD_MAKE_OPTIONAL(pdfBoundaryFactor);

        clazz->setPublic(); // make visible to EAS/spreadsheet
        // get around C++ initialisation issues
        const_cast<string&>(MonteCarlo::DEBUG_GREEK_STD_ERR) = 
            Results::DEBUG_PACKETS_PREFIX + "STD_ERR_";
    }

    static void loadIntoProduct(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER_INTERFACE(MonteCarlo::IIntoProduct, clazz);
        EXTENDS(Model::IModelIntoProduct);
    }

    static IObject* defaultMonteCarlo(){
        return new MonteCarlo();
    }
};

CClassConstSP const MonteCarlo::TYPE = CClass::registerClassLoadMethod(
    "MonteCarlo", typeid(MonteCarlo), MonteCarloHelper::load);

CClassConstSP const MonteCarlo::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("MonteCarlo::IIntoProduct",
                                    typeid(MonteCarlo::IIntoProduct), 
                                    MonteCarloHelper::loadIntoProduct);



/** Override default clone to do shallow copy of 'origPrices' */
IObject* MonteCarlo::clone() const{
    IObjectSP theCopy(CModel::clone());
    MonteCarlo* myCopy = DYNAMIC_CAST(MonteCarlo, theCopy.get());
    myCopy->origPrices = origPrices; // shallow copy
    // the next one is shallow copied - this means that even if the model
    // is cloned in the delta/cross gamma tweaking we are still ok
    myCopy->pricing = pricing; 
    return theCopy.release();
}


/** Determines whether any sort of quick greeks has been requested and 
    verifies that we support it */
void MonteCarlo::quickGreeksChoice(
    Control*                 control,
    IMCQuickGreeks* qckGreeksImpl,
    IMCQuickXGamma* qckXGammaImp) {
    static const string method("MonteCarlo::quickGreeksChoice");
    // to allow us extra control
    if (CString::equalsIgnoreCase(QUICK_GREEK_AND_X_GAMMA, quickGreeks)){
        cachedQckGreeks = IMCProduct::QUICK_GREEKS_BIT | 
            IMCProduct::QUICK_X_GAMMA_BIT;
    } else if (CString::equalsIgnoreCase(QUICK_GREEK_NO_X_GAMMA, quickGreeks)){
        // switch on just QUICK_GREEKS
        cachedQckGreeks = IMCProduct::QUICK_GREEKS_BIT;
    } else if (CString::equalsIgnoreCase(QUICK_X_GAMMA_ONLY, quickGreeks)){
        // switch on just QUICK_X_GAMMA
        cachedQckGreeks = IMCProduct::QUICK_X_GAMMA_BIT;
    } else if (CString::equalsIgnoreCase(QUICK_GREEK_DEFAULT, quickGreeks)){
        /* make it depend upon whether the instrument supports it or not */
        cachedQckGreeks = (qckGreeksImpl? IMCProduct::QUICK_GREEKS_BIT: 0) |
            (qckXGammaImp? IMCProduct::QUICK_X_GAMMA_BIT: 0);
    } else {
        // YES means whatever is available - same as default
        cachedQckGreeks = CBool::fromString(quickGreeks)? 
            ((qckGreeksImpl? IMCProduct::QUICK_GREEKS_BIT: 0) |
             (qckXGammaImp? IMCProduct::QUICK_X_GAMMA_BIT: 0)): 0;
    }
    if ((cachedQckGreeks & IMCProduct::QUICK_GREEKS_BIT) && !qckGreeksImpl){
        throw ModelException(method, "Product doesn't support "
                             "quick greeks");
    }
    if ((cachedQckGreeks & IMCProduct::QUICK_X_GAMMA_BIT) && !qckXGammaImp){
        throw ModelException(method, "Product doesn't support "
                             "quick x gamma");
    }
    // Independent - possible caching at product level
    if (cachePaths.find(MonteCarlo::CACHE_PROD_EXTENSION) != string::npos) {
        cachedQckGreeks |= IMCProduct::CACHE_PRODUCT_BIT;
    }

    if (cachedQckGreeks){
        // now study control and see if we should bother
        SensitivityArrayConstSP sens = control->getSens();
        if (sens->empty()){
            cachedQckGreeks = 0; // don't bother
        } else if (cachedQckGreeks & IMCProduct::QUICK_X_GAMMA_BIT){
            // see if x gamma requested - to do: put in Control.cpp
            for (int i = 0; i < sens->size(); i++){
                if (ISpotCrossDerivative::TYPE->isInstance((*sens)[i].get())){
                    return;
                }
            }
            // switch off x gamma
            cachedQckGreeks &= ~IMCProduct::QUICK_X_GAMMA_BIT;
        }
    }
}

/** "bread & butter" MonteCarlo that can be captured in Pyramid using 
     current IMS */
class MonteCarloLNDefault: public MonteCarlo {
public:
    static CClassConstSP const TYPE;
    
private:
    MonteCarloLNDefault():MonteCarlo(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(MonteCarloLNDefault, clazz);
        SUPERCLASS(MonteCarlo);
        EMPTY_SHELL_METHOD(defaultMonteCarloLNDefault);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultMonteCarloLNDefault(){
        return new MonteCarloLNDefault();
    }
};

CClassConstSP const MonteCarloLNDefault::TYPE = 
CClass::registerClassLoadMethod(
    "MonteCarloLNDefault", typeid(MonteCarloLNDefault), load);

/** "bread & butter" MonteCarlo that can be captured in Pyramid using 
     current IMS */
class MonteCarloImpliedDefault: public MonteCarlo {
public:
    static CClassConstSP const TYPE;
    
private:
    MonteCarloImpliedDefault():MonteCarlo(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(MonteCarloImpliedDefault, clazz);
        SUPERCLASS(MonteCarlo);
        EMPTY_SHELL_METHOD(defaultMonteCarloImpliedDefault);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultMonteCarloImpliedDefault(){
        return new MonteCarloImpliedDefault();
    }
};

CClassConstSP const MonteCarloImpliedDefault::TYPE = 
CClass::registerClassLoadMethod(
    "MonteCarloImpliedDefault", typeid(MonteCarloImpliedDefault), load);


/////////////////////////////////////////////////////////////////////////
// Addin utility to provide list of all dates (for use in PastValues)
/////////////////////////////////////////////////////////////////////////
class SimDates : public CObject {
public:
    static CClassConstSP const TYPE;

    // addin params
    smartPtr<MonteCarlo>  model;
    CInstrumentSP         instrument;
    smartPtr<MarketData>  market;
    DateTimeSP            fromDate;
    
    static IObjectSP getSimDates(SimDates* params) {
        try {
            // copy in case local changes
            CInstrumentSP inst(copy(params->instrument.get()));

            // turn off the "missing future past values exception". 
            // This is a hack but the proper solution requires
            // moving this function so there is visibility of GenericNFBase class while
            // retaining visibility of the MonteCarlo->createProduct. 
            CDataDictionary* ddInst = CDataDictionary::pop2DataDict(inst);
            IObjectSP pastValues = ddInst->get("pastValues");
            CDataDictionary* ddPastValues = CDataDictionary::pop2DataDict(pastValues);
            ddPastValues->put("throwMissingFutureDatesException", false);
            ddInst->put("pastValues", IPastValuesSP::dynamicCast(ddPastValues->pop2Object()));
            inst = CInstrumentSP(CInstrumentSP::dynamicCast(ddInst->pop2Object()));
            delete ddPastValues;
            delete ddInst;

            // do this after the hack since we wish to fill the final 'inst' with market data
            if (params->market.get())
            {
                // Market can be optional since the inst may have been built 
                // complete with market data. If not ...
                params->model->getInstrumentAndModelMarket(params->market.get(), inst.get());
            }

        IMCProductSP prodAP(
                params->model->createProduct(inst.get()));
            const DateTimeArray& simDates = prodAP->getSimSeries()->getAllDates();
            const DateTimeArray& refSimDates = 
                prodAP->getRefLevel()->getAllDates();
            DateTimeArray allDates(DateTime::merge(simDates, refSimDates));
            
            if (!!params->fromDate) {
                return IObjectSP(params->fromDate->getFutureDates(allDates).clone());
            } else {
                return IObjectSP(allDates.clone());
            }
        } catch (ModelException& e){
            // Dump stack trace to error log
            e.printStackTrace();
            return IObjectSP(CString::create("Failed - see error log!"));
        } catch (exception& ) {
            return IObjectSP(CString::create("Failed!"));
        }
    }

    SimDates(): CObject(TYPE) {}

    static IObject* defaultSimDates(){
        return new SimDates();
    }

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(SimDates, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultSimDates);
        FIELD(model,         "MonteCarlo instance");
        FIELD(instrument,    "Instrument");
        FIELD(market,        "Market");
        FIELD_MAKE_OPTIONAL(market);
        FIELD(fromDate,      "Date from which to report simulation dates");
        FIELD_MAKE_OPTIONAL(fromDate);
        Addin::registerInstanceObjectMethod(
            "MC_GET_SIM_DATES",
            Addin::RISK,
            "All dates for Monte Carlo simulation of this instrument",
            TYPE,
            false,
            Addin::expandMulti,
            (Addin::ObjMethod*)getSimDates);
    }
};

CClassConstSP const SimDates::TYPE = 
CClass::registerClassLoadMethod(
    "MonteCarloSimDates", typeid(SimDates),load);

DRLIB_END_NAMESPACE
