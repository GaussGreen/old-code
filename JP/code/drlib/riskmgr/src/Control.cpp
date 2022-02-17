
#include "edginc/config.hpp"
#define QLIB_CONTROL_CPP
#include "edginc/CInstrumentCollection.hpp"
#include "edginc/Delta.hpp"
#include "edginc/CommandLineParams.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/VegaParallel.hpp"
#include "edginc/VegaParallel2Sided.hpp"
#include "edginc/DDeltaDVol.hpp"
#include "edginc/RhoParallel.hpp"
#include "edginc/RhoBorrowParallel.hpp"
#include "edginc/VegaSkewParallel.hpp"
#include "edginc/MuParallel.hpp"
#include "edginc/ImpliedVol.hpp"
#include "edginc/RootTimeVega.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/RiskMgr.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/RiskQuantityEvaluator.hpp"
#include "edginc/SensitiveStrikes.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/Scenario.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/SpreadSheetMode.hpp"
#include "edginc/ModelFilter.hpp"
#include "edginc/Shuffle.hpp"
#include "edginc/MarketData.hpp"

DRLIB_BEGIN_NAMESPACE

 // use theoretical value of assets
const string Control::ASSET_THEO = "ASSET_THEO";
// use mark-to-market value of assets
const string Control::ASSET_MTM = "ASSET_MTM"; 
// do both
const string Control::ASSET_THEO_MTM = "ASSET_THEO_MTM";

Control::~Control(){}

static double computeEstimate(IModel*                 model, 
                              IInstrumentCollectionSP instruments,
                              double                  priceTime) {
    static const string method("Control::computeEstimate");
    try {
        smartPtr<PriceCounter> priceCounter(model->priceCounter());
        CControlSP ctrl(SensitivityFactory::megaControl(false));

        CResultsArraySP dontcare(priceCounter->RunMulti(instruments, ctrl.get()));

        return (priceCounter->pricings()*priceTime);
    } 
    catch (exception& e) {
        throw ModelException(e, method);
    }
}
   
// work out the compute index for Pyramid - a series of lookups based on
// how large the compute estimate is and how many (equity) underlings
int Control::computeIndex(IModel*                 model, 
                          IInstrumentCollectionSP insts,
                          double                  priceTime,
                          Results*                results) {
    static const string method("Control::computeIndex");
    try {

        // In Pyramid, all flavours of SPI are set up with compute index 75
        // to trigger a certain overnight Merlin configuration. To ensure
        // that SPIs are set correctly, hardwire inside QLib rather than
        // rely on some external process.
        CClassConstSP spiType = CClass::forName("SyntheticPortfolioInsurance");

        //see if there's any SPI
        bool spi = false;
        for (int i = 0; i < insts->size() && !spi; i++) {
            spi = spiType->isInstance((*insts)[i]);
        }
        if (spi) {
            return 75;
        }

        // otherwise, go the regular route 
        OutputRequest* request =
            requestsOutput(OutputRequest::COMPUTE_ESTIMATE);
        int cpu;
        if (!request){
            cpu = (int)computeEstimate(model, insts, priceTime);
        } else {
            IObjectConstSP computeEst(results->retrieveRequestResult(
                                          OutputRequest::COMPUTE_ESTIMATE));
            cpu = (int) (dynamic_cast<const 
                         CDouble&>(*computeEst)).doubleValue();
        }
        int cpuidx;

        // figure out if this is 1 or N factor (talking equity here)
        // count number of things tweaked for delta
        Delta delta(Delta::DEFAULT_SHIFT);
        OutputNameArrayConstSP names(delta.names(insts.get()));

        bool isNFactor = names->size() > 1;

        if (!isNFactor) {
            if (cpu <= 60) {
                // less than a minute
                cpuidx = 0;
            }
            else if (cpu <= 300) {
                // less than 5 minutes
                cpuidx = 5;
            }
            else {
                cpuidx = 10;
            }
        }
        else {
            if (cpu <= 60) {
                // less than a minute
                cpuidx = 20;
            }
            else if (cpu <= 600) {
                // less than 5 minutes
                cpuidx = 25;
            }
            else if (cpu <= 3600) {
                // less than an hour
                cpuidx = 30;
            }
            else if (cpu <= 7200) {
                // between 1 and 2 hours
                cpuidx = 50;
            }
            else if (cpu <= 28800) {
                // between 2 and 8 hours
                cpuidx = 60;
            }
        }
        // then catch all for global killers, single factor or otherwise
        if (cpu > 28800) {
            // we're talking "global killers" here > 8 hours
            cpuidx = 70;
        }
    
        return cpuidx;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// here's the control proper

/* Takes copy of sens array. sens and outReqs can be null  */
Control::Control(const SensitivityArrayConstSP&   sens,
                 const OutputRequestArrayConstSP& outReqs,
                 bool                             writeToFile,
                 const string&                    fileName): 
    CObject(TYPE), fileName(fileName), assetPriceSource(ASSET_THEO),
    currentInstrumentIndexInCollection(0), 
    useTheoAssetPrice(true),
    isPricingRun(true), scaleResults(false), writeToFile(writeToFile),
    skipAllShifts(false){
    this->sens = SensitivityArraySP(!sens? 0: sens.clone());
    this->outputRequests = OutputRequestArraySP(!outReqs? 0: outReqs.clone());
    validatePop2Object();
}

/* Takes copy of sens array. sens and outReqs can be null  */
Control::Control(const SensitivityArrayConstSP&   sens,
                 const OutputRequestArrayConstSP& outReqs,
                 bool                             writeToFile,
                 const string&                    fileName,
                 bool                             skipShifts): 
    CObject(TYPE), fileName(fileName), assetPriceSource(ASSET_THEO),
    currentInstrumentIndexInCollection(0), 
    useTheoAssetPrice(true),
    isPricingRun(true), scaleResults(false), writeToFile(writeToFile),
    skipAllShifts(skipShifts){
    this->sens = SensitivityArraySP(!sens? 0: sens.clone());
    this->outputRequests = OutputRequestArraySP(!outReqs? 0: outReqs.clone());
    validatePop2Object();
}

/* Takes copy of sens array. sens and outReqs can be null  */
Control::Control(const SensitivityArrayConstSP&   sens,
                 const OutputRequestArrayConstSP& outReqs,
                 bool                             writeToFile,
                 const string&                    fileName,
                 const string&                    assetPriceSource): 
    CObject(TYPE), fileName(fileName), assetPriceSource(assetPriceSource),
    currentInstrumentIndexInCollection(0), 
    useTheoAssetPrice(true),
    isPricingRun(true), scaleResults(false), writeToFile(writeToFile),
    skipAllShifts(false){
    this->sens = SensitivityArraySP(!sens? 0: sens.clone());
    this->outputRequests = OutputRequestArraySP(!outReqs? 0: outReqs.clone());
    validatePop2Object();
}

void Control::startTiming(){
    startTime = clock();
#ifdef UNIX
    startSec  = time(0);
#endif
}

void Control::endTiming(IModel* model, 
                        IInstrumentCollectionSP instruments,
                        CResultsArraySP resultss){
    // store timings if required
#ifdef UNIX
    time_t  endSec  = time(0);
    // UNIX overflows after 2147 seconds
    double calcTime = difftime(endSec, startSec);
#else
    clock_t endTime = clock();
    double calcTime = (double)(endTime-startTime)/CLOCKS_PER_SEC;
#endif

    if (resultss->size() == 0) return;

    CResultsSP results = (*resultss)[0];

    OutputRequest* request = requestsOutput(OutputRequest::CALC_TIME);
    if (request){
        results->storeRequestResult(request, calcTime);
    }
    OutputNameSP calcTimeName(new OutputName(OutputRequest::CALC_TIME));
    // store it anyway in the debug section
    results->storeScalarGreek(calcTime, Results::DEBUG_PACKET,
                              calcTimeName);
    // typically priceTime will have already been set and this will just
    // be a read 
    double priceTime = recordPriceTime(resultss);
    // must do COMPUTE_ESTIMATE before COMPUTE_INDEX
    if ((request = requestsOutput(OutputRequest::COMPUTE_ESTIMATE))){
        try{
            double cpu = computeEstimate(model, instruments, priceTime);
            results->storeRequestResult(request, cpu);
        } catch (exception& e){
            results->storeRequestResult(request, 
                                        IObjectSP(new Untweakable(e)));
        }
    }
    if ((request = requestsOutput(OutputRequest::COMPUTE_INDEX)) &&
        !request->getHasFinished()){
        try{
            // in case instrument has reported its own index
            int cpuidx = computeIndex(model, instruments, priceTime,
                                      results.get());
            results->storeRequestResult(request, cpuidx);
        } catch (exception& e){
            results->storeRequestResult(request, 
                                        IObjectSP(new Untweakable(e)));
        }
    }
} 

/** calculates price and sensitivities */
void Control::calculate(IModel* model,
                        CInstrument* instrument,
                        Results* results) {
    calculateMulti(model, IInstrumentCollection::singleton(instrument),
                   CResultsArraySP(
                       new CResultsArray(1, ResultsSP::attachToRef(results))));
}

/** calculates prices and sensitivities */
void Control::calculateMulti(IModel* model,
                             IInstrumentCollectionSP instruments,
                             CResultsArraySP resultss) {
    ASSERT(instruments->size() == resultss->size());

    static const string method = "Control::calculateMulti";

    startTiming();

    // then initial pricing
    reset(); // clear what's been calculated 

    tempPackets.resize(resultss->size());
    {for (int i = 0; i < tempPackets.size(); ++i) {
            tempPackets[i] = StringArraySP(new StringArray());
    }}

    try {
        instruments->Price(model, this, resultss);
        // handle all unfulfilled requests. They are all not applicable or
        // may be top level ones still to run like LEGAL_TERMS_FV
        handleUnfulfilledRequests(model, instruments, resultss);
   } catch (exception& e){
        throw ModelException(e, method);
    }
    isPricingRun = false;

    recordPriceTime(resultss); // save price time

    if (!CommandLineParams::hasParameter(CommandLineParams::Pong)) {

        // Do the new-style RiskQuantityFactorySensitivity's first (they are
        // only Sensitivity's for compatibility reasons)

        IRiskQuantityFactoryArraySP greeks(new IRiskQuantityFactoryArray());
        for (int g = 0; g < sens->size(); ++g) {
            IRiskQuantityFactorySP greek(
                dynamic_cast<IRiskQuantityFactory*>((*sens)[g].get()));
            if (!!greek) {
                greeks->push_back(greek);
            }
        }

        // Eval them all together, for speed

        if (!riskQuantityEvaluator) {
            // This is to catch the case where we haven't done getMarket
            riskQuantityEvaluator.reset(new RiskQuantityEvaluator(false));
        }

        riskQuantityEvaluator->storeResults(
            greeks,
            MultiTweakGroup::SP(instruments, IModelSP::attachToRef(model)),
            CControlSP::attachToRef(this),
            resultss);

        // Do remaining "normal" Sensitivity's

        for (int idx = 0; idx < sens->size(); idx++){
            currentSens = (*sens)[idx];
            if (!currentSens){
                throw ModelException(method, "NULL sensitivity control");
            }
            if (IRiskQuantityFactory::TYPE->isInstance(currentSens)) continue;
            try {
                for (int i = 0; i < resultss->size(); ++i) {
                    currentInstrumentIndexInCollection = i;
                    currentSens->calculateSens(model, (*instruments)[i].get(),
                                               this, (*resultss)[i].get());
                }
            } catch (ModelException& e){
                e.addMsg(method + 
                         ": Failed to calculate " +
                         currentSens->getSensOutputName());
                /* note: calculateSens method handles failure to calculate
                   a greek. If an exception propagates to here we must
                   report it up to the caller */
                throw ModelException(e, method);
            }
        }
        
        currentSens.reset();
        // then remove any temporary packets
        for (int i = 0; i < resultss->size(); ++i)
            for (int p = 0; p < tempPackets[i]->size(); p++)
                (*resultss)[i]->removePacket((*tempPackets[i])[p]);
    }


    // flush any caches from model
    model->flush();

    try {
        endTiming(model, instruments, resultss);
    }
    catch (exception& ) {
        // not fatal if these don't work
    }
}

/** calculates the supplied sensivity (the instrument must have already
    been priced). Calls to getCurrentSensitivity will return sens for the
    duration of the call */
void Control::calculateSens(const SensitivitySP&  sens,
                            IModel*               model,
                            CInstrument*          instrument, 
                            Results*              results){
    bool isPricingRunCopy = isPricingRun; // save
    SensitivitySP   currentSensCopy = currentSens; // save
    try{
        isPricingRun = false;
        currentSens = sens;
        sens->calculateSens(model, instrument, this, results);
        isPricingRun = isPricingRunCopy;
        currentSens = currentSensCopy;
    } catch (exception&){
        isPricingRun = isPricingRunCopy;
        currentSens = currentSensCopy;
        throw;
    }
}

/** calculates exposure characteristics for multiple instruments */
void Control::calculateMultiExposures(IModel* model,
        IInstrumentCollectionSP instruments,
        CResultsArraySP resultss) {
    static const string method = "Control::calculateExposures";

    ASSERT(instruments->size() == resultss->size());

    startTiming();

    // then initial pricing
    reset(); // clear what's been calculated 

    tempPackets.resize(resultss->size());
    {for (int i = 0; i < tempPackets.size(); ++i) {
        tempPackets[i] = StringArraySP(new StringArray());
    }}

    try {
        instruments->Price(model, this, resultss);
    } catch (exception& e){
        throw ModelException(e, method);
    }

    isPricingRun = false;

    recordPriceTime(resultss); // save price time

    // Do the new-style RiskQuantityFactorySensitivity's first (they are
    // only Sensitivity's for compatibility reasons)

    IRiskQuantityFactoryArraySP greeks(new IRiskQuantityFactoryArray());
    for (int g = 0; g < sens->size(); ++g) {
        IRiskQuantityFactorySP greek(
            dynamic_cast<IRiskQuantityFactory*>((*sens)[g].get()));
        if (!!greek) {
            greeks->push_back(greek);
        }
    }

    // Eval them all together, for speed
    if (!riskQuantityEvaluator) {
        // This is to catch the case where we haven't done getMarket
        riskQuantityEvaluator.reset(new RiskQuantityEvaluator(false));
    }

    riskQuantityEvaluator->storeExposureResults(
        greeks,
        MultiTweakGroup::SP(instruments, IModelSP::attachToRef(model)),
        CControlSP::attachToRef(this),
        resultss);

    // Do remaining "normal" Sensitivity's
    // Note that we just pretend that we are computing sensitivities in the
    // normal fashion and don't support any special behavior for old-style
    // Sensitivity's.

    for (int idx = 0; idx < sens->size(); idx++){
        currentSens = (*sens)[idx];
        if (!currentSens){
            throw ModelException(method, "NULL sensitivity control");
        }
        if (IRiskQuantityFactory::TYPE->isInstance(currentSens)) continue;
        try {
            for (int i = 0; i < resultss->size(); ++i) {
                currentInstrumentIndexInCollection = i;
                currentSens->calculateSens(model, (*instruments)[i].get(),
                    this, (*resultss)[i].get());
            }
        } catch (ModelException& e){
            e.addMsg(method + 
                ": Failed to calculate " +
                currentSens->getSensOutputName());
            /* note: calculateSens method handles failure to calculate
            a greek. If an exception propagates to here we must
            report it up to the caller */
            throw ModelException(e, method);
        }

        currentSens.reset();
        // then remove any temporary packets
        for (int i = 0; i < resultss->size(); ++i)
            for (int p = 0; p < tempPackets[i]->size(); p++)
                (*resultss)[i]->removePacket((*tempPackets[i])[p]);
    }

    // flush any caches from model
    model->flush();
}

/* returns reference to SensControlArray */
SensitivityArrayConstSP Control::getSens(){
    return sens;
}

void Control::shuffleSens(int seed) { 
    shuffle(*sens, basicRandGen(seed));
}

/** returns reference to OutputRequestArray */
OutputRequestArrayConstSP Control::getOutputRequests() const{
    return outputRequests;
}

/** checks whether a specific output has been requested. Do NOT free the
    pointer returned */
OutputRequest* Control::requestsOutput(const string& output) {
    std::map<string, OutputRequest*>::iterator pos;
    pos = outputRequestSet.find(output);
    if ( pos == outputRequestSet.end() ) {
        return NULL;
    } else {
        return pos->second;
    }
}

/** Deprecated.  */
bool Control::requestsOutput(const string& output, OutputRequest*& request){
    request = requestsOutput(output);
    return request? true: false;
}

/** Returns flag indicating whether a regression test file should 
    be created. */
bool Control::getWriteToFile() const{
    return writeToFile;
}

/** Returns filename for regression file. In particular, if
    validate is true an exception will be thrown if writeToFile is
    true but the file should not be created (eg because a file
    exists already [in a spreadsheet environment]). Also, if
    validate is true throws an exception if the fileName is empty
    and writeToFile is true */
const string& Control::getFileName(bool validate) const{
    static const string method("Control::getFileName");
    if (validate && writeToFile && SpreadSheetMode::isOn()) {
        if (fileName.empty()){
            throw ModelException(method,
                                 "No filename specified for write to file");
        }
        FILE *fp = fopen(fileName.c_str(), "r");
        if (fp) {
            fclose(fp);
            throw ModelException(method,
                                 "file (" + fileName + ") already exists");
        }
    }
    return fileName;
}

/** Switches off the write to file flag */
void Control::switchOffWriteToFile(){
    writeToFile = false;
}

void Control::validatePop2Object(){
    if (writeToFile && fileName == ""){
        throw ModelException("Control::validatePop2Object",
                             "No filename supplied for regression file");
    }
    if (!sens){
        sens = SensitivityArraySP(new SensitivityArray());
    }
    // create cache of requested outputs
    if (!outputRequests){
        outputRequests = OutputRequestArraySP(new OutputRequestArray(0));
    } else {
        for (int i = 0; i < outputRequests->size(); ++i) {
            OutputRequest* request = (*outputRequests)[i].get();
            string requestName = request->getRequestName();
            outputRequestSet[requestName] = request;
        }
    }
}

/** Override clone method to copy our extra data over */
IObject* Control::clone() const{
    // first clone all the registered fields
    IObjectSP  copy = IObjectSP(CObject::clone());
    // copy over non registered fields
    copy->validatePop2Object(); /* cannot just copy maps since they contain 
                                   references to objects inside 
                                   the outputRequests array */
    dynamic_cast<Control&>(*copy).startTime = startTime;
#ifdef UNIX
    dynamic_cast<Control&>(*copy).startSec = startSec;
#endif
    return copy.release();
}


/** hanlde all the requests which have not been calculated yet 
    Flawed design: results cannot be interrogated for existence, yet
    requests are not guaranteed unique. Since there's an addRequest()
    method we can't just make the array unique at construction.
    So, here we make sure we look only at requests which client
    code will see - i.e. those contained in the outputRequestSet. */
void Control::handleUnfulfilledRequests(IModel* model,
                                        IInstrumentCollectionSP instruments,
                                        CResultsArraySP resultss){
    ASSERT(instruments->size() == resultss->size());

    OutputRequestArray& requests = *outputRequests;
    for (int i = 0; i < resultss->size(); ++i) {
        for (int r = 0; r < requests.size(); ++r ){
            string requestName = requests[r]->getRequestName();
            std::map<string, OutputRequest*>::iterator pos = 
                outputRequestSet.find(requestName);
            if (pos == outputRequestSet.end() ||   /* not there - weird, but we
                                                      can react benignly */
                !pos->second->getHasFinished()) {
                //#warning "FIXME this is nasty---need getHasFinished per instrument, and pulling out individual instruments is ugly, in fact need to revisit"
                requests[r]->handleUnfulfilled(model, (*instruments)[i].get(),
                                               (*resultss)[i].get());
            }
        }
    }
}

bool Control::isPricing() const {
    return isPricingRun;
}

/** Returns the first instance of the specified sensitivity in the
    control. The comparision is an exact comparison of types
    (rather than using isAssignableFrom) */
SensitivitySP Control::sensitivityRequested(
    const CClassConstSP& clazz) const{

    //Note it is possible for sens to be NULL: consider VanillaFDProd::postPrice:
    //at the time of writing it creates a 'dummy' Control which has NULL sens.
    if (!!sens)
    {
        for (int i = 0; i < sens->size(); i++){
            if (clazz == (*sens)[i]->getClass()) {
                return (*sens)[i];
            }
        }
    }
    return SensitivitySP(   );
}
    //// Same as above but with clazz = sensToDo->getClass()
SensitivitySP Control::sensitivityRequested(SensitivitySP sensToDo) const{
    return sensitivityRequested(sensToDo->getClass());
}
 
/** Resets the control for a new pricing - clears current sensitivity and
    marks all output requests as not calculated */
void Control::reset(){
    isPricingRun = true;
    currentSens.reset();
    for (int i = 0; i < outputRequests->size(); i++){
        (*outputRequests)[i]->setHasFinished(false);
    }
    tempPackets.clear();
}
 
/** get the sensitivity which is currently being calculated */
SensitivitySP Control::getCurrentSensitivity() const {
    return currentSens;
}

/** populate from market cache */
void Control::getMarket(IModelConstSP model, MarketDataConstSP market,
                        IInstrumentCollectionConstSP instruments) {
    try {
        if (!riskQuantityEvaluator) {
            riskQuantityEvaluator.reset(new RiskQuantityEvaluator(
                model->wantsRiskMapping() == IModel::riskMappingAllowed));
        }

        riskQuantityEvaluator->getMarket(model, market, instruments);

        for (int i = 0; i < sens->size(); i++) {
            if (!!(*sens)[i]) {
                (*sens)[i]->getMarket(model.get(), market.get());
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, "Control::getMarket()");
    }
}

/** get delta shift size of a control.
    returns 0 if delta tweak not requested */
double Control::getDeltaShiftSize() const {
    return scalarShiftSize(Delta::TYPE);
}

double Control::scalarShiftSize(CClassConstSP shiftType) const {
    double shiftSize = 0.0;
    for (int i = 0 ; i < sens->size() ; ++i ) {
        if (shiftType->isInstance((*sens)[i].get())) {
            ScalarShift* scalar = dynamic_cast<ScalarShift*>((*sens)[i].get());
            shiftSize = scalar->getShiftSize();
        }
    }
    return shiftSize;
}

/** add a new sensitivity */
void Control::addSensitivity(SensitivitySP newSens) {
    try {
        sens->push_back(newSens);
    }
    catch (exception& e) {
        throw ModelException(e, "Control::addSensitivity");
    }
}

void Control::removeSensitivity(CClassConstSP sensType) { 

    SensitivityArray &arr = *sens;
    bool done = false;
    while (!done) { 
        int i = 0;
        for (i = 0;i<sens->size();i++) { 
            if (arr[i]->getClass() == sensType) { 
                 arr.erase(arr.begin() + i);
                 break;
            }
        }
        if (i == sens->size()) { 
            done = true;
        }
    }
}

// add a new request */
void Control::addRequest(OutputRequestSP request) {
    try {
        outputRequests->push_back(request);
        string requestName = request->getRequestName();
        outputRequestSet[requestName] = request.get();
    }
    catch (exception& e) {
        throw ModelException(e, "Control::addRequest");
    }
}

//// remove a request */
void Control::removeRequest(const string& requestName){
    for (vector<OutputRequestSP>::iterator iter = outputRequests->begin(); 
         iter != outputRequests->end(); /* inc in loop body */){
        if ((*iter)->getRequestName() == requestName){
            iter = outputRequests->erase(iter);
            outputRequestSet.erase(requestName);
        } else {
            ++iter;
        }
    }
}

/** Processes any command line options which alter the control. 
    The returned ControlSP is either this (and suitably modified) or is 
    a new instance of a Control. */
CControlSP Control::applyCommandLineOptions(){
    // see if anyone wants to add in a greek using its default constructor
    // or go mega
    try {
        CControlSP control(CControlSP::attachToRef(this));
        if (CommandLineParams::hasParameter(CommandLineParams::Mega)) {
            // default everything possible
            control.reset(SensitivityFactory::megaControl());
        }
        if (CommandLineParams::hasParameter(CommandLineParams::AddSens)) {
            IObjectSP param = 
                CommandLineParams::getParameterArgs(CommandLineParams::AddSens);
            const string& s = (dynamic_cast<CString&>(*param)).stringValue();
            SensitivitySP sens(SensitivityFactory::defaultSensitivity(s));
            control->addSensitivity(sens);
        }
        if (CommandLineParams::hasParameter(CommandLineParams::AddReq)) {
            IObjectSP param = 
                CommandLineParams::getParameterArgs(CommandLineParams::AddReq); 
            const string& s = (dynamic_cast<CString&>(*param)).stringValue(); 
            OutputRequestSP req = OutputRequestSP(new OutputRequest(s));
            control->addRequest(req);
        }
        return control;
    }
    catch (exception& e) {
        throw ModelException(e, "Control::applyCommandLineOptions");
    }
}

static SensitivitySP makeSensitivity(const string& name){
    return SensitivitySP(SensitivityFactory::defaultSensitivity(name));
}

static OutputRequestSP makeOutputRequest(const string& name){
    return OutputRequestSP(new OutputRequest(name));
}

/** build up a control from single character flags e.g. D = Delta */
Control* Control::makeFromFlags(const string& flags, double impVolTarget) {
    static string method = "Control::makeFromFlags";
    try {
        SensitivityArraySP   sens(new SensitivityArray(0));
        OutputRequestArraySP req(new OutputRequestArray(0));

        for (int i = 0; i < (int)flags.length(); i++) {
            char iflag = toupper(flags[i]);
            if (iflag == 'D') {
                sens->push_back(makeSensitivity(Delta::NAME));
            } else if (iflag == 'V') {
                sens->push_back(makeSensitivity(VegaParallel::NAME));
            } else if (iflag == 'R') {
                sens->push_back(makeSensitivity(RhoParallel::NAME));
            } else if (iflag == 'B') {
                sens->push_back(makeSensitivity(RhoBorrowParallel::NAME));
            } else if (iflag == 'T') {
                sens->push_back(makeSensitivity(Theta::NAME));
            } else if (iflag == 'S') {
                sens->push_back(makeSensitivity(VegaSkewParallel::NAME));
            } else if (iflag == 'M') {
                sens->push_back(makeSensitivity(MuParallel::NAME));
            } else if (iflag == 'W') {
                sens->push_back(makeSensitivity(VegaParallel2Sided::NAME));
            } else if (iflag == 'E') {
                sens->push_back(makeSensitivity(DDeltaDVol::NAME));
            } else if (iflag == 'X') {
                sens->push_back(makeSensitivity(RootTimeVega::NAME));
            } else if (iflag == 'F') {
                req->push_back(makeOutputRequest(OutputRequest::FWD_AT_MAT));
            } else if (iflag == 'K') {
                req->push_back(makeOutputRequest(OutputRequest::IND_VOL));
            } else if (iflag == 'I') {
                // make a guess - say 16%
                sens->push_back(SensitivitySP(
                                    new ImpliedVol(impVolTarget,0.16,
                                                   0.0001,Results::VALUE)));
            }
        }
        return new Control(sens, req, false, "");
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Returns whether outputs should be scaled or not */
bool Control::scaleOutputs() const
{
    return scaleResults;
}

/** Sets whether outputs should be scaled or not */
void Control::setScaleOutputs(bool scale)
{
    scaleResults = scale;
}


/** Sets whether assets should use theoretical price or mtm */
void Control::setUseTheoAssetPrice(bool useTheoAssetPrice)
{
    useTheoAssetPrice = useTheoAssetPrice;
}

/** Indicates whether assets should use theoretical price or mtm for
    the current pricing */
bool Control::getUseTheoAssetPrice() const{
      return useTheoAssetPrice;
}  

/** Returns the source we should use of asset prices as requested by the
    client. Note the string has multiple values including 'both'. Its
    value should not be used to indicate the rule for the current
    pricing. Instead use getUseTheoAssetPrice */
const string& Control::getAssetPriceSource() const{
    return assetPriceSource;
}

/** Requests that the specified packet is removed from the results at the
    end of the calculate() call */
void Control::removePacketAfterCalc(const string& packetName){
    /* motivation: allow CrossGamma/FXCrossGamma to remove delta packets
       after calculation has finished (if delta not requested). Their
       presence mucks up the output scaling/aggregation. Also we should
       only remove at the end since the delta values are used in both
       CrossGamma and FXCrossGamma */
    tempPackets[currentInstrumentIndexInCollection]->push_back(packetName);
}

/** Records the price time in the results (if requested) assuming that
    the pricing call has now finished. Useful for models that
    calculate their own greeks in the initial pricing call. Returns
    the stored value of priceTime (in seconds) */
double Control::recordPriceTime(CResultsArraySP resultss){
    // store timings if required
    clock_t endPrice = clock();
    bool    storeDbg = false;

    if (resultss->size() == 0)
        return (double)(endPrice-startTime)/CLOCKS_PER_SEC;

    double priceTime;
    CResultsSP results = (*resultss)[0];

    OutputRequest priceTimeRequest(OutputRequest::PRICE_TIME);
    OutputNameSP priceTimeName(new OutputName(OutputRequest::PRICE_TIME));
    // don't overwrite if there (model might do greeks in pricing call)
    // Also, if N/A (bogus value) then calculate
    if ((!results->exists(Results::DEBUG_PACKET, priceTimeName) &&
         !results->exists(&priceTimeRequest)) ||
        results->isNotApplicable(&priceTimeRequest)){
        priceTime = (double)(endPrice-startTime)/CLOCKS_PER_SEC;
        OutputRequest* request = requestsOutput(OutputRequest::PRICE_TIME);
        if (request){
            results->storeRequestResult(request, priceTime);
        }
        storeDbg = true;
    } else {
        if (results->exists(&priceTimeRequest)){
            priceTime = results->retrieveScalarGreek(
                priceTimeRequest.getPacketName(),
                priceTimeName);
            storeDbg = true;
        } else {
            priceTime = results->retrieveScalarGreek(Results::DEBUG_PACKET,
                                                     priceTimeName);
        }
    }
    if (storeDbg){
        // store it anyway in the debug section
        results->storeScalarGreek(priceTime, Results::DEBUG_PACKET,
                                  priceTimeName);
    }
    return priceTime;
}

/** If skipShifts() = true then when a shift to the market is requested
    it can be skipped. This is for a performance gain when assessing
    the total number of calculations an instrument requires */
bool Control::skipShifts() const{
    return skipAllShifts;
}

class ControlHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spread sheet
        REGISTER(Control, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultControl);
        FIELD(sens, "List of sensitivities to calculate");
        FIELD(outputRequests, "List of additional output requests");
        FIELD_MAKE_OPTIONAL(outputRequests);
        FIELD(writeToFile, "Whether to create a regression file");
        FIELD(fileName, "Name for regression file");
        FIELD_MAKE_OPTIONAL(fileName);
        FIELD(scaleResults, "Whether to scale the results");
        FIELD_MAKE_OPTIONAL(scaleResults);
        FIELD(assetPriceSource, "ASSET_THEO (default), ASSET_MTM or "
                     "ASSET_THEO_MTM");
        FIELD_MAKE_OPTIONAL(assetPriceSource);
        FIELD(isPricingRun, "doing price or greeks");
        FIELD_MAKE_TRANSIENT(isPricingRun); // hide from dd interface
        FIELD(currentSens, "sens currently being calculated");
        FIELD_MAKE_TRANSIENT(currentSens); // hide from dd interface
        FIELD(useTheoAssetPrice, "if true use asset theo price"); 
        FIELD_MAKE_TRANSIENT(useTheoAssetPrice); 
        FIELD(currentInstrumentIndexInCollection, "FIXME");
        FIELD_MAKE_TRANSIENT(currentInstrumentIndexInCollection); 
        FIELD(tempPackets, "Packets to remove after pricing");
        FIELD_MAKE_TRANSIENT(tempPackets); 
        FIELD(skipAllShifts, "Skip shifts to market data");
        FIELD_MAKE_TRANSIENT(skipAllShifts); 
        FIELD(riskQuantityEvaluator,
              "Evaluator for 'new-style' greeks, possibly with risk mapping");
        FIELD_MAKE_TRANSIENT(riskQuantityEvaluator);
    }

    static IObject* defaultControl(){
        return new Control();
    }
};

CClassConstSP const Control::TYPE = CClass::registerClassLoadMethod(
    "Control", typeid(Control), ControlHelper::load);

DEFINE_TEMPLATE_TYPE_WITH_NAME("ControlArray", CControlArray);

// for reflection
Control::Control(): CObject(TYPE), 
                    assetPriceSource(ASSET_THEO),
                    currentInstrumentIndexInCollection(0), 
                    useTheoAssetPrice(true), isPricingRun(true), 
                    scaleResults(false), writeToFile(false),
                    skipAllShifts(false) {}


DRLIB_END_NAMESPACE
