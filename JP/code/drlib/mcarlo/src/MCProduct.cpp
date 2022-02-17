#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/HandlePaymentEvents.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/RhoParallel.hpp"
#include "edginc/RhoPointwise.hpp"
#include "edginc/CliquetVolRequest.hpp"
#include "edginc/MCProduct.hpp"
#include "edginc/MCProductIQuicks.hpp"
#include "edginc/MCProductEngineClient.hpp"
#include "edginc/MCPrices.hpp"
#include "edginc/CashSettlePeriod.hpp"

DRLIB_BEGIN_NAMESPACE
// cachedQckGreek is a bitwise variable
const int IMCProduct::QUICK_GREEKS_BIT = 1;
const int IMCProduct::QUICK_X_GAMMA_BIT = 2;
// KKK this is not strictly correct - product caching is NOT a form of quick greeks.
// KKK but it is convenient and doesn't widen the interface when will ultimately probably 
// KKK have just a "cache on/off" driven through "cachePaths".
const int IMCProduct::CACHE_PRODUCT_BIT = 4; 

IMCProduct::~IMCProduct(){}

/** Returns null if not supported. Implemented via dynamic cast. Allows
    per instrument choice */
IMCQuickGreeks* IMCProduct::quickGreeksSupported(){
    return (dynamic_cast<IMCQuickGreeks*>(this));
}

/** Returns null if not supported. Implemented via dynamic cast. Allows
    per instrument choice */
IMCQuickXGamma* IMCProduct::quickXGammaSupported(){
    return (dynamic_cast<IMCQuickXGamma*>(this));
}

/** This method is called every time the path generator is changed
    (which is, at the moment, when the past path generator is
    created, and then when the future path generator is
    created. Default does nothing  */
void IMCProduct::pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){}

void IMCProduct::validate(){
    const static string routine("IMCProduct::validate");
    // some idiot checking
    if (!refLevel){
        throw ModelException(routine, "NULL reference level");
    }
    if (!simSeries){
        throw ModelException(routine, "NULL sim series");
    }
    if (!mcPastValues){
        throw ModelException(routine, "NULL PastValues");
    }
    /* validate numAssets is consistent */
    if (NbAssets != simSeries->getNumAssets()){
        throw ModelException(routine, "MultiFactorsAsset has "+
                             Format::toString(NbAssets)+
                             " but Sim Series has "+
                             Format::toString(simSeries->getNumAssets()));
    }
    if (NbAssets != mcPastValues->getNumAssets()){
        throw ModelException(routine, "MultiFactorsAsset has "+
                             Format::toString(NbAssets)+" assets"
                             " but Past Values has data for"+
                             Format::toString(mcPastValues->getNumAssets())+
                             " assets");
    }
}        

// runs the past without storing results at all - just wraps a call to handlePast()
// may need to keep hold of path generator until finished with to avoid reading freed memeory
MCPathGeneratorSP IMCProduct::runPast(MonteCarlo* mcarlo) {
    IMCPricesSP prices(createOrigPrices(1, 1, 0));
    CControlSP ctrl(new Control(SensitivityArrayConstSP(   ), 
                                OutputRequestArrayConstSP(   ),false,""));
    return handlePast(mcarlo, ctrl.get(), prices, false);
}

MCPathGeneratorSP IMCProduct::handlePast(MonteCarlo*  mcarlo,
                                        Control*     control,
                                        IMCPricesSP&    prices,
                                        bool         needPrices) {
    const static string routine("IMCProduct::handlePast");
        // get path generator for dealing with historic dates
    MCPathGeneratorSP pastPathGenerator(
        mcarlo->getPathConfig()->pastPathGenerator(this));
    try {
        // Handle past early 
        // Tell the product that the generator has changed
        // Do that even if there is no past so that the product gets
        // some state variables e.g. refLevel
        pathGenUpdated(pastPathGenerator.get());

        if (pastPathGenerator->hasPast()) {
            // create our IMCPrices object
            // get hold of the Pricing object
            if (needPrices) {
                MCPricing* pricing = mcarlo->getPricing().get();
                if (!pricing){ // internal error
                    throw ModelException(routine, "Null pricing object");
                }
                prices = IMCPricesSP(pricing->createPastPrices(mcarlo, this));
            }
            pastPathGenerator->generatePath(0); // may well do nothing
            if (mcarlo->inStatelessMode()) {
                MCPastPricing pastPricer;
                pastPricer.runSim(mcarlo, this, control, pastPathGenerator.get(), *prices);
            } else { // use old method 
                payoff(pastPathGenerator.get(), *prices);
            }
        }        
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
    return pastPathGenerator;
}

/** IMCPrices the product. Returns the future path generator or null if there
    was no 'future' to price */
MCPathGeneratorSP IMCProduct::price(MonteCarlo*       mcarlo,
                                   Control*          control, 
                                   Results*          results){
    const static string routine("IMCProduct::price");
    MCPathGeneratorSP futurePathGenerator;
    try{
        bool isPricing = control->isPricing();
		OutputRequest* request;
        // get hold of the Pricing object
        MCPricing* pricing = mcarlo->getPricing().get();
        if (!pricing){ // internal error
            throw ModelException(routine, "Null pricing object");
        }
        
        // deal with historic dates
        IMCPricesSP prices; // note we pass this in empty and it might get filled
        MCPathGeneratorSP pastPathGenerator= handlePast(mcarlo, control, 
                                                         prices, true);
        
        if (hasFuture()){ // sim needed
            // give payoff opportunity to establish future timeline before
            // the PathGen is configured for the simulation
            finaliseSimulationDates();
            // should we request the paths to be cached?
            bool cachePaths = isPricing && mcarlo->getCachedCachePaths() &&
                !control->getSens()->empty();
            /* Some Pricing classes skip paths so the path generator needs
               to provide 'random access' to the paths */
            int cacheMode = 
                (pricing->needRandomAccessToPaths(mcarlo, this, control)?
                 IMCPathConfig::PATH_RANDOM_ACCESS: 0) |
                (cachePaths? IMCPathConfig::PATH_CACHE: 0);
            // now get path generator for future dates
            futurePathGenerator = mcarlo->getPathConfig()->futurePathGenerator(
                cacheMode, mcarlo->getNbIters(true),
                pastPathGenerator, this,
                control, results, pricing->getTimeLine() );
            // switch doingThePast flag
            doingThePast = false;
            MCPathGenerator* pathGen = futurePathGenerator.get(); // for ease
            // tell the product that the generator has changed
            pathGenUpdated(pathGen);

            // where we store the results
            prices = pricing->createFuturePrices(mcarlo, this, 
                                                 control, futurePathGenerator);
            
            // If relevant allow configuration of product cache during tweaks
            if (!isPricing && mcarlo->getCachedCachePaths()) {
                const IMultiFactors* multiFactor = getMultiFactors();
                // then which greek we're doing
                SensitivitySP sens(control->getCurrentSensitivity());
                // getSensitiveAssets() does the real work
                IntArray sensitivePhiAssets(multiFactor->getSensitiveAssets(
                    sens.get(), true)); // include phi etc
                if (sensitivePhiAssets.empty()){
                    IntArray sensitiveAssets(multiFactor->getSensitiveAssets(
                        sens.get(), false)); // exclude phi etc
                    prices->configureCache(sensitiveAssets);
                } else {
                    // mark all path caches as invalid
                    // KKK Would like a nice way to do this...
                    // KKK And preferably inside multiFactor->getSensitiveAssets
                    IntArray allAssets(multiFactor->NbAssets());
                    for (int i = 0; i < allAssets.size(); i++){
                        allAssets[i] = i;
                    }
                    prices->configureCache(allAssets);
                }
            }

            // run loop
            pricing->runSim(mcarlo, this, control, pathGen, *prices);
        }
        // Pricing class is responsible for storing Price
        pricing->postResults(mcarlo, this, control, *prices,
                             discount->getCcy(), results);
        // Allow product to store anything extra
        recordExtraOutput(control, results, *prices);
        if (isPricing){
            // save state of pathConfig if splitting into blocks
            if (mcarlo->getCachedCachePaths()){
                mcarlo->initializePathConfigPostPricing();
            }
            mcarlo->setOrigPrices( prices ); // save for future use
            // do any events that we know how to do
            recordEvents(control, *prices, results);
        }
        if	((isPricing) &&
                 (request = control->requestsOutput(OutputRequest::FWD_AT_MAT)) &&
                 (!request->getHasFinished())) {
            const DateTime& matDate = simSeries->getLastDate();
            if (matDate.isGreaterOrEqual(Today)){
                /* do fwd @ mat by stock */
                 if (mfAsset.get()) {
                      mfAsset->recordFwdAtMat(request, results, matDate);
                 }
            }
        }
    } catch (IPastValues::MissingSampleException& mse){
        mse.addMsg(routine);
        if (mfAsset.get()) {
             mse.setAssetName(mfAsset->getName(mse.getAssetIndex()));
        }
        throw;
    } catch (ModelException& e){
        // check to see if MissingSampleException wrapped in a ModelException
        e.addMsg(routine);
        IPastValues::MissingSampleException* mse = 
            IPastValues::MissingSampleException::getInstance(e);
        if (mse && mfAsset.get()){
            mse->setAssetName(mfAsset->getName(mse->getAssetIndex()));
        }
        throw;
    } catch (exception& e){
        throw ModelException(e, routine);
    }
    return futurePathGenerator;
}

// default implementation to record events (pay dates, cashflows etc)
// will try and use any product's implementation first
void IMCProduct::recordEvents(Control* control,
                             IMCPrices&  prices,
                             Results* results) {
    static const string method("IMCProduct::recordEvents");
    try {
        // see if product handles this first
        IHandlePaymentEvents* handler =
            dynamic_cast<IHandlePaymentEvents*>(this);
  
        if (handler) {
            // this method/interface may need work
            handler->recordEvents(control, results);
        }
        else {
            // built-in implementation for single cashflow payoffs
            OutputRequest* request =
                control->requestsOutput(OutputRequest::PAYMENT_DATES);
            if (request && !request->getHasFinished()) {
                if (paymentDate.empty()) {
                    throw ModelException(method, 
                                         "payment date has not been set");
                }

                DateTimeArray date(1, paymentDate);
                OutputRequestUtil::recordPaymentDates(control,results,&date);
            }

            request =control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
            if (request && !request->getHasFinished() &&
                !simSeries->getLastDate().isGreater(Today) &&
                !settlement->isPhysical()) {
                if (paymentDate.empty()) {
                    throw ModelException(method, 
                                         "payment date has not been set");
                }
                // we know the cashflow
                // This relies upon the happy chance that the prices class
                // holds the non-PVd value.
                double amount;
                double stdErr;
                IMCPricesSimple& simplePrices = dynamic_cast<IMCPricesSimple&>(prices); // FIXME
                simplePrices.getResult(&amount, &stdErr);
            
                CashFlow cf(paymentDate, amount);
                CashFlowArray cfl(1, cf);
                OutputRequestUtil::recordKnownCashflows(control,
                                                        results,
                                                        discount->getCcy(),
                                                        &cfl); 
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** the date to stop tweaking for the supplied tweak. The
    default implementation uses the last simulation date and the 
    payment date. Derived classes overriding pvFromPaymentDate() need to
    override this method too */
DateTime IMCProduct::endDate(const Sensitivity*  sensControl) const{
    const DateTimeArray& simDates = getSimSeries()->getAllDates();
    const DateTimeArray& refSimDates = getRefLevel()->getAllDates();
    DateTimeArray allDates(DateTime::merge(simDates, refSimDates));
    DateTime lastDate = allDates.empty()? DateTime(): allDates.back();
    // short term fix for rho - need better way of determining of deciding
    // value of isRhoSens
    // also to do: take account of stock settlement
    bool isRhoSens = RhoParallel::TYPE->isInstance(sensControl) ||
        RhoPointwise::TYPE->isInstance(sensControl);

    if (isRhoSens && paymentDate.empty()) {
        throw ModelException("IMCProduct::endDate",
                             "payment date has not been set");
    }

    if (isRhoSens && paymentDate.isGreater(lastDate)){
        lastDate = paymentDate;
    }
    return lastDate;
}

/** 'Full' constructor for single factor payoffs */
IMCProduct::IMCProduct(
    const CAsset*               asset,
    const DateTime&             today,
    const YieldCurve*           discount,
    const IRefLevelConstSP&     refLevel,     // how to 'avg in'
    const SimSeriesConstSP&     simSeries,    // sim dates
    const IPastValuesConstSP&   mcPastValues, 
    const InstrumentSettlement* instSettle,   // used for pay date
    const DateTime&             matDate):     // used for pay date
    doingThePast(true), Today(today), discount(discount), 
    mAsset(convertToMultiFactors(asset)), mfAsset(mAsset),
    refLevel(refLevel), simSeries(simSeries), 
    paymentDate(instSettle->settles(matDate, asset)) {
    
    settlement = InstrumentSettlementConstSP(copyIfRef(instSettle));
    // for public use
    NbAssets = mAsset->NbAssets();

    // This does some checking too
    setMCPastValues(mcPastValues);
}

/** 'Full' constructor for multi-factor payoffs */
IMCProduct::IMCProduct(
    const IMultiMarketFactors*  mFactors,
    const DateTime&             today,
    const YieldCurve*           discount,
    const IRefLevelConstSP&     refLevel,     // how to 'avg in'
    const SimSeriesConstSP&     simSeries,    // sim dates
    const IPastValuesConstSP&   mcPastValues, 
    const InstrumentSettlement* instSettle,   // used for pay date
    const DateTime&             matDate):     // used for pay date
    doingThePast(true), Today(today), discount(discount), 
    refLevel(refLevel), simSeries(simSeries), 
    paymentDate(instSettle->settles(matDate, NULL /* asset */)){

    settlement = InstrumentSettlementConstSP(copyIfRef(instSettle));
    mfAsset = IMultiMarketFactorsConstSP::attachToRef(mFactors);

    // for public use
    NbAssets = mfAsset->numFactors(); // terminology is looking a bit skewed

    // This does some checking too
    setMCPastValues(mcPastValues);
}

/** Simpathl constructor */
IMCProduct::IMCProduct(
    const IMultiMarketFactors*  mFactors,
    const DateTime&             today,
    const YieldCurve*           discount):
    Today(today), discount(discount)
{
    mfAsset = IMultiMarketFactorsConstSP::attachToRef(mFactors);
    NbAssets = mfAsset->numFactors(); // terminology is looking a bit skewed
    IPastValuesConstSP pastValues( 
        IPastValues::Util::makeTrivial( 
            today, 
            DoubleArray( mFactors->NbAssets(), 0 ) ) );
    setMCPastValues( pastValues );
}

/** simSeries remains since I cannot see a way to internally build one from SVs
    Ideally simSeries not needed at all, but that is a bigger task */
IMCProduct::IMCProduct(
     const DateTime&             today,
     const YieldCurve*           discount,
     const SimSeriesConstSP&     simSeries): 
     doingThePast(true), Today(today), 
     NbAssets(simSeries->getNumAssets()),
     discount(discount), 
     settlement(new CashSettlePeriod(0)),
     mAsset(0), mfAsset(0),
     refLevel(IRefLevel::Util::makeTrivialAverage(today)),
     simSeries(simSeries),
     paymentDate(simSeries->getAllDates().back())
{
     IPastValuesConstSP pastValues( 
          IPastValues::Util::makeTrivial( 
               today, 
               DoubleArray( NbAssets, 0 ) ) );
     setMCPastValues( pastValues );
}

//// turns Asset into MultiFactors
IMultiFactorsConstSP IMCProduct::convertToMultiFactors(const Asset* asset){
    IMultiFactorsConstSP mAsset(asset->asMultiFactors());
    return mAsset;
}

double IMCProduct::pvFromPaymentDate() const {
    if (paymentDate.empty()) {
        throw ModelException("IMCProduct::pvFromPaymentDate",
                             "payment date has not been set");
    }
    if (paymentDate.isGreater(Today)){
        // Get pv from the settlement
        
        // slight hack for now - as we don't have the matDate
        if (settlement->isMargin()){
            return 1.0;
        }
        return discount->pv(Today, paymentDate);
    }
    return 0.0; // cashflow is in the past
}

//// sets the payment date to given value
void IMCProduct::setPaymentDate(const DateTime& payDate){
    paymentDate = payDate;
}

/** sets the payment date calculated from supplied InstrumentSettlementSP
        and maturity date */
void IMCProduct::setPaymentDate(const InstrumentSettlement* instSettle,
                               const DateTime&             matDate){
    // Meaningful only for cash settlement - physical settlement
    // will barf here since no asset passed, but that's ok
    paymentDate = instSettle->settles(matDate, NULL);
}   
    
/** sets the RefLevel obejct (does not take copy) */
void IMCProduct::setRefLevel(const IRefLevelConstSP& refLevel){
    this->refLevel = refLevel;
}

/** set the SimSeries object (does not take copy) */
void IMCProduct::setSimSeries(const SimSeriesConstSP& simSeries){
    this->simSeries = simSeries;
}

//// Returns the multiFactors
const IMultiFactors* IMCProduct::getMultiFactors() const {
    static const string method("IMCProduct::getMultiFactors");
    if (!mAsset){
        if (!mfAsset.get()) {
             throw ModelException(method,
                                  "Internal error - unexpected call to this method!");
        }
        if (IMultiFactors::TYPE->isInstance(mfAsset)){
            // can just cast up to it
            mAsset = IMultiFactorsConstSP::dynamicCast(mfAsset);
        } else {
            const IAsMultiFactors* asMulti = 
	      dynamic_cast<const IAsMultiFactors*>(mfAsset.get());

            if (!asMulti){
                throw ModelException(method,
                                     "IMultiFactors view not supported by"
                                     " object of type "+
                                     mfAsset->getClass()->getName());
            }
            try {
                mAsset = IMultiFactorsConstSP(asMulti->asMultiFactors());
            } catch (exception& e){
                throw ModelException(e, method);
            }
        }
    }
    return mAsset.get();
}


//// Returns the MultiMarketFactors object
const IMultiMarketFactors* IMCProduct::getMultiMarketFactors() const{
    return mfAsset.get();
}

//// Returns the num of assets (in the payoff view of the world)
int IMCProduct::getNumAssets() const{
    return NbAssets;
}

//// Returns the value date
const DateTime& IMCProduct::getToday() const{
    return Today;
}

//// Returns the SimSeries
const SimSeries* IMCProduct::getSimSeries() const{
    return simSeries.get();
}

//// Returns the RefLevel
const IRefLevel* IMCProduct::getRefLevel() const{
    return refLevel.get();
}

    //// Returns the discount YieldCurve
const YieldCurve* IMCProduct::getDiscount() const{
    return discount;
}


//// Returns the PastValues
const IPastValues* IMCProduct::getMCPastValues() const{
    return mcPastValues.get();
}

/** set the PastValues object (does not take copy) */
void IMCProduct::setMCPastValues(const IPastValuesConstSP& mcPastValues){
    this->mcPastValues = mcPastValues;
}
 
/** Returns an array of 'forward start dates' (ie cliquet start dates)
    which should be used for volatilty interpolation. The default method
    here uses the IMCProductLN interface if the product implements it
    otherwise it fails (so this method would need to be overridden) */
DateTimeArray IMCProduct::getVolStartDates(
    const MCPathGenerator* pathGenerator) const
{
    static const string method("IMCProduct::getVolStartDates");
    const IMCProductLN* mcProductLN = dynamic_cast<const IMCProductLN*>(this);
    if (!mcProductLN){
        throw ModelException(method, "Can't use IMCProductLN to get vol "
                             "start dates. This method needs to be overridden");
    }
    DateTimeArray theDates;
    for (int i = 0; i < NbAssets; i++){
        CVolRequestLNArray requests(mcProductLN->getVolInterp(pathGenerator,
                                                              i));
        if (requests.empty()){
            throw ModelException(method, "No vol requests specified");
        }
        if (CliquetVolRequest::TYPE->isInstance(requests.front())){
            CliquetVolRequestSP cliq(
                CliquetVolRequestSP::dynamicCast(requests.front()));
            DateTimeArray tmpDates(cliq->getCliqStartDates());
            if (i == 0){
                theDates = tmpDates;
            } else if (!DateTime::equals(theDates, tmpDates)){
                throw ModelException(method, "Different cliquet dates per asset"
                                     "not supported");
            }
        } else {
            const DateTime& startDate = requests.front()->getStartDate();
            if (i == 0){
                theDates = DateTimeArray(1, startDate);
            } else if (!startDate.equals(theDates.front())){
                throw ModelException(method, "Different start dates per asset"
                                     "not supported");
            }
        }
    }
    return theDates;
}

//// Returns MAX(simulation start date, today)
const DateTime& IMCProduct::getEffectiveSimStartDate() const{
    return refLevel->getSimStartDate(Today);
}

/** This method is to be retired once we move to 'state variables' for the
    Monte Carlo. This is a hack to get around something that was designed in
    namely the fact that the ref level is responsible for choosing the
    sim start date */
void IMCProduct::setSimStartDateToFirstRefDate(){
    // clone it
    IRefLevelSP newRefLevel(refLevel.clone());
    // and then change it
    newRefLevel->setSimStartDateToFirstRefDate();
    // and then set it
    refLevel = newRefLevel;
}

/** Invoked once the 'past' has been completed (invoked once only
    per pricing)*/
void  IMCProduct::finaliseSimulationDates() {}

/** Default implementation */
IMCPrices* IMCProduct::createOrigPrices(int  nbIter,
                                               int  nbSubSamples,
                                               int  mode) 
{
    return new MCPricesSimple(nbIter, nbSubSamples);
}

/** invoked after final simulated path is run. Default does nothing.
    Allows derived classes to store debug information for example */
void IMCProduct::recordExtraOutput(CControl*     control,
                                  Results*      results,
                                  const IMCPrices& prices) const
{
    // empty
}

// invoked after past is run. Default does nothing.
void IMCProduct::retrieveEvents(EventResults* events) const {}

// go get the events based on this path
void IMCProduct::getEvents(const IMCPathGenerator*  pathGen, 
                          EventResults* events,
                          const DateTime& eventDate) {};


DRLIB_END_NAMESPACE
