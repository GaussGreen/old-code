
#include "edginc/config.hpp"
#define QLIB_MODEL_CPP
#include "edginc/Model.hpp"
#include "edginc/ModelFilter.hpp"
#include "edginc/RiskMgr.hpp"
#include "edginc/Scenario.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/Modifier.hpp"
#include "edginc/Addin.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/Events.hpp"

#include <algorithm>

DRLIB_BEGIN_NAMESPACE

CModel::~CModel(){}

/** Does everything, namely retrieve market data, applies scenario 
    (if non null), calculates price and greeks as well as save to file of
    inputs and outputs. */
CResultsArraySP CModel::go(IInstrumentCollectionSP instruments,
                           ScenarioSP              scenario, // optional
                           CControlSP              control,
                           MarketDataSP            market){
    IModelSP model(IModelSP::attachToRef(this));
    return go(model, instruments, scenario, control, market);
}

/** Same as above go method above but for on a single instrument. The two
    methods are separate to allow eg a different write to file. (This 
    isn't ideal but keeps the ability to create a regression file which
    actually mirrors what was actually done. Alternatives would be a method
    on IInstrumentCollection or some ugly switch on type of instruments */
CResultsSP CModel::go(CInstrumentSP instrument,
                      ScenarioSP    scenario, // optional
                      CControlSP    control,
                      MarketDataSP  market){
    IModelSP model(IModelSP::attachToRef(this));
    return go(model, instrument, scenario, control, market);
}

/** Does everything, namely retrieve market data, applies scenario (if
    non null), calculates price and greeks as well as save to file of
    inputs and outputs. Static method to allow
    other implementations of IModel to use this if they want */
CResultsArraySP CModel::go(IModelSP                model,
                           IInstrumentCollectionSP instruments,
                           ScenarioSP              scenario, // optional
                           CControlSP              control,
                           MarketDataSP            market){
    // do write to file - for backward compatibility need to switch
    // depending on whether scenario is null or not.
    bool writeToFile = control->getWriteToFile();
    if (writeToFile) {
        if (!scenario){
            CRiskMgr::writeInputs(model, instruments, control, market);
        } else {
            scenario->writeInputs(model, instruments, control, market);
        }
    }
    CResultsArraySP results(miniGo(model, instruments,
                                   scenario, control, market));
    if (writeToFile) {
        CRiskMgr::writeOutputs(control.get(), results);
    }
    return results;
}

/** Does everything, namely retrieve market data, applies scenario (if
    non null), calculates price and greeks as well as save to file of
    inputs and outputs. Static method to allow
    other implementations of IModel to use this if they want */
CResultsSP CModel::go(IModelSP      model,
                      CInstrumentSP instrument,
                      ScenarioSP    scenario, // optional
                      CControlSP    control,
                      MarketDataSP  market){
    // do write to file - for backward compatibility need to switch
    // depending on whether scenario is null or not.
    bool writeToFile = control->getWriteToFile();
    if (writeToFile) {
        // Need to make RiskMgr::writeInputs and Scenario::writeInputs public
        if (!scenario){
            CRiskMgr::writeInputs(model, instrument, control, market);
        } else {
            scenario->writeInputs(model, instrument, control, market);
        }
    }
    CResultsArraySP results(miniGo(model, 
                                   IInstrumentCollection::singleton(instrument),
                                   scenario, control, market));
    CResultsSP result(results->front());
    if (writeToFile) {
        CRiskMgr::writeOutputs(control.get(), result);
    }
    return result;
}

/** Same as go except for save to file of inputs and outputs. */
CResultsArraySP CModel::miniGo(
    IModelSP                model,
    IInstrumentCollectionSP instruments,
    ScenarioSP              scenario, // optional
    CControlSP              control,
    MarketDataSP            market){
    // Get market data:
    // work on a copy of the instrument since we're possibly amending data
    // when doing the 0 day theta shift
    IInstrumentCollectionSP insts(instruments.clone());
    // copy model as well since we might be filling it with market data
    IModelSP      mdl(model.clone());
    // copy control in case of state info inside sensitivities
    CControlSP    ctrl(control.clone());

    // For scenario shifts that need to be applied before market data is fetched.
    // Apply scenario if non null.
    if (scenario.get()){
        scenario->preapply(mdl, insts);
    }


    // call model to collect instrument and model market data
    mdl->getInstrumentsAndModelMarket(market, insts);

    // call getMarket to initiate market data selection.
    ctrl->getMarket(mdl, market, insts);

    // apply scenario if non null
    if (scenario.get()){
        scenario->apply(mdl, insts);
    }
    return RiskMgr::calculateTheoAndMtm(mdl, insts, ctrl);
}

/** wrapper round RunMulti */
CResults* CModel::Run(CInstrument* instrument, 
                      CControl*    control) {
    return (*RunMulti(IInstrumentCollection::singleton(instrument),
                      control))[0].release();
}

/** main control - calculates prices and sensitivities. This method is
    implemented and, in general, should not be overridden */
CResultsArraySP CModel::RunMulti(IInstrumentCollectionSP instruments,
                                 CControl* control) {
    CResultsArraySP resultss = instruments->emptyResults();
    instruments->Validate();
    control->calculateMulti(this, instruments, resultss);
    return resultss;
}

void CModel::PriceMulti(IInstrumentCollectionSP instruments, 
                        CControl* control, 
                        CResultsArraySP resultss) {
    ASSERT(instruments->size() == resultss->size());
    for (int i = 0; i < instruments->size(); ++i){
        Price((*instruments)[i].get(), control, (*resultss)[i].get());
    }
}

/** utility wrapper around Price */
double CModel::calcPrice(CInstrument*  instrument, 
                         CControl*     control){
    Results results;
    Price(instrument, control, &results);
    return results.retrievePrice();
}

/** Essentially a wrapper around MarketData::GetData(model, name, type).
    Retrieves the market object with
    the specififed name and type (and retrieves its market data). The
    domesticYCName flags what the current 'domestic' currency is (see
    getDomesticYCName() for more information). */
MarketObjectSP CModel::getMarketObject(const MarketData* market,
                                       const string&     name,
                                       CClassConstSP     type,
                                       const string&     domesticYCName) const{
    /* Check if there is a mdf in place, or create one otherwise */
    checkMDF();

    return mdf->getMarketObject(this, market, name, type, domesticYCName);
}


/** Get the market data with supplied name and type from the given market
    data cache. This gives the model a chance to choose a specific type
    of market data rather than just a general instance. For example,
    the method could request a Black-Scholes Vol rather than just any
    old vol. The default implementation provided by CModel just asks
    the market data for the object of the given type */
MarketObjectSP CModel::GetMarket(const MarketData*    market,
                                 const string&        name,
                                 const CClassConstSP& type) const {
    // make sure valueDate is good - TracerModel contains sub-models that do not populate valueDate
    if( valueDate.empty() )
        const_cast<CModel*>(this)->valueDate = market->GetReferenceDate();

    /* Check if there is a mdf in place, or create one otherwise */
    checkMDF();

    return mdf->fetch(market, name, type, this);
}

/** Method to set the domestic yield curve name in the mdf  */
void CModel::setDomesticYCName(string discountCurveName) const {
    // Check if there is a mdf in place, or create one otherwise
    checkMDF();

    /* Retrieve the discount yield curve name from the instrument and
     * set it on the MDF */
    mdf->setDomesticYCName(discountCurveName);
}

/** Populates valueData field in this object from MarketData - this is
    is done in getInstrumentsAndModelMarket but if you override that method
    this method needs to be called */
void CModel::fetchToday(MarketDataConstSP market){
    valueDate = market->GetReferenceDate();
}

void CModel::getInstrumentsAndModelMarket(MarketDataConstSP market,
                                          IInstrumentCollectionSP insts) {
    try {
        // set valueDate first in case called by any others
        fetchToday(market);

        /* Retrieve the discount yield curve name from the instrument and 
         * set it */
        setDomesticYCName(insts->discountYieldCurveName());

        // Call instrument GetMarket to initiate market data selection.
        insts->GetMarket(this,
                         CMarketDataSP(const_cast<MarketData *>(market.get())));

        // try to get the valueDate if it's still empty, for example
        // in the ported EDG tests
        if (valueDate.empty()){
            fetchToday(market);
        } else {
            market->GetReferenceDate(valueDate);
        }
        // Call model (us) to see if any extra data is required
        this->getMarket(market.get(), insts);
    }
    catch (exception &e) {
        throw ModelException(e,
                             "CModel::getInstrumentAndModelMarket",
                             "failed getting market data");
    }
}

void CModel::getInstrumentAndModelMarket(const MarketData*  market,
                                         CInstrument*       inst) {
    getInstrumentsAndModelMarket(MarketDataConstSP::attachToRef(market),
                                 IInstrumentCollection::singleton(inst));
}

/** Method to check if there is a mdf in place in the current model. Otherwise 
* it creates one */
void CModel::checkMDF() const {
    if (!mdf) {
        mdf.reset(createMDF().get());
    }
}

/** override a control shift (eg for delta on trees)
    returns true if new control is constructed else returns 0 */
SensControl* CModel::AlterControl(const SensControl* currSensControl) const { 
    // Default implementation returns NULL
    return 0;
}

/** called after a series of pricings to indicate that the object
    will not be used again to calculate any
    sensitivities. Basically, state information can be stored
    inside the Model - some of which might be expensive in terms
    of memory. Since the model is constructed by clients, we do
    not control when the Model object is freed. This gives a
    chance for Models to free any expensive caches. The default
    implementation does nothing */
void CModel::flush(){}

/** Having retrieved a market object from the cache (or if presupplied) it
    is necessary to ensure that it has all its market data. Rather than
    call getMarket on the MarketObject, this method should be called 
    instead as it allows the model to have context information 
    regarding any calls to GetMarket. The default implementation,
    of course, just calls getMarket on the MarketObject. */
void CModel::getComponentMarketData(const MarketData*    market,
                                    MarketObjectSP       mo) const{
    /* Check if there is a mdf in place, or create one otherwise */
    checkMDF();

    // delegate to MDF
    mdf->getComponentMarketData(this, market, mo);
}

/** Same as above, but tells the model what the 'domestic' currency is.
    Here domestic means what currency we want things to be in eg the
    instrument's discount currency */
void CModel::getComponentMarketData(const MarketData* market,
                                    MarketObjectSP    mo,
                                    const string&     domesticYCName) const{
    /* Check if there is a mdf in place, or create one otherwise */
    checkMDF();

    mdf->getComponentMarketData(this, market, mo, domesticYCName);
}

 
/** Return the name of the 'domestic' yield curve - only available (at best)
    when fetching market data. See comments in header file */
const string& CModel::getDomesticYCName() const{
    /* Check if there is a mdf in place, or create one otherwise */
    checkMDF();

    return mdf->getDomesticYCName();
}


/** Just route the call to mdf */
CorrelationBaseSP CModel::getCorrelation(const string&     corrName,
                                         CClassConstSP     type1,
                                         CClassConstSP     type2,
                                         CClassConstSP     corrType,
                                         const MarketData* market) const {
    /* Check if there is a mdf in place, or create one otherwise */
    checkMDF();

    return mdf->getCorrelation(this, corrName, type1, type2, corrType, market);
}


/** Just route the call to mdf */
MarketObjectSP CModel::modifyMarketData(
    const MarketData*     market,
    const CClassConstSP&  clazz,      // what type was originally requested
    const MarketObjectSP& mo) const   /* what GetMarket returned or what was
                                       * "inline" already */
{
    /* Check if there is a mdf in place, or create one otherwise */
    checkMDF();

    return mdf->modifyMarketData(this, market, clazz, mo);
}


/** Invoked after instrument has got its market data. Allows model to
    get any extra data required. Default implementation does nothing */
void CModel::getMarket(const MarketData*  market,
                       IInstrumentCollectionSP instruments){}


/** Shifts the object using given shift. */
bool CModel::sensShift(Theta* shift){
    valueDate = shift->rollDate(valueDate);
    return true; // the components can have theta
}


/** get valueDate */
DateTime CModel::getValueDate() const {
    return valueDate;
}

/** Returns the market data fetcher */
MarketDataFetcherSP CModel::getMDF() const {
    /* Check if there is a mdf in place, or create one otherwise */
    checkMDF();

    return mdf;
}

/** Creates a market data fetcher - By default create a brand new mdf */
MarketDataFetcherSP CModel::createMDF() const {
    return MarketDataFetcherSP(new MarketDataFetcher());
}

/** returns a PriceCounter - a "model" that does everything except
actually price, so you can get a count of the times it was asked to price */
PriceCounter* CModel::priceCounter() {
    return new PriceCounter(this);
}

/** returns an ExposureHighlighter - a "model" that does everything except
    actually price, so you get to see what market data it uses 
    Default implementation supplied */
ExposureHighlighter* CModel::exposureHighlighter() {
    return new ExposureHighlighter(this);
}

/** Releases a previously created MDF */
void CModel::releaseMDF() const {
    mdf.reset();
}

CModel::CModel(CClassConstSP clazz): CObject(clazz), 
                                     mdf(NULL){}

//// normally 'today' is populated when retrieving market data, however
//// if this step has been missed then the date needs to be explicitly
//// specified
CModel::CModel(CClassConstSP clazz, const DateTime& today):
    CObject(clazz), valueDate(today){}

class ModelHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CModel, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IModel);
        IMPLEMENTS(Theta::Shift);
        FIELD(valueDate,"");
        FIELD_MAKE_TRANSIENT(valueDate);
    }
};

CClassConstSP const CModel::TYPE = CClass::registerClassLoadMethod(
    "Model", typeid(CModel), ModelHelper::load);

DEFINE_TEMPLATE_TYPE_WITH_NAME("ModelArray", CModelArray);

CClassConstSP const CModel::IModelIntoProduct::TYPE = 
CClass::registerInterfaceLoadMethod(
    "ModelIntoProduct", typeid(CModel::IModelIntoProduct), 0);


void CModel::IProdCreator::addResetDate(DateTime resetDate) {
    DateTimeArray resetDates(1, resetDate);
    addResetDates(resetDates);
}

double CModel::IProdCreator::getValue(DateTime date, CashflowInfo &cfi ) const
{
    cfi.updateAmountType(CashflowInfo::AmountType::UNKNOWN);
    return 0.;
}

void CModel::IProdCreator::setup(const IModel* model, const MarketData* market)
{
    throw ModelException(getClass()->getName()
        +": CModel::IProdCreator::setup(const IModel* model, const MarketData* market) not implemented");
}

void CModel::IProdCreator::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(CModel::IProdCreator, clazz);
    EXTENDS(IObject);
}

typedef CModel::IProdCreator CModelIProdCreator;
CClassConstSP const CModelIProdCreator::TYPE = 
CClass::registerInterfaceLoadMethod(
    "IProdCreator", typeid(CModel::IProdCreator), CModelIProdCreator::load);

typedef CModel::IProdCreatorArray CModelIProdCreatorArray;
DEFINE_TEMPLATE_TYPE_WITH_NAME("IProdCreatorArray", CModelIProdCreatorArray);

class ModelsListAddin: public CObject{
public:
    static CClassConstSP const TYPE;
    // single addin parameter 
    string   instrumentName;

    /** Returns an array of strings which are the model types which can 
        be used with the given product */
    static IObjectSP listModels(ModelsListAddin *params){
        CClassConstSP    clazz = CClass::forName(params->instrumentName);
        CStringArraySP   models(new CStringArray(0));

        // get every class that supports IIntoProduct
        const CClassVec& intoProduct = 
            Model::IModelIntoProduct::TYPE->getAssignableClasses();

        // see if the product supports any of these implementations
        for (unsigned int i = 0; i < intoProduct.size(); i++){
            // hide private classes from spreadsheet
            if (!(Modifier::isPrivate(intoProduct[i]->getModifiers()))) {
                if (clazz != intoProduct[i] && 
                    intoProduct[i]->isAssignableFrom(clazz)) {
                    models->push_back(intoProduct[i]->getName());
                }
            }
        }
        sort(models->begin(), models->end());
        return models;
    }
    ModelsListAddin(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ModelsListAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(instrumentName, "Instrument Name");
        Addin::registerClassObjectMethod("LIST_MODELS",
                                         Addin::UTILITIES,
                                         "Lists models that can be used to"
                                         " price product of supplied type",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)listModels);
    }
    
    static IObject* defaultConstructor(){
        return new ModelsListAddin();
    }
};

CClassConstSP const ModelsListAddin::TYPE = CClass::registerClassLoadMethod(
    "ModelsListAddin", typeid(ModelsListAddin), load);

DRLIB_END_NAMESPACE
