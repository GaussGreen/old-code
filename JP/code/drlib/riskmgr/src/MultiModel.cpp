//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids
//
//   Filename    : MultiModel.cpp
//
//   Description : A model which contains an array of other models
//
//   Author      : Linus Thand
//
//   Date        : 24 May 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MultiModel.hpp"
#include "edginc/Scenario.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/CorrelationBase.hpp"
#include "edginc/MarketData.hpp"

DRLIB_BEGIN_NAMESPACE


MultiModel::~MultiModel() {}

MultiModel::MultiModel() : CObject(TYPE) {} 

IModelSP MultiModel::getModel(int n) const {
    return (*models)[n];
}

/** Does everything, namely retrieve market data, applies scenario 
    (if non null), calculates price and greeks as well as save to file of
    inputs and outputs. */      
CResultsArraySP MultiModel::go(IInstrumentCollectionSP instruments,
                               ScenarioSP              scenario, // optional
                               CControlSP              control,
                               MarketDataSP            market){
    IModelSP model(IModelSP::attachToRef(this));                           
    if (models->size() != instruments->size()) {       
        throw ModelException("MultiModel requires the same number"
                             " of models as instruments");
    }

    bool writeToFile = control->getWriteToFile();
    if (writeToFile) {
        if (!scenario){
            CRiskMgr::writeInputs(model, instruments, control, market);
        } else {
            scenario->writeInputs(model, instruments, control, market);
        }
    }
    CResultsArraySP results(miniGo(instruments,
                                   scenario, control, market));
    if (writeToFile) {
        CRiskMgr::writeOutputs(control.get(), results);
    }
    return results;
}


CResultsSP MultiModel::go(CInstrumentSP instrument,
                          ScenarioSP    scenario, // optional
                          CControlSP    control,
                          MarketDataSP  market){
    IModelSP model(IModelSP::attachToRef(this));
    if (models->size() != 1) {       
        throw ModelException("MultiModel requires the same number"
                             " of models as instruments");
    }
    bool writeToFile = control->getWriteToFile();
    if (writeToFile) {
        // Need to make RiskMgr::writeInputs and Scenario::writeInputs public
        if (!scenario){
            CRiskMgr::writeInputs(model, instrument, control, market);
        } else {
            scenario->writeInputs(model, instrument, control, market);
        }
    }
    CResultsArraySP results(miniGo(IInstrumentCollection::singleton(instrument),
                                   scenario, control, market));
    CResultsSP result(results->front());
    if (writeToFile) {
        CRiskMgr::writeOutputs(control.get(), result);
    }
    return result;
}


/** Same as go except for save to file of inputs and outputs. */
CResultsArraySP MultiModel::miniGo(IInstrumentCollectionSP instruments,
                                   ScenarioSP              scenario, // optional
                                   CControlSP              control,
                                   MarketDataSP            market){
    /** Get market data:
     work on a copy of the instrument since we're possibly amending data
     when doing the 0 day theta shift */
    IInstrumentCollectionSP insts(instruments.clone());
    //// copy model as well since we might be filling it with market data
    IModelSP mdl(IModelSP::attachToRef(this).clone());  
    //// copy control in case of state info inside sensitivities
    CControlSP ctrl(control.clone());

    /** For scenario shifts that need to be applied before.
     Apply scenario if non null. */
    if (scenario.get()) {
        scenario->preapply(mdl, insts);
    }

    for (int i = 0; i < models->size(); ++i) {
        //// Call model to collect instrument and model market data
        (*models)[i]->getInstrumentsAndModelMarket(market, 
                        IInstrumentCollection::singleton((*insts)[i]));     
        //// Call getMarket to initiate market data selection.
        ctrl->getMarket((*models)[i], 
                        market, 
                        IInstrumentCollection::singleton((*insts)[i]));
    } 

    /** For scenario shifts that need to be applied after
     apply scenario if non null */
    if (scenario.get()) {
        scenario->apply(mdl, insts);
    }

    return RiskMgr::calculateTheoAndMtm(IModelSP::attachToRef(this), insts, ctrl);
}

CResultsArraySP MultiModel::RunMulti(IInstrumentCollectionSP instruments,
                                 CControl* control) {
    CResultsArraySP resultss = instruments->emptyResults();
    instruments->Validate();
    control->calculateMulti(this, instruments, resultss);
    return resultss;
}

void MultiModel::PriceMulti(IInstrumentCollectionSP instruments, 
                            CControl* control, 
                           CResultsArraySP resultss) {
    ASSERT(instruments->size() == resultss->size());
    for (int i = 0; i < instruments->size(); ++i){
        (*models)[i]->Price((*instruments)[i].get(), 
                            control, 
                            (*resultss)[i].get());
    }
}

void MultiModel::flush() {
    for (int i = 0; i < models->size(); ++i) {
        (*models)[i]->flush();
    }
}

PriceCounter* MultiModel::priceCounter() {
    throw ModelException("MultiModel::priceCounter() is not implemented.");
}


void MultiModel::Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results) {
    throw ModelException("MultiModel::Price was called unexpectedly.");
}
                                           
IModel::WantsRiskMapping MultiModel::wantsRiskMapping() const {
    throw ModelException(
        "MultiModel::wantsRiskMapping was called unexpectedly.");
}                  

ExposureHighlighter* MultiModel::exposureHighlighter() {
    throw ModelException(
        "MultiModel::exposureHighlighter was called unexpectedly.");
}

CResults* MultiModel::Run(CInstrument* instrument, 
                          CControl*    control) {
    throw ModelException("MultiModel::Run was called unexpectedly.");
}

double MultiModel::calcPrice(CInstrument*  instrument, 
                             CControl*     control) {
    throw ModelException("MultiModel::calcPrice was called unexpectedly.");
}

MarketObjectSP MultiModel::getMarketObject(const MarketData* market,
                                           const string&     name,
                                           CClassConstSP     type,
                                           const string&     domesticYCName) 
                                           const {
    throw ModelException(
        "MultiModel::getMarketObject was called unexpectedly.");
}

MarketObjectSP MultiModel::GetMarket(const MarketData*    market,
                                     const string&        name,
                                     const CClassConstSP& type) const {
    throw ModelException("MultiModel::GetMarket was called unexpectedly.");
}

void MultiModel::setDomesticYCName(string discountCurveName) const {
    throw ModelException(
        "MultiModel::setDomesticYCName was called unexpectedly.");
}

void MultiModel::getInstrumentsAndModelMarket(MarketDataConstSP       market,
                                              IInstrumentCollectionSP insts) {
    throw ModelException(
        "MultiModel::getInstrumentAndModelMarket was called unexpectedly.");
}

void MultiModel::getInstrumentAndModelMarket(const MarketData*  market,
                                             CInstrument*       inst) {
    throw ModelException(
        "MultiModel::getInstrumentAndModelMarket was called unexpectedly.");
}

SensControl* MultiModel::AlterControl(const SensControl* currSensControl) 
    const { 
    throw ModelException("MultiModel::checkMDF was called unexpectedly.");
}

void MultiModel::getComponentMarketData(const MarketData*    market,
                                        MarketObjectSP       mo) const{
    throw ModelException(
        "MultiModel::getComponentMarketData was called unexpectedly.");
}

void MultiModel::getComponentMarketData(const MarketData* market,
                                        MarketObjectSP    mo,
                                        const string&     domesticYCName) const{
    throw ModelException(
        "MultiModel::getComponentMarketData was called unexpectedly.");
}

const string& MultiModel::getDomesticYCName() const{
    throw ModelException(
        "MultiModel::getDomesticYCName was called unexpectedly.");
}

CorrelationBaseSP MultiModel::getCorrelation(const string&     corrName,
                                             CClassConstSP     type1,
                                             CClassConstSP     type2,
                                             CClassConstSP     corrType,
                                             const MarketData* market) const {
    throw ModelException("MultiModel::getCorrelation was called unexpectedly.");
}

MarketObjectSP MultiModel::modifyMarketData(const MarketData*     market,
                                            const CClassConstSP&  clazz,      
                                            const MarketObjectSP& mo) const                   
{
    throw ModelException(
        "MultiModel::modifyMarketData was called unexpectedly.");
}

void MultiModel::getMarket(const MarketData*  market,
                           IInstrumentCollectionSP instruments)
{
    throw ModelException("MultiModel::getMarket was called unexpectedly.");
}

DateTime MultiModel::getValueDate() const {
    throw ModelException("MultiModel::getValueDate was called unexpectedly.");
}

MarketDataFetcherSP MultiModel::getMDF() const {
    throw ModelException("MultiModel::getMDF was called unexpectedly.");
}

MultiModel::MultiModel(CClassConstSP clazz): CObject(clazz) {}

IObject* MultiModel::defaultConstructor(){
    return new MultiModel();
}

/** Invoked when class is 'loaded' */
void MultiModel::load(CClassSP& clazz){
    clazz->setPublic(); //// make visible to EAS/spreadsheet
    REGISTER(MultiModel, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IModel);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(models,"wrapped models");
}

CClassConstSP const MultiModel::TYPE = CClass::registerClassLoadMethod(
    "MultiModel", typeid(MultiModel), load);

bool MultiModelLinkIn ()
{
    return MultiModel::TYPE != NULL;
}


DEFINE_TEMPLATE_TYPE(IModelArray);

DRLIB_END_NAMESPACE
