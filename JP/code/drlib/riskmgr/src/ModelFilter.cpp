//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : ModelFilter.cpp
//
//   Contains    : ModelFilter, PriceCounter, ExposureHighlighter,
//                 IRVegaPointwiseExposureHighlighter
//
//   Description : ModelFilter is an abstract model that doesn't actually price.
//                 How useful is that?  Handy for when you want to find out things
//                 like how may pricings are there (PriceCounter) or what market
//                 objects am I sensitive to (ExposureHighlighter and
//                 IRVegaPointwiseExposureHighlighter) without a full-blown pricing run.
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_MODEL_FILTER_CPP
#include "edginc/ModelFilter.hpp"
#include "edginc/RiskMgr.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/CorrelationBase.hpp"
#include "edginc/MarketData.hpp"

DRLIB_BEGIN_NAMESPACE

// the ModelFilter "non-pricing" model

ModelFilter::ModelFilter(IModel* model) : Model(TYPE), realModel(copy(model)) {}

ModelFilter::ModelFilter(IModel* model, CClassConstSP clazz): 
    Model(clazz), realModel(copy(model)) {}

ModelFilter::ModelFilter(CClassConstSP clazz) : Model(clazz) {}

// IMPORTANT: Does a shallow copy!  Only intended for use by shallowCopy() via clone().
ModelFilter::ModelFilter(const ModelFilter &rhs, CClassConstSP clazz) :
    Model(clazz),
    realModel(rhs.realModel) {}

// all the dull pass-throughs

CResults* ModelFilter::Run(CInstrument* instrument, 
                                   CControl*    control) {
    return CModel::Run(instrument, control);
}

CResultsArraySP ModelFilter::RunMulti(
    IInstrumentCollectionSP instruments, 
    CControl* control) {

    static const string method = "ModelFilter::RunMulti";

    CResultsArraySP resultss = instruments->emptyResults();
    instruments->Validate();

    control->calculateMulti(this, instruments, resultss);

    return resultss;
}

void ModelFilter::PriceMulti(IInstrumentCollectionSP instruments, 
                                     CControl* control, 
                                     CResultsArraySP results) {
    CModel::PriceMulti(instruments, control, results);
}

MarketObjectSP ModelFilter::getMarketObject(
    const MarketData* market,
    const string&     name,
    CClassConstSP     type,
    const string&     domesticYCName) const {
    return realModel->getMarketObject(market,name,type,domesticYCName);
}

void ModelFilter::setDomesticYCName(
    string discountYieldCurveName) const {
    realModel->setDomesticYCName(discountYieldCurveName);
}

void ModelFilter::getInstrumentAndModelMarket(
    const MarketData*  market,
    CInstrument* inst) {
    realModel->getInstrumentAndModelMarket(market, inst);
}

void ModelFilter::getInstrumentsAndModelMarket (
    MarketDataConstSP market,
    IInstrumentCollectionSP insts) {
    realModel->getInstrumentsAndModelMarket(market, insts);
}

MarketObjectSP ModelFilter::GetMarket(
    const MarketData*    market,
    const string&        name,
    const CClassConstSP& type) const {
    return realModel->GetMarket(market, name, type);
}

void ModelFilter::getComponentMarketData(
    const MarketData*    market,
    MarketObjectSP       mo) const {
    realModel->getComponentMarketData(market, mo);
}

void ModelFilter::getComponentMarketData(
    const MarketData* market,
    MarketObjectSP    mo,
    const string&     domesticYCName) const {
    realModel->getComponentMarketData(market, mo, domesticYCName);
}

CorrelationBaseSP ModelFilter::getCorrelation(
    const string&     corrName,
    CClassConstSP     type1,
    CClassConstSP     type2,
    CClassConstSP     corrType,
    const MarketData* market) const {
    return realModel->getCorrelation(corrName,type1,type2,corrType,market);
}

MarketObjectSP ModelFilter::modifyMarketData(
    const MarketData*     market,
    const CClassConstSP&  clazz,     // what type was originally requested
    const MarketObjectSP& mo) const {
    return realModel->modifyMarketData(market,clazz,mo);
}

void ModelFilter::getMarket(
    const MarketData*  market,
    IInstrumentCollectionSP instruments) {
    realModel->getMarket(market, instruments);
}

const string& ModelFilter::getDomesticYCName() const {
    return realModel->getDomesticYCName();
}

DateTime ModelFilter::getValueDate() const {
    return realModel->getValueDate();
}

SensControl* ModelFilter::AlterControl(
    const SensControl* currSensControl) const {
    return realModel->AlterControl(currSensControl);
}

void ModelFilter::flush() {
    realModel->flush();
}

IModel::WantsRiskMapping ModelFilter::wantsRiskMapping() const {
    return realModel->wantsRiskMapping();
}

// Override clone since we likely need to preserve internal state so that
// the PriceCounter's count and ExposureHighlighter's dummyPrice are maintained.
IObject* ModelFilter::clone() const {
    if (getRefCount() == 0)
        // not being accessed via smart pointer - object probably on stack
        return shallowCopy();
    else
        return const_cast<ModelFilter*>(this);
}

// when to stop tweaking (product-provided)
DateTime ModelFilter::endDate(const CInstrument* instrument,
                                      const Sensitivity* sensControl) const
{
    LastProductSensDate *lpsd = dynamic_cast<LastProductSensDate*>(realModel.get());
    if (lpsd)
        return lpsd->endDate(instrument, sensControl);

    throw ModelException("ModelFilter::endDate",
        "Underlying model " + realModel->getClass()->getName() +
        " does not implement product-provided endDate.");
}

class ModelFilterHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ModelFilter, clazz);
        SUPERCLASS(CModel);
        IMPLEMENTS(LastProductSensDate);
        FIELD(realModel, "model that would really be used for pricing");
    }
};

CClassConstSP const ModelFilter::TYPE = CClass::registerClassLoadMethod(
    "ModelFilter", typeid(ModelFilter), ModelFilterHelper::load);

bool ModelFilterLinkIn() {
    return ModelFilter::TYPE != NULL;
}


//////////////////////////////////////////////////////////////////////////
// PriceCounter
//////////////////////////////////////////////////////////////////////////
PriceCounter::PriceCounter(IModel *realModel) :
    ModelFilter(realModel, TYPE), count(0) {}

PriceCounter::PriceCounter(IModel *realModel, CClassConstSP clazz) :
    ModelFilter(realModel, clazz), count(0) {}

PriceCounter::PriceCounter(CClassConstSP clazz) : ModelFilter(clazz), count(0) {}

// See ModelFilter::ModelFilter(const ModelFilter &rhs, CClassConstSP clazz)
PriceCounter::PriceCounter(const PriceCounter &rhs, CClassConstSP clazz): 
    ModelFilter(rhs, clazz), count(rhs.count) {}

/** calculate single price and store result in CResult */
void PriceCounter::Price(CInstrument*  instrument, 
                         CControl*     control, 
                         CResults*     results) {
    static const string method("PriceCounter::Price");
    try {
        // Keep track of number of times we're asked to price and store dummy value.
        ++count;
        results->storePrice(0.0, "");
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

// how often was this model called for pricing?
int PriceCounter::pricings() const {
    return count;
} 

// For clone
PriceCounter* PriceCounter::shallowCopy() const {
    return new PriceCounter(*this, TYPE);
}

class PriceCounterHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to clients
        REGISTER(PriceCounter, clazz);
        SUPERCLASS(ModelFilter);
        EMPTY_SHELL_METHOD(defaultPriceCounter);
        FIELD(count, "");
        FIELD_MAKE_TRANSIENT(count);
    }

    static IObject* defaultPriceCounter(){
        return new PriceCounter(PriceCounter::TYPE);
    }
};

CClassConstSP const PriceCounter::TYPE = CClass::registerClassLoadMethod(
    "PriceCounter", typeid(PriceCounter), PriceCounterHelper::load);


//////////////////////////////////////////////////////////////////////////
// ExposureHighlighter
//////////////////////////////////////////////////////////////////////////
ExposureHighlighter::ExposureHighlighter(IModel *realModel) : 
    ModelFilter(realModel, TYPE), dummyPrice(0.0) {}

ExposureHighlighter::ExposureHighlighter(IModel *realModel, CClassConstSP clazz) : 
    ModelFilter(realModel, clazz), dummyPrice(0.0) {}

ExposureHighlighter::ExposureHighlighter(CClassConstSP clazz): 
    ModelFilter(clazz), dummyPrice(0.0) {}

// See ModelFilter::ModelFilter(const ModelFilter &rhs, CClassConstSP clazz)
ExposureHighlighter::ExposureHighlighter(const ExposureHighlighter &rhs,
                                         CClassConstSP clazz) : 
    ModelFilter(rhs, clazz), dummyPrice(rhs.dummyPrice) {}

/** calculate single price and store result in CResult */
void ExposureHighlighter::Price(CInstrument*  instrument, 
                         CControl*     control, 
                         CResults*     results) {
    static const string method("PriceCounter::Price");
    try {
        // We just need to store an always changing dummy price so that we appear
        // to have sensitivity whenever asked to reprice in a tweaked environment.
        // NOTE: This is only needed for old-style greeks since new-style greeks
        // can determine exposure without actually tweaking and repricing.
        results->storePrice(++dummyPrice, "");
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

CResultsArraySP ExposureHighlighter::RunMulti(
    IInstrumentCollectionSP instruments, 
    CControl* control) {

    static const string method = "ExposureHighlighter::RunMulti";

    CResultsArraySP resultss = instruments->emptyResults();
    instruments->Validate();

    control->calculateMultiExposures(this, instruments, resultss);

    return resultss;
}

// For clone
ExposureHighlighter* ExposureHighlighter::shallowCopy() const {
    return new ExposureHighlighter(*this, TYPE);
}

class ExposureHighlighterHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to clients
        REGISTER(ExposureHighlighter, clazz);
        SUPERCLASS(ModelFilter);
        EMPTY_SHELL_METHOD(defaultExposureHighlighter);
        FIELD(dummyPrice, "");
        FIELD_MAKE_TRANSIENT(dummyPrice);
    }

    static IObject* defaultExposureHighlighter(){
        return new ExposureHighlighter(ExposureHighlighter::TYPE);
    }
};

CClassConstSP const ExposureHighlighter::TYPE = CClass::registerClassLoadMethod(
    "ExposureHighlighter", typeid(ExposureHighlighter), ExposureHighlighterHelper::load);


//////////////////////////////////////////////////////////////////////////
// IRVegaPointwiseExposureHighlighter
//////////////////////////////////////////////////////////////////////////
IRVegaPointwiseExposureHighlighter::IRVegaPointwiseExposureHighlighter(
    IModel *realModel) : ExposureHighlighter(realModel, TYPE) {}

IRVegaPointwiseExposureHighlighter::IRVegaPointwiseExposureHighlighter(
    IModel *realModel, CClassConstSP clazz) : ExposureHighlighter(realModel, clazz) {}

IRVegaPointwiseExposureHighlighter::IRVegaPointwiseExposureHighlighter(CClassConstSP clazz): 
    ExposureHighlighter(clazz) {}

// See ExposureHighlighter::ExposureHighlighter(const ExposureHighlighter &rhs, CClassConstSP clazz)
IRVegaPointwiseExposureHighlighter::IRVegaPointwiseExposureHighlighter(
    const IRVegaPointwiseExposureHighlighter &rhs, CClassConstSP clazz) :
        ExposureHighlighter(rhs, clazz) {}

IRGridPointAbsArraySP IRVegaPointwiseExposureHighlighter::getSensitiveIRVolPoints(
        OutputNameConstSP  outputName,
        const CInstrument* inst) const {
    IRVegaPointwise::ISensitivePoints* sp = 
        dynamic_cast<IRVegaPointwise::ISensitivePoints*>(realModel.get());
    if (!sp)
        throw ModelException("Underlying model (" + realModel->getClass()->getName() +
                             ") of IRVegaPointwiseExposureHighlighter is "
                             "expected to implement ISensitivePoints.");
    return sp->getSensitiveIRVolPoints(outputName, inst);
}

// For clone
IRVegaPointwiseExposureHighlighter* IRVegaPointwiseExposureHighlighter::shallowCopy() const {
    return new IRVegaPointwiseExposureHighlighter(*this, TYPE);
}

class IRVegaPointwiseExposureHighlighterHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to clients
        REGISTER(IRVegaPointwiseExposureHighlighter, clazz);
        SUPERCLASS(ExposureHighlighter);
   //     IMPLEMENTS(IRVegaPointwise::ISensitivePoints);    /* Can't register since IRVegaPointwise::ISensitivePoints is not registered as an interface. */
        EMPTY_SHELL_METHOD(defaultIRVegaPointwiseExposureHighlighter);
    }

    static IObject* defaultIRVegaPointwiseExposureHighlighter(){
        return new IRVegaPointwiseExposureHighlighter(IRVegaPointwiseExposureHighlighter::TYPE);
    }
};

CClassConstSP const IRVegaPointwiseExposureHighlighter::TYPE =
    CClass::registerClassLoadMethod("IRVegaPointwiseExposureHighlighter",
        typeid(IRVegaPointwiseExposureHighlighter),
        IRVegaPointwiseExposureHighlighterHelper::load);

DRLIB_END_NAMESPACE
