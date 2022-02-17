//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids 
//
//   Description : InstrumentCollection with weights and caching
//
//   Author      : Linus Thand 
//
//   Date        : 10 July 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/WeightedInstrumentCollection.hpp"
#include "edginc/MultiModel.hpp"
#include "edginc/ScaleOutputs.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/Results.hpp"
#include "edginc/MarketData.hpp"

DRLIB_BEGIN_NAMESPACE


////WeightedInstrumentPriceCache

CInstrumentSP WeightedInstrumentPriceCache::getInstrument(void) { 
    return inst; 
}

CInstrumentConstSP WeightedInstrumentPriceCache::getInstrument(void) const {
    return CInstrumentConstSP(inst);
}

void WeightedInstrumentPriceCache::markDirty(void) { 
        cacheDirty = true; 
}

WeightedInstrumentPriceCache::WeightedInstrumentPriceCache(CInstrumentSP inst) : 
    CObject(TYPE), inst(inst.clone()), cacheDirty(true)  {}

WeightedInstrumentPriceCache::WeightedInstrumentPriceCache() : 
     CObject(TYPE), cacheDirty(true) {}

WeightedInstrumentPriceCache::~WeightedInstrumentPriceCache() {}

void WeightedInstrumentPriceCache::fieldsUpdated(const CFieldArray& fields) {
    cacheDirty = true;
}

CResultsSP WeightedInstrumentPriceCache::Price(IModelSP model, 
                                               CControlSP control, 
                                               const double weight) 
{
    
    if (control.get() != originalCtrl.get()) {   //// If a new control is used
        localCtrl = CControlSP(control.clone()); //// Keep a local copy
        originalCtrl = control;                  //// Keep reference future comparison
        cacheDirty = true;                       //// Need to calculate new results
    }

    if (model.get() != originalMdl.get()) {     //// If a new model is used
        localMdl = IModelSP(model.clone());     //// Keep a local copy
        originalMdl = model;                    //// Keep reference future comparison
        cacheDirty = true;                      //// Need to calculate new results
    }
  
    CResultsSP results(new CResults());

    if (cacheDirty) { //// Calculate results and populate cache
        localMdl->Price(inst.get(), localCtrl.get(), results.get());
        //// Make sure output requests are calculated correctly.
        if (localCtrl->isPricing()){
            CResultsArraySP myResults = 
                CResultsArraySP(new CResultsArray(1, results));
            localCtrl->handleUnfulfilledRequests(
                localMdl.get(), IInstrumentCollection::singleton(inst), myResults);
        }
        cachedResults.reset(results.clone());
        cacheDirty = false;
    } else { //// Use cached results
        results = CResultsSP(cachedResults.clone());
    }

    results->scale(localCtrl, weight, false); //// Scale by weight

    return results;
}

void WeightedInstrumentPriceCache::load(CClassSP& clazz){
    REGISTER(WeightedInstrumentPriceCache, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(inst, "instrument");
    FIELD(cachedResults, "cached results");
    FIELD(cacheDirty, "is cache dirty?");
    FIELD(originalCtrl, "reference to source of local control");
    FIELD(localCtrl, "local copy of the control");
    FIELD(originalMdl, "reference to source of local model");
    FIELD(localMdl, "local copy of the model");
}

IObject* WeightedInstrumentPriceCache::defaultConstructor(){
    return new WeightedInstrumentPriceCache();
}

CClassConstSP const WeightedInstrumentPriceCache::TYPE = 
    CClass::registerClassLoadMethod("WeightedInstrumentPriceCache", 
                                        typeid(WeightedInstrumentPriceCache), 
                                        load);
/*
template<> CClassConstSP const WeightedInstrumentPriceCacheArray::TYPE =
    CClass::registerClassLoadMethod("WeightedInstrumentPriceCacheArray", 
                                    typeid(WeightedInstrumentPriceCacheArray), 
                                    load);
*/
DEFINE_TEMPLATE_TYPE(WeightedInstrumentPriceCacheArray);

//// WeightedInstrumentCollection

WeightedInstrumentCollection::~WeightedInstrumentCollection() {}

int WeightedInstrumentCollection::size() const {
    return priceCaches->size();
}

InstrumentSP WeightedInstrumentCollection::operator [](int i) {
    return (*priceCaches)[i]->getInstrument();
}

InstrumentConstSP WeightedInstrumentCollection::operator [](int i) const {
    return (*priceCaches)[i]->getInstrument();
}

void WeightedInstrumentCollection::GetMarket(const IModel* model,
                                             const CMarketDataSP market) {
    for (int i = 0; i < size(); ++i) {
        (*this)[i]->GetMarket(model, market);
    }
}

DateTime WeightedInstrumentCollection::getValueDate() const {
    ASSERT(size() > 0);
    return (*this)[0]->getValueDate();
}

string WeightedInstrumentCollection::discountYieldCurveName() const {
    ASSERT(size() > 0);
    return (*this)[0]->discountYieldCurveName();
}

DateTime WeightedInstrumentCollection::endDate(const Sensitivity* sensitivity) 
    const {
    DateTime it = getValueDate();
    DateTime whenever = MaturityPeriod("50Y").toDate(getValueDate()); // FIXME yuk

    for (int i = 0; i < size(); ++i) {
        const LastSensDate *lsd =
            dynamic_cast<const LastSensDate *>((*this)[i].get());
        it = max(it, lsd ? lsd->endDate(sensitivity) : whenever);
    }

    return it;
}

void WeightedInstrumentCollection::scaleOutputs(CControlSP      control,
                                                CResultsArraySP results) {
    ASSERT(results->size() == priceCaches->size());
    for (int i = 0; i < results->size(); ++i) {
        IInstrument *inst = (*this)[i].get();
        if (IScaleOutputs::TYPE->isInstance(inst))
            dynamic_cast<IScaleOutputs *>(inst)->scaleOutputs(control,
                                                              (*results)[i]);
    }
}

TweakOutcome WeightedInstrumentCollection::sensShift(
    const PropertyTweak<WeightedInstrumentTweak>& shift) {
    /** Shifting weight on instrument(s) as specified in the 
        instrumentList */

    IntArrayConstSP insts = shift.tag->instruments;
    
    for (int i = 0; i < insts->size(); ++i) { 
        if ((*insts)[i] > weights.size() - 1) {
            throw ModelException("Instrument to be shifted does"
                                 " not exist in the InstrumentCollection.");
        }
        weights[(*insts)[i]] += shift.coefficient;
    } 
    return TweakOutcome(shift.coefficient, false); 
}

string WeightedInstrumentCollection::sensName(const WeightedInstrumentTweak* tag) 
const {
    return "WeightedInstrumentCollection"; // please don't return "" since that means "don't tweak me"
}

void WeightedInstrumentCollection::sensRestore(
    const PropertyTweak<WeightedInstrumentTweak>& shift){

    IntArrayConstSP insts = shift.tag->instruments;

    /** Unshifting weight on instrument(s) as specified 
        in the instrumentList */

    for (int i = 0; i < insts->size(); ++i) {
        if ((*insts)[i] > weights.size() - 1) {
            throw ModelException("Instrument to be shifted does"
                                 " not exist in the InstrumentCollection.");
        }
        weights[(*insts)[i]] -= shift.coefficient;
    }
}

void WeightedInstrumentCollection::Price(IModel*         model,
                                         CControl*       control,
                                         CResultsArraySP resultss) {
    if (MultiModel::TYPE->isInstance(model)) {
        /** If it's a MultiModel we need to delegate the pricing to the individual
        * models here.
        * Not ideal perhaps, as this functionality already is in MultiModel,
        * but it seems to be the only way of getting the caching to work 
        * in individual instruments. An alternative might be to move the 
        * caching to MultiModel.
        */
        MultiModelSP multiModel =  MultiModelSP::dynamicCast(IModelSP(model));
        for (int i = 0; i < priceCaches->size(); ++i) {
            (*resultss)[i] = (*priceCaches)[i]->Price(multiModel->getModel(i), 
                                                      CControlSP(control), 
                                                      weights[i]);
        }
    } else {
        //// For normal models, don't delegate.
        for (int i = 0; i < priceCaches->size(); ++i) {
            (*resultss)[i] =(*priceCaches)[i]->Price(IModelSP(model), 
                                                     CControlSP(control), 
                                                     weights[i]);
        }
    }
    if (control->isPricing()){
        OutputRequestArrayConstSP requests(control->getOutputRequests());
        for (int i = 0; i < requests->size(); i++){
            (*requests)[i]->setHasFinished(true);
        }
    }
}

WeightedInstrumentCollection::WeightedInstrumentCollection() :
    CInstrumentCollection(TYPE) 
{ }

void WeightedInstrumentCollection::validatePop2Object() {  
    static const string method("WeightedInstrumentCollection::validatePop2Object");
    try {

        //// Create price caches
        priceCaches = WeightedInstrumentPriceCacheArraySP(
            new WeightedInstrumentPriceCacheArray(instruments->size()));

        //// Wrap up each of the instruments in price caches (by cloning)
        for (int i = 0; i < priceCaches->size(); ++i) {
            (*priceCaches)[i] = WeightedInstrumentPriceCacheSP(
                new WeightedInstrumentPriceCache((*instruments)[i]));
        }

        //// Initialise weights
        //// Will create an array of weights of the right length.
        fill_n(weights.back_inserter(),priceCaches->size(),1);

    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

void WeightedInstrumentCollection::Validate() {
    if (instruments->empty()) {
        throw ModelException("WeightedInstrumentCollection::Validate",
                             "'instruments' may not be empty");
    }

    for (int i = 0; i < size(); ++i) {
        (*this)[i]->Validate();
    }
}

/** Invoked when class is 'loaded' */
void WeightedInstrumentCollection::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(WeightedInstrumentCollection, clazz);
    SUPERCLASS(CInstrumentCollection);
    IMPLEMENTS(IRestorableWithRespectTo<WeightedInstrumentTweak>);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(priceCaches, "array of price caches");
    FIELD(instruments, "instruments");
    FIELD(weights, "instrument weights");
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(priceCaches);
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(weights);
    FIELD_MAKE_NONTWEAKABLE(instruments); // Only tweak instruments in caches
}

IObject* WeightedInstrumentCollection::defaultConstructor(){
    return new WeightedInstrumentCollection();
}

CClassConstSP const WeightedInstrumentCollection::TYPE = 
    CClass::registerClassLoadMethod("WeightedInstrumentCollection", 
                                    typeid(WeightedInstrumentCollection), 
                                    load);

bool WeightedInstrumentCollectionLinkIn ()
{
    return WeightedInstrumentCollection::TYPE != NULL;
}

DRLIB_END_NAMESPACE
