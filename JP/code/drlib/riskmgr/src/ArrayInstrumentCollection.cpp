/**
 * @file ArrayInstrumentCollection.cpp
 */

#include "edginc/config.hpp"
#include "edginc/ArrayInstrumentCollection.hpp"
#include "edginc/ScaleOutputs.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/Results.hpp"
#include "edginc/MarketData.hpp"
#include <set>

DRLIB_BEGIN_NAMESPACE

ArrayInstrumentCollection::ArrayInstrumentCollection(CInstrumentArraySP instruments):
    CInstrumentCollection(TYPE),
    instruments(instruments)
{
    if (instruments->empty()) {
        throw ModelException("ArrayInstrumentCollection::ArrayInstrumentCollection",
                             "'instruments' may not be empty");
    }
}

ArrayInstrumentCollection::ArrayInstrumentCollection(): // for reflection
    CInstrumentCollection(TYPE)
{}

ArrayInstrumentCollection::~ArrayInstrumentCollection() {}

int ArrayInstrumentCollection::size() const {
    return instruments->size();
}

InstrumentSP ArrayInstrumentCollection::operator [](int i) {
    return (*instruments)[i];
}

InstrumentConstSP ArrayInstrumentCollection::operator [](int i) const {
    return (*instruments)[i];
}

void ArrayInstrumentCollection::Validate() {
    if (instruments->empty()) {
        throw ModelException("ArrayInstrumentCollection::Validate",
                             "'instruments' may not be empty");
    }

    for (int i = 0; i < size(); ++i) {
        (*this)[i]->Validate();
    }
}

void ArrayInstrumentCollection::GetMarket(const IModel* model,
                                          const CMarketDataSP market) {
    for (int i = 0; i < size(); ++i) 
        (*this)[i]->GetMarket(model, market);
}

DateTime ArrayInstrumentCollection::getValueDate() const {
    ASSERT(size() > 0);
    return (*this)[0]->getValueDate();
}

string ArrayInstrumentCollection::discountYieldCurveName() const {
    ASSERT(size() > 0);
    return (*this)[0]->discountYieldCurveName();
}

DateTime ArrayInstrumentCollection::endDate(
        const Sensitivity* sensitivity) const {
    DateTime it = getValueDate();
    DateTime whenever = MaturityPeriod("50Y").toDate(getValueDate()); // FIXME yuk

    for (int i = 0; i < size(); ++i) {
        const LastSensDate *lsd =
            dynamic_cast<const LastSensDate *>((*instruments)[i].get());
        it = max(it, lsd ? lsd->endDate(sensitivity) : whenever);
    }

    return it;
}

bool ArrayInstrumentCollection::avoidVegaMatrix(const IModel* model) {
    if (size() == 0)
        return false;

    for (int i = 0; i < size(); ++i) {
        IInstrument *inst = (*instruments)[i].get();
        if (!(dynamic_cast<ISensitiveStrikes *>(inst) &&
              dynamic_cast<ISensitiveStrikes *>(inst)->avoidVegaMatrix(model)))
            return false;
    }

    return true;
}

DoubleArraySP ArrayInstrumentCollection::getSensitiveStrikes(
       OutputNameConstSP outputName,
       const IModel* model) {
    set<double> all;

    //FIXME merge strikes if close together -- see VegaMatrix

    for (int i = 0; i < size(); ++i) {
        IInstrument *inst = (*instruments)[i].get();
        if (dynamic_cast<ISensitiveStrikes *>(inst)) {
            DoubleArraySP his = dynamic_cast<ISensitiveStrikes *>(inst)->
                getSensitiveStrikes(outputName, model);
            for (int s = 0; s < his->size(); ++s)
                all.insert((*his)[s]);
        }
    }

    DoubleArraySP them(new DoubleArray());
    for (set<double>::iterator it = all.begin(); it != all.end(); ++it)
        them->push_back(*it);
    return them;
}

void ArrayInstrumentCollection::Price(IModel *model,
                                      CControl *control,
                                      CResultsArraySP resultss) {
    model->PriceMulti(IInstrumentCollectionSP::attachToRef(this),
                      control,
                      resultss);
}

void ArrayInstrumentCollection::scaleOutputs(CControlSP control,
                                             CResultsArraySP results) {
    ASSERT(results->size() == instruments->size());
    for (int i = 0; i < results->size(); ++i) {
        IInstrument *inst = (*instruments)[i].get();
        if (IScaleOutputs::TYPE->isInstance(inst))
            dynamic_cast<IScaleOutputs *>(inst)->scaleOutputs(control,
                                                              (*results)[i]);
    }
}

IObject* ArrayInstrumentCollection::defaultArrayInstrumentCollection() {
    return new ArrayInstrumentCollection();
}

void ArrayInstrumentCollection::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(ArrayInstrumentCollection, clazz);
    SUPERCLASS(CInstrumentCollection);
    IMPLEMENTS(ISensitiveStrikes);
    EMPTY_SHELL_METHOD(defaultArrayInstrumentCollection);
    FIELD(instruments, "instruments");
}

CClassConstSP ArrayInstrumentCollection::TYPE =
    CClass::registerClassLoadMethod(
        "ArrayInstrumentCollection", typeid(ArrayInstrumentCollection),
        load);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
