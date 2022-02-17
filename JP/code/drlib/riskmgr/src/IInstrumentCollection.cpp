/**
 * @file IInstrumentCollection.cpp
 */

#include "edginc/config.hpp"
#define QLIB_IINSTRUMENTCOLLECTION_CPP
#include "edginc/ArrayInstrumentCollection.hpp"

DRLIB_BEGIN_NAMESPACE

IInstrumentCollection::IInstrumentCollection() {}

IInstrumentCollection::~IInstrumentCollection() {}

IInstrumentCollectionSP IInstrumentCollection::singleton(
    InstrumentSP instrument) {

    CInstrumentArraySP justTheOne(new CInstrumentArray(1, instrument));
    return ArrayInstrumentCollectionSP(new ArrayInstrumentCollection(justTheOne));
}

IInstrumentCollectionSP IInstrumentCollection::singleton(
    Instrument* instrument) {

    return singleton(InstrumentSP::attachToRef(instrument));
}

void IInstrumentCollection::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(IInstrumentCollection, clazz);
    EXTENDS(IObject);
}

CClassConstSP IInstrumentCollection::TYPE =
    CClass::registerInterfaceLoadMethod(
        "IInstrumentCollection", typeid(IInstrumentCollection), load);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
