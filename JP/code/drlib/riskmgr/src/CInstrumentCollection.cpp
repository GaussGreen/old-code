/**
 * @file CInstrumentCollection.cpp
 */

#include "edginc/config.hpp"
#include "edginc/CInstrumentCollection.hpp"
#include "edginc/Results.hpp"

DRLIB_BEGIN_NAMESPACE

CInstrumentCollection::CInstrumentCollection(const CClassConstSP& clazz):
    CObject(clazz)
{}

CInstrumentCollection::~CInstrumentCollection() {}

CResultsArraySP CInstrumentCollection::emptyResults() const {
    CResultsArraySP them(new CResultsArray(size()));
    for (int i = 0; i < them->size(); ++i)
        (*them)[i] = CResultsSP(new CResults());
    return them;
}

void CInstrumentCollection::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CInstrumentCollection, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IInstrumentCollection);
}

CClassConstSP CInstrumentCollection::TYPE =
    CClass::registerClassLoadMethod(
        "CInstrumentCollection", typeid(CInstrumentCollection),
        CInstrumentCollection::load);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
