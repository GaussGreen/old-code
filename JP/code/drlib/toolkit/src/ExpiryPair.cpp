/**
 * @file ExpiryPair.cpp
 */

#include "edginc/config.hpp"
#define QLIB_EXPIRYPAIR_CPP
#include "edginc/Expiry.hpp"
#include "edginc/ExpiryPair.hpp"

DRLIB_BEGIN_NAMESPACE

ExpiryPair::ExpiryPair(ExpiryConstSP optExpiry, ExpiryConstSP ulExpiry):
    CObject(TYPE),
    optExpiry(optExpiry), ulExpiry(ulExpiry)
{}

ExpiryPairSP ExpiryPair::SP(ExpiryConstSP optExpiry, ExpiryConstSP ulExpiry) {
    return ExpiryPairSP(new ExpiryPair(optExpiry, ulExpiry));
}

bool ExpiryPair::equals(const ExpiryPair* expiryPair) const {
  return optExpiry->equals(expiryPair->optExpiry.get()) && ulExpiry->equals(expiryPair->ulExpiry.get()); 
}

ExpiryPair::~ExpiryPair() {}

string ExpiryPair::toString() const {
    return optExpiry->toString() + "_" + ulExpiry->toString();
}

static IObject* defaultExpiryPair() {
    return new ExpiryPair(ExpiryConstSP(), ExpiryConstSP());
}

void ExpiryPair::load(CClassSP& clazz) {
    REGISTER(ExpiryPair, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultExpiryPair);
    FIELD(optExpiry, "optional expiry");
    FIELD(ulExpiry, "underlying expiry");
}

CClassConstSP const ExpiryPair::TYPE = CClass::registerClassLoadMethod(
    "ExpiryPair", typeid(ExpiryPair), load);

DEFINE_TEMPLATE_TYPE(ExpiryPairArray);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
