//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : ICDSBootstrappable.cpp
//
//   Description : Interface representing 'bootstrappable' CDSParSpreads
//
//   Date        : 18 August 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_ICDSBOOTSTRAPPABLE_CPP
#include "edginc/ICDSBootstrappable.hpp"
#include "edginc/CreditIndexBasis.hpp"

DRLIB_BEGIN_NAMESPACE


CashFlowArraySP ICDSBootstrappable::asCashFlowArray() const {
    DoubleArrayConstSP spreads = getParSpreads();
    return CashFlow::createCashFlowArray(getExpiryDates(), *(spreads.get()));
}

static void ICDSBootstrappableload(CClassSP& clazz){
    REGISTER_INTERFACE(ICDSBootstrappable, clazz);
    EXTENDS(ICDSParSpreads);
}

CClassConstSP const ICDSBootstrappable::TYPE = CClass::registerInterfaceLoadMethod(
    "ICDSBootstrappable", typeid(ICDSBootstrappable), ICDSBootstrappableload);

DEFINE_TEMPLATE_TYPE(ICDSBootstrappableWrapper);

DRLIB_END_NAMESPACE
