//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 09-Aug-2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ICondLossDistributionsGen.hpp"

DRLIB_BEGIN_NAMESPACE

static void ICondLossDistributionsGenKeyLoad(CClassSP& clazz){
    REGISTER_INTERFACE(ICondLossDistributionsGenKey, clazz);
    EXTENDS(IObject);
}

CClassConstSP const ICondLossDistributionsGenKey::TYPE =
    CClass::registerInterfaceLoadMethod(
        "ICondLossDistributionsGenKey",
        typeid(ICondLossDistributionsGenKey),
        ICondLossDistributionsGenKeyLoad);

DEFINE_TEMPLATE_TYPE(ICondLossDistributionsGenKeyArray);

ICondLossDistributionsGen::ICondLossDistributionsGen()
{}

ICondLossDistributionsGen::~ICondLossDistributionsGen()
{}

DRLIB_END_NAMESPACE
