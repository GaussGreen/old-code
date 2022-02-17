#include "edginc/config.hpp"
#include "edginc/IntegratorND.hpp"

DRLIB_BEGIN_NAMESPACE

static void load(CClassSP& clazz) {
    REGISTER_INTERFACE(IntegratorND, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IntegratorND::TYPE = 
    CClass::registerInterfaceLoadMethod(
        "IntegratorND", 
        typeid(IntegratorND), 
        load);

DEFINE_TEMPLATE_TYPE(IntegratorNDArray);

bool IntegratorNDLoad() {
    return IntegratorND::TYPE != 0;
}

DRLIB_END_NAMESPACE
