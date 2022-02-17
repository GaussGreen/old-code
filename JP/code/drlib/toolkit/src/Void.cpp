/**
 * @file Void.cpp
 */

#include "edginc/config.hpp"
#define QLIB_VOID_CPP
#include "edginc/Void.hpp"

DRLIB_BEGIN_NAMESPACE
Void::~Void(){}
Void::Void(): CObject(TYPE) {}

static IObject* defaultVoid() { return new Void(); }

static void loadVoid(CClassSP &clazz) {
    REGISTER(Void, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultVoid);
};

CClassConstSP const Void::TYPE = CClass::registerClassLoadMethod(
    "Void", typeid(Void), loadVoid);

DEFINE_TEMPLATE_TYPE(VoidArray);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
