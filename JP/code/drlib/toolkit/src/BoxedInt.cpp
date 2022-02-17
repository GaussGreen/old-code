/**
 * @file BoxedInt.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/BoxedInt.hpp"

DRLIB_BEGIN_NAMESPACE

BoxedInt::BoxedInt(int i):
    CObject(TYPE),
    i(i)
{}

BoxedInt::~BoxedInt() {}

BoxedInt* BoxedInt::create(int i) {
    return new BoxedInt(i);
}

IObject* BoxedInt::newOne() {
    return new BoxedInt(0);
}

BoxedIntSP BoxedInt::SP(int i) {
    return BoxedIntSP(new BoxedInt(i));
}

IObject* BoxedInt::clone() const {
    return getRefCount() == 0 ? new BoxedInt(i) : const_cast<BoxedInt*>(this);
}

int BoxedInt::intValue() const {
    return i;
}

string BoxedInt::toString() const {
    return Format::toString(i);
}

void BoxedInt::load(CClassSP& clazz) {
    REGISTER(BoxedInt, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(newOne);
    FIELD(i, "value");
}

CClassConstSP const BoxedInt::TYPE = CClass::registerClassLoadMethod(
    "BoxedInt", typeid(BoxedInt), load);

DEFINE_TEMPLATE_TYPE(BoxedIntArray);

DRLIB_END_NAMESPACE
