/**
 * @file BoxedBool.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/BoxedBool.hpp"

DRLIB_BEGIN_NAMESPACE

BoxedBool::BoxedBool(bool i):
    CObject(TYPE),
    b(b)
{}

BoxedBool::~BoxedBool() {}

BoxedBool* BoxedBool::create(bool b) {
    return new BoxedBool(b);
}

IObject* BoxedBool::newOne() {
    return new BoxedBool(false);
}

BoxedBoolSP BoxedBool::SP(bool b) {
    return BoxedBoolSP(new BoxedBool(b));
}

IObject* BoxedBool::clone() const {
    return getRefCount() == 0 ? new BoxedBool(b) : const_cast<BoxedBool*>(this);
}

bool BoxedBool::boolValue() const {
    return b;
}

string BoxedBool::toString() const {
    return Format::toString(b);
}

void BoxedBool::load(CClassSP& clazz) {
    REGISTER(BoxedBool, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(newOne);
    FIELD(b, "value");
}

CClassConstSP const BoxedBool::TYPE = CClass::registerClassLoadMethod(
    "BoxedBool", typeid(BoxedBool), load);

DEFINE_TEMPLATE_TYPE(BoxedBoolArray);

DRLIB_END_NAMESPACE
