/**
 * @file ExpiryAndStrike.cpp
 */

#include "edginc/config.hpp"
#define QLIB_EXPIRYANDSTRIKE_CPP
#include "edginc/ExpiryAndStrike.hpp"

DRLIB_BEGIN_NAMESPACE

ExpiryAndStrike::ExpiryAndStrike(ExpiryConstSP expiry, double strike):
    CObject(TYPE),
    expiry(expiry),
    strike(strike)
{
    ASSERT(!!expiry);
}

ExpiryAndStrikeConstSP ExpiryAndStrike::SP(ExpiryConstSP expiry, double strike) {
    return ExpiryAndStrikeConstSP(new ExpiryAndStrike(expiry, strike));
}

ExpiryAndStrike::ExpiryAndStrike():
    CObject(TYPE),
    strike(0)
{}

ExpiryAndStrike::~ExpiryAndStrike() {}

ExpiryAndStrikeArrayConstSP ExpiryAndStrike::series(
        ExpiryArrayConstSP exps, const DoubleArray& strikes) {
    ExpiryAndStrikeArraySP them(
        new ExpiryAndStrikeArray(exps->size() * strikes.size()));

    for (int e = 0; e < exps->size(); ++e) {
        for (int s = 0; s < strikes.size(); ++s) {
            (*them)[e * strikes.size() + s].reset(
                new ExpiryAndStrike((*exps)[e], strikes[s]));
        }
    }

    return them;
}

IObject* ExpiryAndStrike::emptyShell() {
    return new ExpiryAndStrike();
}

void ExpiryAndStrike::load(CClassSP& clazz) {
    REGISTER(ExpiryAndStrike, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(emptyShell);
    FIELD(expiry, "expiry");
    FIELD(strike, "strike");
}

CClassConstSP const ExpiryAndStrike::TYPE = CClass::registerClassLoadMethod(
    "ExpiryAndStrike", typeid(ExpiryAndStrike), load);

DEFINE_TEMPLATE_TYPE(ExpiryAndStrikeArray);

DRLIB_END_NAMESPACE
