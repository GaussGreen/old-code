/**
 * @file ExpiryWindow.cpp
 */

#include "edginc/config.hpp"
#define QLIB_EXPIRYWINDOW_CPP
#include "edginc/Expiry.hpp"
#include "edginc/ExpiryWindow.hpp"

DRLIB_BEGIN_NAMESPACE

ExpiryWindow::ExpiryWindow(ExpiryConstSP previous, ExpiryConstSP expiry,
                           ExpiryConstSP next):
    CObject(TYPE),
    previous(previous), expiry(expiry), next(next)
{
    //    ASSERT(!!expiry);
}

ExpiryWindowSP ExpiryWindow::SP(ExpiryConstSP previous, ExpiryConstSP expiry,
                                ExpiryConstSP next) {
    return ExpiryWindowSP(new ExpiryWindow(previous, expiry, next));
}

ExpiryWindowSP ExpiryWindow::around(ExpiryArrayConstSP expiries,
                                    ExpiryConstSP expiry) {
    try {
        ASSERT(!!expiry);
        int e = expiry->search(expiries.get());
        return ExpiryWindowSP(new ExpiryWindow(
            e == 0 ? ExpirySP() : (*expiries)[e-1],
            (*expiries)[e],
            e+1 < expiries->size() ? (*expiries)[e+1] : ExpirySP()));
    }
    catch (exception& e) {
        throw ModelException(e, "ExpiryWindow::around");
    }
}

ExpiryWindowSP ExpiryWindow::find(ExpiryWindowArrayConstSP windows,
                                  ExpiryConstSP expiry) {
    ASSERT(!!expiry);

    for (int e = 0; e < windows->size(); ++e) {
        if (expiry->equals((*windows)[e]->expiry.get())) {
            return (*windows)[e];
        }
    }

    throw ModelException("ExpiryWindow::find",
                         " Expiry ("+expiry->toString()+") not found");
}

ExpiryWindowArraySP ExpiryWindow::series(ExpiryArrayConstSP exps) {
    ExpiryWindowArraySP them(new ExpiryWindowArray(exps->size()));

    for (int e = 0; e < them->size(); ++e) {
        (*them)[e] = ExpiryWindow::SP(
            e == 0 ? ExpirySP() : (*exps)[e-1],
            (*exps)[e],
            e + 1 < exps->size() ? (*exps)[e + 1] : ExpirySP());
    }

    return them;
}

ExpiryArrayConstSP ExpiryWindow::expiries(ExpiryWindowArrayConstSP windows) {
    ExpiryArraySP them(new ExpiryArray(windows->size()));

    for (int w = 0; w < them->size(); ++w)
        (*them)[w] = ExpirySP::constCast((*windows)[w]->expiry);

    return them;
}

ExpiryWindow::~ExpiryWindow() {}

string ExpiryWindow::toString() const {
    return //(!previous ? string() : "[" + previous->toString() + "-]") +
           expiry->toString() /* +
           (!next ? string() : "[-" + next->toString() + "]")*/;
}

static IObject* defaultExpiryWindow() {
    return new ExpiryWindow(ExpiryConstSP(), ExpiryConstSP(), ExpiryConstSP());
}

void ExpiryWindow::load(CClassSP& clazz) {
    REGISTER(ExpiryWindow, clazz);
    clazz->setPublic();
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultExpiryWindow);
    FIELD(previous, "previous");
    FIELD_MAKE_OPTIONAL(previous);
    FIELD(expiry, "expiry");
    FIELD(next, "next");
    FIELD_MAKE_OPTIONAL(next);
}

CClassConstSP const ExpiryWindow::TYPE = CClass::registerClassLoadMethod(
    "ExpiryWindow", typeid(ExpiryWindow), load);

DEFINE_TEMPLATE_TYPE(ExpiryWindowArray);

DRLIB_END_NAMESPACE
