//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StubSimple.cpp
//
//   Description : Simple stub
//
//   Author      : Andrew J Swain
//
//   Date        : 1 March 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/StubSimple.hpp"

DRLIB_BEGIN_NAMESPACE

StubSimple::StubSimple(): Stub(TYPE){
    // empty
}

StubSimple::~StubSimple() {
    // empty
}

/** how big is the stub payment ? */
double StubSimple::payment(const DateTime&           prevCouponDate,
                           const DateTime&           nextCouponDate,
                           const DateTime&           stubStart,
                           const DateTime&           stubEnd,
                           double                    rate,
                           const DayCountConvention* dcc) const {
    static const string method = "StubSimple::payment";
    try {
        double couponYearFrac = dcc->years(prevCouponDate, nextCouponDate);
        double stubFrac = dcc->years(stubStart, stubEnd)/couponYearFrac;
        return (rate * stubFrac);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


string StubSimple::toString() const {
    return "Simple";
}

class StubSimpleHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(StubSimple, clazz);
        SUPERCLASS(Stub);
        EMPTY_SHELL_METHOD(defaultStubSimple);
    }

    static IObject* defaultStubSimple(){
        return new StubSimple();
    }
};

CClassConstSP const StubSimple::TYPE = CClass::registerClassLoadMethod(
    "StubSimple", typeid(StubSimple), StubSimpleHelper::load);

DRLIB_END_NAMESPACE

