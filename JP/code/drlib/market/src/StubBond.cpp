//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StubBond.cpp
//
//   Description : Bond stub
//
//   Author      : Andrew J Swain
//
//   Date        : 1 March 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/StubBond.hpp"

DRLIB_BEGIN_NAMESPACE

StubBond::StubBond(): Stub(TYPE){
    // empty
}

StubBond::~StubBond() {
    // empty
}

/** how big is the stub payment ? */
double StubBond::payment(const DateTime&           prevCouponDate,
                         const DateTime&           nextCouponDate,
                         const DateTime&           stubStart,
                         const DateTime&           stubEnd,
                         double                    rate,
                         const DayCountConvention* dcc) const {
    static const string method = "StubBond::payment";
    try {
        double couponYearFrac = dcc->years(prevCouponDate, nextCouponDate);
        return (rate * couponYearFrac);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


string StubBond::toString() const {
    return "Bond";
}

class StubBondHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(StubBond, clazz);
        SUPERCLASS(Stub);
        EMPTY_SHELL_METHOD(defaultStubBond);
    }

    static IObject* defaultStubBond(){
        return new StubBond();
    }
};

CClassConstSP const StubBond::TYPE = CClass::registerClassLoadMethod(
    "StubBond", typeid(StubBond), StubBondHelper::load);

DRLIB_END_NAMESPACE

