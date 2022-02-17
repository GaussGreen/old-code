//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StubSimple.cpp
//
//   Description : No stub
//
//   Author      : Andrew J Swain
//
//   Date        : 1 March 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/StubNone.hpp"

DRLIB_BEGIN_NAMESPACE

StubNone::StubNone(): Stub(TYPE){
    // empty
}

StubNone::~StubNone() {
    // empty
}

/** how big is the stub payment ? */
double StubNone::payment(const DateTime&           prevCouponDate,
                         const DateTime&           nextCouponDate,
                         const DateTime&           stubStart,
                         const DateTime&           stubEnd,
                         double                    rate,
                         const DayCountConvention* dcc) const {
    static const string method = "StubNone::payment";
    try {
        double couponYearFrac = dcc->years(prevCouponDate, nextCouponDate);
        return (rate * couponYearFrac);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


string StubNone::toString() const {
    return "None";
}

class StubNoneHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(StubNone, clazz);
        SUPERCLASS(Stub);
        EMPTY_SHELL_METHOD(defaultStubNone);
    }

    static IObject* defaultStubNone(){
        return new StubNone();
    }
};

CClassConstSP const StubNone::TYPE = CClass::registerClassLoadMethod(
    "StubNone", typeid(StubNone), StubNoneHelper::load);

DRLIB_END_NAMESPACE

