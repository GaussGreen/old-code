//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Business252.cpp
//
//   Description : Bus/252 day count convention
//
//   Author      : Xiaolan zhang
//
//   Date        : 7 November 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/Business252.hpp"
#include "edginc/Holiday.hpp" 

DRLIB_BEGIN_NAMESPACE

Business252::Business252()
    : DayCountConvention(TYPE, 252), holsSpecified(false), hols(Holiday::weekendsOnly()) {}

Business252::Business252(const string& name)
    : DayCountConvention(TYPE, 252), holsSpecified(true), hols(name) {}

Business252::Business252(const HolidayConstSP& hols)
    : DayCountConvention(TYPE, 252), holsSpecified(true), hols(hols.clone()) {}

Business252::~Business252() {}

int Business252::days(const DateTime &lodate, const DateTime &hidate) const
{
    static const string routine("Business252::days");
    try {
        if( ! hols.get() )
            throw ModelException(routine, "Holidays are not provided for Bus/252 day count convention");

        return hols->businessDaysDiff(lodate, hidate);
    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}

double Business252::years(const DateTime &lodate, const DateTime &hidate) const
{
    return double(days(lodate, hidate))/daysPerYear();
}

string Business252::toString() const
{
    return holsSpecified ?
        Format::toString("Bus/%d(%s)", daysPerYear(), hols->getName().c_str()) :
        Format::toString("Bus/%d", daysPerYear());
}
   
void Business252::setHoliday(const HolidayConstSP& hols) const
{
    holsSpecified = true;
    this->hols.setObject(HolidaySP(hols.clone()));
}

void Business252::getMarket(const IModel* model, const MarketData* market)
{
    hols.getData(model, market);
}

void Business252::validatePop2Object()
{
    holsSpecified = true;
}

class Business252Helper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(Business252, clazz);
        SUPERCLASS(DayCountConvention);
        EMPTY_SHELL_METHOD(defaultBusiness252);
        clazz->enableCloneOptimisations();

        FIELD_NO_DESC(holsSpecified);
        FIELD_MAKE_TRANSIENT(holsSpecified);
        FIELD_NO_DESC(hols);
        FIELD_MAKE_OPTIONAL(hols);
    }

    static IObject* defaultBusiness252(){
        return new Business252();
    }
};

CClassConstSP const Business252::TYPE = CClass::registerClassLoadMethod(
    "Business252", typeid(Business252), Business252Helper::load);

DRLIB_END_NAMESPACE
