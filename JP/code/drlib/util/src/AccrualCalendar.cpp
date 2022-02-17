//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AccrualCalendar.cpp
//
//   Description : Results data holder for Accrual Calendar and Coupon Due output
//                 requests. Stores double, datetime, datetime
//
//   Author      : Andrew McCleery
//
//   Date        : 28 June 2004
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_ACCRUALCALENDAR_CPP
#include "edginc/AccrualCalendar.hpp"

DRLIB_BEGIN_NAMESPACE

// for how many calendar days forward should data be reported
const int AccrualCalendar::ACCA_DATE_WINDOW = 7;
const int AccrualCalendar::CPND_DATE_WINDOW = 7;

AccrualCalendar::AccrualCalendar(const DateTime& calDate, const DateTime& outputDate, double outputAmount): 
    CObject(TYPE), calDate(calDate), outputDate(outputDate), outputAmount(outputAmount) {}

AccrualCalendar::AccrualCalendar() : CObject(TYPE) {}

//// Routes through equalTo. Method added to support
//// instantiating array template
bool AccrualCalendar::operator==(const AccrualCalendar& rhs) const{
    return equalTo(&rhs);
}

class AccrualCalendarHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(AccrualCalendar, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultAccrualCalendar);
        FIELD(calDate, "Calendar Date");
        FIELD(outputDate, "Output Date");
        FIELD(outputAmount, "Output Amount");
    }
    static IObject* defaultAccrualCalendar(){
        return new AccrualCalendar();
    }
};

CClassConstSP const AccrualCalendar::TYPE = CClass::registerClassLoadMethod(
    "AccrualCalendar", typeid(AccrualCalendar), AccrualCalendarHelper::load);

DEFINE_TEMPLATE_TYPE(AccrualCalendarArray);

/** specialisations of arrayObjectCast */
/** Casts array element to an IObject */
IObjectConstSP arrayObjectCast<AccrualCalendar>::toIObject(
                             const AccrualCalendar& value) {
    return IObjectConstSP::attachToRef(&value);
}

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<AccrualCalendar>::toIObject(AccrualCalendar& value) {
    return IObjectSP::attachToRef(&value);
}

/** Turns the IObjectSP into a AccrualCalendar */
AccrualCalendar arrayObjectCast<AccrualCalendar>::fromIObject(IObjectSP& value) {
    AccrualCalendar *dtPtr = dynamic_cast<AccrualCalendar*>(value.get());
    if (!dtPtr){
        throw ModelException("arrayObjectCast::fromIObject", "Object is not"
                             " an AccrualCalendar");
    }
    return *dtPtr;
}

// explicit clone for arrays of accrued data - for performance
IObject* arrayClone<AccrualCalendar>::clone(const CArray* arrayToClone){
    const AccrualCalendarArray& theArray = 
        static_cast<const AccrualCalendarArray&>(*arrayToClone);
    return new AccrualCalendarArray(theArray);
}

DRLIB_END_NAMESPACE

