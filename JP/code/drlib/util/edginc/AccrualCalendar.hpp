//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AccrualCalendar.hpp
//
//   Description : Results data holder for both Accrual Calendar and Coupon Due output
//                 requests. Stores double, datetime, datetime
//
//   Author      : Andrew McCleery
//
//   Date        : 28 June 2004
//
//
//----------------------------------------------------------------------------

#ifndef ACCRUALCALENDAR_HPP
#define ACCRUALCALENDAR_HPP

#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

class UTIL_DLL AccrualCalendar: public CObject {
public:
    static CClassConstSP const TYPE;

    // for how many calendar days forward should data be reported
    static const int ACCA_DATE_WINDOW;
    static const int CPND_DATE_WINDOW;

    AccrualCalendar(const DateTime& calDate, 
                    const DateTime& outputDate, 
                    double outputAmount);

    // AccrualCalendarArray is an array of AccrualCalendar structures, not pointers.
    // So the default constructor must be public and must specialise arrayObjectCase
    // and arrayClone
    AccrualCalendar();

    //// Routes through equalTo. Method added to support
    //// instantiating array template
    bool operator==(const AccrualCalendar& rhs) const;

private:
    friend class AccrualCalendarHelper;

    // fields
    DateTime calDate;
    DateTime outputDate;
    double outputAmount;

};

/** specialisations of arrayObjectCast */
template <> class UTIL_DLL arrayObjectCast<AccrualCalendar>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const AccrualCalendar& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(AccrualCalendar& value);

    /** Turns the IObjectSP into a BarrierLevel */
    static AccrualCalendar fromIObject(IObjectSP& value);
};

/** specialisation of arrayClone */
template <> class UTIL_DLL arrayClone<AccrualCalendar>{
public:
    /** Overridden for performance */
    static IObject* clone(const CArray* arrayToClone);
};

typedef smartPtr<AccrualCalendar> AccrualCalendarSP;
typedef smartConstPtr<AccrualCalendar> AccrualCalendarConstSP;
#ifndef QLIB_ACCRUALCALENDAR_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<AccrualCalendar>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<AccrualCalendar>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<AccrualCalendar>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<AccrualCalendar>);
#endif


typedef array<AccrualCalendar>              AccrualCalendarArray;
typedef smartPtr<AccrualCalendarArray>      AccrualCalendarArraySP;
typedef smartConstPtr<AccrualCalendarArray> AccrualCalendarArrayConstSP;
#ifndef QLIB_ACCRUALCALENDAR_CPP
EXTERN_TEMPLATE(class UTIL_DLL array<AccrualCalendar>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<AccrualCalendarArray>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<AccrualCalendarArray>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL array<AccrualCalendar>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<AccrualCalendarArray>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<AccrualCalendarArray>);
#endif

DRLIB_END_NAMESPACE

#endif
