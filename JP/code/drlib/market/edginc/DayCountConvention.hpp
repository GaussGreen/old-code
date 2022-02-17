//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DayCountConvention.hpp
//
//   Description : Day count convention interface
//
//   Author      : Andrew J Swain
//
//   Date        : 5 November 2000
//
//
//----------------------------------------------------------------------------

#ifndef DAYCOUNTCONVENTION_HPP
#define DAYCOUNTCONVENTION_HPP

#include "edginc/DateTime.hpp"
#include "edginc/config.hpp"
#include <string>

using namespace std;   // string

DRLIB_BEGIN_NAMESPACE
// defines an interface to be implemented by concrete classes
/** A DayCountConvention determines the year fraction between two dates */     
    
class MARKET_DLL DayCountConvention:public CObject {
public:
    static CClassConstSP const TYPE;
    friend class DayCountConventionHelper;

    // coupon types
    static const string REGULAR;
    static const string SHORT_FIRST;    // short first coupon
    static const string LONG_FIRST;     // long first coupon
    static const string SHORT_LAST;     // short last coupon
    static const string LONG_LAST;      // long last coupon
    
    virtual ~DayCountConvention();

    /** how many days in denominator */
    virtual int daysPerYear() const;

    /** how many days between two dates */
    virtual int   days(const DateTime &lodate, 
                       const DateTime &hidate) const = 0;

    /** how many days between two dates that handles end of month 
        adjusted securities */
    virtual int   days(const DateTime &lodate, 
                       const DateTime &hidate,
                       bool  eomAdjSec,
                       bool  eomIgnoreLeapYear) const;

    /** year fraction between two dates */
    virtual double years(const DateTime &lodate, 
                         const DateTime &hidate) const;

    /** year fraction between two dates for accrued purposes*/
    virtual double years(const DateTime &lodate, 
                         const DateTime &hidate,
                         bool eomAdjSec) const;

    /** returns a string description e.g. Act/360 */
    virtual string toString() const = 0;

    // returns year fraction for accrued interest calculations - 
    // regular coupons and irregular first/last coupons
    virtual double accruedFactor(const DateTime &aiDate,
                                 const DateTime &prevDate,
                                 const DateTime &nextDate,
                                 int couponFreq,
                                 const string &couponType,
                                 bool  eomAdjSec,
                                 bool  goneEx) const;

    // returns year fraction for accrued interest calculations - 
    // regular coupons only
    virtual double accruedFactor(const DateTime &aiDate,
                                 const DateTime &prevDate,
                                 const DateTime &nextDate,
                                 bool  eomAdjSec,
                                 bool  goneEx) const;

    /** Used for building instances of day count conventions from
        strings.  requiredType must be a DayCountConvention
        (or derived from a DayCountConvention) */
    static IObject* createFromString(CClassConstSP requiredType,
                                     const string& data);

protected:
    DayCountConvention(CClassConstSP clazz, int denominator);

private:
    static void load(CClassSP& clazz);

    int denominator;
};

typedef smartConstPtr<DayCountConvention>                   DayCountConventionConstSP;
typedef smartPtr<DayCountConvention>                        DayCountConventionSP;
typedef array<DayCountConventionSP, DayCountConvention>     DayCountConventionArray;
typedef smartPtr<DayCountConventionArray>                   DayCountConventionArraySP;
#ifndef QLIB_DAYCOUNTCONVENTION_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<DayCountConvention>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<DayCountConvention>);
EXTERN_TEMPLATE(class MARKET_DLL array<DayCountConventionSP _COMMA_ DayCountConvention>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<DayCountConventionArray>);
EXTERN_TEMPLATE(IObjectSP MARKET_DLL FieldGetSmartPtr<DayCountConventionSP>(
                    DayCountConventionSP* t));
EXTERN_TEMPLATE(void MARKET_DLL FieldSetSmartPtr<DayCountConventionSP>(
                    DayCountConventionSP* t, IObjectSP o));
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<DayCountConvention>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<DayCountConvention>);
INSTANTIATE_TEMPLATE(class MARKET_DLL array<DayCountConventionSP _COMMA_ DayCountConvention>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<DayCountConventionArray>);
INSTANTIATE_TEMPLATE(IObjectSP MARKET_DLL FieldGetSmartPtr<DayCountConventionSP>(
                         DayCountConventionSP* t));
INSTANTIATE_TEMPLATE(void MARKET_DLL FieldSetSmartPtr<DayCountConventionSP>(
                         DayCountConventionSP* t, IObjectSP o));
#endif


DRLIB_END_NAMESPACE
#endif
