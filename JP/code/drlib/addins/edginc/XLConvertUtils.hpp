//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLConvertUtils.hpp
//
//   Description : Provides some utility methods for dealing with 
//                 EDR specific classes (formerly part of XLConvert)
//
//   Author      : Mark A Robson
//
//   Date        : 25 Nov 2003
//
//
//----------------------------------------------------------------------------

#ifndef EDG_XLCONVERT_UTILS_HPP
#define EDG_XLCONVERT_UTILS_HPP

#include "edginc/XLConvert.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Expiry.hpp"
DRLIB_BEGIN_NAMESPACE


/** Provides some utility methods for dealing with EDR specific classes. 
    It's a bit of a lazy approach as classes that deal with EDR specific data
    (eg dates, maturities) derive from this class. This is just so that they
    can see the methods in their namespace (and call protected methods) */
class ADDINS_DLL XLConvertUtils: public XLConvert{
public:
    /** Converts an oper to an integer date */
    static int operToDate(const XL_OPER& oper);

    /** Converts an oper to a 'time of day' */
    static int operToTime(const XL_OPER& oper);

    /** Converts an oper to a DateTime */
    static DateTime opersToDateTime(const XL_OPER& dateOper,
                                    const XL_OPER& timeOper);

    /** construct the xl representation of the date part of a DateTime */
    static void xlDateFromDate(int       date,
                               XL_OPER&  oper);

    /** construct the xl representation of the time part of a DateTime */
    static void xlTimeFromTime(int       time,  /* (I) */
                               XL_OPER&  oper); /* (O) xlType is set too */

    /** construct the xl representation of a DateTime */
    static void xlDateTimeFromDateTime(const DateTime& dateTime,
                                       XL_OPER&        dateOper,
                                       XL_OPER&        timeOper);
protected:
    XLConvertUtils();

    /** Constructs expiry object from supplied opers */
    static ExpirySP opersToExpiry(const XL_OPER& oper1,
                                  const XL_OPER& oper2);

    /** populate the given opers using the supplied expiryObj. Will fail
        if expiryObj is not a recognised expiry */
    static void opersFromExpiry(const IObjectConstSP& expiryObj,
                                XL_OPER&              oper1,
                                XL_OPER&              oper2);

    static bool operIsMaturity(const XL_OPER& oper);


private:
    static int fracOfDayToTime(double  timeFrac);  /* (I) */
};

DRLIB_END_NAMESPACE
#endif
