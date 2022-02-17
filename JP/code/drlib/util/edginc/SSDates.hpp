//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SSDates.hpp
//
//   Description : Code for creation of default Single Stock expiry dates
//
//   Author      : Regis Guichard
//
//   Date        : 20 Nov 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_SS_DATES_HPP
#define EDG_SS_DATES_HPP

#include "edginc/DateTime.hpp"
#include "edginc/Expiry.hpp"

DRLIB_BEGIN_NAMESPACE

class UTIL_DLL SSDates{
public:
    /** Do not alter the array that is returned */
    static DateTimeArrayConstSP getDefault(const DateTime&   valueDate);

private:
    class Cache;
    friend class Cache;
    /** Given the value date, what's the first SS date? */
    static DateTime::MonthDayYear getInitialThirdFriday(const DateTime::MonthDayYear& mdyValueDate,
                                                        int                           monthDivisibility);

    /** Get the third friday in the next month (excluding this month) which 
        is integer divisible by 3 */
    static DateTime::MonthDayYear getNextThirdFriday(const DateTime::MonthDayYear& currentMdyDate,
                                                     int                           monthDivisibility);
};


DRLIB_END_NAMESPACE

#endif
