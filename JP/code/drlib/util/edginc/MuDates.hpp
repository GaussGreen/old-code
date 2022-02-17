//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MuDates.hpp
//
//   Description : Code for creation of default Mu_S and MU_POINTWISE dates
//
//   Author      : Stephen Hope
//
//   Date        : 14 Mar 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_MU_DATES_HPP
#define EDG_MU_DATES_HPP

#include "edginc/DateTime.hpp"
#include "edginc/Expiry.hpp"

DRLIB_BEGIN_NAMESPACE

class UTIL_DLL MuDates{
public:
    /** Public interface to generate the default dates
     */
    static ExpiryArraySP getDefaultMuDates(
        const DateTime&   valueDate);

    // load the symbols
    static bool load();

private:
    /* Given a month and a year, find out which day in this month is the 
       third friday */ 
    static DateTime getThirdFridayInMonth(int month, int year);

    /** Given the value date, what's the first Mu date? */
    static DateTime::MonthDayYear getInitialThirdFriday(const DateTime::MonthDayYear& mdyValueDate);

    /** Get the third friday in the next month (excluding this month) which 
        is integer divisible by 3 */
    static DateTime::MonthDayYear getNextThirdFriday(const DateTime::MonthDayYear& currentMdyDate);
};


DRLIB_END_NAMESPACE

#endif
