//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SSDates.cpp
//
//   Description : Code for creation of default Single Stock expiry dates
//
//   Author      : Regis Guichard
//
//   Date        : 20 Nov 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SSDates.hpp"
#include ext_hash_map

DRLIB_BEGIN_NAMESPACE
class SSDates::Cache{
private:
    /* class needed by hash_map template */
    class HashUtil{
    public:
        bool operator()(const DateTime& date1, 
                        const DateTime& date2) const{
            return (date1.equals(date2));
        }
        size_t operator()(const DateTime& date) const{
            return date.hashCode();
        }
    };
    static hash_map<DateTime, DateTimeArraySP, HashUtil, HashUtil> cache;
    static DateTimeArraySP generateDefault(const DateTime&   valueDate){
        static const string method = "SSDates::Cache::generateDefault";

        static const int nbFreq = 4;
        static const int frequencies[] = 
            {1, 3, 6, 12};   // monthly, quarterly, six-monthly, yearly
        static const int nbDates[] =
            {36, 20, 14, 30};    /* 36 monthly (= 3Y), 
                                    20 quarterly (= 5Y), 
                                    14 six-monthly (= 7Y), 30 yearly (= 30Y) */
        static const int nbWeekly = 4;
        
        try{
            /* Do each frequency in turn and successively merge dates */
            DateTimeArray bucketDates(0);
            DateTime::MonthDayYear mdyValueDate = valueDate.toMDY();
            int iFreq = 0;
            for (; iFreq < nbFreq; ++iFreq){
                DateTimeArray sameFreqDates(nbDates[iFreq]);
                /* For each frequency create dates */
                DateTime::MonthDayYear nextMdyDate = 
                    SSDates::getInitialThirdFriday(mdyValueDate,
                                                   frequencies[iFreq]);
                int iDate = 0;
                for (; iDate < nbDates[iFreq] - 1; ++iDate){
                    sameFreqDates[iDate] = 
                        nextMdyDate.toDateTime(DateTime::END_OF_DAY_TIME); // PM
                    nextMdyDate = 
                        SSDates::getNextThirdFriday(nextMdyDate, 
                                                    frequencies[iFreq]);
                }
                sameFreqDates[iDate] = 
                    nextMdyDate.toDateTime(DateTime::END_OF_DAY_TIME);   // PM
                /* merge */
                bucketDates = DateTime::merge(sameFreqDates, bucketDates);
            }

            /* Add nbWeekly dates (Fridays) to the front */
            DateTime nextFriday(bucketDates[0]);
            do{
                nextFriday = nextFriday.rollDate(-7);
            }
            while(valueDate < nextFriday);
            DateTimeArray weeklyDates(nbWeekly);
            int iWeekly = 0;
            for (; iWeekly < nbWeekly; ++iWeekly){
                nextFriday = nextFriday.rollDate(7);
                weeklyDates[iWeekly] = nextFriday;
            }
            bucketDates = DateTime::merge(weeklyDates, bucketDates);
            return DateTimeArraySP(new DateTimeArray(bucketDates));
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }
public:
    static DateTimeArrayConstSP getDefault(const DateTime&   valueDate){
        DateTimeArraySP& dates = cache[valueDate];
        if (!dates){
            dates = DateTimeArraySP(generateDefault(valueDate));
        }
        return dates;
    }
};
hash_map<DateTime, DateTimeArraySP, SSDates::Cache::HashUtil,
         SSDates::Cache::HashUtil> SSDates::Cache::cache;

/** Do not alter the dates in the array that is returned */
DateTimeArrayConstSP SSDates::getDefault(const DateTime&   valueDate){
    return Cache::getDefault(valueDate);
}

/** Get the third friday in the next month (possibly including this month) which 
    is integer divisible by monthDivisibility */
DateTime::MonthDayYear SSDates::getInitialThirdFriday(const DateTime::MonthDayYear& mdyValueDate,
                                                      int                           monthDivisibility)
{
    int weekday;
    DateTime thirdFriday, valueDate;
    DateTime::MonthDayYear nextMdyDate = mdyValueDate;
    DateTime::MonthDayYear mdyThirdFriday = mdyValueDate;

    /* If month is divisble  by monthDivisibility, this month is a candidate */
    if (mdyThirdFriday.month % monthDivisibility == 0)
    {
        mdyThirdFriday.day = 1;
        /* check which calendar day the first of this month is */
        DateTime temp = mdyThirdFriday.toDateTime(DateTime::START_OF_DAY_TIME);
        weekday = temp.getWeekday();

        /* determine which date is the third friday in the value date's
           months. This will be required to determine whether the value
           date is already past the third friday                        */
        mdyThirdFriday.day = 21 - weekday;

        thirdFriday = mdyThirdFriday.toDateTime(DateTime::START_OF_DAY_TIME);
        valueDate = const_cast<DateTime::MonthDayYear&>(mdyValueDate).toDateTime(DateTime::START_OF_DAY_TIME);

        if (valueDate.isGreaterOrEqual(thirdFriday))
        {
            /* The value date is past the third friday. As we don't want to
               generate past dates, the first mu date will be generated in
               the next month which is integer divisible ny monthDivisibility */
            nextMdyDate = getNextThirdFriday(mdyValueDate, monthDivisibility);
        }
        else
        {
            /* The next third friday which matches the criterium is in the
               same month as the value date */
            nextMdyDate = mdyThirdFriday;
        }
    }
    else
    {
        /* The value date is not in a month which is integer divisible by monthDivisibility
           ==> Roll forward to the next third friday in a month which is
           divisible by monthDivisibility */
        nextMdyDate = getNextThirdFriday(mdyValueDate, monthDivisibility);
    }
    
    return nextMdyDate;
}

/** Get the third friday in the next month (excluding this month) which 
    is integer divisible by monthDivisibility */
DateTime::MonthDayYear SSDates::getNextThirdFriday(const DateTime::MonthDayYear& currentMdyDate,
                                                   int                           monthDivisibility) // 3
{
    int weekday;

    // roll month forward to the next month which is integer divisible ny monthDivisibility 
    DateTime::MonthDayYear nextMdyDate = currentMdyDate;
    nextMdyDate.month = currentMdyDate.month + (monthDivisibility - currentMdyDate.month % monthDivisibility );
    
    // check whether we've rolled over the end of the year
    if ( nextMdyDate.month > 12 ){
        nextMdyDate.month -= 12;
        nextMdyDate.year = currentMdyDate.year + 1;
    }
    else{
        nextMdyDate.year = currentMdyDate.year;
    }

    // determine the calendar day of the first day in the month 
    nextMdyDate.day   = 1;

    // check whether today is after the third friday of this month
    DateTime temp = nextMdyDate.toDateTime(DateTime::START_OF_DAY_TIME);
    weekday = temp.getWeekday();

    /* From the calendar day of the first day of the month, determine the */
    /* date of the third friday in this month */
    nextMdyDate.day = 21 - weekday;

    return nextMdyDate;
}

DRLIB_END_NAMESPACE




