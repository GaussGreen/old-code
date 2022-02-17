//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MuDates.cpp
//
//   Description : Code for creation of default Mu_S and MU_POINTWISE dates
//
//   Author      : Stephen Hope
//
//   Date        : 14 Mar 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MuDates.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/Addin.hpp"
#include "edginc/ClientRunnable.hpp"

DRLIB_BEGIN_NAMESPACE

/** Public interface to generate the default dates
 */
ExpiryArraySP MuDates::getDefaultMuDates(
    const DateTime&   valueDate )
{
    static const string method = "MuDates::getDefaultMuDates";
    const int MU_QUARTERLY_BUCKETS = 8;
    const int MU_ANNUAL_BUCKETS = 9;
    const int MU_ADDITIONAL_BUCKETS = 2;
    const int MU_SERIES_SIZE = MU_QUARTERLY_BUCKETS + MU_ANNUAL_BUCKETS +
            MU_ADDITIONAL_BUCKETS;
    DateTimeArray bucketDates(MU_SERIES_SIZE);
    IntArray muAdditionalDates(MU_ADDITIONAL_BUCKETS);
    muAdditionalDates[0] = 4;
    muAdditionalDates[1] = 20;
    DateTime muDate;

    DateTime::MonthDayYear mdyValueDate = valueDate.toMDY();
    DateTime::MonthDayYear nextMdyDate = getInitialThirdFriday(mdyValueDate);

    int idx = 0;
    int i = 0;

    // Calculate quarterly bucket dates
    for (i = 0; i < MU_QUARTERLY_BUCKETS; ++i)
    {
        bucketDates[i] = nextMdyDate.toDateTime(DateTime::START_OF_DAY_TIME);
        nextMdyDate = getNextThirdFriday(nextMdyDate);
        
        ++idx;
    }
    
    const int referenceMonth = 12; // Month for annual bucketing is December

    // Calculate annual bucket dates
    for (i = 0; i < MU_ANNUAL_BUCKETS; ++i)
    {
        nextMdyDate = bucketDates[idx-1].toMDY();
        muDate = getThirdFridayInMonth(referenceMonth,
                                       (nextMdyDate.month < referenceMonth)?
                                       nextMdyDate.year:nextMdyDate.year+1);
                                                
        bucketDates[idx] = muDate;
        ++idx;
    }

    // Calculate additional bucket dates
    for (i = 0; i < MU_ADDITIONAL_BUCKETS; ++i)
    {
        nextMdyDate = bucketDates[idx-1].toMDY();

        muDate = getThirdFridayInMonth(referenceMonth,
                                       nextMdyDate.year + 
                                       muAdditionalDates[i]);

        bucketDates[idx] = muDate;
        ++idx;
    }

    // convert the bucket dates into an ExpiryArraySP
    ExpiryArraySP expiries(new ExpiryArray(MU_SERIES_SIZE));
    for (i = 0; i < MU_SERIES_SIZE; ++i)
    {
        (*expiries)[i] = ExpirySP(new BenchmarkDate(bucketDates[i]));
    }
    
    return expiries;
}

/** Given the value date, what's the first Mu date? */
DateTime::MonthDayYear MuDates::getInitialThirdFriday(const DateTime::MonthDayYear& mdyValueDate)
{

    int weekday;
    DateTime thirdFriday, valueDate;
    DateTime::MonthDayYear nextMdyDate = mdyValueDate;
    DateTime::MonthDayYear mdyThirdFriday = mdyValueDate;

    /* we're on a potential mu date month - mu dates are generated for */
    /* March, June, September and December                             */
    /* ==> if month is divisble  by 3, this month is a candidate       */
    if (mdyThirdFriday.month % 3 == 0)
    {
        mdyThirdFriday.day = 1;
        // check which calendar day the first of this month is 
        DateTime temp = mdyThirdFriday.toDateTime(DateTime::START_OF_DAY_TIME);
        weekday = temp.getWeekday();

        /* determine which date is the third friday in the value date's */
        /* months. This will be required to determine whether the value */
        /* date is already past the third friday                        */
        mdyThirdFriday.day = 21 - weekday;

        thirdFriday = mdyThirdFriday.toDateTime(DateTime::START_OF_DAY_TIME);
        valueDate = const_cast<DateTime::MonthDayYear&>(mdyValueDate).toDateTime(DateTime::START_OF_DAY_TIME);

        if (valueDate.isGreaterOrEqual(thirdFriday))
        {
            /* The value date is past the third friday. As we don't want to */
            /* generate past dates, the first mu date will be generated in  */
            /* the next month which is integer divisible ny 3               */
            nextMdyDate = getNextThirdFriday(mdyValueDate);
        }
        else
        {
            /* The next third friday which matches the criterium is in the */
            /* same month as the value date */
            nextMdyDate = mdyThirdFriday;
        }
    }
    else
    {
        /* The value date is not in a month which is integer divisible by 3 */
        /* ==> Roll forward to the next third friday in a month which is    */
        /* divisible by 3 */
        nextMdyDate = getNextThirdFriday(mdyValueDate);
    }
    
    return nextMdyDate;
}

/** Get the third friday in the next month (excluding this month) which 
    is integer divisible by 3 */
DateTime::MonthDayYear MuDates::getNextThirdFriday(const DateTime::MonthDayYear& currentMdyDate)
{
    int weekday;

    // roll month forward to the next month which is integer divisible ny 3 
    DateTime::MonthDayYear nextMdyDate = currentMdyDate;
    nextMdyDate.month = currentMdyDate.month + (3 - currentMdyDate.month % 3 );
    
    // check whether we've rolled over the end of the year
    if ( nextMdyDate.month > 12 ) 
    {
        nextMdyDate.month -= 12;
    }

    // determine in which year this date is 
    nextMdyDate.year  = 
        ( ( currentMdyDate.month + ( 3 - currentMdyDate.month % 3 ) ) > 12)?(currentMdyDate.year+1):currentMdyDate.year;

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

/* Given a month and a year, find out which day in this month is the 
   third friday */ 
DateTime MuDates::getThirdFridayInMonth(int month, int year)
{
    DateTime::MonthDayYear mdyDate(1, month, year);
    DateTime muDate;
    int weekday;

    // Check the calendar day of the first day in this month 
    DateTime temp = mdyDate.toDateTime(DateTime::START_OF_DAY_TIME);
    weekday = temp.getWeekday();
    
    // this will be the third friday of the month 
     mdyDate.day = 21 - weekday;

     muDate = mdyDate.toDateTime(DateTime::START_OF_DAY_TIME);

     return muDate;
}

// An addin to generate them - also serves as a mechanism for Pyramid to
// build the dates

class MuDateBuilder: public CObject, public ClientRunnable {
    static CClassConstSP const TYPE;
    friend class MuDates;

    // EdrAction version of addin
    IObjectSP run() {
        return generate(this);
    }

    /** addin takes two parameters - today and the type of dates */
    DateTime today;
    string   mutype;

    /** the 'addin function' - builds array of correct type */
    static IObjectSP generate(MuDateBuilder* params){
        static const string method = "MuDateBuilder::generate";
        try {
            // NB mutype parameter is ignored, but still permitted for
            // compatibility purposes.
            return MuDates::getDefaultMuDates(params->today);
        }
        catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** for reflection */
    MuDateBuilder():  CObject(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(MuDateBuilder, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultMuDateBuilder);
        FIELD(today, "today");
        FIELD(mutype, "deprecated");
        FIELD_MAKE_OPTIONAL(mutype);
        Addin::registerClassObjectMethod("MU_DATES",
                                         Addin::RISK,
                                         "Returns mu bucket dates",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)generate);
    }

    static IObject* defaultMuDateBuilder(){
        return new MuDateBuilder();
    }
    
};

CClassConstSP const MuDateBuilder::TYPE = CClass::registerClassLoadMethod(
    "MuDateBuilder", typeid(MuDateBuilder), load);


bool MuDates::load() {
    return (MuDateBuilder::TYPE && MuDateBuilder::TYPE);
}

DRLIB_END_NAMESPACE




