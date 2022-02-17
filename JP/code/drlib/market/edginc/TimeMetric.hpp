/*****************************************************************************
 *
 *    Group       : Equity Derivatives Research
 *
 *    Description : Defines how time is measured
 *
 *    Author      : Mark A Robson
 *
 *    Date        : 19 Jan 2001
 *
 *
 *****************************************************************************/

#ifndef _EDG_TIMEMETRIC_H
#define _EDG_TIMEMETRIC_H

#include "edginc/DateTime.hpp" 
#include "edginc/Holiday.hpp" 
#include "edginc/GetMarket.hpp"


DRLIB_BEGIN_NAMESPACE

/** Defines how time is measured - used by volatility classes */
class MARKET_DLL TimeMetric: public CObject,
                  virtual public IGetMarket {
public:
    static CClassConstSP const TYPE;
    friend class TimeMetricHelper;

    enum TDirection
    {
        SOY_FORWARD = 0,
        SOY_BACKWARD
    };

    
    /** creates a time metric. Is simple if nonTradingTimeFrac = 1,
        non simple otherwise */
    TimeMetric(double nonTradTimeFrac,
               const  Holiday* marketHols);

    /** overrides default */
    virtual void validatePop2Object();
   
    /** calculates the year fraction between two date times using time metric.
        The formula used is of the form: 
   yearFracTotal = (frac1 + fracRemaining1) + (frac2 + fracRemaining2) + inBetweenYears 
                   -----------------------    -----------------------                   
                           Yn (start)                  Yn (end)                               

        where Yn = Nb + c.Nh   is the normalised year length
        Nb is the number of days in year which are business days
        Nh is the number of days in year which are holidays
        c is the non trading time fraction
        fracRemaining accounts for fractions of a day at the start and end dates
        frac1 = Nb + c.Nh for the interval of time between the start date and end of that year
        frac2 = Nb + c.Nh for the interval of time between the end date and start of that year
        inBetweenYears are the FULL years contained within the date interval 
        (An integer) */
    double yearFrac(
        const DateTime&   dateTime1,          /* (I) */
        const DateTime&   dateTime2) const;   /* (I) */

    /** Calculates yearFracs from fromDate to each toDates.
        yearFracs.size() must be >= toDates.size() */
    void yearFrac(
        const DateTime&      fromDate,
        const DateTimeArray& toDates,
        DoubleArray&         yearFracs) const;

    /** returns effective number of trading (vol) days between [date1, date2] inclusive
        taking into account of non-trading time fraction  */
    double volDays(
        const DateTime&   dateTime1,          /* (I) */
        const DateTime&   dateTime2) const;   /* (I) */
    
    /** calculates the implied dateTime from a year fraction and a date
        time.  i.e start + yearFrac as a date. The actual yearFrac from 
        the start date to the calculated date is also given as the returned 
        value- this may not be exactly the same as yearFrac. */
   DateTime impliedTime(
        const DateTime&  start,         /* (I) */
        double  yearFrac,               /* (I) number of years from start */
        double& actYearFrac) const;     /* (O) actual year fraction */
    
    /** Populate DateTime array with n points with last point at given
        point. Pts are placed equidistant apart (in trading time)
        except that for every day on which a point falls (excluding
        first and last dates) one pt must be on the specified time of
        day (either start or end). Routine only works for case then
        lastPt is before the start point (given by pts[0]) */
    void calcPts(
        int numPts,                   /* (I) */
        DateTime& lastPt,             /* (I) */
        bool  endOfDay,               /* (I) TRUE: pts at end of day */
        DateTimeArray&  pts) const;   /* (M) pts[0] is (I) and contains
                                         start point pts[1] to pts[numPts]
                                         is populated */

    HolidayWrapper getHolidays() const;

    /** sort a date array in ascending order. 
        Duplicated dates removed.
        returns number of dates removed. */
    static int SortDate(bool compTime,
                         DateTimeArray& dateArr);

    /** sort a date array in ascending order. 
        Holidays will be removed if nonTradTimeFrac=0.
        Dates <= date1 or >= date2 are removed.
        Duplicated dates removed.
        returns number of dates removed. */
    static int SortDate(const DateTime& date1, const DateTime& date2, bool compTime,
                         DateTimeArray& dateArr);

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

private:
    /** for reflection */
    TimeMetric();
    TimeMetric(const TimeMetric &rhs);
    TimeMetric& operator=(const TimeMetric& rhs);

    /** Calculates the implied dateTime from a year fraction and a date time.
        The actual yearFrac from the start date to the calculated date is 
        also given - this may not be exactly the same as yearFrac. */
    DateTime impliedTimeOld(const   DateTime& start,     /* (I) */
                            double  yearFrac,            /* (I) */
                            double  normalisedYearFrac,  /* (I) */
                            double& actYearFrac)const;   /* (O) */

    /** general implied time for case where all days in a region have 
        the same weighting. dayWeight must not be 0 */  
    DateTime impliedTime(const   DateTime& start,   /* (I) */
                         double  yearFrac,          /* (I) num years from start */
                         double  dayWeight,         /* (I) weight to apply */
                         double  timeInYear,        /* (I) how many time units in year */
                         double& actYearFrac)const; /* (O) Actual year Frac betwen start and end */


    /** Get the next start of year date from the current date */
    DateTime getNextStartOfYear(const DateTime& currDate, /* (I) */
                                const TDirection  direction)const;    /* (I) */

    /** move a dateTime to the start or end of day */ 
    DateTime movePoint(const bool endOfDay,     /* (I) */
                       const bool fwds,         /* (I) */
                       DateTime&  pt)const;     /* (I) */
    
    /** fills in the non-registered boolean 'simple' field */
    void metricIsSimple();


    double         nonTradTimeFrac;
    HolidayWrapper marketHols;
    // transient
    bool           simple;           /* true - all days are the same */
};

typedef smartConstPtr<TimeMetric> TimeMetricConstSP;
typedef smartPtr<TimeMetric> TimeMetricSP;
typedef array<TimeMetricSP, TimeMetric> TimeMetricArray;
#ifndef QLIB_TIMEMETRIC_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<TimeMetric>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<TimeMetric>);
EXTERN_TEMPLATE(class MARKET_DLL array<TimeMetricSP _COMMA_ TimeMetric>);
EXTERN_TEMPLATE(IObjectSP MARKET_DLL FieldGetSmartPtr<TimeMetricSP>(TimeMetricSP* t));
EXTERN_TEMPLATE(void MARKET_DLL FieldSetSmartPtr<TimeMetricSP>(TimeMetricSP* t, IObjectSP o));
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<TimeMetric>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<TimeMetric>);
INSTANTIATE_TEMPLATE(class MARKET_DLL array<TimeMetricSP _COMMA_ TimeMetric>);
INSTANTIATE_TEMPLATE(IObjectSP MARKET_DLL FieldGetSmartPtr<TimeMetricSP>(TimeMetricSP* t));
INSTANTIATE_TEMPLATE(void MARKET_DLL FieldSetSmartPtr<TimeMetricSP>(TimeMetricSP* t, IObjectSP o));
#endif

DRLIB_END_NAMESPACE

#endif



