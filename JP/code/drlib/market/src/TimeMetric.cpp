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


#include "edginc/config.hpp" 
#define QLIB_TIMEMETRIC_CPP
#include "edginc/TimeMetric.hpp" 
#include "edginc/Maths.hpp"
#include "edginc/Addin.hpp"
#include <algorithm>
 
DRLIB_BEGIN_NAMESPACE

DEFINE_TEMPLATE_TYPE(TimeMetricArray);

#define INT(X) static_cast<int>(X)
#define DBL(X) static_cast<double>(X)
#define MIDDLE_TIME_IN_DAY (INT(0.5 * DBL(DateTime::END_OF_DAY_TIME + \
                            DateTime::START_OF_DAY_TIME)))

/* for reflection */
TimeMetric::TimeMetric(): 
    CObject(TYPE)
{
    // empty
}

/** creates a time metric. Is simple if nonTradingTimeFrac = 1,
    non simple otherwise */
TimeMetric::TimeMetric(double nonTradTimeFrac,
                       const  Holiday* marketHols):
    CObject(TYPE), 
    nonTradTimeFrac(nonTradTimeFrac), marketHols(copy(marketHols)){
    validatePop2Object();
}

/** sort a date array in ascending order. 
    Duplicated dates removed.
    returns number of dates removed. */
int TimeMetric::SortDate(bool compTime,
                         DateTimeArray& dateArr)
{
    int removeCount = 0;

    sort(dateArr.begin(), dateArr.end());

    // remove duplicate
    if (dateArr.size() > 0)
    {
        vector<DateTime>::iterator iter;
        for (iter = dateArr.begin()+1; iter != dateArr.end(); iter++)
        {
            if (iter->equals(*(iter-1), compTime))
            {
                dateArr.erase(iter);
                iter --;
                removeCount ++;
            }
        }
    }
    return removeCount;
}

/** sort a date array in ascending order. 
    Holidays will be removed if nonTradTimeFrac=0. <-- this is removed in version 1.27
    Dates <= date1 or >= date2 are removed.
    Duplicated dates removed.
    returns number of dates removed. */
int TimeMetric::SortDate(const DateTime& date1, const DateTime& date2,
                         bool compTime,
                         DateTimeArray& dateArr)
{
    int removeCount = 0;

    sort(dateArr.begin(), dateArr.end());

    vector<DateTime>::iterator iter = dateArr.begin();
    while(iter != dateArr.end())
    {
        if (date1<(*iter) && (*iter)<date2)
        {
            /* if (nonTradTimeFrac == 0.0)
            {// remove holidays
                if (marketHols->isHoliday(*iter))
                {
                    dateArr.erase(iter);
                    removeCount++;
                }
                else {
                    iter ++;
                }
            }
            else {
                iter ++;
            }*/
            iter++;
    	}
        else // both ends date1 and date2 excluded
        {
            dateArr.erase(iter);
            removeCount ++;
        }
    }
    // remove duplicate
    if (dateArr.size() > 0)
    {
        for (iter = dateArr.begin()+1; iter != dateArr.end(); iter++)
        {
            if (iter->equals(*(iter-1), compTime))
            {
                dateArr.erase(iter);
                iter --;
                removeCount ++;
            }
        }
    }
    return removeCount;
}

/** This gets called after an object is constructed from a data dictionary.
    Not after an object has been copied (see override of clone method below) */
void TimeMetric::validatePop2Object()
{
   static const string method = "TimeMetric::validatePop2Object"; 

   if ((Maths::isNegative(nonTradTimeFrac)) ||
       (nonTradTimeFrac > 1.0))
   {
       throw ModelException(method,
                            "Non-trading time fraction must be "
                            "between 0 and 1");
   }

   metricIsSimple();
   // note that the algorithm means that if all days are business days
   // you get different yearFracs depending on whether 
   // nonTradTimeFrac == 1 or not
#if 0
   // allow trading time on
   if (!Maths::equals(1.0,nonTradTimeFrac))
   {
       throw ModelException(method,
                            "Non-trading time fractions different from 1.0"
                            " have been disabled");
   }
#endif
}

/** fills in the transient boolean 'simple' field */
void TimeMetric::metricIsSimple()
{
    simple = Maths::equals(nonTradTimeFrac, 1.0);
}

/** Calculates yearFracs from fromDate to each toDates.
    yearFracs.size() must be >= toDates.size() */
void TimeMetric::yearFrac(
    const DateTime&      fromDate,
    const DateTimeArray& toDates,
    DoubleArray&         yearFracs) const{
    // this needs optimising
    int numDates = toDates.size();
    if (numDates > 0){
        yearFracs[0] = yearFrac(fromDate, toDates.front());
        for (int i = 1; i < numDates; i++){
            yearFracs[i] = yearFracs[i-1]+yearFrac(toDates[i-1], toDates[i]);
        }
    }
}

/* calculates the year fraction between two date times using time metric.*/
double TimeMetric::yearFrac(
    const DateTime&   date1,          /* (I) */
    const DateTime&   date2) const    /* (I) */
{   
    double yearFracTotal = 0.0;

    // copy the input dates as we need to alter them
    DateTime dateTime1(date1);
    DateTime dateTime2(date2);

    /* Alter year frac calculation so vol is insensitive to
       exact start of day. We need to optimise this to avoid building two
       new datetimes. */
    dateTime1.moveIntoDay();
    dateTime2.moveIntoDay();
    
    if (simple) {
        // just calculate a simple year fraction
        yearFracTotal = dateTime1.yearFrac(dateTime2);
    } else {
        bool   swappedDays = false;
        if (dateTime1.isGreater(dateTime2)) {
            DateTime temp = dateTime1;
            dateTime1 = dateTime2;
            dateTime2 = temp;
            swappedDays = true;
        }
        const Holiday* hols = marketHols.operator->();
        bool dateIsBusDay = hols->isBusinessDay(dateTime1);
        // Calculate the fraction of a day at the start date 
        double fracRemaining1 = (dateIsBusDay? 1.0: nonTradTimeFrac) *
            (DBL(DateTime::END_OF_DAY_TIME - dateTime1.getTime())/
             DBL(DateTime::VOL_TIME_IN_DAY));

        /* Calculate the number of business days within the year of the
           start date and for the year in which the trade ends if not the
           same year. The business day totals for a given year are cached 
           so that the calculation is only done once */
        int dt1Year = dateTime1.getYear();
        int dt2Year = dateTime2.getYear();

        int numBusDays = 0;
        bool isLeapYearStart = false;

        hols->numBusDaysInYear(dt1Year, numBusDays, isLeapYearStart);
        double normalisedYearStart = numBusDays + nonTradTimeFrac *
            ((isLeapYearStart? 366: 365) - numBusDays);
    
        /* If the start and end dates (not times)are the same day, 
           just subtract the start and end day fractions */
        if (dateTime2.equals(dateTime1, false)) {
            // get fraction of day left after time = time2 
            double fracRemaining2 = (dateIsBusDay? 1.0: nonTradTimeFrac) *
                (DBL(DateTime::END_OF_DAY_TIME - dateTime2.getTime())/
                 DBL(DateTime::VOL_TIME_IN_DAY));
            
            // calculate in terms of days to begin with 
            yearFracTotal = fracRemaining1 - fracRemaining2;

            // then convert to years 
            yearFracTotal /= normalisedYearStart;
        } else {
            // Calculate the remaining fraction of the last trading day 
            dateIsBusDay = hols->isBusinessDay(dateTime2);
            double fracRemaining2 = (dateIsBusDay? 1.0: nonTradTimeFrac) *
                (DBL(dateTime2.getTime() - DateTime::START_OF_DAY_TIME)/
                 DBL(DateTime::VOL_TIME_IN_DAY));

            /* Calculate the full years inbetween the start and end dates 
               i.e for every full year within the date range we add 1 
               to the year fraction */
            int inBetweenYears = dt2Year - dt1Year - 1;
            if (inBetweenYears < 1) {
                // There are no 'FULL' years inside the date range 
                inBetweenYears = 0;
            }

            /* If the dates are within the same year simply find the number 
               of business days between them */
            if (dt1Year == dt2Year) {
                if (dateTime2.getDate() == dateTime1.getDate()+1){
                    numBusDays = 0;
                } else {
                    numBusDays = hols->businessDaysDiff(
                        dateTime1,
                        DateTime(dateTime2.getDate()-1, 
                                 dateTime2.getTime()));
                }
                double frac1 = numBusDays + 
                    (nonTradTimeFrac * (dateTime2.getDate()-dateTime1.getDate()
                                        - 1 - numBusDays));
                yearFracTotal = (frac1 + fracRemaining1 + fracRemaining2) /
                                 normalisedYearStart;
            } else {
                bool isLeapYearEnd = false;
                hols->numBusDaysInYear(dt2Year, numBusDays,
                                             isLeapYearEnd);
                double normalisedYearEnd = numBusDays + nonTradTimeFrac *
                    ((isLeapYearEnd? 366: 365) - numBusDays);

                DateTime endOfYearDate(DateTime::endOfYear(dt1Year));
                DateTime startOfYearDateMinus1(dt1Year == dt2Year-1?
                                               endOfYearDate:
                                               DateTime::endOfYear(dt2Year-1));
                numBusDays = hols->businessDaysDiff(dateTime1,
                                                          endOfYearDate); 

                double frac1 = numBusDays + 
                    (nonTradTimeFrac * (endOfYearDate.getDate() -
                                        dateTime1.getDate() - numBusDays));
                numBusDays = hols->businessDaysDiff(
                    startOfYearDateMinus1,
                    DateTime(dateTime2.getDate()-1, 0));
                double frac2 = numBusDays + 
                    (nonTradTimeFrac * 
                     (dateTime2.getDate() - (startOfYearDateMinus1.getDate()+1)
                      - numBusDays));
                
                yearFracTotal = ((frac1 + fracRemaining1)/normalisedYearStart) +
                    ((frac2 + fracRemaining2) / normalisedYearEnd) +
                    inBetweenYears;
            }
        }
        
        // make negative if dates are in reverse chronological order 
        if (swappedDays) {
           yearFracTotal =- yearFracTotal; 
        }
    }
    
    return yearFracTotal;
}

/** returns effective number of trading days between dates */ 
double TimeMetric::volDays(const DateTime&   date1,          /* (I) */
                           const DateTime&   date2) const   /* (I) */
{   
    double days = (double) (date2.getDate() - date1.getDate()) + 1;
    double vDays = days;

    if (!simple)
    {
        vDays = marketHols->businessDaysDiff(date1, date2);
        vDays += nonTradTimeFrac*(days - vDays);
    }
    return vDays;
}
    
/* calculates the implied dateTime from a year fraction and a date
   time.  i.e start + yearFrac as a date. The actual yearFrac from 
   the start date to the calculated date is also given as the returned 
   value- this may not be exactly the same as yearFrac. */
DateTime TimeMetric::impliedTime(
    const DateTime&  start,       /* (I) */
    double  yearFrac,             /* (I) number of years from start */
    double& actYearFrac) const    /* (O) actual year fraction */
{
    DateTime end;
    DateTime nextSoY;
    double fracToSoY = 0.0;

    
    // copy the start date since we may need to modify it
    DateTime startDate = start;
    startDate.moveIntoDay();
    
    if (Maths::isZero(yearFrac))
    {
        end = startDate;
        actYearFrac = 0.0;
    }
    else if (simple)
    {
        int  numDays = INT(yearFrac * 
                           DBL(DateTime::DAYS_PER_YEAR));
        int  numTimeUnits;
        
        // calculate actYearFrac up to start of day 
        actYearFrac = DBL(numDays)/DBL(DateTime::DAYS_PER_YEAR);
        end = DateTime(startDate.getDate()+numDays, end.getTime());
        // Approach here is to truncate
        numTimeUnits = INT(((yearFrac - actYearFrac) *
                            DBL(DateTime::VOL_TIME_IN_YEAR)));
        
        // then add on time portion to yearFrac. 
        actYearFrac += DBL(numTimeUnits)/DBL(DateTime::VOL_TIME_IN_YEAR);
        
        // add on time and correct for 'overflow' in time part 
        int date = end.getDate();
        int time = startDate.getTime() + numTimeUnits;
        if (time > DateTime::END_OF_DAY_TIME)
        {
            time -= DateTime::VOL_TIME_IN_DAY;
            date++;
        }
        else if (time < DateTime::START_OF_DAY_TIME)
        {
            time += DateTime::VOL_TIME_IN_DAY;
            date--;
        }
        end = DateTime(date, time);
    }
    else
    {
        DateTime currPos = startDate;
        DateTime lastPos;
        double   yearFracRemaining = 0.0;
        double   lastYearFrac = 0.0;
        double   yearFracSoFar = 0.0;
        double   normalisedYearFrac = 0.0;
        bool     isLeapYearStart = false;
        int      numBusDaysStartYear = 0;

        /* assuming we're going forward in time 
           move in blocks of bus days/blocks of hols - adding up trading time
           over each block until we've got more than asked for */
        if (yearFrac > 0)
        {
            nextSoY = getNextStartOfYear(currPos,
                                         SOY_FORWARD);

            fracToSoY = this->yearFrac(currPos,
                                       nextSoY);

            yearFracRemaining = yearFrac;

            // check whether we cross a year's boundary at least once #
            if (fracToSoY > yearFrac)
            {
                lastYearFrac = yearFracSoFar;
                lastPos = currPos;

                DateTime::MonthDayYear currPosYMDY = currPos.toMDY();
                
                // now there is only a fraction of one year remaining 
                marketHols->numBusDaysInYear(currPosYMDY.year, 
                                             numBusDaysStartYear, 
                                             isLeapYearStart);
                
                normalisedYearFrac = (numBusDaysStartYear +  
                                      (nonTradTimeFrac *
                                       ((isLeapYearStart? 366: 365) - numBusDaysStartYear)))
                                       * DateTime::VOL_TIME_IN_DAY;
                
                end = impliedTimeOld(currPos,
                                     yearFrac,
                                     normalisedYearFrac,
                                     actYearFrac);
            }
            else
            {
                // at least one year boundary is crossed  
                lastPos = nextSoY;
                yearFracRemaining -= fracToSoY;

                // jump over all 'in-between' years 
                actYearFrac = fracToSoY;
                while ( yearFracRemaining >= 1 )
                {
                    nextSoY = getNextStartOfYear(lastPos,
                                                 SOY_FORWARD);
                    lastPos = nextSoY;
                    yearFracRemaining -= 1.0;
                    actYearFrac++;
                }

                DateTime::MonthDayYear nextSoYMDY = nextSoY.toMDY();

                if (!(Maths::isZero(yearFracRemaining)))
                {
                    // now there is only a fraction of one year remaining 
                    marketHols->numBusDaysInYear(nextSoYMDY.year, 
                                                 numBusDaysStartYear, 
                                                 isLeapYearStart);

                    normalisedYearFrac = (numBusDaysStartYear +  
                                          (nonTradTimeFrac *
                                           ((isLeapYearStart? 366: 365) - 
                                            numBusDaysStartYear)))
                        * DateTime::VOL_TIME_IN_DAY;

                    double extraActYearFrac;
                    end = impliedTimeOld(lastPos,
                                         yearFracRemaining,
                                         normalisedYearFrac,
                                         extraActYearFrac);
                    actYearFrac += extraActYearFrac;
                }
                else
                {
                    end = lastPos;
                }
            }
        }
        else
        {
            // going backwards (yearFrac is negative) 
            
           /* if the current date is a 01/01, we want the previous year, 
              otherwise the current year - subtract one day to determine 
              the correct year    */

            nextSoY = getNextStartOfYear(currPos,
                                         SOY_BACKWARD);

            fracToSoY = this->yearFrac(currPos,
                                       nextSoY);

            yearFracRemaining = yearFrac;

            // check whether we cross a year's boundary at least once 
            if ( yearFrac >= fracToSoY )
            {
                lastPos  = currPos;
                lastYearFrac = yearFracSoFar;

                DateTime::MonthDayYear currPosYMDY = currPos.toMDY();

                if ( currPosYMDY.month == 1 &&
                     currPosYMDY.day   == 1 &&
                     currPos.getTime() == DateTime::START_OF_DAY_TIME )
                {
                    currPosYMDY.year -= 1;
                }

                // now there is only a fraction of one year remaining 
                marketHols->numBusDaysInYear(currPosYMDY.year, 
                                             numBusDaysStartYear, 
                                             isLeapYearStart);

                normalisedYearFrac = (numBusDaysStartYear +  
                                      (nonTradTimeFrac *
                                       ((isLeapYearStart? 366: 365) - 
                                        numBusDaysStartYear)))
                    * DateTime::VOL_TIME_IN_DAY;

                end = impliedTimeOld(currPos,
                                     yearFracRemaining,
                                     normalisedYearFrac,
                                     actYearFrac);
            }
            else
            {
                // at least one year boundary is crossed 
                lastPos = nextSoY;
                yearFracRemaining -= fracToSoY; 

                // jump over all 'in-between' years 
                actYearFrac = fracToSoY;
                while ( yearFracRemaining <= -1.0 )
                {
                    nextSoY = getNextStartOfYear(lastPos,
                                                SOY_BACKWARD);

                    lastPos = nextSoY;
                    yearFracRemaining += 1.0;
                    actYearFrac--;
                }

                DateTime::MonthDayYear nextSoYMDY = nextSoY.toMDY();

                if (!(Maths::isZero(yearFracRemaining)))
                {
                    // now there is only a fraction of one year remaining 
                    marketHols->numBusDaysInYear(nextSoYMDY.year - 1, 
                                                 numBusDaysStartYear, 
                                                 isLeapYearStart);

                    normalisedYearFrac = (numBusDaysStartYear +  
                                          (nonTradTimeFrac *
                                           ((isLeapYearStart? 366: 365) - numBusDaysStartYear)))
                                         * DateTime::VOL_TIME_IN_DAY;
                    double extraActYearFrac;
                    end = impliedTimeOld(lastPos,
                                         yearFracRemaining,
                                         normalisedYearFrac,
                                         extraActYearFrac);
                    actYearFrac += extraActYearFrac;
                }
                else
                {
                   end = lastPos; 
                }
            }
        }
    }

    return end;
}


/** Get the next start of year date from the current date */
    DateTime TimeMetric::getNextStartOfYear(
        const DateTime& currDate,                        /* (I) */
        const TimeMetric::TDirection direction)const     /* (I) */
{
    DateTime::MonthDayYear currDateMDY = currDate.toMDY();
    /* If we're already on a start of the year, we want the start of
       the previous year*/    
    if (direction == SOY_BACKWARD &&
        currDateMDY.month == 1   &&
        currDateMDY.day == 1     &&
        currDate.getTime() == DateTime::START_OF_DAY_TIME)
    {
        currDateMDY.year -= 1;
    }

    DateTime::MonthDayYear soyDateMDY(1, 1,
                                      currDateMDY.year + 
                                      ((direction == SOY_FORWARD)?1:0));

    DateTime soyDate = soyDateMDY.toDateTime();

    return soyDate;
}


/** Calculates the implied dateTime from a year fraction and a date time.
    The actual yearFrac from the start date to the calculated date is 
    also given - this may not be exactly the same as yearFrac. */
DateTime TimeMetric::impliedTimeOld(const   DateTime& start,     /* (I) */
                                    double  yearFrac,            /* (I) */
                                    double  normalisedYearFrac,  /* (I) */
                                    double& actYearFrac) const   /* (O) */
{
    static const string method = "TimeMetric::impliedTimeOld";

    DateTime startDate = start;
    DateTime currPos  = start;
    DateTime end;
    DateTime lastPos;
    DateTime holStart;
    DateTime holEnd;
    double   yearFracSoFar = 0.0;
    double   lastYearFrac = 0.0;
    double   timeInYear = normalisedYearFrac;
    bool     ptInBusDays = false;
    bool     foundHol = true;

    
    startDate.moveIntoDay();
    if (Maths::isZero(yearFrac)) {
        end = start;
        actYearFrac = 0.0;
    } else if(simple) {
        throw ModelException(method,
                             "Simple Trading time is not supported!");
    } else if (yearFrac > 0){
        /* assuming we're going forward in time 
           move in blocks of bus days/blocks of hols - adding up trading time
           over each block until we've got more than asked for */    

        while (yearFracSoFar < yearFrac && foundHol)
        {
            lastPos = currPos;
            lastYearFrac = yearFracSoFar;
            // find start of next block
            marketHols->findNextHol(lastPos,
                                    true,    // go fwds
                                    holStart,
                                    holEnd,
                                    foundHol);

            if (!foundHol)
            {
                ptInBusDays = true; 
            }
            else
            {
               ptInBusDays = false;
               if (currPos.getDate() < holStart.getDate())
               {
                   // add on year frac to start of hol
                   double thisYearFrac = DBL(holStart.getDate() -
                                             currPos.getDate()) *
                       DBL(DateTime::VOL_TIME_IN_DAY);

                   thisYearFrac += DBL(DateTime::START_OF_DAY_TIME -
                                       currPos.getTime());
                   
                   thisYearFrac /= timeInYear;
                   yearFracSoFar += thisYearFrac;
                   // set ptInBusDays to TRUE if have exceeded yearFrac
                   if (yearFracSoFar >= yearFrac)
                   {
                       ptInBusDays = true; 
                   }
                   else
                   {
                       // set lastPos to be at the start of hol
                       lastPos = DateTime(holStart.getDate(),
                                          DateTime::START_OF_DAY_TIME);
                       lastYearFrac = yearFracSoFar;
                   }
               }
               
               if (!ptInBusDays)
               {
                   // update currPos
                   currPos = DateTime(holEnd.getDate() + 1,
                                      DateTime::START_OF_DAY_TIME);
                   // add on year frac over hol - done in stage
                   double  thisYearFrac;
                   thisYearFrac = DBL(currPos.getDate() - lastPos.getDate()) *
                       DBL(DateTime::VOL_TIME_IN_DAY);
                   thisYearFrac += DBL(currPos.getTime() - lastPos.getTime());
                   thisYearFrac *= nonTradTimeFrac;
                   thisYearFrac /= timeInYear;
                   yearFracSoFar += thisYearFrac;
               }
            }
        }
    } else {
        // going backwards (yearFrac is negative)
        while (yearFracSoFar > yearFrac && foundHol)
        {
            lastPos = currPos;
            lastYearFrac = yearFracSoFar;
            // find start of next block
            marketHols->findNextHol(lastPos,
                                    false,    // go backwards
                                    holStart,
                                    holEnd,
                                    foundHol);

           if (!foundHol)
           {
               ptInBusDays = true; 
           } 
           else
           {
               ptInBusDays = false;
               if (currPos.getDate() > holStart.getDate())
               {
                   // add on year frac to start of hol
                   double thisYearFrac = DBL(holStart.getDate() -
                                             currPos.getDate()) *
                       DBL(DateTime::VOL_TIME_IN_DAY);
                   thisYearFrac += DBL(DateTime::END_OF_DAY_TIME -
                                       currPos.getTime());
                   thisYearFrac /= timeInYear;
                   yearFracSoFar += thisYearFrac;
                   // set ptInBusDays to TRUE if have exceeded yearFrac
                   if (yearFracSoFar <= yearFrac)
                   {
                       ptInBusDays = true;
                   }
                   else
                   {
                       // set lastPos to be at the start of hols
                       lastPos = DateTime(holStart.getDate(),
                                          DateTime::END_OF_DAY_TIME);
                       lastYearFrac = yearFracSoFar;
                   }
               }

               if (!ptInBusDays)
               {
                   // update currPos
                   double  thisYearFrac;
                   currPos = DateTime(holEnd.getDate() - 1,
                                      DateTime::END_OF_DAY_TIME);
                   // add on year frac over hol - done in stages
                   thisYearFrac = DBL(currPos.getDate() - lastPos.getDate()) *
                       DBL(DateTime::VOL_TIME_IN_DAY);
                   thisYearFrac += DBL(currPos.getTime() - lastPos.getTime());
                   thisYearFrac *= nonTradTimeFrac;
                   thisYearFrac /= timeInYear;
                   yearFracSoFar += thisYearFrac;
               }
           }
        }
    }

    /* call regular routine with lastPos, reduced yearFrac and scaling 
       factor representing weighting for day (either 1.0 for bus days
       or nonTradTimeFrac for hols */
    if (!ptInBusDays && Maths::isZero(nonTradTimeFrac))
    {
        throw ModelException(method,
                             "Internal error - point apparently lies in "
                             "region with zero time weighting");
    }

    /* recalculate the lastYearFrac, since there can be some numerical 
       errors due to iterated summation of rounded numbers*/
    lastYearFrac = this->yearFrac(startDate,
                                  lastPos);
    end = impliedTime(lastPos,
                      yearFrac - lastYearFrac,
                      ptInBusDays? 1.0: nonTradTimeFrac,
                      timeInYear,
                      actYearFrac);

    actYearFrac += lastYearFrac;

    return end;
}

/** general implied time for case where all days in a region have 
    the same weighting. dayWeight must not be 0 */  
DateTime TimeMetric::impliedTime(const   DateTime& start,   /* (I) */
                                 double  yearFrac,          /* (I) */
                                 double  dayWeight,         /* (I) */
                                 double  timeInYear,        /* (I) */
                                 double& actYearFrac)const  /* (O) */
{
    int endDate;
    int endTime;
    int timeInDay;
    // Approach here is to round to nearest
    double  halfSecond     = 0.5 / timeInYear;
    int     numTimeUnits = (int)(yearFrac * timeInYear / dayWeight +
                                 (yearFrac > 0? halfSecond: -halfSecond));
    int     numDays       = (int) 
        (DBL(numTimeUnits)/DBL(DateTime::VOL_TIME_IN_DAY));

    // calculate actYearFrac
    endDate = start.getDate() + numDays;
    timeInDay = numTimeUnits - numDays * DateTime::VOL_TIME_IN_DAY;

    actYearFrac = (DBL(numTimeUnits)/timeInYear) * dayWeight;

    //  add on time and correct for 'overflow' in time part
    endTime = start.getTime() + timeInDay;
    if (endTime > DateTime::END_OF_DAY_TIME)
    {
        endTime -= DateTime::VOL_TIME_IN_DAY;
        endDate++;
    }
    else if (endTime < DateTime::START_OF_DAY_TIME)
    {
        endTime += DateTime::VOL_TIME_IN_DAY;
        endDate--;
    }

    DateTime end(endDate, endTime);
    
    return end;
}

/** move a dateTime to the start or end of day */ 
DateTime TimeMetric::movePoint(const bool endOfDay,
                               const bool fwds,
                               DateTime&  pt)const
{
    int date = pt.getDate();
    int time;

    if (endOfDay) {
        time = DateTime::END_OF_DAY_TIME;
        if (!fwds) {
            date--;
        }
    } else {
        time = DateTime::START_OF_DAY_TIME;
        if (fwds) {
            date++;
        }
    }

    return (DateTime(date, time));
}

HolidayWrapper TimeMetric::getHolidays() const {
    return marketHols;
}

/** Populate DateTime array with n points with last point at given
    point. Pts are placed equidistant apart (in trading time)
    except that for every day on which a point falls (excluding
    first and last dates) one pt must be on the specified time of
    day (either start or end). Routine only works for case then
    lastPt is before the start point (given by pts[0]) */
void TimeMetric::calcPts(
    int             numPts,  
    DateTime&       lastPt,  
    bool            endOfDay,
    DateTimeArray&  pts) const
{
    static const string method = "TimeMetric::calcPts";
    bool goingBackwards = pts[0].isGreater(lastPt);


    if (numPts < 1)
    {
        throw ModelException(method,
                             "Need at least one point");

    }
    if (!goingBackwards)
    {
        throw ModelException(method,
                             "Forwards algorithm not yet implemented");
    }

    pts[numPts] = lastPt;
    if (numPts > 1)
    {
        double   totalLen;
        double   lenPerPt;
        int      idx;

        lastPt.moveIntoDay();
        // first of all, position dates equally
        totalLen = this->yearFrac(pts[0], lastPt);
        lenPerPt = totalLen/numPts;

        for (idx = numPts -1; idx > 0; idx--)   
        {
            double actYearFrac;
            pts[idx] = this->impliedTime(pts[idx+1],
                                         -lenPerPt,
                                         actYearFrac);
        }

        /* now adjust dates so hit end of day/start of day. Strategy:
           (A) if less than one point per day then just move each point
           to nearest boundary. This criterium may not work, if the non-
           trading-time fraction is less than 100%. In this case use the 
           algorithm which counts the number of nodes for each day. */
        
        // if (-lenPerPt > 1.0/DateTime::DAYS_PER_YEAR)
        if (simple && (numPts <= pts[0].getDate() - lastPt.getDate()))
        {
            for (idx = numPts -1; idx > 0; idx--)
            {
                this->movePoint(endOfDay,
                                // fwds or backwards
                                pts[idx].getTime() > MIDDLE_TIME_IN_DAY,
                                pts[idx]);
            }
        }
        else
        {
            /* (B): 
               1. Count number of points on each day.
               2. if eod move last point to end else move first point to end
               3. smooth out position of remaining days */
            for (idx = numPts; idx > 0;) //in inner loop 
            {
                int startOfDay;
                int timeInDay;
                int jdx = idx;

                // find how many pts lie on current date
                for (idx--; idx >= 0 && pts[idx].getDate() == pts[idx+1].getDate();
                     idx--); // empty loop 
                
                if (jdx == numPts)
                {
                    // special case of last day
                    startOfDay = lastPt.getTime();
                    timeInDay = DateTime::END_OF_DAY_TIME - lastPt.getTime();
                    if (endOfDay && jdx - idx > 1)
                    {
                        pts[idx+1] = DateTime(pts[idx+1].getDate(),
                                              DateTime::END_OF_DAY_TIME);
                        jdx--; // don't adjust lastPt 
                    }
                }
                else if (idx == -1)
                {
                    // special case of first day
                    startOfDay = DateTime::START_OF_DAY_TIME;
                    timeInDay = pts[0].getTime() - DateTime::START_OF_DAY_TIME;
                    if (!endOfDay && jdx - idx > 1)
                    {
                        pts[jdx] = DateTime(pts[jdx].getDate(),
                                            DateTime::START_OF_DAY_TIME);
                        idx++; // don't adjust startPt 
                    }
                }
                else
                {
                    startOfDay = DateTime::START_OF_DAY_TIME;
                    timeInDay = DateTime::VOL_TIME_IN_DAY; 
                    // adjust first/last point in day
                    if (endOfDay)
                    {
                        pts[idx+1] = DateTime(pts[idx+1].getDate(),
                                              DateTime::END_OF_DAY_TIME);
                    }
                    else
                    {
                        pts[jdx] = DateTime(pts[jdx].getDate(),
                                            DateTime::START_OF_DAY_TIME);
                    }
                }
                if (endOfDay)
                {
                    int kdx; 
                    // space all except last point equally
                    for (kdx = jdx; kdx > idx + 1; kdx--)
                    {
                        int time = startOfDay +
                            INT(DBL(timeInDay * (jdx - kdx + 1))/
                                DBL(jdx - idx));
                        
                        pts[kdx] = DateTime(pts[kdx].getDate(), time);
                    }
                }
                else
                {
                    int kdx;
                    // space all except first point equally
                    for (kdx = jdx -1; kdx > idx; kdx--)
                    {
                        int time = startOfDay +
                            INT(DBL(timeInDay * (jdx - kdx))/
                                DBL(jdx - idx));

                        pts[kdx] = DateTime(pts[kdx].getDate(), time);
                    }
                }
            }
        }
    }
    return;
}
                                  
/** populate from market cache */
void TimeMetric::getMarket(const IModel* model, const MarketData* market) {
    marketHols.getData(model, market);
}


class TimeMetricHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(TimeMetric, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultTimeMetric);
        FIELD(nonTradTimeFrac, "non trad time frac");
        FIELD(marketHols, "market holidays");
        FIELD(simple, "cached simple flag");
        FIELD_MAKE_TRANSIENT(simple); // hide from dd interface
    }

    static IObject* defaultTimeMetric(){
        return new TimeMetric();
    }
};

CClassConstSP const TimeMetric::TYPE = CClass::registerClassLoadMethod(
    "TimeMetric", typeid(TimeMetric), TimeMetricHelper::load);

class TimeMetricCreateAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // addin parameters
    double    nonTradTimeFrac;        // non trading time fraction
    HolidaySP marketHols;             // payment dates

    // create a Time Metric
    static IObjectSP create(TimeMetricCreateAddin *params)
        {
            static const string method = "TimeMetricCreateAddin::create";

            // do some validation on the input parameters
            if ((Maths::isNegative(params->nonTradTimeFrac)) ||
                (params->nonTradTimeFrac > 1.0))
            {
                throw ModelException(method,
                                     "Non-trading time fraction must be "
                                     "between 0 and 1");
            }
            
            TimeMetricSP newTimeMetric(new TimeMetric(params->nonTradTimeFrac,
                                                      params->marketHols.get()));

            return newTimeMetric;
        }


    TimeMetricCreateAddin(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(TimeMetricCreateAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultTimeMetricCreateAddin);
        // order of registration effects order of parameters in addin function
        FIELD(nonTradTimeFrac, "non-trading time fraction");
        FIELD(marketHols, "holidays");
        Addin::registerClassObjectMethod("TIME_METRIC",
                                         Addin::MARKET,
                                         "Creates a handle to a time metric ",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)create);
    }
    
    static IObject* defaultTimeMetricCreateAddin(){
        return new TimeMetricCreateAddin();
    }
 
};

CClassConstSP const TimeMetricCreateAddin::TYPE = CClass::registerClassLoadMethod(
    "TimeMetricCreateAddin", typeid(TimeMetricCreateAddin), load);


class TimeMetricYearFrac: public CObject{
public:
    static CClassConstSP const TYPE;

    // addin parameters
    TimeMetricSP        metric;      // the time metric
    DateTimeArray       startDate;   // start date
    DateTimeArray       endDate;     // end date

    static IObjectSP yearFracStatic(TimeMetricYearFrac* params){
        return params->yearFrac();
    }

    /** set an object in a data dictionary */
    IObjectSP yearFrac(){
        int numStart = startDate.size();
        int numEnd = endDate.size();
        if (numStart == 0 || numEnd == 0){
            throw ModelException("Must supply at least one start and one end"
                                 " date");
        }
        int numDates = Maths::max(numStart, numEnd);
        DoubleArraySP results(new DoubleArray(numDates));
        for (int i =0; i < numDates; i++){
            const DateTime& start = numStart == 1? startDate[0]: startDate[i];
            const DateTime& end = numEnd == 1? endDate[0]: endDate[i];

            (*results)[i] = metric->yearFrac(start, end);
        }
        return results;
    }

    TimeMetricYearFrac(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(TimeMetricYearFrac, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultTimeMetricYearFrac);
        FIELD(metric, "time metric");
        FIELD(startDate, "start dates");
        FIELD(endDate, "end dates");
        Addin::registerInstanceObjectMethod(
            "TIME_METRIC_YEAR_FRAC",
            Addin::MARKET,
            "Calculates the year fraction between two sets of dates",
            TYPE,
            false,
            Addin::expandSimple,
            (Addin::ObjMethod*)yearFracStatic);
    }

    static IObject* defaultTimeMetricYearFrac(){
        return new TimeMetricYearFrac();
    }
};

CClassConstSP const TimeMetricYearFrac::TYPE = CClass::registerClassLoadMethod(
    "TimeMetricYearFrac", typeid(TimeMetricYearFrac), load);

class ImpliedTimeAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // addin parameters
    TimeMetricSP        metric;      // the time metric
    DateTimeArray       startDate;   // start date
    CDoubleArray        yearFrac;    // yearFracs

    static IObjectSP impliedTimeStatic(ImpliedTimeAddin* params){
        return params->impliedTime();
    }

    /** set an object in a data dictionary */
    IObjectSP impliedTime(){
        int numStart = startDate.size();
        int numEnd = yearFrac.size();
        if (numStart == 0 || numEnd == 0){
            throw ModelException("Must supply at least one start and one year"
                                 " frac");
        }
        int numDates = Maths::max(numStart, numEnd);
        DateTimeArraySP dates(new DateTimeArray(numDates));
        DoubleArraySP   actYearFracs(new DoubleArray(numDates));
        for (int i =0; i < numDates; i++){
            const DateTime& start = numStart == 1? startDate[0]: startDate[i];
            double yf = numEnd == 1? yearFrac[0]: yearFrac[i];
            (*dates)[i] = metric->impliedTime(start, yf, (*actYearFracs)[i]);
        }
        ObjectArraySP output(new ObjectArray(2));
        (*output)[0] = dates;
        (*output)[1] = actYearFracs;
        return output;
    }

    ImpliedTimeAddin(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ImpliedTimeAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultImpliedTimeAddin);
        FIELD(metric, "time metric");
        FIELD(startDate, "Start Dates");
        FIELD(yearFrac, "Year Fractions");
        Addin::registerInstanceObjectMethod(
            "IMPLIED_TIME",
            Addin::MARKET,
            "Calculates the date which is the year fraction away "
            "from the start date using the given time metric",
            TYPE,
            false,
            Addin::expandMulti,
            (Addin::ObjMethod*)impliedTimeStatic);
    }

    static IObject* defaultImpliedTimeAddin(){
        return new ImpliedTimeAddin();
    }
};

CClassConstSP const ImpliedTimeAddin::TYPE = CClass::registerClassLoadMethod(
    "ImpliedTimeAddin", typeid(ImpliedTimeAddin), load);


DRLIB_END_NAMESPACE


