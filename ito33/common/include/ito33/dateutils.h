/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/dateutils.h
// Purpose:     Some useful functions related to Date
// Created:     2004/05/13
// RCS-ID:      $Id: dateutils.h,v 1.9 2006/06/09 15:59:10 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/dateutils.h
    @brief Some useful functions related to Date
 */

#ifndef _ITO33_DATEUTILS_H_
#define _ITO33_DATEUTILS_H_

#include "ito33/constants.h"
#include "ito33/date.h"
#include "ito33/list.h"

namespace ito33
{

/**
   Helper function to transform a date into a double

   @param date the date to be transformed.

   @return a double corresponding to the given date
 */
inline double GetDoubleFrom(const Date& date)
{
  return date.GetExcel() * ONEDAY;
}

/**
   Helper function to transform a double time into a date(biggest but 
   still smaller than the time). 

   Note that no verification is done.

   @param dTime the time to be transformed.

   @return a date corresponding to the given date(smaller)
 */
inline Date GetDateFrom(double dTime)
{
   return Date( (unsigned long)( (dTime + TIMETOLERANCE) * INVERSEONEDAY ) );
}

/**
   Helper function to add months to a date with an adjustment for the end of
   month.

   @param date the date at which we want to add months
   @param iNbMonths number of months to add to the given date

   @return the new date
 */
inline Date AddMonthsAdjustedForEndOfMonth(Date date, int iNbMonths)
{
  bool 
    bIsEndOfMonth 
      = (    Date::GetNumOfDaysInMonth(date.GetYear(), date.GetMonth())
          == date.GetDay() );

  date.AddMonths(iNbMonths);

  if (bIsEndOfMonth)
    return Date( date.GetYear(), date.GetMonth(), 
                 Date::GetNumOfDaysInMonth(date.GetYear(), date.GetMonth()) );

  return date; 
}

/**
   Helper function to generate a serie of regular dates from a first date.

   @param startDate the date start at which we generate
   @param endDate the end date(will be included if it is a regular date) 
   @param iInterval number of month between two succesive regular dates
   @param bIsEOM flag for special case when start date is at end of month

   @return a vector of generated dates, first date not included
 */
inline std::list<Date> 
GenerateRegularDates
(Date startDate, Date endDate, int iInterval, bool bIsEOM)
{
  std::list<Date> dates;

  Date date = startDate;
  int i = 0;
  if  ( bIsEOM )
    while ( date <= endDate )
    {
      i += iInterval;
      date = AddMonthsAdjustedForEndOfMonth(startDate, i);
      
      if ( date <= endDate)
        dates.push_back(date);
    }
  else // Non end of month case
  {
    while ( date <= endDate )
    {
      date = startDate;
      i += iInterval;
      
      date.AddMonths(i);

      if ( date <= endDate )
        dates.push_back(date);
    }
  }

  return dates; 
}

} // namespace ito33

#endif // #ifndef _ITO33_DATEUTILS_H_
