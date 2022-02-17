/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/timefunction.h
// Purpose:     TimeFunction class
// Created:     2006/07/05
// RCS-ID:      $Id: timefunction.h,v 1.1 2006/07/21 12:43:16 zhang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/timefunction.h
    @brief Declaration of TimeFunction class
 */

#ifndef _ITO33_FINANCE_TIMEFUNCTION_H_
#define _ITO33_FINANCE_TIMEFUNCTION_H_

#include "ito33/common.h"
#include "ito33/date.h"
#include "ito33/vector.h"

namespace ito33
{

namespace finance
{

/**
    TimeFunction class defines a function of time.

    Actually, it corresponds to a piecewise constant function of time. 
 */
class ITO33_DLLDECL TimeFunction
{
public:
  /**
      Creates TimeFunction object.

      @param dates the points that split the time into n+1 intervals, where
                   n is the size of the dates vector.
      @param values the values of the time function at each interval.
   */
  TimeFunction(const std::vector<Date>& dates,
               const std::vector<double>& values);

  /**
      Gets the time array.

      @return the time array of the time function.
   */
  const std::vector<Date>& GetDates() const
  {
    return m_pDates;
  }

  /**
      Gets the value array of the time function.

      @return the value array of the time function.
   */
  const std::vector<double>& GetValues() const
  {
    return m_pdValues;
  }

  /**
      Changes the values of the time function.

      @param values the new value arrays of the time function.
   */
  void ChangeValues(const std::vector<double>& values);

private:
  /// Date array of the time function.
  std::vector<Date> m_pDates;
  
  /// Value array of the time function.
  std::vector<double> m_pdValues;

};

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_TIMEFUNCTION_H_
