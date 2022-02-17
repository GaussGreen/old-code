/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/timefunction.cpp
// Purpose:     Implemenation of the TimeFunction class
// Created:     2006/03/17
// RCS-ID:      $Id: timefunction.cpp,v 1.3 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"

#include "ito33/finance/timefunction.h"
#include "ito33/finance/error.h"

extern const ito33::finance::Error
  ITO33_TIMEFUNCTION_EMPTY,
  ITO33_TIMEFUNCTION_INVALID_SIZE,
  ITO33_TIMEFUNCTION_INVALID_DATES;


namespace ito33
{

namespace finance
{

TimeFunction::TimeFunction(const std::vector<Date>& dates,
                           const std::vector<double>& values)
  : m_pDates(dates),
    m_pdValues(values)
{
  CHECK_COND( !values.empty(), ITO33_TIMEFUNCTION_EMPTY);

  CHECK_COND( m_pDates.size() == values.size(), 
             ITO33_TIMEFUNCTION_INVALID_SIZE);

  for (size_t n = 1;
       n < m_pDates.size() && m_pDates[n] > m_pDates[n - 1];
       n++)
    ;

  CHECK_COND( n == m_pDates.size(), ITO33_TIMEFUNCTION_INVALID_DATES);
}

void TimeFunction::ChangeValues(const std::vector<double>& values)
{
  CHECK_COND( m_pDates.size() == values.size(), 
              ITO33_TIMEFUNCTION_INVALID_SIZE);

  m_pdValues = values;
}

} // namespace finance

} // namespace ito33
