/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/timeunit.h
// Purpose:     TimeUnit enum
// Author:      ZHANG Yunzhi
// Created:     2005/04/11
// RCS-ID:      $Id: timeunit.h,v 1.1 2005/04/18 13:53:56 zhang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/timeunit.h
    @brief TimeUnit enum.
 */

#ifndef ITO33_FINANCE_TIMEUNIT_H_
#define ITO33_FINANCE_TIMEUNIT_H_

namespace ito33
{

namespace finance
{

enum TimeUnit
{
  TimeUnit_Day = 0,

  TimeUnit_Month = 1,

  TimeUnit_Year = 2
  
#ifndef __CPP2ANY__
  ,
  // @noexport
  TimeUnit_Max = 3
#endif
};


} // namespace finance

} // namespace ito33


#endif // #ifndef ITO33_FINANCE_TIMEUNIT_H_

