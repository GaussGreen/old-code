/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/predicatetime.h
// Purpose:     Predicate function for time
// Author:      Nabil
// Created:     2004/04/21
// RCS-ID:      $Id: predicatetime.h,v 1.5 2004/10/05 09:13:38 pedro Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_NUMERIC_PREDICATETIME_H_
#define _ITO33_NUMERIC_PREDICATETIME_H_

#include "ito33/beforestd.h"
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/constants.h"

namespace ito33
{

namespace numeric
{


inline bool AreTimesEqual(double dVal1, double dVal2)
{
  return fabs(dVal1 - dVal2) < TIMETOLERANCE;
}

inline bool IsBefore(double dTime1, double dTime2)
{
  return dTime1 < dTime2 - TIMETOLERANCE;
}

inline bool IsAfter(double dTime1, double dTime2)
{
  return dTime1 > dTime2 + TIMETOLERANCE;
}

inline bool IsEqualOrBefore(double dTime1, double dTime2)
{
  return dTime1 < dTime2 + TIMETOLERANCE;
}

inline bool IsEqualOrAfter(double dTime1, double dTime2)
{
  return dTime1 > dTime2 - TIMETOLERANCE;
}

} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_PREDICATETIME_H_
