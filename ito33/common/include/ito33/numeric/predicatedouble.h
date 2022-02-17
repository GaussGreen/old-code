/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/predicatedouble.h
// Purpose:     Predicate function for double
// Author:      Nabil
// Created:     2004/03/22
// RCS-ID:      $Id: predicatedouble.h,v 1.8 2005/10/31 20:01:17 dave Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_NUMERIC_PREDICATEDOUBLE_H_
#define _ITO33_NUMERIC_PREDICATEDOUBLE_H_

#include "ito33/beforestd.h"
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/constants.h"

namespace ito33
{

namespace numeric
{


inline bool IsEqual(double dVal1, double dVal2)
{
  double
    dMax = fabs(dVal1) > fabs(dVal2) ? fabs(dVal1) : fabs(dVal2);
  bool
    bResult = fabs(dVal1 - dVal2) <= DOUBLETOLERANCE * dMax;
  return bResult;
}

/**
   Helper function to compare spot with a trigger (call or conversion).
 */
inline bool LessThanTrigger(double dSpot, double dTrigger)
{
  return dSpot < dTrigger * (1. - DOUBLETOLERANCE);
}

/**
  Helper function to compare two double
*/
inline bool IsEqualOrLess(double dVal1, double dVal2)
{
  return dVal1 < dVal2 + DOUBLETOLERANCE;
}

/**
  Helper function to compare two double
*/
inline bool IsEqualOrGreater(double dVal1, double dVal2)
{
  return dVal1 > dVal2 - DOUBLETOLERANCE;
}

} //namespace numeric

} //namespace ito33

#endif //_ITO33_NUMERIC_PREDICATEDOUBLE_H_
