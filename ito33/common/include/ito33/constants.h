/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/constants.h
// Purpose:     some constants
// Author:      ITO 33
// Created:     12.02.04
// RCS-ID:      $Id: constants.h,v 1.17 2006/06/27 15:38:08 yann Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/constants.h
    @brief constant definitions
 */

#ifndef _ITO33_CONSTANTS_H_
#define _ITO33_CONSTANTS_H_
  
#include <cstddef>

namespace ito33
{

/// The max value of type size_t. It is used to check the validity an index
const size_t INVALIDINDEX = size_t(-1);

/// @name time concerned
//@{

/// float value of one day
const double ONEDAY = 1. / 365.; 

/// number of days in one year, as an int
const int INVERSEONEDAY = 365;

/**
   Tolerance used to compare two times. 
   
   Typical time is around one hundred. 
   
   This number should take care of the machine precision, so that round error 
   can be safely ignored. 

   The minimum time step should be limited to a number bigger than this. We 
   don't yet explicitly require it. So, don't refine too much the time mesh. 

   For the moment, we choose a more or less reasonable number.
 */
const double TIMETOLERANCE = 1.e-10;

//@}

/**
   Tolerance used to comprate two doubles
   For the moment, we choose a more or less reasonable number.
 */
const double DOUBLETOLERANCE = 1.e-12;

/**
   This value can never be reached when pricing any contracts. 
   This variables is currently used to 
   avoid outputing the market price of the derivative.
   Since by default it is not set.
 */
const double INVALIDPRICE = 1.e99;

/// Shift
const double SHIFT = 1.e-6; 


} // namespace ito33

#endif // #ifndef _ITO33_CONSTANTS_H_
