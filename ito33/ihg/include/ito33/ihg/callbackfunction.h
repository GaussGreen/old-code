/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/callbackfunction.h
// Purpose:     common callback function for vol and hazard rate
// Author:      Wang
// Created:     2004/03/17
// RCS-ID:      $Id: callbackfunction.h,v 1.6 2004/10/04 18:04:05 pedro Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/ihg/callbackfunction.h
   @brief definition of common callback function for vol and hazard rate 

   @todo the parameter iUserData is sufficient for now but maybe not enough
         for general use.
 */
  
#ifndef _ITO33_IHG_CALLBACKFUNTION_H_
#define _ITO33_IHG_CALLBACKFUNTION_H_

#include "ito33/cwrap.h"    // for ITO33CALLCONV

namespace ito33
{

namespace ihg
{

/**
    Prototype for the callback function used by VolatilityCallBack and
    HazardRateCallBack classes.

    @param dTime the time for which pdValues should be filled
    @param pdStocks the "spot" values 
    @param pdValues the output parameter
    @param nNbPoints the number of points in pdStocks
    @param iUserData client data
 */
typedef
void (ITO33CALLCONV *CallBackFunction)(double dTime, 
                                       const double *pdStocks,
                                       double *pdValues,
                                       size_t nNbPoints,   
                                       int iUserData);

} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_CALLBACKFUNTION_H_
