///////////////////////////////////////////////////////////////////////////////
// File:             ito33/finance/putcallparity.h  
// Purpose:          Functions supporting put-call parity calculations
// Author:           ITO 33
// Created:          2005/05/12
// RCS-ID:           $Id: putcallparity.h,v 1.2 2005/07/11 18:34:36 yann Exp $                                                    
// Copyright:        (C) 2005- Trilemma LLP
///////////////////////////////////////////////////////////////////////////////

/**
  
  @file ito33/finance/putcallparity.h 

  @brief Functions supporting put-call parity calculations

 */

#ifndef _ITO33_FINANCE_PUTCALLPARITY_H_
#define _ITO33_FINANCE_PUTCALLPARITY_H_

#include "ito33/sharedptr.h"

#include "ito33/finance/option.h"

#include "ito33/finance/theoreticalmodel.h"

namespace ito33
{

namespace finance
{


/**
  Check Put call parity

  Put = Call - Call(K = 0) + K exp( -r * T )

  @param option option contract
  @param tm Theoretical model
  @param dTolerance tolerance of put call parity check

  @return true/false
*/
bool IsPutCallParityValid(const Option &option, const TheoreticalModel& tm, 
                          double dTolerance = 1.e-2);
                        
} //namespace finance

} //namespace ito33

#endif // #ifndef _ITO33_FINANCE_PUTCALLPARITY_H_
