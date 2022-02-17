///////////////////////////////////////////////////////////////////////////////
// File:             ito33/finance/exoticoption/onetouchutils.h                                
// Purpose:          Description of One touch using market convention                                                
// Created:          2005/12/08                                            
// RCS-ID:           $Id: onetouchutils.h,v 1.4 2006/08/19 19:11:50 wang Exp $
// Copyright         (c) 2005 -  Trilemma LLP                                        //
///////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/finance/exoticoption/onetouchutils.h
 */

#ifndef _ITO33_FINANCE_EXOTICOPTION_ONETOUCHUTILS_H_
#define _ITO33_FINANCE_EXOTICOPTION_ONETOUCHUTILS_H_

#include "ito33/sharedptr.h"
#include "ito33/date.h"

namespace ito33
{

namespace finance
{

class ITO33_DLLDECL SessionData;

/**
    Converts from market quotation to definition that is more easier to use.

    @param pSessionData The underlying of the one touch
    @param maturityDate The maturity date of the one touch
    @param dBSBarrier The value that infers the barrier of the one touch
    @param barrierType The type of the barrier (up/down).
    @param dVol The reference volatility of the one touch

    @return The barrier of the one touch or 0 if not found.
 */
double 
GetBarrierFromDelta
(const shared_ptr<SessionData>& pSessionData, 
 Date maturityDate, double dBSBarrier, BarrierType barrierType, double dVol);

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_EXOTICOPTION_ONETOUCHUTILS_H_
