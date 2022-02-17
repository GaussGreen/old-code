/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/rebatetype.h
// Purpose:     enum for different rebates
// Created:     2005/07/05
// RCS-ID:      $Id: rebatetype.h,v 1.1 2005/07/05 10:02:21 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/finance/rebatetype.h
   @brief Enumeration of the type of rebate (immediate/at maturity)
 */

#ifndef _ITO33_FINANCE_REBATETYPE_H_
#define _ITO33_FINANCE_REBATETYPE_H_

namespace ito33
{

namespace finance
{


/// type of the rebate: immediate or at maturity
enum RebateType
{
  Rebate_Immediate,
  
  Rebate_AtMaturity

  #ifndef __CPP2ANY__
  , 

  /// noexport
  Rebate_Max

  #endif

}; // enum RebateType


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_REBATETYPE_H_
