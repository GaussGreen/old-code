/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/cocotype.h
// Purpose:     type of contingent conversion
// Author:      ZHANG Yunzhi
// Created:     2004 july 5
// RCS-ID:      $Id: cocotype.h,v 1.17 2006/03/23 09:27:14 yann Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/cocotype.h
    @brief type of contingent conversion
    
 */

#ifndef _ITO33_FINANCE_BONDLIKE_COCOTYPE_H_
#define _ITO33_FINANCE_BONDLIKE_COCOTYPE_H_

namespace ito33
{

namespace finance
{


///   Different contingent conversions.
enum CoCoType
{
  /**
     Quarterly monitoring-- Should the CoCo condition be met on a check date 
     (ie end of any quarter), then the conversion is allowed as of this check 
     date.
   */
  CoCoType_CheckQuarterlyAndConvertAsOfCheckDate,

  /**
     Quarterly monitoring-- Should the CoCo condition be met on a check date 
     (ie end of any quarter), then the conversion is allowed during the 
     following quarter. 
   */
  CoCoType_CheckQuarterlyAndConvertDuringNextQuarter,

  /**
     Daily monitoring. Each check date is any date of the conversion period
     Should the CoCo condition be met on a check date, then the conversion
     is allowed as of the check date.
   */
  CoCoType_CheckAnyTimeAndConvertAsOfCheckDate,

  /**
     Daily monitoring. Each check date is any date of the conversion period-- 
     Should the CoCo condition be met on a check date, then the conversion
     is allowed on this check date.
   */
  CoCoType_CheckAnyTimeAndConvertOnCheckDate

  #ifndef __CPP2ANY__
  ,
  /// @noexport
  CoCoType_Max
  #endif

}; // enum CoCoType


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_COCOTYPE_H_
