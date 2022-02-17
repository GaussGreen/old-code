/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/barriertype.h
// Purpose:     enum for different barriers
// Created:     2005/07/05
// RCS-ID:      $Id: barriertype.h,v 1.2 2005/07/05 18:34:56 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/finance/barriertype.h
   @brief Enumeration of the type of barrier (up/down)
 */

#ifndef _ITO33_FINANCE_BARRIERTYPE_H_
#define _ITO33_FINANCE_BARRIERTYPE_H_

namespace ito33
{

namespace finance
{


/// type of the barrier: up or down
enum BarrierType
{
  Barrier_UpAndOut,
  
  Barrier_DownAndOut

  #ifndef __CPP2ANY__
  , 

  /// noexport
  Barrier_Max

  #endif

}; // enum BarrierType


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BARRIERTYPE_H_
