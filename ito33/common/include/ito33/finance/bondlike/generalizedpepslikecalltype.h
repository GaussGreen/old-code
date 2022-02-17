 /////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/generalizedpepslikecalltype.h
// Purpose:     Generalized PEPS-like call type enum
// Author:      ZHANG Yunzhi
// Created:     2005/03/12
// RCS-ID:      $Id: generalizedpepslikecalltype.h,v 1.3 2006/03/23 09:27:14 yann Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/generalizedpepslikecalltype.h
    @brief declaration of the enumeration
 */

#ifndef _ITO33_FINANCE_BONDLIKE_GENERALIZED_PEPSLIKE_CALLTYPE_H_
#define _ITO33_FINANCE_BONDLIKE_GENERALIZED_PEPSLIKE_CALLTYPE_H_


namespace ito33
{

namespace finance
{

/// Declaration of the different enum type of generalizedPEPSLikeCallType
enum GeneralizedPEPSLikeCallType
{
  GeneralizedPEPSLikeCallType_FixedShare = 0,

  GeneralizedPEPSLikeCallType_VariableShare = 1

  #ifndef __CPP2ANY__
  ,
  /// @noexport
  GeneralizedPEPSLikeCallType_Max
  #endif

};

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_GENERALIZED_PEPSLIKE_CALLTYPE_H_
