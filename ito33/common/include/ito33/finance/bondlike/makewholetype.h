/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/makewholetype.h
// Purpose:     Make-Whole type
// Author:      ZHANG Yunzhi
// Created:     2004 july 5
// RCS-ID:      $Id: makewholetype.h,v 1.9 2006/01/09 11:13:25 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/makewholetype.h
    @brief Type of make whole.
 */

#ifndef _ITO33_FINANCE_BONDLIKE_MAKEWHOLETYPE_H_
#define _ITO33_FINANCE_BONDLIKE_MAKEWHOLETYPE_H_

namespace ito33
{

namespace finance
{


/// Make-Whole type
enum MakeWholeType
{
  /// Specifies that the make-whole provision is premium make-whole.
  MakeWholeType_Premium,

  /// Specifies that the make-whole provision is coupon make-whole.
  MakeWholeType_Coupon
  
  #ifndef __CPP2ANY__
  ,
  /// @noexport
  MakeWholeType_Max
  #endif
};


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_MAKEWHOLETYPE_H_
