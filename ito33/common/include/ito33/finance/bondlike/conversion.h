/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/conversion.h
// Purpose:     common terms for all conversions 
// Created:     2005/01/06
// RCS-ID:      $Id: conversion.h,v 1.3 2005/04/22 10:08:17 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/conversion.h
    @brief common terms for all conversions
    
 */

#ifndef _ITO33_FINANCE_BONDLIKE_CONVERSION_H_
#define _ITO33_FINANCE_BONDLIKE_CONVERSION_H_

#include "ito33/common.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{


/**
   Class for holding common terms of all conversions.

   @noexport COM
   @nocreate
 */
class ITO33_DLLDECL Conversion
{
public:

  /**
     Constructs an empty conversion schedule. 

     KeepAccrued is set to false and ForfeitCoupon is set to false.
   */
  Conversion() : m_bKeepAccrued(false), m_bForfeitCoupon(false)
  { }

  /// virtual desctructor for base class
  virtual ~Conversion() { }

  ///@name methods for initializing Conversion
  //@{

  /**
     The keep Accrued flag, true if security holder receives accrued interest
     upon conversion.

     @param bKeepAccrued Keep accrued flag.
   */
  void SetKeepAccrued(bool bKeepAccrued) { m_bKeepAccrued = bKeepAccrued; }
  
  /**
     The forfeit coupon flag, true if the coupon is not paid should the 
     conversion coincide with a coupon payment date.

     @param bForfeitCoupon Forfeit coupon flag.
   */
  void SetForfeitCoupon(bool bForfeitCoupon)
  {
    m_bForfeitCoupon = bForfeitCoupon;
  }

  //@}  // name methods for initializing Conversion


  ///@name methods for accessing ConversionSchedule
  //@{

  /**
     The keep Accrued flag, true if security holder receives accrued interest
     upon conversion.

     @return Keep accrued flag.
   */
  bool GetKeepAccrued() const { return m_bKeepAccrued; }
  
  /**
     The forfeit coupon flag, true if the coupon is not paid should the 
     conversion coincide with a coupon payment date.

     @return Forfeit coupon flag
   */
  bool GetForfeitCoupon() const { return m_bForfeitCoupon; }

  //@}  // name methods for accessing Conversion


protected:

  /**
     Writes myself to tag root of specific Conversion.

     @param tagRoot CallSchedule etc.
   */
  void DumpMe(XML::Tag& tagRoot) const;

  /// Keep accrued flag.
  bool m_bKeepAccrued;

  /// Forfeit coupon flag.
  bool m_bForfeitCoupon;

}; // class Conversion


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_CONVERSION_H_
