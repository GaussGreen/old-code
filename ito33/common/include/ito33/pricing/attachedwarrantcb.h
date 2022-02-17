/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/attachedwarrantcb.h
// Purpose:     contract class for attached warrant convertible bond
// Author:      Ito33
// Created:     2005/01/13
// RCS-ID:      $Id: attachedwarrantcb.h,v 1.7 2006/08/19 22:01:27 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/attachedwarrantcb.h
    @brief The declaration of the  contract class.

    The base class for attached warrant convertible bond contracts.  
 */

#ifndef _ITO33_PRICING_ATTACHEDWARRANTCB_H_
#define _ITO33_PRICING_ATTACHEDWARRANTCB_H_

#include "ito33/pricing/cb.h"
#include "ito33/pricing/sharedependentconversion.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL AttachedWarrantConvertibleBond;
}

namespace pricing
{
// Note:
// cb warrant is derived from cb however, there are
// currently two conversions object inside AttachedWarrantConvertibleBond
// there is the classic cbconversion and the sharedependentconversion.
// Both are needed since the conversion --object-- is stored in the class
// object is stored to allow copy for path dependent and call notice.
// so to ensure that the right conversion is called inside this
// class the function getconversion is made virtual in cb.h 
// and a new setconversion with the right type information is provided.

/// The declaration of the (backward) warrant contract class.
class AttachedWarrantConvertibleBond : public CB 
{

public:
  /**
     Creates a warrant by financial sharedependent object.

     @param warrant financial warrant
   */
  AttachedWarrantConvertibleBond
  (const finance::AttachedWarrantConvertibleBond& warrant);

  /// virtual dtor for base class
  virtual ~AttachedWarrantConvertibleBond() { }
 
  /**
      Gets the reset time at which the conversion ratio is fixed.

      @return the reset time
   */
  double GetResetTime() const 
  { 
    return m_conversionsShareDependent.GetResetTime();
  }

  /**
      Indicates whether or not a reset time has been set.

      @return true/false if a reset time has been set/not set
   */
  bool HasResetTime() const
  {
    return m_conversionsShareDependent.HasResetTime();
  }

  /**
      Gets the base conversion ratio.

      @return The base conversion ratio
   */
  double GetBaseRatio() const 
  { 
    return m_conversionsShareDependent.GetBaseRatio(); 
  }

  /**
      Gets the maximum conversion ratio.

      @return the maximum conversion ratio
   */
  double GetCapRatio() const 
  { 
    return m_conversionsShareDependent.GetCapRatio(); 
  }
  
  /**
      Gets the current conversion ratio.

      @return the current conversion ratio.
   */
  double GetCurrentConversionRatio() const
  {
    return m_dCurrentConversionRatio;
  }
  
  /**
      Get Conversions.
   */
  ConversionProvisions* GetConversions() 
  { 
    return &m_conversionsShareDependent; 
  }
  
  /**
      Sets the conversions.
   */
  void SetConversions(const ShareDependentConversion& conversionsShD) 
  { 
    m_conversionsShareDependent = conversionsShD;
  }

protected:
 
  ///sharedependentconversion
  ShareDependentConversion m_conversionsShareDependent;

  ///current conversion ratio
  double m_dCurrentConversionRatio;

}; // class AttachedWarrantConvertibleBond;


} // namespace pricing

} // namespace ito33

#endif // _ITO33_PRICING_ATTACHEDWARRANTCB_H_
