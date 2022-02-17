/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/convertiblebond.h
// Purpose:     convertible bond class
// Author:      ZHANG Yunzhi
// Created:     2004 may 3
// RCS-ID:      $Id: convertiblebond.h,v 1.46 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/convertiblebond.h
    @brief declaration of the financial convertible bond class.  
 */

#ifndef _ITO33_FINANCE_BONDLIKE_CONVERTIBLEBOND_H_
#define _ITO33_FINANCE_BONDLIKE_CONVERTIBLEBOND_H_

#include "ito33/finance/bondlike/cb_base.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{

class ITO33_DLLDECL BondTerms;
class ITO33_DLLDECL ConversionSchedule;
class ITO33_DLLDECL SessionData;


/**
   A Convertible bond.
 */
class ITO33_DLLDECL ConvertibleBond : public CBBase
{
public:
  /**
     Creates a convertible bond using BondTerms and ConversionSchedule.

     @param pBondTerms the principle characteristics of the bond.
     @param pConversionSchedule conversion provisions
   */
  ConvertibleBond(const shared_ptr<BondTerms>& pBondTerms,
                  const shared_ptr<ConversionSchedule>& pConversionSchedule);

  /// virtual dtor for base class
  virtual ~ConvertibleBond() { }
  

  /**
      @name Methods for accessing convertible bond specific date.
   */
  //@{

  /**
     Gets the conversions data.

     @return conversions
   */
  const shared_ptr<ConversionSchedule>& GetConversionSchedule() const
  {
    return m_pConversionSchedule;
  }

  /**
     @internal
     @biref Sets the ConversionSchedule. Added temporarily for testing stuff,
            should be removed later.

     @noexport
   */
  void SetConversionSchedule(const shared_ptr<ConversionSchedule>& pConversions)
  {
     m_pConversionSchedule = pConversions;
  }

  //@}

  virtual void Validate() const;
  virtual void ValidateWith(const SessionData& sessionData) const;

  virtual void Visit(DerivativeVisitor& visitor) const;
  virtual void Visit(DerivativeModifyingVisitor& visitor);
  virtual XML::Tag Dump(XML::Tag& tagParent) const;


protected:
  
  /// conversion schedule
  shared_ptr<ConversionSchedule> m_pConversionSchedule;

}; // class ConvertibleBond


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_CONVERTIBLEBOND_H_
