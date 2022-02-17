/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/attachedwarrantconvertiblebond.h
// Purpose:     convertible bond with share dependent conversion
// Author:      ITO 33
// Created:     2004/12/31
// RCS-ID:      $Id: attachedwarrantconvertiblebond.h,v 1.9 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/attachedwarrantconvertiblebond.h
    @brief declaration of the financial convertible bond class with
           share dependent conversion (equivalently, convertible bond
           with attached warrant)
    
 */

#ifndef _ITO33_FINANCE_BONDLIKE_ATTACHEDWARRANTCONVERTIBLEBOND_H_
#define _ITO33_FINANCE_BONDLIKE_ATTACHEDWARRANTCONVERTIBLEBOND_H_

#include "ito33/finance/bondlike/cb_base.h"

namespace ito33
{

namespace finance
{
/// forward declaration
class ITO33_DLLDECL ShareDependentConversion;

/**
   A convertible bond with share dependent conversion.
 */
class ITO33_DLLDECL AttachedWarrantConvertibleBond : public CBBase
{
public:
 
  /**
     Creates a convertible bond with attached warrant.  Equivalent to a 
     convertible bond with a share dependent conversion ratio.

     @param pBondTerms the principle characteristics of the bond
     @param pShareDepConversion share dependent conversion specification
   */
  AttachedWarrantConvertibleBond(const shared_ptr<BondTerms>& pBondTerms,
        const shared_ptr<ShareDependentConversion>& pShareDepConversion);

  /// virtual dtor for base class
  virtual ~AttachedWarrantConvertibleBond() { }

  /**
      @name Methods for accessing convertible bond specific data.
   */
  //@{

  /**
     Get the share dependent conversion specification.

     @return share dependent conversion specification
   */
  const shared_ptr<ShareDependentConversion>& 
  GetShareDependentConversion() const
  {
    return m_pShareDepConversion;
  }

  //@}

  virtual void ValidateWith(const SessionData& sessionData) const;

  virtual void Visit(DerivativeVisitor& visitor) const;
  virtual void Visit(DerivativeModifyingVisitor& visitor);
  virtual XML::Tag Dump(XML::Tag& tagParent) const;


protected:

  /// share dependent conversion specification (share factor, strike, etc)
  shared_ptr<ShareDependentConversion> m_pShareDepConversion;

}; // class AttachedWarrantConvertibleBond


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_ATTACHEDWARRANTCONVERTIBLEBOND_H_
