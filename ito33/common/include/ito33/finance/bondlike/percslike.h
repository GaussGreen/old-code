/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/percslike.h
// Purpose:     PERCS-like financial class
// Author:      Wang
// Created:     2004/08/16
// RCS-ID:      $Id: percslike.h,v 1.22 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/percslike.h
    @brief declaration of the financial PERCS-like class   
 */

#ifndef _ITO33_FINANCE_BONDLIKE_PERCSLIKE_H_
#define _ITO33_FINANCE_BONDLIKE_PERCSLIKE_H_

#include "ito33/finance/bondlike/convertiblelike.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{

class ITO33_DLLDECL BondLikeTerms;
class ITO33_DLLDECL CallSchedule; 


/**
   The financial PERCS-like class.
 */
class ITO33_DLLDECL PERCSLike : public ConvertibleLike
{
public:
  
  /**
     Creates a PERCS-like using BondLikeTerms, cap price and 
     conversion ratio at maturity.

     @param pBondLikeTerms the terms of the PERC-like security 
							that are common to a bond.
     @param dCapPrice the maximum payoff of the PERCS.
     @param dMaxConversionRatio the max conversion ratio at maturity
   */
  PERCSLike(const shared_ptr<BondLikeTerms>& pBondLikeTerms,
            double dCapPrice,
            double dMaxConversionRatio = 1.);

  // Default copy ctor and dtor are ok

  /**
     The call schedule of the PERCS-like security.

     @param pCallSchedule shared pointer to the CallSchedule
   */
  void SetCallSchedule(const shared_ptr<CallSchedule>& pCallSchedule);

  /**
      @name Methods for accessing the PERCS-like security.
   */
  //@{

  /**
     Gets the cap price, the maximum possible payoff of the PERCS-like.

     @return the cap price
   */
  double GetCapPrice() const { return m_dCapPrice; }

  /**
     Gets the maximum conversion ratio. 

     @return the maximum conversion ratio
   */
  double GetMaxConversionRatio() const
  {
    return m_dMaxConversionRatio;
  }

  /**
     The call schedule.
    
     @return call provision
   */
  const shared_ptr<CallSchedule>& GetCallSchedule() const { return m_pCall; }

  //@}

  void Visit(DerivativeVisitor& visitor) const;

  virtual XML::Tag Dump(XML::Tag& tagParent) const;


private:
  
  /// the cap price, the maximum possible payoff of the PERCS-like security.
  double m_dCapPrice;

  /// the (max) conversion ratio
  double m_dMaxConversionRatio;

  /// call provisions
  shared_ptr<CallSchedule> m_pCall;

}; // class PERCSLike


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_PERCSLIKE_H_
