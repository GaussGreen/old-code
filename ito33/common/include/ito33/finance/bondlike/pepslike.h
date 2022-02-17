 /////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/pepslike.h
// Purpose:     PEPS-like financial class
// Author:      Wang
// Created:     2004/08/17
// RCS-ID:      $Id: pepslike.h,v 1.27 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/pepslike.h
    @brief declaration of the financial PEPS-like class   
 */

#ifndef _ITO33_FINANCE_BONDLIKE_PEPSLIKE_H_
#define _ITO33_FINANCE_BONDLIKE_PEPSLIKE_H_

#include "ito33/finance/bondlike/convertiblelike.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{

class ITO33_DLLDECL BondLikeTerms;
class ITO33_DLLDECL CallSchedule; 
class ITO33_DLLDECL CallFixedShare;

/**
   The financial PEPS-like class.
 */
class ITO33_DLLDECL PEPSLike : public ConvertibleLike
{
public:
  
  /**
     Creates a PEPS-like by BondLikeTerms and the max and min conversion
     ratios at maturity.

     @param pBondLikeTerms the terms of the PEPS that are common 
                           to a the bond terms.
     @param dMaxConversionRatio the maximum conversion ratio at maturity 
     @param dMinConversionRatio the minimum conversion ratio at maturity
   */
  PEPSLike(const shared_ptr<BondLikeTerms>& pBondLikeTerms,
           double dMaxConversionRatio,
           double dMinConversionRatio);

  // Default copy ctor and dtor are ok

  /*
    For the moment, only fixed share conversion is implemented.  Support may 
    be added later for variable share optional conversion.
  */
  
  /**
     The security has optional conversion.
   */
  void EnableOptionalConversion()
  {
    m_bHasOptionalConversion = true;
  }

  /**
     The fixed cash call of PEPS-like instruments.

     @param pCallSchedule shared pointer to the fixed cash call.
     
     @method
   */
  void SetCallFixedCash(const shared_ptr<CallSchedule>& pCallSchedule);

  /**
     The fixed share call provision of PEPS-like instruments.

     Note that this is indeed the optional conversion at issuer's option.

     @param pCallFixedShare shared pointer to the fixed share call

     @method
   */
  void SetCallFixedShare(const shared_ptr<CallFixedShare>& pCallFixedShare);

  /**
      @name Methods for accessing the PEPS like instrument.
   */
  //@{

  /**
     Gets the maximum conversion ratio at maturity.

     @return the maximum conversion ratio at maturity 
   */
  double GetMaxConversionRatio() const
  {
    return m_dMaxConversionRatio;
  }

  /**
     Gets the minimum conversion ratio at maturity.

     @return the minimum conversion ratio at maturity 
   */
  double GetMinConversionRatio() const
  {
    return m_dMinConversionRatio;
  }

  /**
     Whether the PEPS-like security has optional conversion (fixed share).

     @return true if optional conversion at holder's option exists.
   */
  bool HasOptionalConversion() const
  {
    return m_bHasOptionalConversion;
  }

  /**
     The fixed cash call of PEPS-like instruments.
    
     @return call provision
   */
  const shared_ptr<CallSchedule>& GetCallFixedCash() const
  {
    return m_pCallSchedule; 
  }

  /**
     The fixed share call provision of PEPS-like instruments.
     
     Note that this is indeed the optional conversion at issuer's option.

     @return fixed share call provision
   */
  const shared_ptr<CallFixedShare>& GetCallFixedShare() const
  {
    return m_pCallFixedShare;
  }

  //@}

  void Visit(DerivativeVisitor& visitor) const;

  virtual XML::Tag Dump(XML::Tag& tagParent) const;


private:
  
  /// the maximum conversion ratio at maturity
  double m_dMaxConversionRatio;

  /// the minimum conversion ratio at maturity
  double m_dMinConversionRatio;

  /// Is the optional conversion enabled?
  bool m_bHasOptionalConversion;

  /// fixed cash call provision
  shared_ptr<CallSchedule> m_pCallSchedule;

  /// fixed share call provision
  shared_ptr<CallFixedShare> m_pCallFixedShare;

  /// variable share call provision
  // shared_ptr<CallVariableShare> m_pCallVariableShare;

  /// whether or not the call is set
  bool m_bCallSet;
}; // class PEPSLike


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_PEPSLIKE_H_
