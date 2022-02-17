// Name:        ito33/finance/derivativevisitors/bondlikevisitor.h
// Purpose:     Visitor for bond like classes
// Author:      ZHANG Yunzhi
// Created:     2004-09-09
// RCS-ID:      $Id: bondlikevisitor.h,v 1.10 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////
 
/**
   @file ito33/finance/derivativevisitors/bondlikevisitor.h
 */
 
#ifndef _ITO33_FINANCE_DERIVATIVEVISITORS_BONDLIKEVISITOR_H_
#define _ITO33_FINANCE_DERIVATIVEVISITORS_BONDLIKEVISITOR_H_

#include "ito33/finance/derivative_visitor.h"

#include "ito33/finance/bondlike/bond.h"
#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/generalizedpepslike.h"
#include "ito33/finance/bondlike/pepslike.h"
#include "ito33/finance/bondlike/percslike.h"
#include "ito33/finance/bondlike/reset.h"
#include "ito33/finance/bondlike/attachedwarrantconvertiblebond.h"
#include "ito33/finance/bondlike/cboption.h"

#include "ito33/dlldecl.h"


namespace ito33
{

namespace finance
{

/**
    bondlike visitor

    An object of the class must be passed by the user code
    to any functions which may access heterogeneous collections of Derivatives,
    e.g. XML::Reader in IHG currently.
 */
class ITO33_DLLDECL BondLikeVisitor: public DerivativeVisitor
{
public:

  /// Gets visited Bond, otherwise returns 0
  shared_ptr<Bond> GetBond() { return m_pBond; }

   
  /// Gets visited convertibleBond, otherwise returns 0
  shared_ptr<ConvertibleBond> GetConvertibleBond() 
  { 
    return m_pConvertibleBond;
  } 

  /// Gets visited Reset, otherwise returns 0
  shared_ptr<Reset> GetReset() 
  { 
    return m_pReset;
  } 

  /// Gets visited GeneralizedPEPSLike, otherwise returns 0
  shared_ptr<GeneralizedPEPSLike> GetGeneralizedPEPSLike() 
  { 
    return m_pGeneralizedPEPSLike;
  } 

  /// Gets visited PEPSLike, otherwise returns 0
  shared_ptr<PEPSLike> GetPEPSLike() 
  { 
    return m_pPEPSLike;
  } 

  /// Gets visited PERCSLike, otherwise returns 0
  shared_ptr<PERCSLike> GetPERCSLike() 
  { 
    return m_pPERCSLike;
  } 

  /// Gets visited Warrant, otherwise returns 0
  shared_ptr<finance::AttachedWarrantConvertibleBond> 
    GetAttachedWarrantConvertibleBond()
  {
    return m_pWarrant;
  }

  /// Gets visited CBOption, otherwise returns 0
  shared_ptr<CBOption> GetCBOption() 
  { 
    return m_pCBOption;
  } 

public:

  /// Called for a Bond
  virtual void OnBond(const Bond& contract)  
  { 
    m_pBond = shared_ptr<Bond>(new Bond(contract));
  } 

  /// Called for a ConvertibleBond
  virtual void OnConvertibleBond(const ConvertibleBond& contract)
  { 
    m_pConvertibleBond
      = shared_ptr<ConvertibleBond>(new ConvertibleBond(contract));
  } 

  /// Called for a Reset
  virtual void OnReset(const Reset& contract)  
  {   
    m_pReset = shared_ptr<Reset>(new Reset(contract));
  } 

  virtual void OnGeneralizedPEPSLike(const GeneralizedPEPSLike& contract)
  {
    m_pGeneralizedPEPSLike
        = shared_ptr<GeneralizedPEPSLike>(new GeneralizedPEPSLike(contract));
  }

  /// Called for a PEPS
  virtual void OnPEPSLike(const PEPSLike& contract) 
  { 
    m_pPEPSLike = shared_ptr<PEPSLike>(new PEPSLike(contract));
  } 

  /// Called for a PERCS
  virtual void OnPERCSLike(const PERCSLike& contract)  
  {   
    m_pPERCSLike = shared_ptr<PERCSLike>(new PERCSLike(contract));
  } 

  /// Called for warrant
  virtual void OnAttachedWarrantConvertibleBond(
                const AttachedWarrantConvertibleBond &contract)
  {
    m_pWarrant = shared_ptr<AttachedWarrantConvertibleBond>
      ( new AttachedWarrantConvertibleBond(contract) );
  }

  /// Called for a CBOption
  virtual void OnCBOption(const CBOption& contract)
  { 
    m_pCBOption
      = shared_ptr<CBOption>(new CBOption(contract));
  } 

private:

  shared_ptr<GeneralizedPEPSLike> m_pGeneralizedPEPSLike;

  shared_ptr<PEPSLike> m_pPEPSLike;

  shared_ptr<PERCSLike> m_pPERCSLike;

  shared_ptr<ConvertibleBond> m_pConvertibleBond;

  shared_ptr<Bond> m_pBond;
  
  shared_ptr<Reset> m_pReset;

  shared_ptr<AttachedWarrantConvertibleBond> m_pWarrant;

  shared_ptr<CBOption> m_pCBOption;

};

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_DERIVATIVEVISITORS_BONDLIKEVISITOR_H_
