/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/derivativevisitors/derivative_visitor_goodtype.h
// Purpose:     Visitor for Derivative-derived classes
// Author:      Vadim Zeitlin
// Created:     2004-05-11
// RCS-ID:      $Id: derivative_visitor_goodtype.h,v 1.19 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/derivativevisitors/derivative_visitor_goodtype.h
    @brief Visitor class accepting different Derivative kinds.

    This is the central part of the implementation of the visitor pattern for
    the derivatives. To do something different depending on the exact type of
    a Derivative object you have to define a new class inheriting from this one
    and do whatever is required in its methods.
 */

#ifndef _ITO33_FINANCE_DERIVATIVEVISITORS_DERIVATIVE_VISITOR_GOODTYPE_H_
#define _ITO33_FINANCE_DERIVATIVEVISITORS_DERIVATIVE_VISITOR_GOODTYPE_H_

#include "ito33/finance/derivative_visitor.h"

#include "ito33/finance/option.h"
#include "ito33/finance/exoticoption/asianoption.h"
#include "ito33/finance/referencecds.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/parbond.h"
#include "ito33/finance/eds.h"
#include "ito33/finance/varianceswap.h"
#include "ito33/finance/bondlike/bond.h"
#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/cboption.h"
#include "ito33/finance/bondlike/reset.h"
#include "ito33/finance/bondlike/pepslike.h"
#include "ito33/finance/bondlike/generalizedpepslike.h"
#include "ito33/finance/bondlike/percslike.h"
#include "ito33/finance/bondlike/attachedwarrantconvertiblebond.h"

namespace ito33
{

namespace finance
{

/**
    Derivatives visitor to get shared ptr to the right type of the object.

    An object of the class must be passed by the user code
    to any functions which may access heterogeneous collections of Derivatives,
    e.g. XML::Reader in IHG currently.
 */
class DerivativeVisitorGoodType : public DerivativeVisitor
{
public:
  
  /**
      Gets a bond like (if any). We are supposing implicitly that there is
      at most one bond like, otherwise, the behavior is not well defined.

      Currently, mandatory is ignored.
   */
  shared_ptr<ConvertibleLike> GetBondLike() const
  {
    if ( m_pConvertibleBond )
      return m_pConvertibleBond;
    
    if ( m_pAttachedWarrantCB )
      return m_pAttachedWarrantCB;

    if ( m_pReset )
      return m_pReset;

    return shared_ptr<ConvertibleLike>();
  }

  /// Gets shared_ptr to option if option has been visited, otherwise returns 0
  shared_ptr<Option> GetOption() const { return m_pOption; }

  /// Gets shared_ptr to 2nd option, otherwise returns 0
  shared_ptr<Option> GetOption2() const { return m_pOption2; }

  /// Gets visited CDS, otherwise returns 0
  shared_ptr<CDS> GetCDS() const { return m_pCDS; }

  /// Gets visited ReferenceCDS, otherwise returns 0
  shared_ptr<ReferenceCDS> GetReferenceCDS() const { return m_pReferenceCDS; }

  /// Gets visited ParBond, otherwise returns 0
  shared_ptr<ParBond> GetParBond() const { return m_pParBond; }

  /// Gets visited Bond, otherwise returns 0
  shared_ptr<Bond> GetBond() const { return m_pBond; }

  /// Gets visited attached warrant, otherwise returns 0
  shared_ptr<AttachedWarrantConvertibleBond> GetAttachedWarrantCB() const
  {
    return m_pAttachedWarrantCB;
  }   

  /// Gets visited convertibleBond, otherwise returns 0
  shared_ptr<ConvertibleBond> GetConvertibleBond() const
  { 
    return m_pConvertibleBond;
  } 

  /// Gets visited cb option, otherwise returns 0
  shared_ptr<CBOption> GetCBOption() const
  { 
    return m_pCBOption;
  }  

  /// Gets visited reset, otherwise returns 0
  shared_ptr<Reset> GetReset() const
  { 
    return m_pReset;
  } 

  /// Gets visited GeneralizedPEPSLike, otherwise returns 0
  shared_ptr<GeneralizedPEPSLike> GetGeneralizedPEPSLike() const
  { 
    return m_pGeneralizedPEPSLike;
  } 

  /// Gets visited PEPSLike, otherwise returns 0
  shared_ptr<PEPSLike> GetPEPSLike() const
  { 
    return m_pPEPSLike;
  } 

  /// Gets visited PERCSLike, otherwise returns 0
  shared_ptr<PERCSLike> GetPERCSLike() const
  { 
    return m_pPERCSLike;
  } 
 
  /// Gets visited EDS, otherwise returns 0
  shared_ptr<EDS> GetEDS() const
  {
    return m_pEDS;
  }

  /// Gets visited Asian Option, otherwise returns 0
  shared_ptr<AsianOption> GetAsianOption() const
  {
    return m_pAsianOption;
  }

  /// Gets visited variance swap, otherwise returns 0
  shared_ptr<VarianceSwap> GetVarianceSwap() const
  {
    return m_pVarianceSwap;
  }

public:
  
  /// Called for a CDS
  virtual void OnCDS(const CDS& contract)  
  { 
    m_pCDS = shared_ptr<CDS>(new CDS(contract));
  } 

  /// Called for a CDS
  virtual void OnReferenceCDS(const ReferenceCDS& contract)  
  { 
    m_pReferenceCDS = shared_ptr<ReferenceCDS>(new ReferenceCDS(contract));
  } 

  /// Called for a ParBond
  virtual void OnParBond(const ParBond& contract)  
  { 
    m_pParBond = shared_ptr<ParBond>(new ParBond(contract));
  } 

  /// Called for an EDS
  void OnEDS(const EDS &contract)  
  { 
    m_pEDS = shared_ptr<EDS>(new EDS(contract) );
  }

  /// Called for an Option
  virtual void OnOption(const Option& contract)  
  { 
    // The volpower calibration requires two options. Hence, two options
    // can appear in an xml file, and need to be stored. For now, just
    // store them individually, and have two Get functions.  If other
    // parametrizations (or whatever) require more than one derivative
    // of another type, and they are not stored in term structures, a more 
    // flexible storage solution may be required
    if (m_pOption)
      m_pOption2 = shared_ptr<Option>(new Option(contract));
    else
      m_pOption = shared_ptr<Option>(new Option(contract));
  } 
  
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

  /// Called for a cb option
  virtual void OnCBOption(const CBOption& contract)
  { 
    m_pCBOption
      = shared_ptr<CBOption>(new CBOption(contract));
  }  

  /// Called for a resettable convertible bond
  virtual void OnReset(const Reset& contract)  
  {   
    m_pReset = shared_ptr<Reset>(new Reset(contract));
  } 

  /// Called for a GeneralizedPEPSLike
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

  /// Called for attachedwarrant cb
  virtual void OnAttachedWarrantConvertibleBond(const 
    AttachedWarrantConvertibleBond  &contract)
  {
    m_pAttachedWarrantCB = shared_ptr<AttachedWarrantConvertibleBond>(
          new AttachedWarrantConvertibleBond(contract) );
  }

  /// Called for Asian option
  virtual void OnAsianOption(const AsianOption &contract)
  {
    m_pAsianOption = shared_ptr<AsianOption>(new AsianOption(contract));
  }

  /// Called for variance swap
  virtual void OnVarianceSwap(const VarianceSwap &contract)
  {
    m_pVarianceSwap = shared_ptr<VarianceSwap>(new VarianceSwap(contract));
  }


private:

  shared_ptr<AsianOption> m_pAsianOption;

  shared_ptr<CDS> m_pCDS;

  shared_ptr<ReferenceCDS> m_pReferenceCDS;

  shared_ptr<ParBond> m_pParBond;

  shared_ptr<EDS> m_pEDS;

  shared_ptr<Option> m_pOption;

  shared_ptr<Option> m_pOption2;

  shared_ptr<GeneralizedPEPSLike> m_pGeneralizedPEPSLike;

  shared_ptr<PEPSLike> m_pPEPSLike;

  shared_ptr<PERCSLike> m_pPERCSLike;

  shared_ptr<AttachedWarrantConvertibleBond> m_pAttachedWarrantCB;

  shared_ptr<ConvertibleBond> m_pConvertibleBond;

  shared_ptr<CBOption> m_pCBOption;

  shared_ptr<Reset> m_pReset;

  shared_ptr<Bond> m_pBond;

  shared_ptr<VarianceSwap> m_pVarianceSwap;

};


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_DERIVATIVEVISITORS_DERIVATIVE_VISITOR_GOODTYPE_H_
