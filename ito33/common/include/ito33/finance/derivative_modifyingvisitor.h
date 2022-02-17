/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/derivative_modifyingvisitor.h
// Purpose:     Visitor for Derivative-derived where the visitor pattern can
//              modify the object visited
// Author:      Zhang
// Created:     October 7, 2005
// RCS-ID:      $Id: derivative_modifyingvisitor.h,v 1.7 2006/05/26 20:24:04 dave Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/derivative_modifyingvisitor.h
    @brief Visitor class accepting different Derivative kinds where
    the visitor can modify the object being visited.
 */

#ifndef _ITO33_FINANCE_DERIVATIVE_MODIFYING_VISITOR_H_
#define _ITO33_FINANCE_DERIVATIVE_MODIFYING_VISITOR_H_

#include "ito33/debug.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{

/// @name Forward declaration
//@{
class ITO33_DLLDECL AsianOption;
class ITO33_DLLDECL Option;
class ITO33_DLLDECL OneTouch;
class ITO33_DLLDECL CDS;
class ITO33_DLLDECL ReferenceCDS;
class ITO33_DLLDECL EDS;
class ITO33_DLLDECL ParBond;
class ITO33_DLLDECL VarianceSwap;

class ITO33_DLLDECL Bond;
class ITO33_DLLDECL AttachedWarrantConvertibleBond;
class ITO33_DLLDECL ConvertibleBond;
class ITO33_DLLDECL GeneralizedPEPSLike;
class ITO33_DLLDECL PEPSLike;
class ITO33_DLLDECL PERCSLike;
class ITO33_DLLDECL Reset;
//@}

/// Derivatives modifying visitor.
class DerivativeModifyingVisitor
{
public:
  
  /*
    Normally, this should be a pure virtual class. However, every time we 
    add a new derivative, all specific vistors would be broken. That is why
    we implement, by default, all OnXXX() functions by FAIL.
  */
  
  /// Called for an Option
  virtual void OnOption(Option&)  { FAIL("Not implemented."); } 
  
  /// Called for a OneTouch
  virtual void OnOneTouch(OneTouch&)  { FAIL("Not implemented."); } 

  /// Called for a CDS
  virtual void OnCDS(CDS&)  { FAIL("Not implemented."); }

  /// Called for a ReferenceCDS
  virtual void OnReferenceCDS(ReferenceCDS&)  { FAIL("Not implemented."); }

  /// Called for a ParBond
  virtual void OnParBond(ParBond&)  { FAIL("Not implemented."); }

  /// Called for a EDS
  virtual void OnEDS(EDS&)  { FAIL("Not implemented."); }
 
  /// Called for a Variance Swap
  virtual void OnVarianceSwap(VarianceSwap&) { FAIL("Not implemented."); }

  /// Called for an Asian option
  virtual void OnAsianOption(AsianOption&) { FAIL("Not implemented."); }

  /// Called for a Bond
  virtual void OnBond(Bond&)  { FAIL("Not implemented."); }
   
  /// Called for a ConvertibleBond
  virtual void OnConvertibleBond(ConvertibleBond&) {FAIL("Not implemented."); } 

  /// Called for a CB with attached warrrant
  virtual 
  void OnAttachedWarrantConvertibleBond(AttachedWarrantConvertibleBond&)
  { 
    FAIL("Not implemented."); 
  }

  /// Called for a PEPS
  virtual void OnPEPSLike(PEPSLike&)  { FAIL("Not implemented."); } 

  /// Called for a GeneralizedPEPSLike
  virtual void OnGeneralizedPEPSLike(GeneralizedPEPSLike&)
  {
    FAIL("Not implemented."); 
  }

  /// Called for a PERCS
  virtual void OnPERCSLike(PERCSLike&)  { FAIL("Not implemented."); }

  /// Called for a Reset
  virtual void OnReset(Reset&)  { FAIL("Not implemented."); }

  /// Virtual dtor for any base class
  virtual ~DerivativeModifyingVisitor() { }
};

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_DERIVATIVE_MODIFYING_VISITOR_H_
