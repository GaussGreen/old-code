/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/basket_visitor.h
// Purpose:     Visitor for Derivative-derived classes used in calibration
// Created:     2006/08/07
// RCS-ID:      $Id: basket_visitor.h,v 1.5 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/basket_visitor.h
    @brief Visitor class accepting different Derivatives used in calibration.
 */

#ifndef _ITO33_FINANCE_BASKET_VISITOR_H_
#define _ITO33_FINANCE_BASKET_VISITOR_H_

#include "ito33/sharedptr.h"

#include "ito33/finance/derivative_visitor.h"

namespace ito33
{

namespace finance
{

/// Forward declaration of the calibration basket
class ITO33_DLLDECL Derivatives;

/// Forward declaration for collection of European options
class ForwardOption;

/// Forward declaration of generic derivative, no special handling
class ITO33_DLLDECL Derivative;

/// Forward declaration of instrument that will be handled specially
class ITO33_DLLDECL Option;
class ITO33_DLLDECL OneTouch;
class ITO33_DLLDECL FXOneTouch;
class ITO33_DLLDECL CDS;
class ITO33_DLLDECL ReferenceCDS;
class ITO33_DLLDECL EDS;
class ITO33_DLLDECL ParBond;
class ITO33_DLLDECL VarianceSwap;
class ITO33_DLLDECL VarianceSwaption;

/**
    Derivatives visitor for calibration.

    This class can be used for preprocessing of the calibration basket. It
    will try to use forward equation whenever possible, or one backward
    equation for multiple instruments; it will set up the default weights
    that will be multiplied to the absolute error on price.
 */
class BasketVisitor : protected DerivativeVisitor
{
public:

  /**
      Default ctor takes the calibration basket that will be handled.

      By default, forward equation will be used for European options.

      @param derivatives The calibration basket
   */
  BasketVisitor(const Derivatives& derivatives);

  /**
      Enables/disables forward equation for European options.

      @param bUseForward Use or not forward equation for European option
   */
  void EnableForwardForOption(bool bUseForward = true)
  { 
    m_bUseForward = bUseForward;
  }

  /**
      Enables/disables advanced target eventually different than market price.

      @param bUseAdvancedTarget Use or not adavanced target
   */
  void EnableAdvancedTarget(bool bUseAdvancedTarget = true)
  {
    m_bUseAdvancedTarget = bUseAdvancedTarget;
  }

  /**
      Enables/disables relative error on target.

      @param bUseRelativeError Use or not relative error on target
   */
  void EnableRelativeError(bool bUseRelativeError = true)
  {
    m_bUseRelativeError = bUseRelativeError;
  }

  /**
      Runs the basket visitor to preprocess the calibration basket.
   */
  void Run();
     
  /**
      Returns the collection of European options if any.

      @return The collection of European options or null
   */
  shared_ptr<ForwardOption> GetForwardOption() const
  {
    return m_pForwardOption;
  }

  /**
      Returns the generic instruments if any.

      @return The collection of generic instruments or null
   */
  shared_ptr<Derivatives> GetGenericDerivatives() const
  {
    return m_pGenericDerivatives;
  }

  virtual void OnOption(const Option&);

  virtual void OnOneTouch(const OneTouch&);

private:
  
  /// The original basket
  const Derivatives& m_derivatives;

  /// Called for any derivative that is not handled specially
  void OnGeneric();

  /// If forward equation will be used for European option
  bool m_bUseForward;

  /// If advanced target that might be different than market price will be used
  bool m_bUseAdvancedTarget;

  /// If relative error on target will be used
  bool m_bUseRelativeError;

  /// Current derivative in the calibration basket
  shared_ptr<Derivative> m_pDerivative;

  /// Current weight
  double m_dWeight;

  /// European options can be handled by forward equation
  shared_ptr<ForwardOption> m_pForwardOption;

  /// Sets of generic derivatives
  shared_ptr<Derivatives> m_pGenericDerivatives;

  NO_COPY_CLASS(BasketVisitor);

}; // class BasketVisitor

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BASKET_VISITOR_H_
