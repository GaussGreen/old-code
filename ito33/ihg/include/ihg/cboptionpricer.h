/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/cboptionpricer.h
// Purpose:     cb option pricer class
// Author:      Nabil
// Created:     2005/10/18
// RCS-ID:      $Id: cboptionpricer.h,v 1.3 2006/03/20 14:54:11 yann Exp $
// Copyright:   (c) 2004-2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/cboptionpricer.h
    @brief cb option pricer class

    Implementation of the pricer class for cb options.
*/

#ifndef _IHG_CBOPTIONPRICER_H_ 
#define _IHG_CBOPTIONPRICER_H_

#include "ito33/autoptr.h"
#include "ito33/vector.h"
#include "ito33/dlldecl.h"

#include "ito33/pricing/specialengine.h"


namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL ComputationalFlags;
}

namespace pricing
{
  class CBOptionParams;
  class CBMeshManager;
}


namespace ihg
{

class Model;
class CBOptionInstData;
class CBOptionStepper;
class CBOptionNumOutput;

typedef pricing::SpecialEngine<pricing::CBOptionParams, pricing::CBMeshManager,
  CBOptionInstData, CBOptionStepper, CBOptionNumOutput> CBOptionEngine;

/**
   CB Option Pricer class

   This class simply constructs the relevant classes for cb option pricing,
   and calls engine to do the pricing.
 */
class CBOptionPricer
{

public:

  /// Construct from cb option parameters
  CBOptionPricer(pricing::CBOptionParams& cboptionparams, 
                 Model &model, 
                 const finance::ComputationalFlags& flags);

  /** 
      Price the cb option
      
      @return CBOptionNumOutput object containing all price data
  */
  AutoPtr<CBOptionNumOutput> Price()
  {     
    if ( IsPathDependent() )
      return PricePathDep();
    else
      return PriceNormal();
  }


protected:  
  
  /**
    Determine if the convertible bond is path dependent

    @return true if path dependent, false otherwise
  */
  bool IsPathDependent();

  /**
    Price a cb option with a normal underlying CB contract 
    (without path dependent features).

    @return CBOptionNumOutput autoptr containing all (numerical) pricing data
  */
  AutoPtr<CBOptionNumOutput> PriceNormal();

  /**
    Price a CBOption contract with path dependent features

    @return CBOptionNumOutput autoptr containing all (numerical) pricing data
  */
  AutoPtr<CBOptionNumOutput> PricePathDep();

  /// The cb option params of the cb option to price
  pricing::CBOptionParams& m_cboptionparams;
  
  /// The ihg model parameters
  Model& m_model;

  /// The computational flags (ie. what data to compute)
  const finance::ComputationalFlags& m_flags;


private:

  NO_COPY_CLASS(CBOptionPricer);

}; // class CBOptionPricer


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_CBOPTIONPRICER_H_
