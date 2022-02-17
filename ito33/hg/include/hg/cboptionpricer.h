/////////////////////////////////////////////////////////////////////////////
// Name:        hg/cboptionpricer.h
// Purpose:     cb option pricer class
// Created:     2006/01/19
// RCS-ID:      $Id: cboptionpricer.h,v 1.2 2006/03/20 14:54:09 yann Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/cboptionpricer.h
    @brief cb option pricer class

    Implementation of the pricer class for cb options.
 */

#ifndef _HG_CBOPTIONPRICER_H_ 
#define _HG_CBOPTIONPRICER_H_

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

namespace hg
{

class Model;
class CBOptionInstData;
class CBOptionStepper;
class CBOptionNumOutput;


typedef pricing::SpecialEngine<pricing::CBOptionParams, pricing::CBMeshManager,
  CBOptionInstData, CBOptionStepper, CBOptionNumOutput> CBOptionEngine;

/**
    Class simply constructs the relevant classes for cb option pricing,
    and calls engine to do the pricing.
 */
class CBOptionPricer
{

public:

  /// Construct from cb option parameters
  CBOptionPricer(pricing::CBOptionParams& cboptionparams, 
                 Model& model, 
                 const finance::ComputationalFlags& flags);

  /** 
      Prices the cb option.
      
      @return CBOptionNumOutput object containing all price data
   */
  AutoPtr<CBOptionNumOutput> Price();

protected:  
  
  /**
      Determines if the convertible bond is path dependent.

      @return true if path dependent, false otherwise
   */
  bool IsPathDependent();

  /**
      Prices a cb option with a normal underlying CB contract 
      (without path dependent features).

      @return CBOptionNumOutput autoptr containing all (numerical) pricing data
   */
  AutoPtr<CBOptionNumOutput> PriceNormal();

  /**
      Prices a CBOption contract with path dependent features

      @return CBOptionNumOutput autoptr containing all (numerical) pricing data
   */
  AutoPtr<CBOptionNumOutput> PricePathDep();

  /// The cb option params of the cb option to price
  pricing::CBOptionParams& m_cboptionparams;
  
  /// The hg model parameters
  Model& m_model;

  /// The computational flags (ie. what data to compute)
  const finance::ComputationalFlags& m_flags;


private:

  NO_COPY_CLASS(CBOptionPricer);

}; // class CBOptionPricer


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_CBOPTIONPRICER_H_
