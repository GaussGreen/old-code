/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/cbpricer.h
// Purpose:     cb pricer class
// Author:      Nabil
// Created:     2004/04/14
// RCS-ID:      $Id: cbpricer.h,v 1.30 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/cbpricer.h
    @brief cb pricer class

    Implementation of the pricer class for convertible bonds.
*/

#ifndef _IHG_CBPRICER_H_ 
#define _IHG_CBPRICER_H_

#include "ito33/autoptr.h"
#include "ito33/vector.h"

#include "ito33/pricing/specialengine.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{
  class Payoff;
  class ITO33_DLLDECL ComputationalFlags;
}

namespace pricing
{
  class CBLikeParams;
  class CBMeshManager;
  class PathDepEvent;
}


namespace ihg
{

class Model;
class CBInstData;
class CBStepper;
class CBNumOutput;


typedef pricing::SpecialEngine<pricing::CBLikeParams, pricing::CBMeshManager,
  CBInstData, CBStepper, CBNumOutput> CBEngine;


/**
   CB Pricer class

   This class simply constructs the relevant classes for cb pricing,
   and calls engine to do the pricing.
 */
class CBPricer
{
public:

  /// Construct from cblike parameters
  CBPricer(pricing::CBLikeParams& cbparams, 
           Model &model, 
           const finance::ComputationalFlags& flags);

  /**
      Sets an external initial condition/payoff.

      @param pPayoff The external initial condition
  */
  void SetInitialValue(const shared_ptr<finance::Payoff>& pPayoff)
  {
    m_pPayoff = pPayoff;
  }

  /** 
      Price the cb
      
      @return CBNumOutput object containing all price data
  */
  AutoPtr<CBNumOutput> Price();


protected:

  /**
    Price a normal CB contract without path dependent features

    @return CBNumOutput autoptr containing all (numerical) pricing data
  */
  AutoPtr<CBNumOutput> PriceNormal();

  /**
    Price a CB contract with path dependent features

    @return CBNumOutput autoptr containing all (numerical) pricing data
  */
  AutoPtr<CBNumOutput> PricePathDep();

  /// The cb contract to price
  pricing::CBLikeParams& m_cbparams;
  
  /// The ihg model parameters
  Model& m_model;

  /// The computational flags (ie. what data to compute)
  const finance::ComputationalFlags& m_flags;

  /// External payoff (eg. used for call notice)
  shared_ptr<finance::Payoff> m_pPayoff;


private:

  NO_COPY_CLASS(CBPricer);

}; // class CBPricer


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_CBPRICER_H_

