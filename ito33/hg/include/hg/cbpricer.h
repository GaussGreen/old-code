/////////////////////////////////////////////////////////////////////////////
// Name:        hg/cbpricer.h
// Purpose:     cb pricer class
// Created:     2005/04/11
// RCS-ID:      $Id: cbpricer.h,v 1.5 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/cbpricer.h
    @brief cb pricer class

    Implementation of the pricer class for convertible bonds.
*/

#ifndef _HG_CBPRICER_H_ 
#define _HG_CBPRICER_H_

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
  class CBLikeParams;
  class CBMeshManager;
  class PathDepEvent;
}

namespace hg
{

class Model;
class CBInstData;
class CBStepper;
class CBNumOutput;
class Payoff;

typedef pricing::SpecialEngine< pricing::CBLikeParams, 
                                pricing::CBMeshManager,
                                CBInstData, 
                                CBStepper,
                                CBNumOutput
                              > CBEngine;

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
           Model& model, 
           const finance::ComputationalFlags& flags);

  /**
     Sets the initial values defined from outside
   */
  void SetInitialValue(const shared_ptr<Payoff>& pPayoff)
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
     Determine if the convertible bond is path dependent

     @return true if path dependent, false otherwise
   */
  bool IsPathDependent();

  /**
     Price a normal CB contract without path dependent features.

     @return CBNumOutput autoptr containing all (numerical) pricing data
   */
  AutoPtr<CBNumOutput> PriceNormal();

  /**
     Price a CB contract with path dependent features.

     @return CBNumOutput autoptr containing all (numerical) pricing data
   */
  AutoPtr<CBNumOutput> PricePathDep();

  /** 
     Create a path dependent pricing event list
   */
  void 
  CreatePathDepEvents
  (std::list< shared_ptr<pricing::PathDepEvent> >& pathDepEvents,
   std::vector< AutoPtr<pricing::CBLikeParams> >& cblikeParams);

  /// The cb contract to price
  pricing::CBLikeParams& m_cbparams;
  
  /// The HG model parameters
  Model& m_model;

  /// The computational flags (ie. what data to compute)
  const finance::ComputationalFlags& m_flags;

  shared_ptr<Payoff> m_pPayoff;


private:

  NO_COPY_CLASS(CBPricer);

}; // class CBPricer


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_CBPRICER_H_
