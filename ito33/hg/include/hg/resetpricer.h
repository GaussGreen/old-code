/////////////////////////////////////////////////////////////////////////////
// Name:        hg/resetpricer.h
// Purpose:     resettable cb pricer class
// Created:     2006/04/17
// RCS-ID:      $Id: resetpricer.h,v 1.2 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/resetpricer.h
    @brief resettable cb pricer class

    Implementation of the reset pricer class.
*/

#ifndef _HG_RESETPRICER_H_ 
#define _HG_RESETPRICER_H_

#include "ito33/autoptr.h"
#include "ito33/common.h"
#include "ito33/dateutils.h"
#include "ito33/dlldecl.h"

#include "ito33/pricing/specialengine.h"
#include "ito33/pricing/resetparams.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL ComputationalFlags;
}

namespace pricing
{
  class CBMeshManager;
}


namespace hg
{

class CBInstData;
class CBStepper;
class ResetNumOutput;
class Payoff;


typedef pricing::SpecialEngine<pricing::CBLikeParams, pricing::CBMeshManager,
  CBInstData, CBStepper, ResetNumOutput> ResetEngine;


/**
   Resettable CB Pricer class.

   This class simply constructs the relevant pricing and event classes 
   for reset pricing, and calls the appropriate engine to do the pricing.
 */
class ResetPricer
{
public:

  /// Construct from reset parameters
  ResetPricer(pricing::ResetParams& resetparams, 
              Model &model, 
              const finance::ComputationalFlags& flags)
    : m_resetparams(resetparams), 
      m_model(model),
      m_flags(flags)
  {
    if ( flags.GetAnalysisDate().IsValid() )
      m_resetparams.SetAnalysisTime( GetDoubleFrom(flags.GetAnalysisDate()) );
  }
  
  /**
      Sets an external initial condition/payoff.

      @param pPayoff The external initial condition
  */
  void SetInitialValue(const shared_ptr<Payoff>& pPayoff)
  {
    m_pPayoff = pPayoff;
  }

  /** 
      Price the reset.
      
      @return CBNumOutput object containing all price data
  */
  AutoPtr<CBNumOutput> Price();


protected:

  /// The reset params being used for pricing
  pricing::ResetParams& m_resetparams;
  
  /// The hg model parameters
  Model& m_model;

  /// The computational flags (ie. what data to compute)
  const finance::ComputationalFlags& m_flags;

  /// External payoff (eg. used for call notice)
  shared_ptr<Payoff> m_pPayoff;


  /// check whether or not a path dependent or normal solve is needed
  bool IsPathDependent() const;

  /** 
      Price the reset using the path dependent framework
      
      @return CBNumOutput object containing all price data
   */
  AutoPtr<CBNumOutput> PricePathDependent();

  /** 
      Price the reset after transforming to a one dimensional problem.
      Limited to
        - yield dividend
        - previous ratio reset
        - initial ratio reset but with only one date
       
      @return CBNumOutput object containing all price data
  */
  AutoPtr<CBNumOutput> PriceOneD();


private:

  NO_COPY_CLASS(ResetPricer);

}; // class ResetPricer


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_RESETPRICER_H_
