/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/attachedwarrantcbpricer.h
// Purpose:     Attached warrant cb pricer class
// Author:      Ito33
// Created:     2005/01/13
// RCS-ID:      $Id: attachedwarrantcbpricer.h,v 1.8 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/attachedwarrantcbpricer.h
    @brief Attached warrant cb pricer class

    Implementation of the pricer class for Attached warrant.
*/

#ifndef _IHG_ATTACHEDWARRANTCBPRICER_H_ 
#define _IHG_ATTACHEDWARRANTCBPRICER_H_

#include "ito33/autoptr.h"
#include "ito33/common.h"
#include "ito33/dateutils.h"
#include "ito33/dlldecl.h"

#include "ito33/pricing/specialengine.h"
#include "ito33/pricing/attachedwarrantcbparams.h"

#include "ito33/finance/computationalflags.h"

namespace ito33
{

namespace finance
{
  class Payoff;
}

namespace pricing
{
  class CBMeshManager;
}


namespace ihg
{

class CBInstData;
class CBStepper;
class CBNumOutput;

typedef pricing::SpecialEngine<pricing::CBLikeParams, pricing::CBMeshManager,
                               CBInstData, CBStepper, CBNumOutput> 
        AttachedWarrantConvertibleBondEngine;


/**
   Warrant Pricer class

   This class simply constructs the relevant pricing end event classes 
   for cb pricing, and calls the appropriate engine to do the pricing.
 */
class AttachedWarrantConvertibleBondPricer
{
public:

  /// Construct from  parameters
  AttachedWarrantConvertibleBondPricer(
              pricing::AttachedWarrantConvertibleBondParams& warrantparams, 
              Model &model, 
              const finance::ComputationalFlags& flags)
    : m_warrantparams(warrantparams), m_model(model),m_flags(flags)
  {
    m_pFlagsPricer = make_ptr( new finance::ComputationalFlags(flags) );

    if ( flags.GetAnalysisDate().IsValid() )
      m_warrantparams.SetAnalysisTime( GetDoubleFrom(flags.GetAnalysisDate()) );
  }
  
  /**
      Sets an external initial condition/payoff.

      @param pPayoff The external initial condition
  */
  void SetInitialValue(const shared_ptr<finance::Payoff>& pPayoff)
  {
    m_pPayoff = pPayoff;
  }

  /** 
      Price the instrument
      
      @return CBNumOutput object containing all price data
  */
  AutoPtr<CBNumOutput> Price();


  /** 
      Price the instrument without computing the greeks
      
      @return CBNumOutput object containing all price data
  */
  AutoPtr<CBNumOutput> PriceOnly();


  /**
    indicate whether or not the greeks can be computed
    by greek we mean, rho and vega. These greeks can
    sometime be not computed when a mix of pricer is used.
  */
  bool CanCalculatePriceAndGreeks();

protected:

  /// 1D pricing when put and call occur on same date
  AutoPtr<CBNumOutput> Price1D();

  /// Normal 2D pricing when pricing extends past the reset date
  AutoPtr<CBNumOutput> Price2D();

  /// Create payoff for "1D" pricing
  shared_ptr<finance::Payoff> ConstructPutCallPayoff();

  /// The params being used for pricing
  pricing::AttachedWarrantConvertibleBondParams& m_warrantparams;
  
  /// The ihg model parameters
  Model& m_model;

  /// The computational flags (ie. what data to compute)
  const finance::ComputationalFlags& m_flags;

  ///copy of the user flags, this way the
  //user flags are never changed. This is particularly useful
  //when solving the acwb when turning off the
  //flags to compute the greeks correclty.
  shared_ptr<finance::ComputationalFlags> m_pFlagsPricer;

  /// External payoff (eg. used for call notice)
  shared_ptr<finance::Payoff> m_pPayoff;

private:

  NO_COPY_CLASS(AttachedWarrantConvertibleBondPricer);

}; // class AttachedWarrantConvertibleBondPricer


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_ATTACHEDWARRANTCBPRICER_H_
