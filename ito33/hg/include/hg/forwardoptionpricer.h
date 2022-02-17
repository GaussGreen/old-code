/////////////////////////////////////////////////////////////////////////////
// Name:        hg/forwardoptionpricer.h
// Purpose:     implementation of pricer class for forward options
// Created:     2005/05/05
// RCS-ID:      $Id: forwardoptionpricer.h,v 1.3 2006/04/17 17:58:53 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/forwardoptionpricer.h
   @brief implementation of pricer class for forward options
 */

#ifndef _HG_FORWARDOPTIONPRICER_H_ 
#define _HG_FORWARDOPTIONPRICER_H_

#include "ito33/autoptr.h"
#include "ito33/common.h"

#include "ito33/pricing/forwardoptionmeshmanager.h"

#include "hg/forwardoptioninstdata.h"

namespace ito33
{

namespace finance
{
  class ComputationalFlags;
}

namespace pricing
{
  class ForwardOptionParams;
}

namespace hg
{
  class Model;
  class ForwardOptionNumOutput;

/**
   ForwardOption Pricer class

   This class simply constructs the relevant classes for forward option 
   pricing, and calls engine to do the pricing.  It then calculates theta 
   values (if requested by the user) using finite differencing.
 */
class ForwardOptionPricer
{
public:

  /// Construct from forward option parameters
  ForwardOptionPricer(pricing::ForwardOptionParams& params, 
                      Model& model,
                      const finance::ComputationalFlags& flags);
  
  /** 
      Price the forward option
      
      @return ForwardOptionNumOutput object containing all price data
   */
  AutoPtr<ForwardOptionNumOutput> Price();


protected:

  /// The contract params
  pricing::ForwardOptionParams& m_params;

  /// The underlying hg model
  Model& m_model;

  /// The computational flags (i.e. what Greeks to compute, surfaces, etc)
  const finance::ComputationalFlags& m_flags;

  /// The mesh 
  pricing::ForwardOptionMeshManager m_meshes;
  
  /// The instantaneous data
  ForwardOptionInstData m_instdata;
  
  /// The internal numerical output
  AutoPtr<ForwardOptionNumOutput> m_pNumOutput;
  

private:

  NO_COPY_CLASS(ForwardOptionPricer);

}; // class ForwardOptionPricer


} // namespace hg

} // namespace ito33

#endif // #ifndef _FORWARDHG_OPTIONPRICER_H_
