/////////////////////////////////////////////////////////////////////////////
// Name:        hg/optionpricer.h
// Purpose:     implementation of pricer class for regular options
// Created:     2005/01/13
// RCS-ID:      $Id: optionpricer.h,v 1.3 2006/04/17 17:58:53 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/optionpricer.h
   @brief implementation of pricer class for regular options
 */

#ifndef _HG_OPTIONPRICER_H_ 
#define _HG_OPTIONPRICER_H_

#include "ito33/autoptr.h"
#include "ito33/common.h"

#include "ito33/pricing/optionmeshmanager.h"

#include "hg/optioninstdata.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL ComputationalFlags;
}

namespace pricing
{
  class OptionParams;
}

namespace hg
{
  class Model;
  class OptionNumOutput;

/**
   Option Pricer class

   This class simply constructs the relevant classes for option pricing,
   and calls engine to do the pricing.  It then calculates theta values
   (if requested by the user) using finite differencing.
 */
class OptionPricer
{
public:

  /// Construct from option parameters
  OptionPricer(pricing::OptionParams& params, 
               Model& model,
               const finance::ComputationalFlags& flags);
  
  /** 
      Price the option
      
      @return optionNumOutput object containing all price data
   */
  AutoPtr<OptionNumOutput> Price();


protected:

  pricing::OptionParams& m_params;

  Model& m_model;

  const finance::ComputationalFlags& m_flags;

  pricing::OptionMeshManager m_meshes;
  
  OptionInstData m_instdata;
  
  AutoPtr<OptionNumOutput> m_pNumOutput;
  

private:

  NO_COPY_CLASS(OptionPricer);

}; // class OptionPricer


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_OPTIONPRICER_H_

