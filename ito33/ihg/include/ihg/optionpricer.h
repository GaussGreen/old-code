/****************************************************************************
 * Name:        ihg/optionpricer.h
 * Purpose:     option pricer class
 * Author:      David
 * Created:     2003/12/10
 * RCS-ID:      $Id: optionpricer.h,v 1.21 2006/03/20 14:54:11 yann Exp $
 * Copyright:   (c) 2003-2003 Trilemma LLP
 ****************************************************************************/

/**
    @file ihg/optionpricer.h
    @brief Option pricer class

    Implementation of the pricer class for options.
*/

#ifndef _IHG_OPTIONPRICER_H_ 
#define _IHG_OPTIONPRICER_H_

#include "ito33/autoptr.h"
#include "ito33/common.h"
#include "ito33/dlldecl.h"

#include "ito33/pricing/optionmeshmanager.h"

#include "ihg/optioninstdata.h"
#include "ihg/optionstepper.h"

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

namespace ihg
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
  
  OptionStepper m_stepper;
  
  AutoPtr<OptionNumOutput> m_pNumOutput;
  

private:

  NO_COPY_CLASS(OptionPricer);

}; // class OptionPricer


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_OPTIONPRICER_H_

