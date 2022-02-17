/****************************************************************************
 * Name:        ihg/asianoptionpricer.h
 * Purpose:     asian option pricer class
 * Author:      April 6, 2005
 * Created:     ITO 33 Canada
 * RCS-ID:      $Id: asianoptionpricer.h,v 1.4 2006/05/23 15:13:45 yann Exp $
 * Copyright:   (c) 2005 Trilemma LLP
 ****************************************************************************/

/**
    @file ihg/asianoptionpricer.h
    @brief Option pricer class

    Implementation of the pricer class for asian options.
*/

#ifndef _IHG_ASIANOPTIONPRICER_H_ 
#define _IHG_ASIANOPTIONPRICER_H_

#include "ito33/beforestd.h"
#include <vector>
#include <list>
#include "ito33/afterstd.h"
#include "ito33/dlldecl.h"

#include "ito33/sharedptr.h"
#include "ito33/autoptr.h"
#include "ito33/common.h"

#include "ito33/pricing/optionmeshmanager.h"
#include "ito33/pricing/pathdepevent.h"

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
  class AsianOptionParams;
}

namespace ihg
{
  class Model;
  class OptionNumOutput;

/**
   Asian Option Pricer class

   This class simply constructs the relevant classes for option pricing,
   and calls engine to do the pricing.  It then calculates theta values
   (if requested by the user) using finite differencing.
 */
class AsianOptionPricer
{
public:

  /// Construct from option parameters
  AsianOptionPricer(pricing::AsianOptionParams & params, 
    Model& model, const finance::ComputationalFlags& flags);
  
  /** 
      Price the option
      
      @return optionNumOutput object containing all price data
   */
  AutoPtr<OptionNumOutput> Price();


protected:


  pricing::AsianOptionParams& m_params;

  Model& m_model;

  const finance::ComputationalFlags& m_flags;
  
  AutoPtr<OptionNumOutput> m_pNumOutput;
  

private:

  /**
    Used to determined whether or not we can do a similarity transform
  */
  bool IsPathDependent();


  NO_COPY_CLASS(AsianOptionPricer);

}; // class AsianOptionPricer


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_ASIANOPTIONPRICER_H_

