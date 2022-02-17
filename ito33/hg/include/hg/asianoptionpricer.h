/////////////////////////////////////////////////////////////////////////////
// Name:        hg/asianoptionpricer.h
// Purpose:     asian option pricer class
// Created:     2006/03/02
// RCS-ID:      $Id: asianoptionpricer.h,v 1.3 2006/05/23 15:13:45 yann Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/asianoptionpricer.h
    @brief Asian option pricer class

    Implementation of the pricer class for asian options.
 */

#ifndef _HG_ASIANOPTIONPRICER_H_ 
#define _HG_ASIANOPTIONPRICER_H_

#include "ito33/vector.h"
#include "ito33/list.h"

#include "ito33/sharedptr.h"
#include "ito33/autoptr.h"
#include "ito33/common.h"

#include "ito33/pricing/optionmeshmanager.h"
#include "ito33/pricing/pathdepevent.h"

#include "hg/optioninstdata.h"
#include "hg/optionstepper.h"

namespace ito33
{
 
namespace finance { class ITO33_DLLDECL ComputationalFlags; }

namespace pricing
{
  class AsianOptionParams;
}

namespace hg
{
  class Model;
  class ComputationalFlags;
  class OptionNumOutput;

/**
    Asian Option Pricer class.

    This class simply constructs the relevant classes for Asian option 
    pricing, and calls engine to do the pricing.  
 */
class AsianOptionPricer
{
public:

  /// Constructs from Asian option parameters
  AsianOptionPricer(pricing::AsianOptionParams& params, 
                    Model& model, const finance::ComputationalFlags& flags);
  
  /** 
      Prices the Asian option.
      
      @return optionNumOutput object containing all price data
   */
  AutoPtr<OptionNumOutput> Price();

protected:

  pricing::AsianOptionParams& m_params;

  Model& m_model;

  const finance::ComputationalFlags& m_flags;
  
  AutoPtr<OptionNumOutput> m_pNumOutput;
  
private:
 

  NO_COPY_CLASS(AsianOptionPricer);

}; // class AsianOptionPricer


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_ASIANOPTIONPRICER_H_
