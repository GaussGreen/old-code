/////////////////////////////////////////////////////////////////////////////
// Name:        hg/varianceswap_closedform.h
// Purpose:     HG variance swap pricer for vanilla vs
// Created:     2007/07/06
// RCS-ID:      $Id: varianceswap_closedform.h,v 1.4 2006/08/17 10:58:06 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/varianceswap_closedform.h
    @brief HG variance swap pricer class
 */

#ifndef _HG_VARIANCESWAP_CLOSEDFORM_H_
#define _HG_VARIANCESWAP_CLOSEDFORM_H_

#include "ito33/dlldecl.h"
#include "ito33/autoptr.h"

namespace ito33
{

namespace finance 
{ 
  class ITO33_DLLDECL ComputationalFlags; 
}

namespace pricing
{
  class VarianceSwapParams;
}

namespace hg
{

  class Model;
  class BackwardNumOutput;

/**
    Vanilla variance swap pricer class.

    Vanilla variance swap prices can be computed analytically.
 */
class VarianceSwapClosedForm
{
public:

  /**
      Constructor.

      @param params reference to variance swap pricing object
      @param model the HG model to use
      @param flags the computational flags specifying what to compute
   */
  VarianceSwapClosedForm(pricing::VarianceSwapParams& params, 
                         Model& model, 
                         const finance::ComputationalFlags& flags);
  
  /** 
      Prices the variance swap.
      
      @return object containing all price data
   */
  AutoPtr<BackwardNumOutput> Price();

  /** 
      Prices the variance swap by using the log contract.
      
      @return object containing all price data
   */
  AutoPtr<BackwardNumOutput> PriceByLog();

protected:

  /** 
      Prices the gamma variance swap.
      
      @return object containing all price data
   */
  AutoPtr<BackwardNumOutput> PriceGammaSwap();

  /// The variance swap params being used for pricing
  pricing::VarianceSwapParams& m_params;

  /// The HG model parameters
  Model& m_model;

  /// The computational flags (ie. what data to compute)
  const finance::ComputationalFlags& m_flags;
  
private:
  
  NO_COPY_CLASS(VarianceSwapClosedForm);

}; // class VarianceSwapClosedForm


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_VARIANCESWAP_CLOSEDFORM_H_
