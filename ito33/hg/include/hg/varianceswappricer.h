/////////////////////////////////////////////////////////////////////////////
// Name:        hg/varianceswappricer.h
// Purpose:     HG variance swap pricer class 
// Created:     2006/03/05
// RCS-ID:      $Id: varianceswappricer.h,v 1.6 2006/07/25 21:43:46 dave Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/varianceswappricer.h
    @brief HG variance swap pricer class
 */

#ifndef _HG_VARIANCESWAPPRICER_H_
#define _HG_VARIANCESWAPPRICER_H_

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
  class VarianceSwapNumOutput;

/**
    Variance swap pricer class.

    In general, variance swaps can be priced by introducing two new state
    variables: previous spot and current average of the squared returns.
    Hence, it is a 3D problem.  A similarity reduction can sometimes be
    used to reduce the dimension by one.  The path dependent pricer is used
    regardless.
 */
class VarianceSwapPricer
{
public:

  /**
      Constructor.

      @param params reference to variance swap pricing object
      @param model the HG model to use
      @param flags the computational flags specifying what to compute
   */
  VarianceSwapPricer(pricing::VarianceSwapParams& params, 
                     Model& model, 
                     const finance::ComputationalFlags& flags);
  
  /** 
      Prices the variance swap.
      
      @return VarianceSwapNumOutput object containing all price data
   */
  AutoPtr<VarianceSwapNumOutput> Price();

protected:

  /**
      Checks if a similarity reduction can be used.

      A similarity reduction can be used to reduce the problem by
      one dimension if:
      - there are no cash dividends between valuation and maturity
      - the variance swap does not have any exotic features

      Other checks may also be done (such as the "use similarity reduction"
      computational flag).

      @return true if a similarity reduction is possible, false otherwise
   */
  bool HasSimilarityReduction();

  /**
      Price the fixed part of the variance swap.

      Primarily used for conditional variance swaps.  The fixed part requires
      an additional state variable.  It is more efficient to solve the two
      parts of the payoff separately (two low order path dependent problems)
      instead of solving the full problem in one pass (one higher order
      path dependent problem).

      @return VarianceSwapNumOutput object containing all price data
   */
  AutoPtr<VarianceSwapNumOutput> PriceFixedLeg();


  /// The variance swap params being used for pricing
  pricing::VarianceSwapParams& m_params;

  /// The HG model parameters
  Model& m_model;

  /// The computational flags (ie. what data to compute)
  const finance::ComputationalFlags& m_flags;
  
private:
  
  NO_COPY_CLASS(VarianceSwapPricer);

}; // class VarianceSwapPricer


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_VARIANCESWAPPRICER_H_
