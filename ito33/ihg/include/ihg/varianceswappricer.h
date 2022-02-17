/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/varianceswappricer.h
// Purpose:     variance swap pricer class
// Created:     2006/02/21
// RCS-ID:      $Id: varianceswappricer.h,v 1.4 2006/07/25 21:43:46 dave Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/varianceswappricer.h
    @brief Variance swap pricer class

    Implementation of the pricer class for variance swaps.
 */

#ifndef _IHG_VARIANCESWAPPRICER_H_ 
#define _IHG_VARIANCESWAPPRICER_H_

#include "ito33/vector.h"
#include "ito33/autoptr.h"

namespace ito33
{

namespace finance { class ITO33_DLLDECL ComputationalFlags; }

namespace pricing
{
  class VarianceSwapParams;
  class PathDepEvent;
}

namespace ihg
{
  class Model;
  class VarianceSwapNumOutput;

/**
    Variance swap pricer class.

    In general, variance swaps can be priced by introducing two new state
    variables: previous spot and current average of the squared returns.
    Hence, it is a 3D problem.  A similarity reduction can sometimes be
    used to reduce the dimension by one.  
 */
class VarianceSwapPricer
{
public:

  /**
      Constructor.

      @param params reference to variance swap pricing object
      @param model the ihg model to use
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
      - the model is homogeneous (vol and hr are stock independent)
      - the variance swap does not have any exotic features

      Other checks may also be done (such as the "use similarity reduction"
      computational flag).

      @return true if a similarity reduction is possible, false otherwise
   */
  bool HasSimilarityReduction();

  /// The variance swap params being used for pricing
  pricing::VarianceSwapParams& m_params;

  /// The ihg model parameters
  Model& m_model;

  /// The computational flags (ie. what data to compute)
  const finance::ComputationalFlags& m_flags;
  
 
private:
  
  NO_COPY_CLASS(VarianceSwapPricer);

}; // class VarianceSwapPricer


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_VARIANCESWAPPRICER_H_
