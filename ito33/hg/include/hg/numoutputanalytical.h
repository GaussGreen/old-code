/////////////////////////////////////////////////////////////////////////////
// Name:      hg/numoutputanalytical.h
// Purpose:   implementation of NumOutput class using analytical formula
// Created:   2006/08/02
// RCS-ID:    $Id: numoutputanalytical.h,v 1.4 2006/08/19 23:46:52 wang Exp $
// Copyright: (c) 2006  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/numoutputanalytical.h
    @brief Numerical output class for HG model using analytical formula
 */

#ifndef _HG_NUMOUTPUTANALYTICAL_H_
#define _HG_NUMOUTPUTANALYTICAL_H_

#include "ito33/vector.h"
#include "ito33/sharedptr.h"

#include "hg/backwardnumoutput.h"

namespace ito33
{

namespace hg
{

  class ModelOutput;

/**
    Numerical output class for HG model using analytical formula.

    Can be used for path dep problem where prices can only be obtained
    at valuation date, so no surface output is supported, and analysis date
    is supported only when it equals to valuation date.
 */
class NumOutputAnalytical : public BackwardNumOutput
{
public:
  
  NumOutputAnalytical(pricing::Params& params);

  // Default dtor is ok

  /**
      Sets the final values, including the final spots.
      It defines also the analysis data at valuation date if required.

      @param pdS the space mesh at maturity
      @param pdPrices the final prices 
   */
  void SetFinalValues(const std::vector<double>& pdS,
                      const std::vector<double>& pdPrices);

  /**
      Sets the price at spot share price.

      @param dPrice price computed by analytical formula
   */
  void SetPrice(double dPrice)
  {
    m_dPrice = dPrice;
  }

  /**
      Sets the delta at spot share price.

      @param dDelta delta computed by analytical formula
   */
  void SetDelta(double dDelta)
  {
    m_dDelta = dDelta;
  }

  /**
      Sets the gamma at spot share price.

      @param dGamma gamma computed by analytical formula
   */
  void SetGamma(double dGamma)
  {
    m_dGamma = dGamma;
  }

  virtual shared_ptr<ModelOutput> GetModelOutput();

private:

  NO_COPY_CLASS(NumOutputAnalytical);

}; // class NumOutputAnalytical


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_NUMOUTPUTANALYTICAL_H_
