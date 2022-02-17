/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/varianceswap.cpp
// Purpose:     contracts class for variance swap pricing (backward)
// Created:     2006/02/21
// RCS-ID:      $Id: varianceswap.cpp,v 1.24 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file pricing/varianceswap.cpp
    @brief implementation of the contract for variance swaps
 */

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"
#include "ito33/dateutils.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/varianceswapterms.h"
#include "ito33/finance/varianceswaplike.h"
#include "ito33/finance/varianceswap.h"
#include "ito33/finance/gammavarianceswap.h"
#include "ito33/finance/conditionalvarianceswap.h"
#include "ito33/finance/optionvarianceswap.h"
#include "ito33/finance/swappayofftype.h"

#include "ito33/pricing/varianceswap.h"
                      

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(pricing::VarianceSwap);
  
namespace pricing
{

VarianceSwap::VarianceSwap(const finance::VarianceSwap& varianceSwap)
: Contract(varianceSwap)
{
  InitCommonTerms(varianceSwap);

  // Get terms specific to a standard payoff
  m_swapPayoffType = finance::SwapPayoff_Standard;
}


VarianceSwap::VarianceSwap(const finance::GammaVarianceSwap& varianceSwap)
: Contract(varianceSwap)
{
  InitCommonTerms(varianceSwap);

  // Get terms specific to a gamma payoff
  m_swapPayoffType = finance::SwapPayoff_Gamma;

  m_dStartSharePrice = varianceSwap.GetStartSharePrice();

  Date valuationDate = varianceSwap.GetSessionData()->GetValuationDate();

  Date startOfSampling = varianceSwap.GetTerms()->GetStartOfSamplingPeriod();

  if ( valuationDate == startOfSampling )
  {
    // start share price must equal previous share price
    m_dStartSharePrice = m_dPreviousSharePrice;
  }

}


VarianceSwap::VarianceSwap(const finance::ConditionalVarianceSwap& varianceSwap)
: Contract(varianceSwap)
{
  InitCommonTerms(varianceSwap);

  // Get terms specific to a conditional payoff
  m_swapPayoffType = finance::SwapPayoff_Conditional;

  Date valuationDate = varianceSwap.GetSessionData()->GetValuationDate();

  Date startOfSampling = varianceSwap.GetTerms()->GetStartOfSamplingPeriod();

  m_iCurrentConditionalCount = 0;
  if ( startOfSampling < valuationDate )
    m_iCurrentConditionalCount = varianceSwap.GetCurrentConditionalCount();
}


VarianceSwap::VarianceSwap(const finance::OptionVarianceSwap& varianceSwap)
: Contract(varianceSwap)
{
  InitCommonTerms(varianceSwap);

  // Get terms specific to an option payoff
  if ( varianceSwap.GetOptionType() == finance::Option_Call )
    m_swapPayoffType = finance::SwapPayoff_Call;
  else
    m_swapPayoffType = finance::SwapPayoff_Put;
}


void VarianceSwap::InitCommonTerms(const finance::VarianceSwapLike& varianceSwap)
{
  // Copy values from the financial object
  m_dVolatilityStrike = varianceSwap.GetVolatilityStrike();
  
  const finance::VarianceSwapTerms& terms( *varianceSwap.GetTerms() );

  m_swapType = terms.GetSwapType();

  m_returnType = terms.GetReturnType();

  Date startOfSamplingPeriod = terms.GetStartOfSamplingPeriod();

  m_dStartTimeOfSamplingPeriod = GetDoubleFrom( startOfSamplingPeriod );

  m_nNbSamplingReturns = terms.GetNbSamplingReturns();  
  
  m_dCapMultiplier = terms.GetCapMultiplier();

  m_dUpCorridorBarrier = terms.GetUpCorridorBarrier();

  m_dDownCorridorBarrier = terms.GetDownCorridorBarrier();

  m_nAnnualReturnFrequency = terms.GetAnnualReturnFrequency();

  m_dPreviousSharePrice = 
    varianceSwap.GetSessionData()->GetPreviousSharePrice();

  // Check for data that depends if the variance swap if forward starting 
  // or not.
  Date valuationDate = varianceSwap.GetSessionData()->GetValuationDate();

  if ( valuationDate < startOfSamplingPeriod )
  {
    // Set the previous spot to S
    m_dPreviousSharePrice = varianceSwap.GetSessionData()->GetSpotSharePrice();
  }

  // Initialize assuming that valuation is before the start of the
  // sampling period.  Update the variables if otherwise.
  m_dCurrentAvgSqrReturn = 0.0;
  m_nNbSamplesUsed = 0;

  if ( startOfSamplingPeriod < valuationDate )
  {
    m_dCurrentVolatility = varianceSwap.GetCurrentVolatility();

    m_nNbSamplesUsed = varianceSwap.GetNbSamplesUsed();

    // Convert from volatility to the average of the squared returns,
    // which is the state variable for pricing
    m_dCurrentAvgSqrReturn = m_dCurrentVolatility * m_dCurrentVolatility;    

    // Remove the annualization factor which appears in the volatility
    // definition
    m_dCurrentAvgSqrReturn *= 1. / double(m_nAnnualReturnFrequency);    
  }

}

double VarianceSwap::GetPayoffValue(double dZ) const
{

  double dValue = 0.0;
  double dVariance = double(m_nAnnualReturnFrequency) * dZ;

  // For conditional swaps the fixed leg is valued separately, so should not 
  // be included in the payoff
  double dFixed = m_dVolatilityStrike;
  if ( m_swapPayoffType == finance::SwapPayoff_Conditional )
    dFixed = 0.0;

  switch ( m_swapType )
  {
  case finance::Swap_Variance:
    {
      if ( m_dCapMultiplier > 0 )
      {
        double dCap = m_dCapMultiplier * m_dCapMultiplier 
                    * m_dVolatilityStrike * m_dVolatilityStrike;
        dVariance = std::min( dCap, dVariance );
      }
      
      dValue = dVariance - dFixed * dFixed;
      break;
    }
  case finance::Swap_Volatility:  
    {
      double dVol = sqrt(dVariance);
      
      if ( m_dCapMultiplier > 0 )
      {
        double dCap = m_dCapMultiplier * m_dVolatilityStrike;
        dVol = std::min( dCap, dVol );
      }

      dValue = dVol - dFixed;
      break;
    }
        
  default: 
    FAIL("Payoff cannot be made for the given swap type.");   
  }

  return dValue;

} //GetPayoffValue(double dX)

} // namespace pricing

} // namespace ito33
