/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/asianoption.cpp
// Author:      ITO 33 Canada
// Created:     March 31, 2005
// RCS-ID:      $Id: asianoption.cpp,v 1.8 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file pricing/asianoption.cpp
   @brief implementation of the contract for asian option
 */

#include "ito33/autoptr.h"
#include "ito33/useexception.h"
#include "ito33/dateutils.h"

#include "ito33/finance/sessiondata.h"

#include "ito33/pricing/asianoption.h"

extern const ito33::Error ITO33_UNEXPECTED;

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(pricing::AsianOption); 

namespace pricing
{
  

AsianOption::AsianOption(const finance::AsianOption & asianOption):
  Option(asianOption)
{
  double dValuationTime = 
    GetDoubleFrom( asianOption.GetSessionData()->GetValuationDate() );
  
  m_dAverageEndTime = GetDoubleFrom( asianOption.GetAverageEndDate() );

  m_dAverageStartTime = GetDoubleFrom( asianOption.GetAverageStartDate() );

  m_nNumberOfSamplingAverages  = asianOption.GetNbSamplingAverages();

  m_nNbSamplesUsed = asianOption.GetNbSamplesUsed();

  m_dCurrentAverage = asianOption.GetCurrentAverage();
 
  m_dStrike = asianOption.GetFixedStrike();

  if ( dValuationTime <=  m_dAverageStartTime )
  {
   // Current average can be set However
   // if valuation date is before the average start time
   // it does not make sense for the current average to be different
   // from the stock price
    m_dCurrentAverage = asianOption.GetSessionData()->GetSpotSharePrice();

    // set the number of sample used to be zero.
    // Otherwise the event are going to be wrong.
    m_nNbSamplesUsed = 0;

  }
  else
  {
    if ( m_dCurrentAverage == 0. ) 
      throw EXCEPTION_MSG
        (
          ITO33_UNEXPECTED,
          "Current Average must be specified."
        );
  }

  switch ( asianOption.GetOptionType() )
  { 
  case finance::Option_Call:
    if ( HasFixedStrike() )
      m_optionType = pricing::AsianOption_FixedStrikeCall;
    else
      m_optionType = pricing::AsianOption_FloatingStrikeCall;  
    break;

  case finance::Option_Put:
    if ( HasFixedStrike() )  
      m_optionType = pricing::AsianOption_FixedStrikePut;
    else
      m_optionType = pricing::AsianOption_FloatingStrikePut;
    break;

  default:
    FAIL("No such Asian option type.");      
  } //end switch statement

}


} // namespace pricing

} // namespace ito33
