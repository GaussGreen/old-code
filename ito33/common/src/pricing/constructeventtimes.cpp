/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/constructeventtimes.cpp
// Purpose:     Helper fnction to construct event times
// Author:      Ito33 team Canada
// Created:     May 22, 2006
// RCS-ID:      $Id: constructeventtimes.cpp,v 1.3 2006/07/04 18:23:36 dave Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/numeric/predicatetime.h"

#include "ito33/pricing/constructeventtimes.h"


namespace ito33
{

namespace pricing
{

 
void ConstructSamplingEventTimes(double dValuationTime, 
                             double dStartSamplingTime,
                             double dStopSamplingTime,
                             size_t nNbSampling,
                             size_t nNbSamplesUsed,
                             std::vector<double> &pdEventTimes,
                             std::vector<size_t> &pnSamplingNumbers)
{
  // Observations at: t_i = T_start + i/N (T_stop - T_start)
  // for i = 1..N.  T_start is the start of the sampling period,
  // and is not a return calculation date
  size_t nSamplingNumber = 1;
  size_t nNbSamplesRemaining = nNbSampling;

  // Take into account when the valuation date is
  // past the stop sampling date
  if ( numeric::IsEqualOrAfter(dValuationTime, dStopSamplingTime) )
    return;

  // Adjust if already past the sampling start date
  if ( numeric::IsAfter(dValuationTime, dStartSamplingTime) )
  {
    dStartSamplingTime  = dValuationTime;
    nNbSamplesRemaining = nNbSampling - nNbSamplesUsed;
    nSamplingNumber     = nNbSamplesUsed + 1;
  }

  // Calculate what will be a "trading day"
  double dStep = (dStopSamplingTime-dStartSamplingTime) / nNbSamplesRemaining;

  // Create nNbSamplesRemaining events, starting at the trading day past
  // the start time
  double dEventTime = dStartSamplingTime + dStep;  

  // Add the last event manually to avoid cumulative round-off error
  while ( nSamplingNumber < nNbSampling )
  {
    pdEventTimes.push_back( dEventTime );
    pnSamplingNumbers.push_back( nSamplingNumber );

    // Update for next event
    nSamplingNumber++;
    dEventTime += dStep;
  } // end while loop creating events

  pdEventTimes.push_back( dStopSamplingTime );
  pnSamplingNumbers.push_back( nNbSampling );

}
} // namespace pricing

} // namespace ito33
