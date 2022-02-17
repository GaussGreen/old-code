/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/constructeventtimes.h
// Purpose:     Helper functions to create event times
// Author:      Ito33 Canada
// Created:     May 22, 2006
// RCS-ID:      $Id: constructeventtimes.h,v 1.1 2006/05/23 17:20:49 yann Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file  ito33/pricing/constructeventtimes.h
    @brief Helper to declare function to create event times.
    It is possible that the same event schedule can be used
    for multiple contracts. For example the variance swap
    and averaging event construction are the same. This
    prevents also code duplication.
 */

#ifndef _ITO33_PRICING_CONSTRUCT_EVENT_TIMES_H_
#define _ITO33_PRICING_CONSTRUCT_EVENT_TIMES_H_

#include "ito33/vector.h"
#include "ito33/date.h"
 
namespace ito33
{

namespace pricing
{

/**  
    Construct the sampling event times.

    @param dValuationTime valuation time
    @param dStartSamplingTime start of sampling period
    @param dStopSamplingTime stop of sampling period
    @param nNbSampling number of sampling("observation") points
    @param nNbSamplesUsed number of sampling already used

    @param pdEventTimes array of event times
    @param pnSamplingNumbers array of observation numbers
*/
void ConstructSamplingEventTimes(double dValuationTime, 
                             double dStartSamplingTime,
                             double dStopSamplingTime,
                             size_t nNbSampling,
                             size_t nNbSamplesUsed,
                             std::vector<double> &pdEventTimes,
                             std::vector<size_t> &pnSamplingNumbers);


} // namespace pricing

} // namespace ito33 

#endif // #ifndef _ITO33_PRICING_CONSTRUCT_EVENT_TIMES_H_
