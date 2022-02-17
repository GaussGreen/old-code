/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/conditionalswapevent.h
// Purpose:     conditional variance swap event
// Created:     2006/07/14
// RCS-ID:      $Id: conditionalswapevent.h,v 1.2 2006/08/03 21:25:18 dave Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file  ito33/pricing/conditionalswapevent.h
    @brief conditional variance swap event.
 */

#ifndef _ITO33_PRICING_CONDITIONALSWAPEVENT_H_
#define _ITO33_PRICING_CONDITIONALSWAPEVENT_H_

#include "ito33/vector.h"

#include "ito33/pricing/pathdepevent.h"
#include "ito33/pricing/varianceswap.h"
 
namespace ito33
{

namespace pricing
{

class PathDepStructure;

class ConditionalSwapEvent : public PathDepEvent
{

public:

  /** 
      Constructor.

      @param dTime time of the event
      @param bIsLastEvent indicate if this is the last event
      @param varianceSwap the underlying variance swap from which we get
             corridors, payoff type, etc
   */
  ConditionalSwapEvent(double dTime, bool bIsLastEvent, VarianceSwap& varianceSwap)
    : PathDepEvent(dTime, dTime, dTime), 
      m_bIsLastEvent(bIsLastEvent),
      m_varianceSwap(varianceSwap)
  {    
  }

  virtual ~ConditionalSwapEvent () {}
    

protected:

  virtual void ApplyAtStartTime(PathDepStructure& ) const {}

  virtual void ApplyAtEndTime(PathDepStructure& ) const {}

  virtual void ApplyAtTime(PathDepStructure& pathDepStruct) const;


protected:

  /// Indicates is this is the last event
  bool m_bIsLastEvent;

  /// Indicates if a recursive call to ApplyAtTime has been made
  mutable bool m_bRecurse;

  /// Reference to underlying variance swap contract
  VarianceSwap& m_varianceSwap;

private:

  NO_COPY_CLASS(ConditionalSwapEvent);

};

} // namespace pricing

} // namespace ito33 

#endif // #ifndef _ITO33_PRICING_CONDITIONALSWAPEVENT_H_
