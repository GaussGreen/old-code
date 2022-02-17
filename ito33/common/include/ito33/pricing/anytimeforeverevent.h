/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/anytimeforever.h
// Purpose:     any time forever event
// Author:      Yann and David
// Created:     18/08/2004
// RCS-ID:      $Id: anytimeforeverevent.h,v 1.7 2004/10/05 09:13:38 pedro Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/anytimeforeverevent.h
    @brief anytime forever event
*/

#ifndef _ITO33_PRICING_ANYTIMEFOREVEREVENT_H_
#define _ITO33_PRICING_ANYTIMEFOREVEREVENT_H_

#include "ito33/pricing/pathdepevent.h"
 
namespace ito33
{


namespace pricing
{

 class PathDepStructure;
 class ConversionProvisions;


class AnyTimeForeverEvent : public PathDepEvent
{

public:

  /** Constructor

      @param dTime time of the payment
      @param dCutoff The cutoff value
  */
  AnyTimeForeverEvent (
        const double dStartTime,
        const double dEndTime,
        const double dTime, 
        const double dTriggerRate,
        const ConversionProvisions* pConversionProvisions,
        const size_t nCallDays = 0)
    : PathDepEvent(dStartTime,dEndTime,dTime), 
      m_dTriggerRate(dTriggerRate),
      m_pConversionProvisions(pConversionProvisions),
      m_nCallPaths(nCallDays + 1)
  {    
  }


  ~AnyTimeForeverEvent () {};
    
  /**
    Get the trigger rate

    @return trigger rate
   */
  double GetTriggerRate()const
  { 
    return m_dTriggerRate; 
  }

protected:

  /// The trigger rate.
  double m_dTriggerRate;

  /// The conversion provisions used to figure out the actual trigger level
  const ConversionProvisions* m_pConversionProvisions;

  /// The number of call paths ("days above plus 1")
  size_t m_nCallPaths;

  /**
    Function to apply event at the start time

    @param pathDepStruct the path dependent structure containing prices, grids
  */
  void ApplyAtStartTime(PathDepStructure& pathDepStruct) const;

  /**
    Function to apply event End Time

    @param pathDepStruct the path dependent structure containing prices, grids
  */
  void ApplyAtEndTime(PathDepStructure& pathDepStruct) const;

  /**
    Function to apply event at time m_dTime

    @param pathDepStruct the path dependent structure containing prices, grids
  */
  void ApplyAtTime(PathDepStructure& pathDepStruct) const;

};

} // namespace pricing

} // namespace ito33 

#endif // #ifndef _ITO33_PRICING_ANYTIMEFOREVEREVENT_H_
