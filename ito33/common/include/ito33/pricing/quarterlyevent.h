/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/quarterlyevent.h
// Purpose:     quartely discrete event
// Author:      Yann and David
// Created:     18/08/2004
// RCS-ID:      $Id: quarterlyevent.h,v 1.10 2006/08/03 21:25:18 dave Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file  ito33/pricing/quarterlyevent.h
    @brief Coco end-of-quarter forever.
*/

#ifndef _ITO33_PRICING_QUARTERLYEVENT_H_
#define _ITO33_PRICING_QUARTERLYEVENT_H_

#include "ito33/pricing/pathdepevent.h"
 
namespace ito33
{


namespace pricing
{

 class PathDepStructure;
 class ConversionProvisions;


class QuarterlyEvent : public PathDepEvent
{

public:

  /** Constructor

      @param dTime time of the event
      @param bInitializeOnly indicates if the event is initialized only
      @param dTriggerRate is the trigger rate (multipled by conversion strike)
      @param pConversionProvisions is needed to get the conversion strike
      @param nCallDays is the number of call window period days
  */
  QuarterlyEvent (
        const double dTime,
        const bool bInitializeOnly,
        const double dTriggerRate,
        const ConversionProvisions* pConversionProvisions,
        const size_t nCallDays = 0)
    : PathDepEvent(dTime, dTime, dTime), 
      m_bInitializeOnly(bInitializeOnly),
      m_dTriggerRate(dTriggerRate),
      m_pConversionProvisions(pConversionProvisions),
      m_nCallPaths(nCallDays + 1)
  {    
  }

  ~QuarterlyEvent () {};
    
  /**
    Get the trigger rate.

    @return trigger rate
   */
  double GetTriggerRate()const
  { 
    return m_dTriggerRate; 
  }


protected:

  /// Indicate if the event should just be initialized, and not applied
  // In this case, copy from base paths to quaterly paths
  bool m_bInitializeOnly;

  /// The trigger rate
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

  /**
    Function to initialize the paths by copying from base to quarterly paths

    @param pathDepStruct the path dependent structure containing prices, grids
  */
  void Initialize(PathDepStructure& pathDepStruct) const;

  /**
    Function to turn off the paths

    @param pathDepStruct the path dependent structure containing prices, grids
  */
  void TurnOffPaths(PathDepStructure& pathDepStruct) const;




};

} // namespace pricing

} // namespace ito33 

#endif // #ifndef _ITO33_PRICING_QUARTELYEVENT_H_
