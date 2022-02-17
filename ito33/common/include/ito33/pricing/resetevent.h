/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/resetevent.h
// Purpose:     Reset event
// Author:      Ito33 team
// Created:     05/10/2004
// RCS-ID:      $Id: resetevent.h,v 1.8 2006/08/03 21:25:18 dave Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file  ito33/pricing/resetevent.h
    @brief reset events.
*/

#ifndef _ITO33_PRICING_RESETEVENT_H_
#define _ITO33_PRICING_RESETEVENT_H_

#include "ito33/pricing/pathdepevent.h"
 
namespace ito33
{


namespace pricing
{

 class PathDepStructure;


/**
   Path dependent reset event class.
*/
class ResetEvent : public PathDepEvent
{

public:

  /** Constructor

      @param dTime time of the event
      @param dCapRate cap on the conversion ratio
      @param dFloorRate floor on the conversion ratio
      @param dMultiplier 
      @param dInitialConvPrice Initial conversion price
      @param dCurrentConvPrice Current conversion price
      @param bIsResetFlooredBy_Prevailing indicate whether the new conversion
                                        ratio is computed based on the previous
                                        one or not
      @param bIsLastReset indicate if this event is the last reset event
  */
  ResetEvent (double dTime,
              double dCapRate,
              double dFloorRate,
              double dMultiplier,
              double dInitialConvPrice,
              double dCurrentConvPrice,
              bool   bIsResetFlooredBy_Prevailing,
              bool   bIsLastReset)
    : PathDepEvent(dTime, dTime, dTime),
    m_dCapRate(dCapRate),
    m_dFloorRate(dFloorRate),
    m_dMultiplier(dMultiplier),
    m_dInitialConvPrice(dInitialConvPrice),
    m_dCurrentConvPrice(dCurrentConvPrice),
    m_bIsResetFlooredBy_Prevailing(bIsResetFlooredBy_Prevailing),
    m_bIsLastReset(bIsLastReset)
  {    
  }



  ~ResetEvent () {};
    

protected:

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
    Function to turn off the paths

    @param pathDepStruct the path dependent structure containing prices, grids
  */
  void TurnOffPaths(PathDepStructure& pathDepStruct) const;


  /**
    Function to apply event at time m_dTime, using second method.

    @param pathDepStruct the path dependent structure containing prices, grids
  */
  void ApplySecondMethod(PathDepStructure& pathDepStruct) const;
  
private:

  ///cap for conversion rate
  double m_dCapRate;
 
  ///conversion rate floor
  double m_dFloorRate;

  ///mutliplier
  double m_dMultiplier;

  ///store a copy of the initial conversion price (issue date)
  double m_dInitialConvPrice;

  ///store a copy of the current conversion price (valuation date)
  double m_dCurrentConvPrice;

  ///use the previous or initial conversion ratio for estimating the new ratio
  bool m_bIsResetFlooredBy_Prevailing;

  ///indicate if this is the last reset date
  bool m_bIsLastReset;

};

} // namespace pricing

} // namespace ito33 

#endif // #ifndef _ITO33_PRICING_RESETEVENT_H_
