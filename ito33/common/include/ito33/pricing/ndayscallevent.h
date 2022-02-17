/////////////////////////////////////////////////////////////////////////////
// Name:        common/include/ito33/pricing/ndayscallevent.h
// Purpose:     instrument can be called at time 
//              t only when the stock has been traded 
//              above trigger during N days prior to t.
//              Copy the solution from above if S is above trigger
// Author:      Yann and David
// Created:     08/09/2004
// RCS-ID:      $Id: ndayscallevent.h,v 1.8 2006/03/22 13:10:23 yann Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ndayscallevent.h
    @brief instrument can be called at time 
           t only when the stock has been traded 
           above trigger during N days prior to t.
           Copy the solution from above if S is above trigger
*/

#ifndef _ITO33_PRICING_NDAYSCALLEVENT_H_
#define _ITO33_PRICING_NDAYSCALLEVENT_H_

#include "ito33/pricing/pathdepevent.h"
 
namespace ito33
{

namespace pricing
{
 class PathDepStructure;
 class ConversionProvisions;

/// Declaration of the NDaysCallEvent class
class NDaysCallEvent : public PathDepEvent
{

public:

  /** Constructor
      @param dStartTime beginning of the event (forward in time)
      @param dEndTime   end of event (forward in time) 
      @param dTime    time of the payment
      @param dTrigger trigger value
  */
  NDaysCallEvent(
        const double dStartTime,
        const double dEndTime,
        const double dTime, 
        const double dTriggerRate,
        const ConversionProvisions *pConversionProvisions,
        const size_t nCallDays,
        const size_t nPathToSave = 0)
    : PathDepEvent(dStartTime,dEndTime,dTime),
      m_dTriggerRate(dTriggerRate),
      m_pConversionProvisions(pConversionProvisions),
      m_nPathToSave(nPathToSave),
      m_nCallPaths(nCallDays + 1)
  {  }


  ~NDaysCallEvent () {}
    
  /**
    Function to get the trigger value.

    @return trigger value
   */
  double GetTriggerRate()const
  { 
    return m_dTriggerRate; 
  }

  /**
    This event is partially continuous and partially discrete.
    It is continuous because we want to do something special at start
    and end times.
    It is discrete since event application is discrete (end of each day).

    Indicate this by returning false from this function, but
    maintaining different start and end times.  Event application
    looks at start and end times and calls functions appropriately.
   */
  bool IsContinuous() const
  {
     return false;
  }


protected:

  /// Path to save, usually the trigger history
  size_t m_nPathToSave;

  /// The trigger value.
  double m_dTriggerRate;

  /// The number of days above trigger plus 1
  size_t m_nCallPaths;

  /// The call provisions used to figure out the actual trigger level
  const ConversionProvisions* m_pConversionProvisions;

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

#endif // #ifndef _ITO33_PRICING_NDAYSCALLEVENT_H_
