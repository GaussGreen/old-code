///////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/common/pathdependent/discreteevent.h
// Purpose:     discrete  event
// Author:      ITO 33 
// Created:     04/02/2005
// RCS-ID:      $Id: discreteevent.h,v 1.2 2005/06/27 16:23:26 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
///////////////////////////////////////////////////////////////////////////////

/**
    @file continuous call event
    @brief
*/

#ifndef _ITO33_PRICING_DISCRETEEVENT_H_
#define _ITO33_PRICING_DISCRETEEVENT_H_

#include "ito33/pricing/pathdepevent.h"
 
namespace ito33
{


namespace pricing
{

 class PathDepStructure;


class DiscreteEvent : public PathDepEvent
{

public:

  /** Constructor

      @param dTime time of the payment
      @param dCutoff The cutoff value
  */
  DiscreteEvent(
        const double dStartTime,
        const double dEndTime,
        const double dTime)
    : PathDepEvent(dStartTime,dEndTime,dTime)
  {    
  }


  ~DiscreteEvent() {};
    

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

};

} // namespace pricing

} // namespace ito33 

#endif // #ifndef _ITO33_PRICING_DISCRETEEVENT_H_

