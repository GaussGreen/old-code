/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/sharedependentconversionevent.h
// Purpose:     Share dependent conversion event
// Author:      Ito33 team
// Created:     2005/01/14
// RCS-ID:      $Id: sharedependentconversionevent.h,v 1.6 2006/08/03 21:25:18 dave Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file  ito33/pricing/sharedependentconversionevent.h
    @brief reset events.
*/

#ifndef _ITO33_PRICING_SHAREDEPENDENTCONVERSIONEVENT_H_
#define _ITO33_PRICING_SHAREDEPENDENTCONVERSIONEVENT_H_

#include "ito33/pricing/pathdepevent.h"
#include "ito33/pricing/sharedependentconversion.h"

namespace ito33
{

namespace pricing
{

 class PathDepStructure;


/**
   Path dependent event class for share dependent conversion contracts.
*/
class ShareDependentConversionEvent : public PathDepEvent
{

public:

  /** Constructor

      @param dTime time of the event
      @param pShareDependentConversion is needed to get the conversion price

  */
  ShareDependentConversionEvent(
              double dTime,
              const ShareDependentConversion* pShareDependentConversion)
    : PathDepEvent(dTime, dTime, dTime),
     m_pShareDependentConversion(pShareDependentConversion)
  {    
    m_dBaseRatio = m_pShareDependentConversion->GetBaseRatio();
    m_dCapRatio  = m_pShareDependentConversion->GetCapRatio();
    m_dIncrementalShareFactor = 
      m_pShareDependentConversion->GetIncrementalShareFactor();
  }

  //dtor
  ~ShareDependentConversionEvent () {};
    
protected:

  /**
    Function to apply event at the start time

    @param pathDepStruct the path dependent structure containing prices, grids
  */
  void ApplyAtStartTime(PathDepStructure& ) const { }

  /**
    Function to apply event End Time

    @param pathDepStruct the path dependent structure containing prices, grids
  */
  void ApplyAtEndTime(PathDepStructure& ) const { }

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
  
private:

  /// base ratio
  double m_dBaseRatio;

  /// cap ratio
  double m_dCapRatio;

  /// incremental share factor
  double m_dIncrementalShareFactor;

  /// The conversion provisions used to figure out the conversion price
  const ShareDependentConversion* m_pShareDependentConversion;

};

} // namespace pricing

} // namespace ito33 

#endif // _ITO33_PRICING_SHAREDEPENDENTCONVERSIONEVENT_H_
