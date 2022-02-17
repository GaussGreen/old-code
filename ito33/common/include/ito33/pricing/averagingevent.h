/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/averagingevent.h
// Purpose:     Averaging Event
// Author:      Ito33 Canada
// Created:     April 11, 2005
// RCS-ID:      $Id: averagingevent.h,v 1.9 2006/08/03 21:25:18 dave Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file  ito33/pricing/averagingevent.h
    @brief Averaging event.
 */

#ifndef _ITO33_PRICING_AVERAGINGEVENT_H_
#define _ITO33_PRICING_AVERAGINGEVENT_H_

#include "ito33/pricing/pathdepevent.h"
 
namespace ito33
{

namespace pricing
{

 class PathDepStructure;

/**
    Path dependent averaging event class.

    Used by Asian options, for example.
 */
class AveragingEvent : public PathDepEvent
{

public:

  /** 
      Constructor.

      @param dTime time of the event
      @param nObservation Observation number
      @param bIsLastEvent indicate if this is the last event
      @param bIsSimilarityReduction indicate if a similarity reduction is
             possible
   */
  AveragingEvent (double dTime, size_t nObservation, bool bIsLastEvent,
                  bool bHasSimilarityReduction=false)
    : PathDepEvent(dTime, dTime, dTime), m_nObservation(nObservation),
    m_bIsLastEvent(bIsLastEvent),
    m_bHasSimilarityReduction(bHasSimilarityReduction),
    m_bRecurse(false)
  {    
  }

  virtual ~AveragingEvent () {}
    
protected:

  ///  Function to apply event at the start time
  void ApplyAtStartTime(PathDepStructure& ) const {}

  /// Function to apply event End Time
  void ApplyAtEndTime(PathDepStructure& ) const {}

  /// Function to apply event at time m_dTime
  void ApplyAtTime(PathDepStructure& pathDepStruct) const;


  /**
      Given current spot and average, compute new average.

      @param dA current average
      @param dS current spot

      @return new computed average
   */
  virtual double GetNewAverage(double dA, double dS) const; 

  /**
      Interpolates on the upstream point. 
      
      The interpolation can be either linear or quadratic.

      @param dY the new average
      @param dS the share price
      @param pathDepStruct the path dependent structure
      @param the average grid
   */
  double Interpolate(double dY, double dS, 
    PathDepStructure& pathDepStruct, std::vector<double> &pdGridY) const;

protected:

  /// Observation number
  size_t m_nObservation;

  /// Indicate is this is the last event
  bool m_bIsLastEvent;

  /// Indicate if a similarity reduction is possible
  bool m_bHasSimilarityReduction;

  /// Indicate if a recursive call to ApplyAtTime has been made
  mutable bool m_bRecurse;

};

} // namespace pricing

} // namespace ito33 

#endif // #ifndef _ITO33_PRICING_ASIANEVENT_H_
