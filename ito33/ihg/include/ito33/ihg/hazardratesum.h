/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/hazardratesum.h
// Purpose:     hazard rate which is the sum of time comp. and space comp.
// Author:      (z)
// Created:     03/11/04
// RCS-ID:      $Id: hazardratesum.h,v 1.15 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////


#ifndef _ITO33_IHG_HAZARDRATESUM_H_
#define _ITO33_IHG_HAZARDRATESUM_H_

#include "ito33/sharedptr.h"

#include "ito33/ihg/hazardratewithtimecomponent.h"
#include "ito33/ihg/spotcomponent.h"

namespace ito33
{

namespace ihg
{

/**
  hazard rate which is the sum of time component and space component
  */
class HazardRateSum : public HazardRateWithTimeComponent
{
public:
  
  /**
    constructor. The time component is defined by two arrays.

    @param pSpotcomponent space component
    @param pDates the times that split the time into nNbTimes intervals
    @param pdValues the values at each intervals
    @param nNbTimes the number of time intervals
   */
  HazardRateSum(shared_ptr<SpotComponent> &pSpotcomponent,
                const Date* pDates, const double *pdValues, size_t nNbTimes)
              : HazardRateWithTimeComponent(pDates, pdValues, nNbTimes),
                m_pSpotComponent(pSpotcomponent) { }

  void GetHazardRates(double dTime, const double *pdS, double *pdValues, 
                      size_t nNumber) const;

private:

  /// space component
  shared_ptr<SpotComponent> m_pSpotComponent;

}; // class HazardRateSum


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_HAZARDRATESUM_H_
