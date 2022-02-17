/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/hazardratemultiply.h
// Purpose:     hazard rate which is the multipl. of time comp. and space comp.
// Author:      Wang
// Created:     2004/06/04
// RCS-ID:      $Id: hazardratemultiply.h,v 1.6 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////


#ifndef _ITO33_IHG_HAZARDRATEMULTIPLY_H_
#define _ITO33_IHG_HAZARDRATEMULTIPLY_H_

#include "ito33/sharedptr.h"

#include "ito33/ihg/hazardratewithtimecomponent.h"
#include "ito33/ihg/spotcomponent.h"


namespace ito33
{

namespace ihg
{


/**
  hazard rate which is the multiplication of time component
  and space component.
  */
class HazardRateMultiply : public HazardRateWithTimeComponent
{
public:
  
  /**
    constructor. The time component is defined by two arrays.

    @param pSpotcomponent space component
    @param pDates the times that split the time into nNbTimes intervals
    @param pdValues the values at each intervals
    @param nNbTimes the number of time intervals
   */
  HazardRateMultiply(shared_ptr<SpotComponent> &pSpotcomponent,
                     const Date* pDates, 
                     const double *pdValues, 
                     size_t nNbTimes);

  void GetHazardRates(double dTime, const double *pdS, double *pdValues, 
                      size_t nNumber) const;

private:

  shared_ptr<SpotComponent> m_pSpotComponent;

}; // class HazardRateMultiply


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_HAZARDRATEMULTIPLY_H_
