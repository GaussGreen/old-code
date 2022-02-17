/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/hazardratedecay.h
// Purpose:     exponentially decaying (space only) hazard rate class
// Author:      David (based on idea in calldef project)
// Created:     03/12/19
// RCS-ID:      $Id: hazardratedecay.h,v 1.20 2005/08/23 12:00:15 zhang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
  @file ito33/ihg/hazardratedecay.h
  @brief exponentially decaying (space only) hazard rate class
  */


#ifndef _ITO33_IHG_HAZARDRATEDECAY_H_
#define _ITO33_IHG_HAZARDRATEDECAY_H_

#include "ito33/ihg/hazardrate.h"


namespace ito33
{

namespace ihg
{
  class HazardRateVisitor;

/**
    @internal 
    @brief Exponentially decaying (space only) hazard rate class 

    The analytic formula is h(S) = alpha * exp(-S/S0)
    ie. It has a value of alpha at S=0, and the value alpha/exp(1) at S=S0
    Typically, S0 will be the current spot, or something close

    This file is meant for testing, and does not neccessarily reflect a 
    market hazard rate function.

    @noexport
*/
class HazardRateDecay : public HazardRate
{
public:
  /**
    usuel constructor 

    @param dAlpha alpha value
    @param dS0 s0 value
    */
  HazardRateDecay(double dAlpha, double dS0)
    : m_dAlpha(dAlpha), m_dS0(dS0) { }

  /// dtor
  virtual ~HazardRateDecay() { }

  
  void GetHazardRates(double dTime, const double *pdS, double *pdValues, 
                      size_t nNumber) const;

  virtual void Dump(XML::Tag& tagParent) const;

  virtual void Visit(HazardRateVisitor& visitor) const;


private:
  
  double m_dAlpha;
  double m_dS0;


};

} // namespace ihg

} // namespace ito33


#endif // #ifndef _ITO33_IHG_HAZARDRATEDECAY_H_

