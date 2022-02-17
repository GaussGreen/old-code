/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/hazardratetimeonly.h
// Purpose:     time only hazard rate class
// Author:      ZHANG Yunzhi
// Created:     03/11/14
// RCS-ID:      $Id: hazardratetimeonly.h,v 1.43 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2003 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/hazardratetimeonly.h
    @brief time only hazard rate class
 */

#ifndef _ITO33_IHG_HAZARDRATETIMEONLY_H_
#define _ITO33_IHG_HAZARDRATETIMEONLY_H_

#include "ito33/date.h"

#include "ito33/ihg/hazardratewithtimecomponent.h"

namespace ito33
{

namespace ihg
{

  class HazardRateVisitor;

/// Class for a time only, piecewise constant hazard rate
class ITO33_IHG_DLLDECL HazardRateTimeOnly : public HazardRateWithTimeComponent
{
  
public:

  /**
      @internal
      @brief Dummy ctor is useful for calibration

      @noexport
   */
  HazardRateTimeOnly() : HazardRateWithTimeComponent() {}

  /**
      @internal 
      @brief Ctor constructs from two arrays
     
      @param pDates the times that split the tike into nNbTimes intervals
      @param pdValues the values at each intervals
      @param nNbTimes the number of time intervals

      @noexport
   */
  HazardRateTimeOnly(const Date* pDates, const double *pdValues, 
                     size_t nNbTimes);

  HazardRateTimeOnly(const std::vector<Date>& dates,
                     const std::vector<double>& values);
  
  // Default dtor is ok 
  
  /**
      Gets all the values.

      @return values for all the time intervals
   */
  const std::vector<double>& GetValues() const 
  {
    return GetTimeComponentValues();
  } 
  
  /**
      @internal
      @brief Gets the the hazard rate at a given time.

      @param dTime the time at which we the hazard rate value is evaluated
      @return the evaluated hazard rate at given time

      @noexport
   */
  virtual double GetValueAtTime(double dTime) const;

  /** 
      Returns a perturbed hazard rate

      To calculate the sensibilities with respect to the "spread" curves by finite 
      difference, we must be able to shift the curve by some 
      (small) amount. This method must be implemented by the derived classes
      to do this.

      Note that it doesn't modify this object at all but rather returns a new,
      perturbed, copy of it.

      @param dShift the perturbation parameter, its meaning varies for
                    different derived classes

      @return a new perturbed hazard rate

      @noexport
   */
  virtual shared_ptr<HazardRate> Perturb(double dShift);

  //

  virtual bool IsTimeOnly() const { return true; }
  
  virtual void 
  GetHazardRates(double dTime, const double *pdS, double *pdValues, 
                 size_t nNbS) const;


  virtual void Dump(ito33::XML::Tag& tagParent) const;

  virtual void Visit(ito33::ihg::HazardRateVisitor& visitor) const;

  virtual void
  GetModelParameters(finance::ModelParametersConsumer& visitor) const;

}; // class HazardRateTimeOnly


} // namespace ihg

} // namespace ito33


#endif // #ifndef _ITO33_IHG_HAZARDRATETIMEONLY_H_

