/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/hazardratewithtimecomponent.h
// Purpose:     base hazard rate class
// Author:      (z)
// Created:     03/11/04
// RCS-ID:      $Id: hazardratewithtimecomponent.h,v 1.31 2006/08/20 09:47:41 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/ihg/hazardratewithtimecomponent.h
   @brief Base class for hazard rate with a time component
          (HazardRateTimeOnly, HazardRateCombo)
 */

#ifndef _ITO33_IHG_HAZARDRATEWITHTIMECOMPONENT_H_
#define _ITO33_IHG_HAZARDRATEWITHTIMECOMPONENT_H_

#include "ito33/sharedptr.h"
#include "ito33/vector.h"
#include "ito33/date.h"

#include "ito33/ihg/hazardrate.h"

namespace ito33
{

namespace numeric
{
  class PiecewiseConstantFunctor;
}

namespace ihg
{

  class HazardRateVisitor;

/**
    Base class for hazard rate with a time component
    (HazardRateTimeOnly, HazardRateCombo)

    @nocreate
 */
class ITO33_IHG_DLLDECL HazardRateWithTimeComponent : public HazardRate
{

public:

  /**
     @internal
     @brief Dummy ctor is useful for calibration

     @noexport
   */
  HazardRateWithTimeComponent();

  /**
     @internal
     @brief Ctor constructs the time component by using two arrays

     Note that the element pdTimes[nNbTimes - 1] is not needed given the 
     convention used.

     @param pDates the times that split the time into nNbTimes intervals
     @param pdValues the values at each intervals
     @param nNbTimes the number of time intervals

     @noexport
   */
  HazardRateWithTimeComponent(const Date* pDates, const double *pdValues, 
                              size_t nNbTimes);

  /**
     Ctor constructs the time component by using two arrays

     Note that the element dates[nNbTimes - 1] has no need with the actual 
     convention.

     @param dates the times that split the time into intervals
     @param values the values at each intervals
   */
  HazardRateWithTimeComponent(const std::vector<Date>& dates,
                              const std::vector<double>& values);

  /// Dummy virtual dtor for base class
  virtual ~HazardRateWithTimeComponent() { }

  /**
     Gets the time points.

     @return the time points
   */
  const std::vector<Date>& GetDates() const { return m_pDates; }

  /**
     Gets the values of the time component.

     @return the values of the time component
   */
  const std::vector<double>& GetTimeComponentValues() const;

  /**
     @internal
     @brief Resets the time component.

     Helps calibration on the time component. 

     @noexport
   */
  void ResetTimeComponent(const Date* pDates, const double *pdValues, 
                          size_t nNbTimes);

  void GetSpecialTimes(numeric::mesh::SpecialTimes& specialTimes) const;


  /*
    This is an intermediate class. So the following virtual functions
    don't need to be implemented.

    GetHazardRates 
    Dump
    Visit 
    GetModelParameters
   */


protected:
  
  /// time component
  shared_ptr<numeric::PiecewiseConstantFunctor> m_pTimeComponent;

  /**
     Date vector
  
     @internal We need a seperate date vector (of type date) here as the
     the time points is saved to double type in PiecewiseConstantFunctor.
     In ctor, we just store the input date values in this vector. With 
     this vector, implementation of Dump() and GetTime() function is easier.
   */
  std::vector<Date> m_pDates;

}; // class HazardRateWithTimeComponent


} // namespace ihg

} // namespace ito33

#endif // #ifndef _ITO33_IHG_HAZARDRATEWITHTIMECOMPONENT_H_

