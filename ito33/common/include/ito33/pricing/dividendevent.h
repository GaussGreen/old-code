/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/dividendevent.h
// Purpose:     Base dividend event class
// Created:     2005/06/02
// RCS-ID:      $Id: dividendevent.h,v 1.5 2005/06/29 20:02:35 dave Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/pricing/dividendevent.h
   @brief The declaration of base dividend event class.
 */

#ifndef _ITO33_PRICING_DIVIDENDEVENT_H_
#define _ITO33_PRICING_DIVIDENDEVENT_H_

#include "ito33/pricing/event.h"
#include "ito33/numeric/extrapolationmode.h"
#include "ito33/numeric/interpolation.h"
 
namespace ito33
{

  namespace numeric
  {
    class InterpolationMatrix;
  }

namespace pricing
{


/// base class for dividend events
class DividendEvent : public Event
{

public:

  /** 
     Default ctor sets the event type to ET_Dividend

     @param dTime time of the dividend
   */
  DividendEvent(double dTime) 
  : Event(dTime),
    m_emLeft(numeric::ExtrapolationMode_Linear),
    m_emRight(numeric::ExtrapolationMode_Linear),
    m_interpolationMethod(numeric::InterpolationMethod_Quadratic)
  {
    m_eventType = ET_Dividend;
  }

  virtual void 
  ApplyToSpots(const double* pdS, double* pdNewS, size_t nNbS) const = 0;

  /** 
      Apply the dividend to the prices.
  
      Default implementation is for backward.

      @param pdS Pointer to the grid
      @param pdValues Pointer to the current colution values
      @param nNbS Number of grid points/solution values
   */
  virtual void 
  ApplyToPrice(const double *pdS, double* pdValues, size_t nNbS) const;

  /** 
      Apply the proportional dividend to Greek values

      @param pdS Pointer to the grid
      @param pdValues Pointer to the current colution values
      @param nNbS Number of grid points/solution values
   */
  void 
  ApplyToGreek(const double *pdS, double* pdValues, size_t nNbS) const
  {
    ApplyToPrice(pdS, pdValues, nNbS);
  }

  /**
     Gets the interpolation matrix. Dividend can be translated into an
     interpolation.

     Default inplementation is for backward
   */
  virtual numeric::InterpolationMatrix* 
  GetInterpolationMatrix
  (const double* pdS, size_t nNbS, size_t nNbSubSystem) const;

  /**
    Set the left extrapolation mode.

    @param emLeft extrapolation method to the left
   */
  void SetEMLeft(numeric::ExtrapolationMode emLeft)
  {
    m_emLeft = emLeft;
  }

  /**
    Set the right extrapolation mode.

    @param emRight extrapolation method to the right
   */
  void SetEMRight(numeric::ExtrapolationMode emRight)
  {
    m_emRight = emRight;
  }

  /**
    Set the interpolation method.

    @param interpolationMethod interpolation method
   */
  void SetInterpolationMethod(numeric::InterpolationMethod interpolationMethod)
  {
    m_interpolationMethod = interpolationMethod;
  }


protected:

  /// Method used to extrapolate to the left (default to linear)
  numeric::ExtrapolationMode m_emLeft;

  /// Method used to extrapolate to the right (default to linear)
  numeric::ExtrapolationMode m_emRight;

  /// Method used to interpolate (default to quadratic)
  numeric::InterpolationMethod m_interpolationMethod;

};


} // namespace pricing

} // namespace ito33 

#endif // #ifndef _ITO33_PRICING_DIVIDENDEVENT_H_
