/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/reseteventoned.h
// Purpose:     the dividend events 
// Author:      David Pooley
// Created:     2004/10/24
// RCS-ID:      $Id: reseteventoned.h,v 1.3 2005/03/30 11:23:36 wang Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/reseteventoned.h
    @brief The declaration of reset classes.

*/

#ifndef _ITO33_PRICING_RESETEVENT1D_H_
#define _ITO33_PRICING_RESETEVENT1D_H_

#include "ito33/pricing/event.h"
 
namespace ito33
{

namespace pricing
{


/** 
    Reset class

*/
class ResetEvent1D : public Event
{

public:

  /** Constructor

      @param dEventTime time of the event
      @param dConversionRatioCapRate cap on the conversion ration
      @param dConversionRatioFloorRate floor on the conversion ratio
      @param dMultiplier ????
      @param dNominal nominal of the bond
   */
  ResetEvent1D(
             double dEventTime, 
             double dConversionRatioCapRate,
             double dConversionRatioFloorRate,
             double dMultiplier,
             double dNominal) :
     Event(dEventTime),
     m_dConversionRatioCapRate(dConversionRatioCapRate),
     m_dConversionRatioFloorRate(dConversionRatioFloorRate),
     m_dMultiplier(dMultiplier),
     m_dNominal(dNominal)
  {
     m_eventType = ET_Reset;

     /// reset event must be applied at t-
     m_IsAppliedAfterConstraints = true;
  }

  
  /** 
      Apply event to the prices

      @param pdS Pointer to the grid
      @param pdValues Pointer to the current colution values
      @param nNbS number of grid points/solution values
  */
  void ApplyToPrice(const double *pdS, double* pdValues, size_t nNbS) const;

  /** 
      Apply the event to Greek values

      @param pdS Pointer to the grid
      @param pdValues Pointer to the current colution values
      @param nNbS number of grid points/solution values
  */
  void ApplyToGreek(const double *pdS, double* pdValues, size_t nNbS) const
  {
    ApplyToPrice(pdS, pdValues, nNbS);
  }

protected:

  ///cap for conversion rate
  double m_dConversionRatioCapRate;
 
  ///covnersion rate floor
  double m_dConversionRatioFloorRate;

  ///mutliplier
  double m_dMultiplier;

  ///nominal   
  double m_dNominal;

};

}// namespace pricing

}//namespace ito33

#endif // #ifndef _ITO33_PRICING_RESETEVENT1D_H_


