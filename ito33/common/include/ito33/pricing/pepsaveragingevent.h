/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/pepsaveragingevent.h
// Purpose:     PEPS Averaging Event
// Author:      Ito33 Canada
// Created:     June 27, 2005
// RCS-ID:      $Id: pepsaveragingevent.h,v 1.3 2006/08/03 21:25:18 dave Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file  ito33/pricing/pepsaveragingevent.h
    @brief Averaging event.
*/

#ifndef _ITO33_PRICING_PEPS_AVERAGINGEVENT_H_
#define _ITO33_PRICING_PEPS_AVERAGINGEVENT_H_

#include "ito33/pricing/averagingevent.h"
#include "ito33/pricing/mandatoryconversion.h"
 
namespace ito33
{


namespace pricing
{

 class PathDepStructure;


 class PEPSAveragingEvent : public AveragingEvent
{

public:

  PEPSAveragingEvent(double dTime, size_t nObservation, bool bIsLastEvent,
    const MandatoryConversion *pConversion, 
    bool bIsStockAveraging)
    :AveragingEvent(dTime, nObservation, bIsLastEvent), 
    m_bIsStockAveraging(bIsStockAveraging),
    m_pConversion(pConversion)
  {
  }

  ~PEPSAveragingEvent () {};
   
protected:
  
  double GetNewAverage(double dA, double dS) const; 
  
 
private:

  ///local pointer of the mandatory conversion to determine 
  ///the current conversion ratio
   const MandatoryConversion  *m_pConversion;

  ///whether stock averaging or conversion ratio is used
  bool m_bIsStockAveraging;

};

} // namespace pricing

} // namespace ito33 

#endif // #ifndef _ITO33_PRICING_PEPS_AVERAGINGEVENT_H_
