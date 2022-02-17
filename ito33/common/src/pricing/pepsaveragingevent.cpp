/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/pepsaveragingevent.cpp
// Purpose:     PEPS Averaging Event
// Author:      Ito33 team Canada
// Created:     June 27, 2005
// RCS-ID:      $Id: pepsaveragingevent.cpp,v 1.6 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"

#include "ito33/pricing/mandatoryconversion.h"

#include "ito33/pricing/pepsaveragingevent.h"

namespace ito33
{

namespace pricing
{

 
double PEPSAveragingEvent::GetNewAverage(double dA, double dS) const
{
  double dANew = -1.0;

  if ( m_bIsStockAveraging )  
    dANew = dA + ( dS - dA ) / double( m_nObservation );
  else
  {
    double dConversion = m_pConversion->GetConversionRatio(dS);
    
    dANew = dA + ( dConversion - dA ) / double(m_nObservation);
  }
  
  return dANew;

} // PEPSAveragingEvent::GetNewAverage(double dA, double dS) const

 
} // namespace pricing

} // namespace ito33
