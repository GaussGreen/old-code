/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/hrtimecomponentcalibrator.cpp
// Purpose:     calibrate time component of HazardRateWithtimeComponent
// Author:      Wang
// Created:     2004/06/11
// RCS-ID:      $Id: hrtimecomponentcalibrator.cpp,v 1.13 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ihg/src/calibration/hrtimecomponentcalibrator.cpp
   @brief calibration of time component of HazardRateWithtimeComponent 

   @todo imrpove of exception handling?
 */

#include "ito33/numeric/nonlinearsolver.h"

#include "ito33/finance/modeloutput.h"

#include "ito33/ihg/volatility.h"
#include "ito33/ihg/hazardratewithtimecomponent.h"
#include "ito33/ihg/hazardratetimeonly.h"

#include "ihg/hrtimecomponentcalibrator.h"

namespace ito33
{

namespace ihg
{


double HazardRateTimeComponentCalibrator::operator()(double dHazardRate)
{
  m_pdValues[m_nIdx] = dHazardRate;

  m_pHazardRate->ResetTimeComponent
                 (&m_pMaturityDates[0], &m_pdValues[0], m_nIdx + 1);

  m_theoreticalModel.SetHazardRate(m_pHazardRate);

  shared_ptr<finance::ModelOutput> 
    output = m_theoreticalModel.Compute( *m_pDerivative);

  double dScale = fabs( m_dMarketPrice );
  if (dScale < 1.e-6)
    dScale = 1.0;

  return (output->GetPrice() - m_dMarketPrice)/dScale;    
}

shared_ptr<HazardRateWithTimeComponent> 
HazardRateTimeComponentCalibrator::Calibrate
    (
      const finance::TermStructureDerivative& tsDerivative,
      const shared_ptr<Volatility>& pVolatility, 
      shared_ptr<HazardRateWithTimeComponent> pHazardRate
    )
{ 
  m_theoreticalModel.SetVolatility(pVolatility);

  const finance::TermStructureDerivative::Elements& 
    derivs = tsDerivative.GetAll(); 
  
  m_pMaturityDates.clear();
 
  finance::TermStructureDerivative::Elements::const_iterator 
    ppDerivative = derivs.begin();

  for ( ; ppDerivative != derivs.end(); ++ppDerivative)
    m_pMaturityDates.push_back( (*ppDerivative)->GetMaturityDate() );

  m_pdValues.resize( derivs.size() );

  if (pHazardRate)
    m_pHazardRate = pHazardRate;
  else
    m_pHazardRate = make_ptr( new HazardRateTimeOnly
                                  ( &m_pMaturityDates[0], &m_pdValues[0],
                                    m_pdValues.size() ) );

  for (ppDerivative = derivs.begin(), m_nIdx = 0; 
       ppDerivative != derivs.end(); 
       ++ppDerivative, ++m_nIdx)
  {
    m_pDerivative = *ppDerivative;

    // Save values needed by operator (), which is called by the root finder
    m_dMarketPrice = m_pDerivative->GetMarketPrice();

    numeric::RegulaFalsi solver(1.e-4, 200);
      
    m_pdValues[m_nIdx] = solver(*this, 0, 10);

  }

  m_pHazardRate->ResetTimeComponent
                 ( &m_pMaturityDates[0], &m_pdValues[0], m_pdValues.size() );

  return m_pHazardRate;
}


} // namespace ihg

} // namespace ito33

