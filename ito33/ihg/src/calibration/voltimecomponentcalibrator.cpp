/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/voltimecomponentcalibrator.cpp
// Purpose:     calibrate time component of VolatilityWithtimeComponent
// Created:     2005/03/04
// RCS-ID:      $Id: voltimecomponentcalibrator.cpp,v 1.6 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/numeric/nonlinearsolver.h"

#include "ito33/finance/modeloutput.h"

#include "ito33/ihg/volatility.h"
#include "ito33/ihg/volatilitywithtimecomponent.h"
#include "ito33/ihg/volatilitytimeonly.h"

#include "ihg/voltimecomponentcalibrator.h"

namespace ito33
{

namespace ihg
{


double VolTimeComponentCalibrator::operator()(double dVol)
{
  m_pdValues[m_nIdx] = dVol;

  m_pVolatility->ResetTimeComponent
                 (&m_pMaturityDates[0], &m_pdValues[0], m_nIdx + 1);

  m_theoreticalModel.SetVolatility(m_pVolatility);

  shared_ptr<finance::ModelOutput> 
    output = m_theoreticalModel.Compute(*m_pDerivative);

  return (output->GetPrice() - m_dMarketPrice) * m_dScale;    
}

shared_ptr<VolatilityWithTimeComponent> 
VolTimeComponentCalibrator::Calibrate
    (
      const finance::TermStructureDerivative& tsDerivative,
      shared_ptr<VolatilityWithTimeComponent> pVolatility
    )
{ 
  // Get the maturity dates from the term structure.  These dates will
  // be used in the output VolatilityWithTimeComponent
  const finance::TermStructureDerivative::Elements& 
    derivs = tsDerivative.GetAll(); 
  
  m_pMaturityDates.clear();
 
  finance::TermStructureDerivative::Elements::const_iterator 
    ppDerivative = derivs.begin();

  for ( ; ppDerivative != derivs.end(); ++ppDerivative)
    m_pMaturityDates.push_back( (*ppDerivative)->GetMaturityDate() );

  m_pdValues.resize( derivs.size() );

  // Create the vol if needed
  if (pVolatility)
    m_pVolatility = pVolatility;
  else
    m_pVolatility = make_ptr( new VolatilityTimeOnly
                                  ( &m_pMaturityDates[0], &m_pdValues[0],
                                    m_pdValues.size() ) );

  // Loop through the term structure, building the volatility along the way
  for (ppDerivative = derivs.begin(), m_nIdx = 0; 
       ppDerivative != derivs.end(); 
       ++ppDerivative, ++m_nIdx)
  {
    m_pDerivative = *ppDerivative;

    // Save values needed by operator (), which is called by the root finder
    m_dMarketPrice = m_pDerivative->GetMarketPrice();

    // Also the scale
    m_dScale = fabs( m_dMarketPrice );
    if (m_dScale < 1.e-6)
      m_dScale = 1.0;
    
    m_dScale = 1. / m_dScale;

    // Copy the logic from the ImpliedVol class
    numeric::RegulaFalsi solver(1.e-4, 300);

    double dVol = 0;
    try
    {
      // don't begin with 5: too high to be realistic and often problematic
      dVol = solver(*this, 0., 1.);
    }
    catch(const ito33::numeric::Exception& /* e */)
    {
      // Try again with higher starting values
      dVol = solver(*this, 1., 4.99);
    }

    m_pdValues[m_nIdx] = dVol;
  }

  m_pVolatility->ResetTimeComponent
                 ( &m_pMaturityDates[0], &m_pdValues[0], m_pdValues.size() );

  return m_pVolatility;
}


} // namespace ihg

} // namespace ito33

