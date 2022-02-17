/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/voltimeonlyhrtimecomponentcalibrator.cpp
// Purpose:     calibrator on a time only vol and a hr with time component
// Created:     2005/07/15
// RCS-ID:      $Id: voltimeonlyhrtimecomponentcalibrator.cpp,v 1.4 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/list.h"
#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/numeric/exception.h"
#include "ito33/numeric/newton2d.h"
#include "ito33/numeric/numericerror.h"

#include "ito33/finance/error.h"
#include "ito33/finance/derivative.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/ihg/hazardratewithtimecomponent.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/volatilitytimeonly.h"

#include "ihg/voltimeonlyhrtimecomponentcalibrator.h"

extern const ito33::Error ITO33_MAX_ITER;
extern const ito33::finance::Error ITO33_CALIBRATION_FAIL;

namespace ito33
{

namespace ihg
{

using ito33::numeric::Exception;


void VolTimeOnlyHRTimeComponentCalibrator::operator () 
   (double dVol, double dHR, double &dPrice1, double &dPrice2, double &dF)
{

  // Set the current guess
  m_pdVolValues[m_nIdx] = dVol;
  m_pdHRValues[m_nIdx] = dHR;

  m_pVolatility->ResetTimeComponent
                 (&m_pMaturityDates[0], &m_pdVolValues[0], m_nIdx + 1);

  m_pHazardRate->ResetTimeComponent
                 (&m_pMaturityDates[0], &m_pdHRValues[0], m_nIdx + 1);

  m_theoreticalModel.SetHazardRate(m_pHazardRate);
  m_theoreticalModel.SetVolatility(m_pVolatility);

  dPrice1 = m_theoreticalModel.Compute(*m_pDeriv1)->GetPrice();
  dPrice2 = m_theoreticalModel.Compute(*m_pDeriv2)->GetPrice();

  dPrice1 -= m_dMarketPrice1;
  dPrice2 -= m_dMarketPrice2;

  dF = .5*( dPrice1*dPrice1*m_dScale1*m_dScale1 
          + dPrice2*dPrice2*m_dScale2*m_dScale2 );

} 


void VolTimeOnlyHRTimeComponentCalibrator::Calibrate
     (
     const std::list< shared_ptr<finance::Derivative> >& derivatives, 
     shared_ptr<HazardRateWithTimeComponent> pHazardRate
     )
{ 

  // Initialize the dates and values in the time components. It is
  // assumed that the derivative list has already been sorted, and
  // has an even number of entries
  size_t nNbDerivatives = derivatives.size();

  m_pMaturityDates = std::vector<Date>(nNbDerivatives / 2);
  m_pdVolValues = std::vector<double>(nNbDerivatives / 2, 0.2);
  m_pdHRValues = std::vector<double>(nNbDerivatives / 2, 0.0);

  std::list< shared_ptr<finance::Derivative> >::const_iterator iter;
  iter = derivatives.begin();  
  size_t nCounter = 0;
  while (iter != derivatives.end() )
  {
    ++iter;
    m_pMaturityDates[nCounter++] = (*iter)->GetMaturityDate();
    ++iter;
  }

  // Init output objects to default values in case calibration fails
  m_pVolatility = make_ptr( new VolatilityTimeOnly
                                ( m_pMaturityDates, m_pdVolValues ) );

  if (pHazardRate)
    m_pHazardRate = pHazardRate;
  else
    m_pHazardRate = make_ptr( new HazardRateTimeOnly
                                  (m_pMaturityDates, m_pdHRValues) );

  // Newton solver used during each period
  double dVolLowerBound = 0.0;
  double dVolUpperBound = 4.99;
  double dHRLowerBound = 0.0;
  double dHRUpperBound = 1.99;
  // with a hr spot component, the hr parameters have a different meaning
  if (pHazardRate)
    dHRUpperBound = 10.0;

  double dTol = 1.e-8;
  size_t nIterMax = 80;

  numeric::Newton2D solver(dVolLowerBound, dVolUpperBound,
                           dHRLowerBound, dHRUpperBound,
                           dTol, nIterMax);

  // Loop over the derivatives in pairs.  The second derivative in
  // each pair defines the end of the period.  For each period,
  // find the vol and hr time components to fit the pair of derivatives.
  m_nIdx = 0;
  double dLastVol = 0.2;
  double dLastHR = 0.05;

  iter = derivatives.begin();
  while (iter != derivatives.end() )
  {
    // Move through the list in pairs
    m_pDeriv1 = *iter;
    ++iter;
    m_pDeriv2 = *iter;

    m_dMarketPrice1 = m_pDeriv1->GetMarketPrice();
    m_dMarketPrice2 = m_pDeriv2->GetMarketPrice();

    m_dScale1 = fabs(m_dMarketPrice1);
    if (m_dScale1 > 1.e-6)
      m_dScale1 = 1.0/m_dScale1;
    else
      m_dScale1 = 1.0;

    m_dScale2 = fabs(m_dMarketPrice2);
    if (m_dScale2 > 1.e-6)
      m_dScale2 = 1.0/m_dScale2;
    else
      m_dScale2 = 1.0;

    // Initial guess is the last calibrated value
    double dVol = dLastVol;
    double dHR = dLastHR;

    // Actually calibrate.  This calls the () operator
    numeric::NumericError err = solver(*this, dVol, dHR);

    if ( err == numeric::ITO33_NOT_CONVERGED )
      throw EXCEPTION(ITO33_CALIBRATION_FAIL);

    if ( err == numeric::ITO33_TOO_MANY_ITERATION )
      throw EXCEPTION(ITO33_MAX_ITER);

    m_pdVolValues[m_nIdx] = dVol;
    m_pdHRValues[m_nIdx] = dHR;

    // Move to next index and next derivative
    m_nIdx++;
    ++iter;
  }

  // Save the final values
  m_pHazardRate->ResetTimeComponent
                 (&m_pMaturityDates[0], &m_pdHRValues[0], m_pdHRValues.size());

  m_pVolatility->ResetTimeComponent
               (&m_pMaturityDates[0], &m_pdVolValues[0], m_pdVolValues.size());

}


} // namespace ihg

} // namespace ito33
