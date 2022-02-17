/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/translator.cpp
// Purpose:     implementation of IHG helper class Translator for calibration
// RCS-ID:      $Id: translator.cpp,v 1.9 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/*
 @TODO  Add functions to set initial guess.
*/

#include <cmath>

#include "ito33/date.h"
#include "ito33/vector.h"
#include "ito33/list.h"
#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/derivatives.h"

#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/volatilitypower.h"
#include "ito33/ihg/volatilitytanh.h"

#include "ito33/ihg/hazardrateflat.h"
#include "ito33/ihg/hazardratepower.h"
#include "ito33/ihg/hazardratecombo.h"
#include "ito33/ihg/hrspotcomponentpower.h"

#include "ihg/translator.h"

#include "ito33/ihg/theoreticalmodel.h"

namespace ito33
{

  ITO33_IMPLEMENT_AUTOPTR(ihg::Translator);

} // namespace ito33

extern const ito33::Error ITO33_BAD_PARAM;

namespace ito33
{

namespace ihg
{


Translator::Translator(VolType volType, HRType hrType, double dSpot)
: m_VolType(volType), m_HRType(hrType), m_dSpot(dSpot)
{
  
  // Get the number of vol and hr params so we can initialize arrays to
  // the correct size
  // Get the total number of vol params first.  
  size_t nNbVol = 0;
  switch (m_VolType)
  {
  case VolType_flat:
    nNbVol = 1;
    break;

  case VolType_power:
    nNbVol = 2;
    break;

  case VolType_tanh:
    nNbVol = 3;
    break;

  default:
    throw EXCEPTION_MSG
          (
            ITO33_BAD_PARAM,
            TRANS("Unknown volatility type in IHG translator.")
          );
  }

  // Now get number of hr params
  size_t nNbHR = 0;

  switch (m_HRType)
  {
  case HRType_flat:
    nNbHR = 1;
    break;

  case HRType_power:
    nNbHR = 2;
    break;

  default:
    throw EXCEPTION_MSG
          (
            ITO33_BAD_PARAM,
            TRANS("Unknown hazard rate type in IHG translator.")
          );
  }

  // Allocate storage and setup default param values and bounds
  m_nNbFlags = nNbVol + nNbHR;
  m_flags.resize(m_nNbFlags);
  m_pdParams.resize(m_nNbFlags);
  m_pdLowerBounds.resize(m_nNbFlags);
  m_pdUpperBounds.resize(m_nNbFlags);

  // initialize class arrays
  InitVolParams(0);
  InitHRParams(nNbVol);
 
  // init the flags
  for (size_t n = 0; n < m_nNbFlags; n++)
    m_flags[n] = true;

  // Clear the dates, so it can act like a flag
  m_TimeComponentDates.clear();
}


Translator::Translator(VolType volType, 
                       SpotType spotType, 
                       finance::Derivatives derivatives,
                       double dSpot)
  : m_VolType(volType), m_spotType(spotType), m_dSpot(dSpot)
{
  // Get the number of vol and hr params so we can initialize arrays to
  // the correct size
  // Get the total number of vol params first.  
  size_t nNbVol = 0;
  switch (m_VolType)
  {
  case VolType_flat:
    nNbVol = 1;
    break;

  case VolType_power:
    nNbVol = 2;
    break;

  case VolType_tanh:
    nNbVol = 3;

  default:
    throw EXCEPTION_MSG
          (
            ITO33_BAD_PARAM,
            TRANS("Unknown volatility type in IHG translator.")
          );
  }

  // Now get number of spot hr params
  size_t nNbSpotHR = 0;

  switch (m_spotType)
  {
  case SpotType_power:
    nNbSpotHR = 1;
    break;

  default:
    throw EXCEPTION_MSG
          (
            ITO33_BAD_PARAM,
            TRANS("Unknown spot hazard rate type in IHG translator.")
          );
  }

  // Now get the number of time components
  m_TimeComponentDates = GetMaturityDates(derivatives);
  size_t nNbMaturities = m_TimeComponentDates.size();

  // Allocate storage and setup default param values and bounds
  m_nNbFlags = nNbVol + nNbSpotHR + nNbMaturities;
  m_flags.resize(m_nNbFlags);
  m_pdParams.resize(m_nNbFlags);
  m_pdLowerBounds.resize(m_nNbFlags);
  m_pdUpperBounds.resize(m_nNbFlags);

  // initialize class arrays
  InitVolParams(0);
  InitSpotHRParams(nNbVol);
  InitTimeComponentParams(nNbVol + nNbSpotHR, nNbMaturities);

  // init the flags
  for (size_t n = 0; n < m_nNbFlags; n++)
    m_flags[n] = true;

  // Since the power spot component does not have an "alpha" parameter, 
  // scaling is determined by the hazard rate time components. Otherwise
  // the first time component should be frozen (flag set to false).
  
}


shared_ptr<ihg::TheoreticalModel> 
Translator::operator()(const double* pdXOriginal)
{
  // Apply the bounds (note that the NAG optimizer can sometimes make a guess
  // that violates the bounds, especially the lower bounds)
  const std::vector<double> pdXNew = ApplyBounds(pdXOriginal);
  const double* pdX = &pdXNew[0];

  // counters to know where we are in the param arrays
  size_t nFlagCounter = 0;
  size_t nXCounter = 0;

  shared_ptr<Volatility> pVolatility;
  pVolatility = RestoreVolatility(pdX, nFlagCounter, nXCounter);

  shared_ptr<HazardRate> pHazardRate;
  if (m_TimeComponentDates.size() == 0)
  {
    pHazardRate = RestoreHazardRate(pdX, nFlagCounter, nXCounter);
  }
  else
  {
    pHazardRate = RestoreComboHR(pdX, nFlagCounter, nXCounter);
  }

  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel());
  pModel->SetVolatility(pVolatility);
  pModel->SetHazardRate(pHazardRate);

  return pModel;
}


shared_ptr<ihg::TheoreticalModel> 
Translator::operator()(const std::vector<double>& pdX)
{
  return operator()(&pdX[0]);
}


void Translator::GetParameters(double* pdX) const
{
  size_t nXCounter = 0;
  size_t nFlagCounter = 0;
 
  // loop over the params
  for (size_t nIdx = 0; nIdx < m_pdParams.size(); nIdx++)
  {
    if ( m_flags[nFlagCounter] )
      pdX[nXCounter++] = m_pdParams[nFlagCounter];

    nFlagCounter++;
  }

}


std::vector<double> Translator::GetLowerBounds() const
{
  std::vector<double> pdBounds( GetNbParameters() );

  size_t nXCounter = 0;
 
  // loop over the lower bounds
  for (size_t nIdx = 0; nIdx < m_nNbFlags; nIdx++)
  {
    if ( m_flags[nIdx] )
      pdBounds[nXCounter++] = m_pdLowerBounds[nIdx];
  }

  return pdBounds;
}


std::vector<double> Translator::GetUpperBounds() const
{
  std::vector<double> pdBounds( GetNbParameters() );

  size_t nXCounter = 0;
 
  // loop over the upper bounds
  for (size_t nIdx = 0; nIdx < m_nNbFlags; nIdx++)
  {
    if ( m_flags[nIdx] )
      pdBounds[nXCounter++] = m_pdUpperBounds[nIdx];
  }

  return pdBounds;
}


void Translator::InitVolParams(size_t nCounter)
{
  
  // Assume vol type has been initialized correctly
  switch (m_VolType)
  {
  case VolType_flat:
    m_pdParams[nCounter] = 0.2;
    m_pdLowerBounds[nCounter] = 0.0;
    m_pdUpperBounds[nCounter] = 4.99;
    nCounter++;
    break;

  case VolType_power:    
    // alpha
    m_pdParams[nCounter] = 0.2;  
    m_pdLowerBounds[nCounter] = 0.0;
    m_pdUpperBounds[nCounter] = 2.0;
    nCounter++;

    // beta
    m_pdParams[nCounter] = 0.0;
    m_pdLowerBounds[nCounter] = -1.99;
    m_pdUpperBounds[nCounter] = 1.99;
    nCounter++;

    break;

  case VolType_tanh:
    // value at -ve infinity
    m_pdParams[nCounter] = 0.3;
    m_pdLowerBounds[nCounter] = 1.e-12;
    m_pdUpperBounds[nCounter] = 2.0;
    nCounter++;

    // value at +ve infinity
    m_pdParams[nCounter] = 0.1;
    m_pdLowerBounds[nCounter] = 1.e-12;
    m_pdUpperBounds[nCounter] = 2.0;
    nCounter++;

    // scale factor
    m_pdParams[nCounter] = 1.0;
    m_pdLowerBounds[nCounter] = 1.e-12;
    m_pdUpperBounds[nCounter] = 100.0;
    nCounter++;

  default:
    break;
  }

}


shared_ptr<Volatility> 
Translator::RestoreVolatility(const double* pdX, size_t& nFlagCounter,
                              size_t& nXCounter) const
{
  shared_ptr<Volatility> pVolatility;

  switch (m_VolType)
  {
  case VolType_flat:
    {
    if (m_flags[nFlagCounter])
      pVolatility = make_ptr( new VolatilityFlat(pdX[nXCounter++]) );
    else
      pVolatility = make_ptr( new VolatilityFlat(m_pdParams[nFlagCounter]) );

    nFlagCounter++;

    break;
    }
  case VolType_power:
    {
    double dAlpha, dBeta;
    
    if (m_flags[nFlagCounter])
      dAlpha = pdX[nXCounter++];
    else 
      dAlpha = m_pdParams[nFlagCounter];
    nFlagCounter++;

    if (m_flags[nFlagCounter])
      dBeta = pdX[nXCounter++];
    else 
      dBeta = m_pdParams[nFlagCounter];
    nFlagCounter++;

    pVolatility = make_ptr( new VolatilityPower(dAlpha, dBeta, m_dSpot) );

    break;
    }
  case VolType_tanh:
    {
    double dLeft, dRight, dScale;

    if (m_flags[nFlagCounter])
      dLeft = pdX[nXCounter++];
    else 
      dLeft = m_pdParams[nFlagCounter];
    nFlagCounter++;

    if (m_flags[nFlagCounter])
      dRight = pdX[nXCounter++];
    else 
      dRight = m_pdParams[nFlagCounter];
    nFlagCounter++;

    if (m_flags[nFlagCounter])
      dScale = pdX[nXCounter++];
    else 
      dScale = m_pdParams[nFlagCounter];
    nFlagCounter++;

    pVolatility = make_ptr( new VolatilityTanh
                                (dLeft, dRight, dScale, m_dSpot) );

    break;
    }
  default:
    break;
  }

  return pVolatility;
}


void Translator::InitHRParams(size_t nCounter)
{
  // Assume m_HRType is valid
  switch (m_HRType)
  {
  case HRType_flat:
    m_pdParams[nCounter] = 0.2;
    m_pdLowerBounds[nCounter] = 0.0;
    m_pdUpperBounds[nCounter] = 4.99;
    nCounter++;
    break;

  case HRType_power:
    // alpha
    m_pdParams[nCounter] = 0.01;  
    m_pdLowerBounds[nCounter] = 0.0;
    m_pdUpperBounds[nCounter] = 2.0;
    nCounter++;

    // beta
    m_pdParams[nCounter] = 0.01;
    m_pdLowerBounds[nCounter] = 0.0;
    m_pdUpperBounds[nCounter] = 1.99;
    nCounter++;

    break;

  default:
    break;
  }
}


shared_ptr<HazardRate> 
Translator::RestoreHazardRate(const double* pdX, size_t& nFlagCounter, 
                              size_t& nXCounter) const
{
  shared_ptr<HazardRate> pHazardRate;

  switch (m_HRType)
  {
  case HRType_flat:
    if (m_flags[nFlagCounter])
      pHazardRate = make_ptr( new HazardRateFlat(pdX[nXCounter++]) );
    else
      pHazardRate = make_ptr( new HazardRateFlat(m_pdParams[nFlagCounter]) );

    nFlagCounter++;
    break;

  case HRType_power:
    double dAlpha, dBeta;
    
    if (m_flags[nFlagCounter])
      dAlpha = pdX[nXCounter++];
    else 
      dAlpha = m_pdParams[nFlagCounter];
    nFlagCounter++;

    if (m_flags[nFlagCounter])
      dBeta = pdX[nXCounter++];
    else 
      dBeta = m_pdParams[nFlagCounter];
    nFlagCounter++;

    // Bug in nag can set beta to tiny, negative value
    if ( fabs(dBeta) < 1.e-16 )
      dBeta = 0.0;

    pHazardRate = make_ptr( new HazardRatePower(dAlpha, dBeta, m_dSpot) );
    break;

  default:
    break;
  }

  return pHazardRate;
}


void Translator::InitSpotHRParams(size_t nCounter)
{
  // Assume m_spotType is valid
  switch (m_spotType)
  {
  case SpotType_power:
    // beta
    m_pdParams[nCounter] = 0.01;
    m_pdLowerBounds[nCounter] = 0.0;
    m_pdUpperBounds[nCounter] = 1.99;
    nCounter++;

    break;

  default:
    break;
  }
}


shared_ptr<HazardRate> 
Translator::RestoreComboHR(const double* pdX, size_t& nFlagCounter, 
                           size_t& nXCounter) const
{
  shared_ptr<HazardRate> pHazardRate;
  shared_ptr<SpotComponent> pSpotComponent;

  switch (m_spotType)
  {
  case SpotType_power:
    double dBeta;

    if (m_flags[nFlagCounter])
      dBeta = pdX[nXCounter++];
    else 
      dBeta = m_pdParams[nFlagCounter];
    nFlagCounter++;

    pSpotComponent = make_ptr( new HRSpotComponentPower(dBeta, m_dSpot) );
    break;

  default:
    break;
  }

  size_t nNbTimes = m_TimeComponentDates.size();
  std::vector<double> pdValues(nNbTimes);

  ASSERT_MSG(nNbTimes == m_flags.size() - nFlagCounter,
             "Mismatch in number of parameters");

  for (size_t nIdx = 0; 
       nFlagCounter < m_flags.size(); 
       nFlagCounter++, nIdx++)
  {
    if (m_flags[nFlagCounter])
      pdValues[nIdx] = pdX[nXCounter++];
    else
      pdValues[nIdx] = m_pdParams[nFlagCounter];
  }
  
  pHazardRate = make_ptr( new HazardRateCombo
                              (pSpotComponent, 
                               m_TimeComponentDates, pdValues) );

  return pHazardRate;
}


void Translator::InitTimeComponentParams(size_t nCounter, 
                                         size_t nNbParams)
{

  // Only multiplication is allowed
  double dValue = 1.0;
  double dLowerBound = 0.0;
  double dUpperBound = 100.0;

  // initialize the arrays to the default value
  for (size_t nIdx = nCounter; nIdx < nCounter + nNbParams; nIdx++)
  {
    m_pdParams[nIdx] = dValue;
    m_pdLowerBounds[nIdx] = dLowerBound;
    m_pdUpperBounds[nIdx] = dUpperBound;  
  }

}


std::vector<Date> 
Translator::GetMaturityDates(finance::Derivatives derivatives)
{

  // Extract all maturity dates into a giant list
  std::list<Date> dateListTmp;
  finance::Derivatives::Elements::const_iterator iter;
  finance::Derivatives::Elements elements = derivatives.GetAll();

  for (iter = elements.begin();
       iter != elements.end();
       ++iter)
  {
    dateListTmp.push_back( iter->first->GetMaturityDate() );
  }

  // Sort and remove duplicates
  dateListTmp.sort();

  std::list<Date> dateList;
  dateList.push_back(dateListTmp.front());
  std::list<Date>::const_iterator iterDates;
  for (iterDates = dateListTmp.begin(); 
       iterDates != dateListTmp.end(); 
       ++iterDates)
  {
    if ( *iterDates > dateList.back() )
      dateList.push_back( *iterDates );
  }

  // Convert to vector and return
  std::vector<Date> pdDates(dateList.size());
  size_t nCounter = 0;
  std::list<Date>::const_iterator iterList;
  for (iterList = dateList.begin();
       iterList != dateList.end();
       ++iterList)
  {
    pdDates[nCounter++] = *iterList;
  }

  return pdDates;
}

} // namespace ihg

} // namespace ito33
