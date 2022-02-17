/////////////////////////////////////////////////////////////////////////////
// Name:        testcdp.cpp
// Purpose:     Acceptance test for cumulative default probabilities
// Created:     13/01/2006
// RCS-ID:      $Id: testcdp.cpp,v 1.4 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------
#include "ito33/beforestd.h"
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/link.h"
#include "ito33/sharedptr.h"
#include "ito33/date.h"
#include "ito33/vector.h"
#include "ito33/dateutils.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/numeraire.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardrateflat.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratepower.h"

#include "testcdp.h"

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::ihg;

namespace ito33
{

namespace ihg
{

namespace CDPTEST
{

// Called by the constructor
shared_ptr<finance::SessionData> CDPTest::MakeSessionData()
{
  Date valuationDate(2003, Date::Feb, 1);

  shared_ptr<Numeraire> pNumeraire(new Numeraire("EUR"));

  shared_ptr<Equity> pEquity(new Equity(50.0, pNumeraire));

  shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.0));
  
  pEquity->SetBorrowCurve(pyf);

  shared_ptr<YieldCurve> pyc( new YieldCurveFlat(0.0) );

  shared_ptr<RateData> pRateData(new RateData);

  pRateData->SetYieldCurve(pNumeraire, pyc);
  
  shared_ptr<SessionData> pSessionData( new 
    SessionData(pRateData, pEquity, valuationDate) );

  return pSessionData;

} // MakeSessionData


std::vector<Date> CDPTest::MakeDates()
{
  // Create the dates at which the probabilities are desired
  std::vector<Date> pDates( 3, m_pSessionData->GetValuationDate() );
  for (size_t nIdxDate = 0; nIdxDate < pDates.size(); nIdxDate++)
    pDates[nIdxDate].AddYears(nIdxDate+1);

  return pDates;
}


double CDPTest::CDPForTimeOnlyHR(shared_ptr<HazardRateTimeOnly> pHR,
                                   double dStartTime,
                                   double dEndTime)
{

  // Formula is: cdp = 1 - exp( - \int hr(t) dt), where the integral is 
  // from the start time to the endtime
  // Start by computing the integral.
  const std::vector<double> pdValues = pHR->GetValues();
  const std::vector<Date> pdDates = pHR->GetDates();
  size_t nNbTimes = pdDates.size();

  // Find the first hr time past the start time
  size_t nIdx = 0;  
  while (nIdx < nNbTimes - 1 && GetDoubleFrom(pdDates[nIdx]) < dStartTime)
    nIdx++;

  // compute cumulative sum for each interval
  double dSum = 0.0;
  double dIntervalStart = dStartTime;
  double dIntervalEnd = dStartTime;
  double dCurrentHR = pdValues[nIdx];
  while ( dIntervalEnd < dEndTime && nIdx < nNbTimes)
  {
    // Check if we went passed the end time, or if end time is passed
    // the last hazard rate time
    dIntervalEnd = GetDoubleFrom(pdDates[nIdx]);
    if (dIntervalEnd > dEndTime || nIdx == nNbTimes - 1)
      dIntervalEnd = dEndTime;

    dSum += (dIntervalEnd - dIntervalStart) * dCurrentHR;

    // Get readly for next iterval. Use hazard rate at right side.
    dIntervalStart = dIntervalEnd;
    nIdx++;
    if (nIdx < nNbTimes)
      dCurrentHR = pdValues[nIdx];
  }

  // Formula is: cdp = 1 - exp( - \int hr(t) dt), where the integral is 
  // from the start time to the endtime
  double dCDP = exp(-dSum);
  dCDP = 1.0 - dCDP;

  return dCDP;

} // CDPForTimeOnlyHR


void CDPTest::ZeroHR() 
{
  // Create the model with no hazard rate
  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  shared_ptr<ihg::VolatilityFlat> pVol( new ihg::VolatilityFlat(0.2) ); 
  pModel->SetVolatility( pVol );

  shared_ptr<ihg::HazardRateFlat> pHR( new ihg::HazardRateFlat(0.0) );
  pModel->SetHazardRate( pHR );

  // Do the computation
  std::vector<double> pdProbs = 
    pModel->ComputeCumulativeDefaultProbability(*m_pSessionData, m_pDates);

  // with no hazard rate, default probability should be zero
  for (size_t nIdx = 0; nIdx < pdProbs.size(); nIdx++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(pdProbs[nIdx], 0.0, 1.e-8);  

} // ZeroHR()


void CDPTest::ConstHR() 
{

  // Create the model with constant hazard rate
  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  shared_ptr<ihg::VolatilityFlat> pVol( new ihg::VolatilityFlat(0.2) ); 
  pModel->SetVolatility( pVol );

  double dHR = 0.05;
  shared_ptr<ihg::HazardRateFlat> pHR( new ihg::HazardRateFlat(dHR) );
  pModel->SetHazardRate( pHR );

  // Do the computation
  std::vector<double> pdProbs = 
    pModel->ComputeCumulativeDefaultProbability(*m_pSessionData, m_pDates);

  // compare to analytic formula  
  for (size_t nIdx = 0; nIdx < pdProbs.size(); nIdx++)
  {
    double dValuationTime = GetDoubleFrom(m_pSessionData->GetValuationDate());
    double dDateTime = GetDoubleFrom(m_pDates[nIdx]);
    double dTime = dDateTime - dValuationTime;

    double dAnalytic = 1.0 - exp(-dHR*dTime);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dAnalytic, pdProbs[nIdx], 1.e-4);  
  }

} // ConstHR()


void CDPTest::TimeOnlyHR1() 
{

  // Test tolerance
  double dTol = 2.e-2;

  // Create the model with time only hazard rate
  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  shared_ptr<ihg::VolatilityFlat> pVol( new ihg::VolatilityFlat(0.2) ); 
  pModel->SetVolatility( pVol );

  std::vector<Date> pDates(1);
  pDates[0] = m_pSessionData->GetValuationDate();
  std::vector<double> pdValues(1);
  pdValues[0] = 0.05;
  
  shared_ptr<ihg::HazardRateTimeOnly> pHR( 
    new ihg::HazardRateTimeOnly(pDates, pdValues) );
  pModel->SetHazardRate( pHR );

  // Do the computation
  std::vector<double> pdProbs = 
    pModel->ComputeCumulativeDefaultProbability(*m_pSessionData, m_pDates);

  // compare to analytic formula  
  for (size_t nIdx = 0; nIdx < pdProbs.size(); nIdx++)
  {
    double dValuationTime = GetDoubleFrom(m_pSessionData->GetValuationDate());
    double dTime = GetDoubleFrom( m_pDates[nIdx] );

    double dAnalytic = CDPForTimeOnlyHR(pHR, dValuationTime, dTime);
    double dError = fabs(pdProbs[nIdx] - dAnalytic) / dAnalytic;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, dError, dTol);  
  }

  // Repeat with HR date before valuation date
  pDates[0].AddYears(-1);
  pHR = make_ptr( new ihg::HazardRateTimeOnly(pDates, pdValues) );
  pModel->SetHazardRate( pHR );

  pdProbs = 
    pModel->ComputeCumulativeDefaultProbability(*m_pSessionData, m_pDates);

  // compare to analytic formula  
  for (size_t nIdx = 0; nIdx < pdProbs.size(); nIdx++)
  {
    double dValuationTime = GetDoubleFrom(m_pSessionData->GetValuationDate());
    double dTime = GetDoubleFrom( m_pDates[nIdx] );

    double dAnalytic = CDPForTimeOnlyHR(pHR, dValuationTime, dTime);
    double dError = fabs(pdProbs[nIdx] - dAnalytic) / dAnalytic;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, dError, dTol);  
  }

  // Repeat with HR date after last date
  pDates[0] = m_pDates[m_pDates.size() - 1];
  pDates[0].AddYears(1);
  pHR = make_ptr( new ihg::HazardRateTimeOnly(pDates, pdValues) );
  pModel->SetHazardRate( pHR );

  pdProbs = pModel->ComputeCumulativeDefaultProbability(*m_pSessionData, m_pDates);

  // compare to analytic formula  
  for (size_t nIdx = 0; nIdx < pdProbs.size(); nIdx++)
  {
    double dValuationTime = GetDoubleFrom(m_pSessionData->GetValuationDate());
    double dTime = GetDoubleFrom( m_pDates[nIdx] );

    double dAnalytic = CDPForTimeOnlyHR(pHR, dValuationTime, dTime);
    double dError = fabs(pdProbs[nIdx] - dAnalytic) / dAnalytic;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, dError, 1.e-3);  
  }

} // TimeOnlyHR1()


void CDPTest::TimeOnlyHR2() 
{

  // Test tolerance
  double dTol = 2.e-2;

  // Create the model with time only hazard rate
  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  shared_ptr<ihg::VolatilityFlat> pVol( new ihg::VolatilityFlat(0.2) ); 
  pModel->SetVolatility( pVol );

  std::vector<Date> pDates(4, m_pDates[0]);
  pDates[0].AddMonths(-15);
  pDates[1].AddMonths(-8);
  pDates[2].AddMonths(3);
  pDates[3].AddMonths(10);
  std::vector<double> pdValues(4);
  pdValues[0] = 0.03;
  pdValues[1] = 0.05;
  pdValues[2] = 0.1;
  pdValues[3] = 0.07;
  
  shared_ptr<ihg::HazardRateTimeOnly> pHR( 
    new ihg::HazardRateTimeOnly(pDates, pdValues) );
  pModel->SetHazardRate( pHR );

  // Do the computation
  std::vector<double> pdProbs = 
    pModel->ComputeCumulativeDefaultProbability(*m_pSessionData, m_pDates);

  // compare to analytic formula  
  for (size_t nIdx = 0; nIdx < pdProbs.size(); nIdx++)
  {
    double dValuationTime = GetDoubleFrom(m_pSessionData->GetValuationDate());
    double dTime = GetDoubleFrom( m_pDates[nIdx] );

    double dAnalytic = CDPForTimeOnlyHR(pHR, dValuationTime, dTime);
    double dError = fabs(pdProbs[nIdx] - dAnalytic) / dAnalytic;
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, dError, dTol);  
  }

  // Repeat with HR dates that go past last test date
  pDates[0].AddMonths(5);
  pDates[1].AddMonths(1);
  pDates[2].AddMonths(10);
  pDates[3].AddMonths(25);

  pHR = make_ptr( new ihg::HazardRateTimeOnly(pDates, pdValues) );
  pModel->SetHazardRate( pHR );

  pdProbs = pModel->ComputeCumulativeDefaultProbability(*m_pSessionData, m_pDates);

  // compare to analytic formula  
  for (size_t nIdx = 0; nIdx < pdProbs.size(); nIdx++)
  {
    double dValuationTime = GetDoubleFrom(m_pSessionData->GetValuationDate());
    double dTime = GetDoubleFrom( m_pDates[nIdx] );

    double dAnalytic = CDPForTimeOnlyHR(pHR, dValuationTime, dTime);
    double dError = fabs(pdProbs[nIdx] - dAnalytic) / dAnalytic;
    
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, dError, dTol);  
  }

} // TimeOnlyHRHR2()


void CDPTest::IncreasingHR() 
{

  // Create the model with time only hazard rate
  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  shared_ptr<ihg::VolatilityFlat> pVol( new ihg::VolatilityFlat(0.2) ); 
  pModel->SetVolatility( pVol );

  // Set initial probability values to zero
  std::vector<double> pdProbsOld(m_pDates.size(), 0.0);

  for (double dAlpha = 0.05; dAlpha < 0.25; dAlpha += 0.02)
  {

    // Set the hazard rate. Increases each iteration.
    double dBeta = 0.6;
    double dS0 = m_pSessionData->GetSpotSharePrice();

    shared_ptr<ihg::HazardRatePower> pHR( 
      new ihg::HazardRatePower(dAlpha, dBeta, dS0) );
    pModel->SetHazardRate( pHR );

    // Do the computation
    std::vector<double> pdProbs = 
      pModel->ComputeCumulativeDefaultProbability(*m_pSessionData, m_pDates);

    // Values should increase each iteration
    for (size_t nIdx = 0; nIdx < pdProbs.size(); nIdx++)
      CPPUNIT_ASSERT(pdProbs[nIdx] > pdProbsOld[nIdx]);  

    // Get ready for next iteration
    pdProbsOld = pdProbs;

  } // loop over increasing alpha values

} // IncreasingHR


} // namespae CDPTest

} // namespace ihg

} // namespace ito33

