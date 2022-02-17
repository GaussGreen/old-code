/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/tests/onetouch.cpp
// Purpose:     Implementation of tests on homogeneity applied to one touch
// Created:     2006/02/20
// RCS-ID:      $Id: onetouch.cpp,v 1.7 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include <cmath>

#include "ito33/date.h"

#include "ito33/finance/dividends.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/exoticoption/onetouch.h"
#include "ito33/finance/exoticoption/onetouches.h"

#include "ito33/numeric/interpolation.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/multioutput.h"

#include "hg/numoutput.h"
#include "hg/priceonetouches.h"
#include "hg/tests/onetouch.h"

namespace ito33
{

namespace hg
{

namespace OneTouchTest
{

// The homogeneity is well respected so we are using a high precision.
static const double dPrecision = 5 * 1.e-4;

// Global object to work around CppUnit limitation
shared_ptr<TheoreticalModel> m_pModel;
shared_ptr<finance::OneTouch> m_pOneTouch;

void Setup(shared_ptr<TheoreticalModel> pModel,
           shared_ptr<finance::OneTouch> pOneTouch)
{
  m_pModel = make_ptr( pModel->Clone() );
  m_pOneTouch = pOneTouch;
}

bool OneTouchTest::HomogeneityCanApply()
{
  shared_ptr<finance::SessionData>
    pSessionData = m_pOneTouch->GetSessionData();

  shared_ptr<finance::Dividends>
    pDividends = pSessionData->GetEquity()->GetDividends();

  if (    pDividends 
       && pDividends->HasCashBetween
              (pSessionData->GetValuationDate(), m_pOneTouch->GetMaturityDate())
     )
  {
    return false;
  }

  if (    m_pOneTouch->GetBarrierType() == finance::Barrier_UpAndOut 
       && pSessionData->GetSpotSharePrice() >= m_pOneTouch->GetBarrier())
  {
    return false;
  }
  else if (   m_pOneTouch->GetBarrierType() == finance::Barrier_DownAndOut
           && pSessionData->GetSpotSharePrice() <= m_pOneTouch->GetBarrier())
  {
    return false;
  }

  return true;
}

void OneTouchTest::Homogeneity()
{
  if ( !HomogeneityCanApply() )
    return;

  shared_ptr<finance::SessionData>
    pSessionData = m_pOneTouch->GetSessionData();

  // Determine a new barrier
  double dShift = 0.1;
  double dBarrier2;
  if ( m_pOneTouch->GetBarrierType() == finance::Barrier_DownAndOut )
    dBarrier2 = m_pOneTouch->GetBarrier() * (1 - dShift);
  else
    dBarrier2 = m_pOneTouch->GetBarrier() * (1 + dShift);

  // Create a new one touch from the original one:
  shared_ptr<finance::OneTouch>
    pOneTouch2(new finance::OneTouch
                   ( m_pOneTouch->GetMaturityDate(),
                     dBarrier2, m_pOneTouch->GetBarrierType(),
                     m_pOneTouch->GetRebateType() 
              ));

  pOneTouch2->SetSessionData(m_pOneTouch->GetSessionData());

  // Price this new one touch
  shared_ptr<finance::ModelOutput> mo2(m_pModel->Compute(*pOneTouch2));

  // Price at first given OneTouch, setting analysis date to valuation date
  shared_ptr<finance::ComputationalFlags> pFlags(new finance::ComputationalFlags);
  pFlags->SetAnalysisDate(m_pOneTouch->GetSessionData()->GetValuationDate());

  m_pOneTouch->SetComputationalFlags(pFlags);
  
  shared_ptr<finance::ModelOutput> mo(m_pModel->Compute(*m_pOneTouch));
  
  // Compare the price with the interpolated price
  double dSpot2 = m_pOneTouch->GetBarrier() / pOneTouch2->GetBarrier()
                * m_pOneTouch->GetSessionData()->GetSpotSharePrice();

  double dPrice2;

  // TODO: add a function in MO to get price for a given spot at analysis date
  numeric::Interpolate(&mo->GetSpotsAtAnalysisDate()[0],
                       &mo->GetPricesAtAnalysisDate()[0],
                       mo->GetSpotsAtAnalysisDate().size(),
                       &dSpot2, &dPrice2, 1);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(dPrice2, mo2->GetPrice(), dPrecision);
}

void OneTouchTest::HomogeneityImplementation()
{
  if ( !HomogeneityCanApply() )
    return;

  const size_t nNb = 10;

  double dBarrier0 = m_pOneTouch->GetBarrier();
  double dSpot0 = m_pOneTouch->GetSessionData()->GetSpotSharePrice();
  double dShift = 0.8 * fabs(  dSpot0 - m_pOneTouch->GetBarrier() );
  double dStep = dShift / nNb;
  double dSpot = dSpot0 + dShift;
  double dBarrier;

  // create a one touch list sorted with barriers decreasing
  finance::OneTouches::Elements oneTouchList;

  shared_ptr<finance::ComputationalFlags> 
    pFlags(new finance::ComputationalFlags);
  
  pFlags->ActivateAllSensitivities(true);

  std::vector<double> pdSpots;
  std::vector<double> pdPrices;
  std::vector<double> pdSensitivities;
  for (size_t nIdx = 0; dSpot >= dSpot0 - dShift; nIdx++)
  {
    pdSpots.push_back(dSpot);
    dBarrier = dSpot / dSpot0 * dBarrier0;

    shared_ptr<finance::OneTouch>
      pOneTouch(new finance::OneTouch
                             ( m_pOneTouch->GetMaturityDate(),
                               dBarrier, 
                               m_pOneTouch->GetBarrierType(),
                               m_pOneTouch->GetRebateType()
                             )
               );

    // just set some prices
    double dMarketPrice = 0.5 * (nIdx + 1) / (2 * nNb);
    pOneTouch->SetMarketPrice(dMarketPrice);

    pOneTouch->SetSessionData( m_pOneTouch->GetSessionData() );

    pOneTouch->SetComputationalFlags(pFlags);

    shared_ptr<finance::ModelOutput>
      pOutput( m_pModel->Compute(*pOneTouch) );

    double dPrice = pOutput->GetPrice();

    pdPrices.push_back(dPrice);

    oneTouchList.push_back(pOneTouch);

    dSpot -= dStep;

    shared_ptr<NumOutput>
      pNumOutput( static_pointer_cast<NumOutput>( pOutput->GetNumOutput() ) );

    CPPUNIT_ASSERT( pNumOutput->HasSensitivities() );

    const std::vector<double>& 
      pdSensitivitiesTmp(pNumOutput->GetSensitivities());

    if ( pdSensitivities.empty() )
      pdSensitivities.resize(pdSensitivitiesTmp.size());
    else
      CPPUNIT_ASSERT( pdSensitivities.size() == pdSensitivitiesTmp.size() );

    for (size_t n = 0; n < pdSensitivities.size(); n++)
      pdSensitivities[n] += (dPrice - dMarketPrice) * pdSensitivitiesTmp[n];
  }

  finance::OneTouches oneTouches(oneTouchList);

  shared_ptr<MultiOutput>
    pOutput( PriceOneTouches(*m_pModel, oneTouches) );

  const std::vector<double>& pdNewPrices(pOutput->GetPrices());

  CPPUNIT_ASSERT( pdNewPrices.size() == pdPrices.size() );

  for (size_t n = 0; n < pdPrices.size(); n++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(pdNewPrices[n], pdPrices[n], dPrecision);

  CPPUNIT_ASSERT( pOutput->IsSensitivityOnObjectif() ); 
    
  shared_ptr<NumOutput>
    pNumOutput( static_pointer_cast<NumOutput>( pOutput->GetNumOutput() ) );

  const std::vector<double>&
    pdNewSensitivities(pNumOutput->GetSensitivities());
  CPPUNIT_ASSERT( pdNewSensitivities.size() == pdSensitivities.size() );
  
  for (size_t n = 0; n < pdSensitivities.size(); n++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(pdNewSensitivities[n], pdSensitivities[n], dPrecision);
}

} // namespace OneTouchTest

} // namespace ito33

} // namespace hg
