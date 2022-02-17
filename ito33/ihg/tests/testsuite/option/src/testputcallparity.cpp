/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/common/putcallparity/testputcallparity.cpp
// Purpose:     Acceptance test for put call parity helper functions
// Author:      ITO 33
// Created:     13/05/2005
// RCS-ID:      $Id: testputcallparity.cpp,v 1.13 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/sharedptr.h"

#include "ito33/finance/option.h"
#include "ito33/finance/putcallparity.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/yieldcurve_annuallycompounded.h"
#include "ito33/finance/computationalflags.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilitypower.h"
#include "ito33/ihg/hazardratepower.h"

#include "testputcallparity.h"

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::ihg;


shared_ptr<ihg::TheoreticalModel> 
PutCallParityTest::ConstructModel()
{
  // Construct a model with a non-trivial hazard and volatility to make
  // testing harder
  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);
    
  shared_ptr<ihg::Volatility> pVol( new ihg::VolatilityPower(0.2, 0.2, 100.0) );

  pModel->SetVolatility( pVol );

  shared_ptr<ihg::HazardRate> pHR( new HazardRatePower(0.15, 1.2, 100.0) );
  pModel->SetHazardRate( pHR );

  return pModel;

}

shared_ptr<finance::SessionData> 
PutCallParityTest::ConstructSessionData(int iDivType)
{
  
  Date valuationDate(2003, Date::Feb, 1);

  // Use non-trivial yield curves to make a harder test
  shared_ptr<YieldCurveAnnuallyCompounded> pyc( new
    YieldCurveAnnuallyCompounded( Date(2002, Date::Jan, 1), 10) );

  pyc->AddLeg(90, 0.02);
  pyc->AddLeg(180, 0.0230);
  pyc->AddLeg(270, 0.0250);
  pyc->AddLeg(365, 0.0269);
  pyc->AddLeg(545, 0.0287);
  pyc->AddLeg(2*365, 0.0303);
  pyc->AddLeg(3*365, 0.0318);
  pyc->AddLeg(5*365, 0.0331);
  pyc->AddLeg(7*365, 0.0343);
  pyc->AddLeg(10*365, 0.0353);

  shared_ptr<Numeraire> pNumeraire(new Numeraire("EUR"));

  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pNumeraire, pyc);  
  
  shared_ptr<Equity> pEquity(new Equity(100.0, pNumeraire));

  shared_ptr<YieldCurveAnnuallyCompounded> pyf( new
    YieldCurveAnnuallyCompounded( Date(2002, Date::Jan, 1), 10));

  pyf->AddLeg(90, 0.01);
  pyf->AddLeg(180, 0.0105);
  pyf->AddLeg(270, 0.013);
  pyf->AddLeg(365, 0.0135);
  pyf->AddLeg(545, 0.016);
  pyf->AddLeg(2*365, 0.017);
  pyf->AddLeg(3*365, 0.0185);
  pyf->AddLeg(5*365, 0.02);
  pyf->AddLeg(7*365, 0.0205);
  pyf->AddLeg(10*365, 0.022);

/*
  pyf->AddLeg(90, 0.0);
  pyf->AddLeg(180, 0.0);
  pyf->AddLeg(270, 0.0);
  pyf->AddLeg(365, 0.0);
  pyf->AddLeg(545, 0.0);
  pyf->AddLeg(2*365, 0.0);
  pyf->AddLeg(3*365, 0.0);
  pyf->AddLeg(5*365, 0.0);
  pyf->AddLeg(7*365, 0.0);
  pyf->AddLeg(10*365, 0.0);
*/
  pEquity->SetBorrowCurve(pyf);
  
  switch (iDivType)
  {
  case 1:
    {
    shared_ptr<Dividends> pDividends( new Dividends());
    pDividends->AddCash( Date(2003, Date::Feb, 15), 5.0);
    pDividends->AddCash( Date(2004, Date::Feb, 15), 4.0);
    pDividends->AddCash( Date(2004, Date::Dec, 15), 8.0);
    pEquity->SetDividends(pDividends);
    }
    break;

  default:
    // No dividends
    break;
  }

  shared_ptr<SessionData> 
    pSessionData( new SessionData(pRateData, pEquity, valuationDate) );

  return pSessionData;
}


shared_ptr<finance::Option> 
PutCallParityTest::ConstructPutOption(shared_ptr<finance::SessionData> pSessionData,
                                      shared_ptr<ihg::TheoreticalModel> pModel)
{
  shared_ptr<Option> pOpt(new Option(96., Date(2005, Date::Jan, 1),
                                    Option_Put, ExerciseType_European));
  pOpt->SetSessionData(pSessionData);

  // Price the option and set the market price
  shared_ptr<finance::ModelOutput> pOutput = pModel->Compute(*pOpt);
  double dComputedPrice = pOutput->GetPrice();

  pOpt->SetMarketPrice(dComputedPrice);

  return pOpt;

}


void PutCallParityTest::TestConvertPutToCall1()
{

  // Construct with no dividends
  shared_ptr<SessionData> pSessionData = ConstructSessionData(0);
  shared_ptr<ihg::TheoreticalModel> pModel = ConstructModel();
  shared_ptr<Option> pPutOption = ConstructPutOption(pSessionData, pModel);

  bool bResult = IsPutCallParityValid(*pPutOption,*pModel);

  CPPUNIT_ASSERT( bResult );

} //PutCallParityTest::TestConvertPutToCall1()


void PutCallParityTest::TestConvertPutToCall2()
{
  
  // Construct with cash dividend
  shared_ptr<SessionData> pSessionData = ConstructSessionData(1);
  shared_ptr<ihg::TheoreticalModel> pModel = ConstructModel();
  shared_ptr<Option> pPutOption = ConstructPutOption(pSessionData, pModel);

  bool bResult = IsPutCallParityValid(*pPutOption,*pModel);

  CPPUNIT_ASSERT( bResult );

} //PutCallParityTest::TestConvertPutToCall2()


void PutCallParityTest::TestConvertPutToCall3()
{

  // Construct with dividends
  shared_ptr<SessionData> pSessionData = ConstructSessionData(1);
  //construct with hazard rate
  shared_ptr<ihg::TheoreticalModel> pModel = ConstructModel();
  shared_ptr<Option> pPutOption = ConstructPutOption(pSessionData, pModel);

  bool bResult = IsPutCallParityValid(*pPutOption,*pModel);

  CPPUNIT_ASSERT( bResult );

} //PutCallParityTest::TestConvertPutToCall3()