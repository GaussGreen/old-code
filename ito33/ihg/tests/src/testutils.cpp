/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/src/testutils.cpp
// Purpose:     Implenation of helper functions for coreinterface tests
// Created:     2006/03/23
// RCS-ID:      $Id: testutils.cpp,v 1.3 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/date.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/numeraire.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/issuer.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/yieldcurve_flat.h"

namespace ito33
{

namespace ihg
{

shared_ptr<finance::SessionData> MakeSessionData(double dSpot)
{

  // SessionData needs RateData, Equity, and valuation date

  // The valuation date
  Date valuationDate(2003, Date::Feb, 1);

  // The equity
  Date fiscalYearStart(2003, Date::Jul, 1);
  shared_ptr<finance::Issuer> pIssuer(new finance::Issuer);
  pIssuer->SetFiscalYearStartDate(fiscalYearStart);

  shared_ptr<finance::Numeraire> pCurrencyEUR(new finance::Numeraire("EUR") );

  shared_ptr<finance::Equity> pEquity(new finance::Equity(dSpot, pCurrencyEUR));

  pEquity->SetIssuer(pIssuer);

  shared_ptr<finance::YieldCurve> pyf(new finance::YieldCurveFlat(0.01));    
  pEquity->SetBorrowCurve(pyf);

  // The rate data  
  shared_ptr<finance::RateData> pRateData(new finance::RateData);

  shared_ptr<finance::YieldCurve> pycEUR( new finance::YieldCurveFlat(0.05) );
  pRateData->SetYieldCurve(pCurrencyEUR, pycEUR);

  // For cross-currency examples
  shared_ptr<finance::SpotFXRates> pSpotFXRates(new finance::SpotFXRates);

  double dSpotFXRate = 0.8; 
  shared_ptr<finance::Numeraire> pCurrencyUSD(new finance::Numeraire("USD"));
  
  pSpotFXRates->SetFXRate(pCurrencyEUR, pCurrencyUSD, dSpotFXRate);

  shared_ptr<finance::YieldCurve> pycUSD( new finance::YieldCurveFlat(0.02) );
  pRateData->SetYieldCurve(pCurrencyUSD, pycUSD);

  pRateData->SetSpotFXRates(pSpotFXRates);

  // Create the SessionData and return
  shared_ptr<finance::SessionData> 
    pSessionData(new finance::SessionData(pRateData, pEquity, valuationDate));

  return pSessionData;

}

} // namespace ihg

} // namespace ito33
