/////////////////////////////////////////////////////////////////////////////
// Name:        tests/bondlike/testnewshare.cpp
// Purpose:     test file for the new share
// Author:      Nabil
// Created:     2005/04/11
// RCS-ID:      $Id: testnewshare.cpp,v 1.7 2006/08/19 23:22:40 wang Exp $
// Copyright:   (c) 2005- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/common.h"
#include "ito33/debug.h"
#include "ito33/exception.h"
#include "ito33/date.h"
#include "ito33/dateutils.h"
#include "ito33/array.h"
#include "ito33/cppunit.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/bondlike/bond.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/conversionschedule.h"

#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/issuer.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/yieldcurve_flat.h"

#include "ito33/tests/testnewshare.h"

using namespace ito33;

using namespace ito33::finance;

shared_ptr<SessionData> NewShareTest::InitSessionData()
{
   Date valuationDate = Date(2003, Date::Jan, 14);

   shared_ptr<Numeraire> pCurrency( new Numeraire("EUR") );

  // Setup the equity, and attach to session
  shared_ptr<Equity> pEquity(new Equity(95., pCurrency));
  
  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>(new YieldCurveFlat(0.0)) );
  
  // Setup the issuer, and attach to the session  
  shared_ptr<RateData> pRateData(new RateData() );
  pRateData->SetYieldCurve(pCurrency, 
    shared_ptr<YieldCurve>(new YieldCurveFlat(0.04) ) );
  
  shared_ptr<SessionData> 
    pSessionData( new SessionData(pRateData, pEquity, valuationDate) );

  return pSessionData;
}

void NewShareTest::NewShareAndExchangeable()
{
  
  Date
    issueDate = Date(2002, Date::Jan, 14),
    maturityDate = Date(2053, Date::Jan, 14);

  double
    dIssuePrice = 1,
    dNominal = 110,
    dRedemptionRate = 1,
    dRecoveryRate = 0.;

  // Create the bondterms ============================================

  shared_ptr<BondTerms> 
    bc( new BondTerms(issueDate, dIssuePrice,
                      maturityDate, dNominal, dRedemptionRate,
                      dRecoveryRate) );
  
  // Create the session =====================================================
    
  shared_ptr<SessionData> pSessionData = InitSessionData();

  // Create the CB ========================================================
  
  shared_ptr<ConversionSchedule> conv( new ConversionSchedule() );
  conv->SetKeepAccrued(false);
  
  conv->AddConversionPeriod
        (shared_ptr<ConversionPeriod>
         (new ConversionPeriod(issueDate, maturityDate, 1.)));

  shared_ptr<ConvertibleBond> 
    pCB( new ConvertibleBond(bc, conv) );

  pCB->SetSessionData(pSessionData);

  pCB->SetConvertIntoNewShare(true);

  pCB->SetExchangeable(false);
}

void NewShareTest::ExchangeableAndNewShare()
{
  
  Date
    issueDate = Date(2002, Date::Jan, 14),
    maturityDate = Date(2053, Date::Jan, 14);

  double
    dIssuePrice = 1,
    dNominal = 110,
    dRedemptionRate = 1,
    dRecoveryRate = 0.;

  // Create the bondterms ============================================

  shared_ptr<BondTerms> 
    bc( new BondTerms(issueDate, dIssuePrice,
                      maturityDate, dNominal, dRedemptionRate,
                      dRecoveryRate) );
  
  // Create the session =====================================================
    
  shared_ptr<SessionData> pSessionData = InitSessionData();

  // Create the CB ========================================================
  
  shared_ptr<ConversionSchedule> conv( new ConversionSchedule() );
  conv->SetKeepAccrued(false);
  
  conv->AddConversionPeriod
        (shared_ptr<ConversionPeriod>
         ( new ConversionPeriod(issueDate, maturityDate, 1.)));

  shared_ptr<ConvertibleBond> pCB( new ConvertibleBond(bc, conv) );

  pCB->SetSessionData(pSessionData);

  pCB->SetExchangeable(false);

  pCB->SetConvertIntoNewShare(true);
}
