/////////////////////////////////////////////////////////////////////////////
// Name:        tests/bondlike/testyieldtomaturity.cpp
// Purpose:     test file for the yield-to-maturity
// Author:      Nabil
// Created:     2005/02/17
// RCS-ID:      $Id: testyieldtomaturity.cpp,v 1.16 2006/08/19 23:22:40 wang Exp $
// Copyright:   (c) 2005- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------
#include "ito33/beforestd.h"
#include <iostream>
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/common.h"
#include "ito33/debug.h"
#include "ito33/exception.h"
#include "ito33/date.h"
#include "ito33/dateutils.h"
#include "ito33/list.h"
#include "ito33/vector.h"
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
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/numeraire.h"

#include "ito33/tests/testyieldtomaturity.h"

//#define TEST_YTMDEBUG

const double YTM_RELATIVE_ERROR = 1.e-4;

using namespace ito33;

using namespace ito33::finance;

shared_ptr<SessionData> InitSessionData()
{
   Date valuationDate = Date(2003, Date::Jan, 14);

   shared_ptr<Numeraire> pCurrency( new Numeraire("USD"));
  // Setup the equity, and attach to session
  shared_ptr<Equity> pEquity(new Equity(95., pCurrency));
  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>( new YieldCurveFlat(0.0) ) );
  
  // Setup the issuer, and attach to the session  
  shared_ptr<RateData> pRateData( new RateData() );
  pRateData->SetYieldCurve( pCurrency, shared_ptr<YieldCurve>( new YieldCurveFlat(0.04) ) );
  
  shared_ptr<SessionData> 
    pSessionData( new SessionData(pRateData, pEquity, valuationDate) );

  return pSessionData;
}

void YTMAtCouponDates(Frequency frequency)
{
  Date
    issueDate = Date(2002, Date::Jan, 14),
    firstPaymentDate = Date(2003, Date::Jan, 14),
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
                       dRecoveryRate)
      );
  
  // Add cash Flow Stream ============================================

  double dCouponRate = 0.02;

  shared_ptr<CashFlowStream>
    pInterests( new CashFlowStreamUniform
                    (
                      issueDate, firstPaymentDate, maturityDate, dCouponRate,
                      Date::DayCountConvention_30360, frequency 
                    )
              );
    
  bc->SetCashDistribution(pInterests);
  
  // Create the session =====================================================
    
  shared_ptr<SessionData> pSessionData = InitSessionData();

  // Create the bond ========================================================
  
  shared_ptr<finance::Bond> 
    pBond( new finance::Bond(bc) );

  pBond->SetSessionData(pSessionData);

  // Test ==================================================================

  double dYTM;

  dYTM = pBond->ComputeYieldToMaturity
                ( dNominal, frequency, Date::DayCountConvention_30360 );
  
  CPPUNIT_ASSERT( fabs(dYTM - dCouponRate) 
                  < YTM_RELATIVE_ERROR * dCouponRate );
  
# ifdef TEST_YTMDEBUG
  std::cout << std::endl << pSession->GetValuationDate().Format("%Y/%m/%d") 
      << std::endl;  
    std::cout << "=>" << fabs(dYTM - dCouponRate) << std::endl;    
# endif
  
  CashFlowStream::Elements::const_iterator iter = pInterests->begin();
/*      
  for(iter = pInterests->begin(); 
      iter != pInterests->end() && iter->first < maturityDate; 
      ++iter)
 */
  {
    pSessionData->SetValuationDate(iter->first);

#   ifdef TEST_YTMDEBUG
    std::cout << std::endl << (iter->first).Format("%Y/%m/%d") << std::endl;    
#   endif
    
    dYTM = pBond->ComputeYieldToMaturity
                  ( dNominal, frequency, Date::DayCountConvention_30360 );
 
#   ifdef TEST_YTMDEBUG   
    std::cout << "=>" << fabs(dYTM - dCouponRate) << std::endl;
#   endif

    CPPUNIT_ASSERT( fabs(dYTM - dCouponRate) 
                    < YTM_RELATIVE_ERROR * dCouponRate );  
    
    double dYTM = pBond->ComputeYieldToMaturity
                  ( 20 * dNominal, frequency, Date::DayCountConvention_30360 );

    CPPUNIT_ASSERT( dYTM < 0 );
  }

}

void YieldToMaturityTest::YTMAtAnnualCouponDates()
{
  YTMAtCouponDates(Frequency_Annual);
}

void YieldToMaturityTest::YTMAtSemiAnnualCouponDates()
{
  YTMAtCouponDates(Frequency_SemiAnnual);
}

void YieldToMaturityTest::YTMAtQuarterCouponDates()
{
  YTMAtCouponDates(Frequency_Quarterly);
}

void YieldToMaturityTest::YTMAtBiMonthCouponDates()
{
  YTMAtCouponDates(Frequency_BiMonthly);
}

void YieldToMaturityTest::YTMAtMonthCouponDates()
{
  YTMAtCouponDates(Frequency_Monthly);
}

void YieldToMaturityTest::YTMForAZeroCouponBond()
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
                       dRecoveryRate)
      );

  Frequency cmpFreq = Frequency_Annual;

  bc->SetYieldCompoundingFrequency(cmpFreq);

  bc->SetYieldDayCountConvention(Date::DayCountConvention_ActAct);
  
  // Create the session =====================================================
    
  shared_ptr<SessionData> pSessionData = InitSessionData();

  // Create the bond ========================================================
  
  shared_ptr<finance::Bond> 
    pBond( new finance::Bond(bc) );

  pBond->SetSessionData(pSessionData);

  // Test ==================================================================

  double dYTM;

  dYTM = pBond->ComputeYieldToMaturity
                ( dNominal, cmpFreq, Date::DayCountConvention_ActAct );
  
  CPPUNIT_ASSERT( fabs(dYTM) < YTM_RELATIVE_ERROR );
  
# ifdef TEST_YTMDEBUG
  std::cout << std::endl << pSession->GetValuationDate().Format("%Y/%m/%d") 
      << std::endl;  
    std::cout << "=>" << fabs(dYTM) << std::endl;    
# endif

}

void YieldToMaturityTest::YTMForAConvertibleBond()
{
  
  Date
    issueDate = Date(2002, Date::Jan, 14),
    firstPaymentDate = Date(2003, Date::Jan, 14),
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
                       dRecoveryRate)
      );
  
  // Add cash Flow Stream ============================================

  double dCouponRate = 0.02;

  Frequency cmpFreq = Frequency_Annual;

  shared_ptr<CashFlowStream>
    pInterests( new CashFlowStreamUniform
                    (
                      issueDate, firstPaymentDate, maturityDate, dCouponRate,
                      Date::DayCountConvention_30360, cmpFreq 
                    )
              );
    
  bc->SetCashDistribution(pInterests);
  
  // Create the session =====================================================
    
  shared_ptr<SessionData> pSessionData = InitSessionData();

  // Create the CB ========================================================
  
  shared_ptr<ConversionSchedule> conv( new ConversionSchedule() );
  conv->SetKeepAccrued(false);
  
  conv->AddConversionPeriod(shared_ptr<ConversionPeriod>
                            ( new ConversionPeriod(issueDate, maturityDate, 1.) ));

  shared_ptr<finance::ConvertibleBond> 
    pCB( new finance::ConvertibleBond(bc, conv) );

  pCB->SetSessionData(pSessionData);

  // Test ==================================================================

  double dYTM;

  dYTM = pCB->ComputeYieldToMaturity
              ( dNominal, cmpFreq, Date::DayCountConvention_30360 );
  
  CPPUNIT_ASSERT( fabs(dYTM - dCouponRate) 
                  < YTM_RELATIVE_ERROR * dCouponRate );
  
# ifdef TEST_YTMDEBUG
    std::cout << std::endl << pSession->GetValuationDate().Format("%Y/%m/%d") 
        << std::endl;  
    std::cout << "=>" << fabs(dYTM - dCouponRate) << std::endl;    
# endif
  
  CashFlowStream::Elements::const_iterator iter;
      
  for(iter = pInterests->begin(); 
      iter != pInterests->end() && iter->first < maturityDate; 
      ++iter)
  {
    pSessionData->SetValuationDate(iter->first);

#   ifdef TEST_YTMDEBUG
    std::cout << std::endl << (iter->first).Format("%Y/%m/%d") << std::endl;    
#   endif
    
    dYTM = pCB->ComputeYieldToMaturity
                ( dNominal, cmpFreq, Date::DayCountConvention_30360 );
 
#   ifdef TEST_YTMDEBUG   
    std::cout << "=>" << fabs(dYTM - dCouponRate) << std::endl;
#   endif

    CPPUNIT_ASSERT( fabs(dYTM - dCouponRate) 
                    < YTM_RELATIVE_ERROR * dCouponRate );  
  }

}
