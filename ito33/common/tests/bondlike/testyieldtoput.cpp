/////////////////////////////////////////////////////////////////////////////
// Name:        tests/bondlike/testyieldtoput.cpp
// Purpose:     test file for the yield-to-put
// Author:      Nabil
// Created:     2005/02/23
// RCS-ID:      $Id: testyieldtoput.cpp,v 1.13 2006/08/19 23:22:40 wang Exp $
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
#include "ito33/finance/bondlike/putschedule.h"

#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/numeraire.h"

#include "ito33/tests/testyieldtoput.h"

//#define TEST_YTPDEBUG

const double YTP_RELATIVE_ERROR = 1.e-4;

using namespace ito33;

using namespace ito33::finance;

shared_ptr<SessionData> YieldToPutTest::InitSessionData()
{ 
  Date valuationDate = Date(2003, Date::Jan, 14);
  
  shared_ptr<Numeraire> pCurrency( new Numeraire("USD") );

  // Setup the equity, and attach to session
  shared_ptr<Equity> pEquity(new Equity(95., pCurrency) );
  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>( new YieldCurveFlat(0.0) ) );

  shared_ptr<RateData> pRateData( new RateData() ) ;
  pRateData->SetYieldCurve( pCurrency, 
    shared_ptr<YieldCurve>(new YieldCurveFlat(0.04) ) );
  
  shared_ptr<SessionData> 
    pSessionData( new SessionData(pRateData, pEquity, valuationDate) );

  return pSessionData;
}

void YieldToPutTest::YTPAtCouponDates(Frequency frequency)
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

  //  Add just one put at maturity
  shared_ptr<PutSchedule> puts ( new PutSchedule() );
  puts->AddPutWithStrike(maturityDate, 1.);

  pBond->SetPutSchedule(puts);

  // Test ==================================================================

  double dYTP;

  dYTP = pBond->ComputeYieldToPut
                ( dNominal, frequency, Date::DayCountConvention_30360 );
  
  CPPUNIT_ASSERT( fabs(dYTP - dCouponRate) 
                  < YTP_RELATIVE_ERROR * dCouponRate );
  
# ifdef TEST_YTPDEBUG
  std::cout << std::endl << pSession->GetValuationDate().Format("%Y/%m/%d") 
      << std::endl;  
    std::cout << "=>" << fabs(dYTP - dCouponRate) << std::endl;    
# endif
  
  CashFlowStream::Elements::const_iterator iter = pInterests->begin();
  /*    
  for(iter = pInterests->begin(); 
      iter != pInterests->end() && iter->first < maturityDate; 
      ++iter)
      */
  {
    pSessionData->SetValuationDate(iter->first);

#   ifdef TEST_YTPDEBUG
    std::cout << std::endl << (iter->first).Format("%Y/%m/%d") << std::endl;    
#   endif
    
    dYTP = pBond->ComputeYieldToPut
                  ( dNominal, frequency, Date::DayCountConvention_30360 );
 
#   ifdef TEST_YTPDEBUG   
    std::cout << "=>" << fabs(dYTP - dCouponRate) << std::endl;
#   endif

    CPPUNIT_ASSERT( fabs(dYTP - dCouponRate) 
                    < YTP_RELATIVE_ERROR * dCouponRate );  

    dYTP = pBond->ComputeYieldToPut
                  ( 10 * dNominal, frequency, Date::DayCountConvention_30360 );

    CPPUNIT_ASSERT( dYTP < 0);
  }

}

void YieldToPutTest::YTPAtAnnualCouponDates()
{
  YTPAtCouponDates(Frequency_Annual);
}

void YieldToPutTest::YTPAtSemiAnnualCouponDates()
{
  YTPAtCouponDates(Frequency_SemiAnnual);
}

void YieldToPutTest::YTPAtQuarterCouponDates()
{
  YTPAtCouponDates(Frequency_Quarterly);
}

void YieldToPutTest::YTPAtBiMonthCouponDates()
{
  YTPAtCouponDates(Frequency_BiMonthly);
}

void YieldToPutTest::YTPAtMonthCouponDates()
{
  YTPAtCouponDates(Frequency_Monthly);
}

void YieldToPutTest::YTPForAZeroCouponBond()
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
  
  //  Add just one put at maturity
  shared_ptr<PutSchedule> puts ( new PutSchedule() );
  puts->AddPutWithStrike(maturityDate, 1.);

  pBond->SetPutSchedule(puts);

  // Test ==================================================================

  double dYTP;

  dYTP = pBond->ComputeYieldToPut
                ( dNominal, cmpFreq, Date::DayCountConvention_ActAct );
  
  CPPUNIT_ASSERT( fabs(dYTP) < YTP_RELATIVE_ERROR );
  
# ifdef TEST_YTPDEBUG
  std::cout << std::endl << pSession->GetValuationDate().Format("%Y/%m/%d") 
      << std::endl;  
    std::cout << "=>" << fabs(dYTP) << std::endl;    
# endif

}

void YieldToPutTest::YTPForAConvertibleBond()
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
  
  conv->AddConversionPeriod
        ( shared_ptr<ConversionPeriod>
          ( new ConversionPeriod(issueDate, maturityDate, 1.) ) );

  shared_ptr<finance::ConvertibleBond> 
    pCB( new finance::ConvertibleBond(bc, conv) );

  pCB->SetSessionData(pSessionData);
  
  //  Add just one put at maturity
  shared_ptr<PutSchedule> puts ( new PutSchedule() );
  puts->AddPutWithStrike(maturityDate, 1.);

  pCB->SetPutSchedule(puts);

  // Test ==================================================================

  double dYTP;

  dYTP = pCB->ComputeYieldToPut
              ( dNominal, cmpFreq, Date::DayCountConvention_30360 );
  
  CPPUNIT_ASSERT( fabs(dYTP - dCouponRate) 
                  < YTP_RELATIVE_ERROR * dCouponRate );
  
# ifdef TEST_YTPDEBUG
    std::cout << std::endl << pSession->GetValuationDate().Format("%Y/%m/%d") 
        << std::endl;  
    std::cout << "=>" << fabs(dYTP - dCouponRate) << std::endl;    
# endif
  
  CashFlowStream::Elements::const_iterator iter;
      
  for(iter = pInterests->begin(); 
      iter != pInterests->end() && iter->first < maturityDate; 
      ++iter)
  {
    pSessionData->SetValuationDate(iter->first);

#   ifdef TEST_YTPDEBUG
    std::cout << std::endl << (iter->first).Format("%Y/%m/%d") << std::endl;    
#   endif
    
    dYTP = pCB->ComputeYieldToPut
                ( dNominal, cmpFreq, Date::DayCountConvention_30360 );
 
#   ifdef TEST_YTPDEBUG   
    std::cout << "=>" << fabs(dYTP - dCouponRate) << std::endl;
#   endif

    CPPUNIT_ASSERT( fabs(dYTP - dCouponRate) 
                    < YTP_RELATIVE_ERROR * dCouponRate );  
  }

}

void YieldToPutTest::YTPForBondWithoutPutSchedule()
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

  double dYield = 0.02;
  bc->SetAccretingBond (dYield);
  
  // Create the session =====================================================
    
  shared_ptr<SessionData> pSessionData = InitSessionData();

  // Create the bond ========================================================
  
  shared_ptr<finance::Bond> 
    pBond( new finance::Bond(bc) );

  pBond->SetSessionData(pSessionData);

  // Test ==================================================================

  double dYTP;

  dYTP = pBond->ComputeYieldToPut
                ( dNominal, cmpFreq, Date::DayCountConvention_ActAct ); 
}

void YieldToPutTest::YTPForConvertibleBondWithoutPutSchedule()
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
  
  conv->AddConversionPeriod
        ( shared_ptr<ConversionPeriod>
          ( new ConversionPeriod(issueDate, maturityDate, 1.) ) );

  shared_ptr<finance::ConvertibleBond> 
    pCB( new finance::ConvertibleBond(bc, conv) );

  pCB->SetSessionData(pSessionData);

  // Test ==================================================================

  double dYTP;

  dYTP = pCB->ComputeYieldToPut
              ( dNominal, cmpFreq, Date::DayCountConvention_30360 );

}
