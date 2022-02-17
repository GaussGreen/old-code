/////////////////////////////////////////////////////////////////////////////
// Name:        testcds.cpp
// Purpose:     Acceptance test for CDS 
// Author:      Yann d'Halluin
// Created:     17/11/2004
// RCS-ID:      $Id: testcds.cpp,v 1.14 2006/08/19 23:22:40 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------


#include "ito33/sharedptr.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/referencecds.h"
#include "ito33/finance/cashflowstream_uniform.h"

#include "ito33/tests/utilexml.h"

#include "ito33/tests/testcds.h"

#include "ito33/xml/write.h"

using namespace ito33;
using namespace ito33::finance;

	void CDSTest::NegativeRecoveryRate()
{
 
  double dRecoveryRate = -0.1;
  double dAnnualAmt    = .25;

  shared_ptr<finance::CashFlowStreamUniform>
	  pSpreadStream( new finance::CashFlowStreamUniform 
                          (
                            Date(2003, Date::Feb, 1),
                            Date(2003, Date::May, 1),
                            Date(2011, Date::May, 1),
                            dAnnualAmt,
                            Date::DayCountConvention_Act365,
                            finance::Frequency_BiMonthly
                          )
                 );

  
   CDS cds(dRecoveryRate, pSpreadStream);


} // end negative recovery rate

void CDSTest::RecoveryRateTooLarge()
{
 
  double dRecoveryRate = 2.;
  double dAnnualAmt    = .25;

  shared_ptr<finance::CashFlowStreamUniform>
	  pSpreadStream( new finance::CashFlowStreamUniform 
                          (
                            Date(2003, Date::Feb, 1),
                            Date(2003, Date::May, 1),
                            Date(2011, Date::May, 1),
                            dAnnualAmt,
                            Date::DayCountConvention_Act365,
							              finance::Frequency_BiMonthly
                          )
                 );


    CDS  cds(dRecoveryRate, pSpreadStream);

  
} // end large recovery rate


void CDSTest::NoCashFlowStream()
{
  double dRecoveryRate = 2.;
  
  CDS  cds(dRecoveryRate, shared_ptr<CashFlowStreamUniform>());
 
} // CDSTest::NoCashFlowStream
  
void CDSTest::Dump()
{
 
  double dRecoveryRate = 0.1;
  double dAnnualAmt     = .25;

  shared_ptr<finance::CashFlowStreamUniform>
	  pSpreadStream( new finance::CashFlowStreamUniform 
                          (
                            Date(2003, Date::Feb, 1),
                            Date(2003, Date::May, 1),
                            Date(2011, Date::May, 1),
                            dAnnualAmt,
                            Date::DayCountConvention_Act365,
                            finance::Frequency_BiMonthly
                          )
                 );

  
  CDS cds(dRecoveryRate, pSpreadStream);

  std::ostringstream oss;

  ExpectedXML expected(oss,
    "<?xml version=\"1.0\"?>"
    "<root>\n"
    "<cds>\n"
    "<spread_stream>\n"
    "<cash_flow_stream_uniform>\n"
    "<contracting_date>2003-02-01</contracting_date>\n"
    "<day_count_convention>act365</day_count_convention>\n"
    "<payment_frequency>bimonthly</payment_frequency>\n"
    "<first_date>2003-05-01</first_date>\n"
    "<last_date>2011-05-01</last_date>\n"
    "<annual_amount>0.25</annual_amount>\n"    
    "</cash_flow_stream_uniform>\n"
    "</spread_stream>\n"
    "<recovery_rate>0.1</recovery_rate>\n"
    "</cds>\n"
    "</root>\n");

  ito33::XML::RootTag root("root",oss);

  cds.Dump(root);

} //CDSTest::Dump

void CDSTest::SlidingCDS()
{
  // Today 
  Date today(2005, Date::Jul, 18);
  
  Date firstPaymentDate(2005, Date::Sep, 20);

  Date maturityDate(2006, Date::Sep, 20);

  // The frequency is usually quarterly
  Frequency freq = Frequency_Quarterly;

  // Use whatever day count convention
  Date::DayCountConvention dcc = Date::DayCountConvention_30360;

  // Create a one year sliding CDS
  ReferenceCDS cds(12, freq, dcc, 0.5);

  cds.SetSpread(0.10);

  // create underlying session data
  shared_ptr<Equity> pEquity(new Equity(shared_ptr<Numeraire>(new Numeraire("EUR"))));

  shared_ptr<SessionData> pSessionData(new SessionData(shared_ptr<RateData>(new RateData()),
                                                      pEquity,
                                                      today));
  cds.SetSessionData(pSessionData);

  CPPUNIT_ASSERT( cds.GetSpreadStream()->GetFirstPaymentDate() == firstPaymentDate );
  CPPUNIT_ASSERT( cds.GetMaturityDate() == maturityDate );

  // Also check the case where today is already after the last regular date
  // of this year
  today = Date(2005, Date::Dec, 23);
  firstPaymentDate = Date(2006, Date::Mar, 20);
  maturityDate = Date(2007, Date::Mar, 20);

  // Create a one year sliding CDS
  ReferenceCDS cds2(12, freq, dcc, 0.5);
  cds2.SetSpread(0.10);

  pSessionData = shared_ptr<SessionData>
                                     (new SessionData(shared_ptr<RateData>(new RateData()),
                                                      pEquity,
                                                      today));
  cds2.SetSessionData(pSessionData);

  CPPUNIT_ASSERT( cds2.GetSpreadStream()->GetFirstPaymentDate() == firstPaymentDate );
  CPPUNIT_ASSERT( cds2.GetMaturityDate() == maturityDate );
}
