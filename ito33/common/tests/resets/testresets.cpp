/////////////////////////////////////////////////////////////////////////////
// Name:        testresets.cpp
// Purpose:     Acceptance test for resets  
// Author:      Yann d'Halluin
// Created:     17/11/2004
// RCS-ID:      $Id: testresets.cpp,v 1.19 2006/08/19 23:22:40 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------
#include "ito33/beforestd.h"
#include <vector>
#include "ito33/afterstd.h"
#include "ito33/date.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/bondlike/conversionpricereset.h"
#include "ito33/finance/bondlike/resetconversionschedule.h"
#include "ito33/finance/bondlike/resetflooredby.h"
#include "ito33/finance/bondlike/reset.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/callperiod.h"

#include "ito33/tests/utilexml.h"

#include "ito33/tests/testresets.h"

#include "ito33/xml/write.h"

using namespace ito33;
using namespace ito33::finance;

// ----------------------------------------------------------------------------
// Conversion schedule reset
// ----------------------------------------------------------------------------


void ResetTest::ConversionScheduleExists()
{
  Date issueDate(2000,Date::Jan,1);
  Date maturityDate(2003,Date::Jan,1);

  double dIssuePrice     = 1;
  double dParValue       = 110;
  double dRedemptionRate = 1;
  double dRecoveryRate   = 0.;  

  shared_ptr<BondTerms> 
      pBondTerms( new BondTerms(issueDate, dIssuePrice,
                         maturityDate, dParValue, dRedemptionRate,
                         dRecoveryRate)
                );
 
  Reset reset(pBondTerms, shared_ptr<ResetConversionSchedule>());

}//ResetTest::ConversionScheduleExists()

void ResetTest::ContingentCallNotSupported()
{
  Date startDate(2001,Date::Jan,1);
  Date endDate(2002,Date::Jan,1);
  Date issueDate(2000,Date::Jan,1);
  Date maturityDate(2003,Date::Jan,1);

  ResetFlooredBy flooredBy = ResetFlooredBy_PrevailingConversionPrice;

  double dInitialConvPrice = 1.0;
  double dCurrentConvPrice = 1.0;

  shared_ptr<ResetConversionSchedule> 
    pResetConvSch( new ResetConversionSchedule
                   ( startDate, endDate, dInitialConvPrice, 
                    dCurrentConvPrice, flooredBy ) 
                 );

    double dFloorRate = 1.0;
  Date resetDate(2001, Date::Mar, 1);

  shared_ptr<ConversionPriceReset> 
    pConvPriceReset( new ConversionPriceReset(resetDate, dFloorRate) );

  pResetConvSch->AddConversionPriceReset(pConvPriceReset);

   double dIssuePrice     = 1;
   double dParValue       = 110;
   double dRedemptionRate = 1;
   double dRecoveryRate   = 0.;  

    shared_ptr<BondTerms> 
      pBondTerms( new BondTerms(issueDate, dIssuePrice,
                         maturityDate, dParValue, dRedemptionRate,
                         dRecoveryRate)
                );

    Reset reset(pBondTerms,pResetConvSch);
  
    shared_ptr<CallSchedule> pCallSchedule( new CallSchedule() );
    pCallSchedule->SetTriggerCheckPeriod(20, 5);
    shared_ptr<CallPeriod> 
      pCallPeriod( CallPeriod::CreateWithStrike
                               (issueDate, maturityDate, 100.0) );
    pCallSchedule->AddCallPeriod(pCallPeriod);

    reset.SetCallSchedule(pCallSchedule);

}//ResetTest::ContigentCallNotSupported()


void ResetTest::ConvSchStartDateAfterBondIssueDate()
{

  Date startDate(1999,Date::Jan,1);
  Date endDate(2002,Date::Jan,1);
  Date issueDate(2000,Date::Jan,1);
  Date maturityDate(2003,Date::Jan,1);

  ResetFlooredBy flooredBy = ResetFlooredBy_PrevailingConversionPrice;

  double dInitialConvPrice = 1.0;
  double dCurrentConvPrice = 1.0;

  shared_ptr<ResetConversionSchedule> 
    pResetConvSch( new ResetConversionSchedule
                   ( startDate, endDate, dInitialConvPrice, 
                    dCurrentConvPrice, flooredBy ) 
                 );

  double dFloorRate = 1.0;
  Date resetDate(2001, Date::Mar, 1);

  shared_ptr<ConversionPriceReset> 
    pConvPriceReset( new ConversionPriceReset(resetDate, dFloorRate) );

  pResetConvSch->AddConversionPriceReset(pConvPriceReset);

   double dIssuePrice     = 1;
   double dParValue       = 110;
   double dRedemptionRate = 1;
   double dRecoveryRate   = 0.;  

    shared_ptr<BondTerms> 
      pBondTerms( new BondTerms(issueDate, dIssuePrice,
                         maturityDate, dParValue, dRedemptionRate,
                         dRecoveryRate) );

    Reset reset(pBondTerms,pResetConvSch);


}//ResetTest::ConvSchStartDateAfterBondIssueDate()


void ResetTest::ConvSchEndDateBeforeBondMaturityDate()
{

  Date startDate(2001,Date::Jan,1);
  Date endDate(2004,Date::Jan,1);
  Date issueDate(2000,Date::Jan,1);
  Date maturityDate(2003,Date::Jan,1);

  ResetFlooredBy flooredBy = ResetFlooredBy_PrevailingConversionPrice;

  double dInitialConvPrice = 1.0;
  double dCurrentConvPrice = 1.0;

  shared_ptr<ResetConversionSchedule> 
    pResetConvSch( new ResetConversionSchedule
                   ( startDate, endDate, dInitialConvPrice, 
                    dCurrentConvPrice, flooredBy ) 
                 );

  double dFloorRate = 1.0;
  Date resetDate(2001, Date::Mar, 1);

  shared_ptr<ConversionPriceReset> 
    pConvPriceReset( new ConversionPriceReset(resetDate, dFloorRate) );

  pResetConvSch->AddConversionPriceReset(pConvPriceReset);


   double dIssuePrice     = 1;
   double dParValue       = 110;
   double dRedemptionRate = 1;
   double dRecoveryRate   = 0.;  

    shared_ptr<BondTerms> 
      pBondTerms( new BondTerms(issueDate, dIssuePrice,
                         maturityDate, dParValue, dRedemptionRate,
                         dRecoveryRate)
                );

    Reset reset(pBondTerms,pResetConvSch);

}//ResetTest::ConvSchEndDateBeforeBondMaturityDate()

void ResetTest::ConvSchContainsAtLeastOneResetDate()
{
  Date startDate(2001,Date::Jan,1);
  Date endDate(2002,Date::Jan,1);
  Date issueDate(2000,Date::Jan,1);
  Date maturityDate(2003,Date::Jan,1);

  ResetFlooredBy flooredBy = ResetFlooredBy_PrevailingConversionPrice;

  double dInitialConvPrice = 1.0;
  double dCurrentConvPrice = 1.0;

  shared_ptr<ResetConversionSchedule> 
    pResetConvSch( new ResetConversionSchedule
                   ( startDate, endDate, dInitialConvPrice, 
                    dCurrentConvPrice, flooredBy ) 
                 );

  double dIssuePrice     = 1;
  double dParValue       = 110;
  double dRedemptionRate = 1;
  double dRecoveryRate   = 0.;  

  shared_ptr<BondTerms> 
    pBondTerms( new BondTerms(issueDate, dIssuePrice,
                        maturityDate, dParValue, dRedemptionRate,
                        dRecoveryRate)
              );

  Reset reset(pBondTerms,pResetConvSch);

  reset.Validate();

}//ResetTest::ConvSchContainsAtLeastOneResetDate()


void ResetTest::LastResetDateDifferentFromBondMaturity()
{
  Date startDate(2001,Date::Jan,1);
  Date endDate(2002,Date::Jan,1);
  Date issueDate(2000,Date::Jan,1);
  Date maturityDate(2003,Date::Jan,1);

  ResetFlooredBy flooredBy = ResetFlooredBy_PrevailingConversionPrice;

  double dInitialConvPrice = 1.0;
  double dCurrentConvPrice = 1.0;

  shared_ptr<ResetConversionSchedule> 
    pResetConvSch( new ResetConversionSchedule
                   ( startDate, endDate, dInitialConvPrice, 
                    dCurrentConvPrice, flooredBy ) 
                 );

  double dFloorRate = 1.0;
  Date resetDate(2003, Date::Mar, 1);

  shared_ptr<ConversionPriceReset> 
    pConvPriceReset( new ConversionPriceReset(resetDate, dFloorRate) );

  pResetConvSch->AddConversionPriceReset(pConvPriceReset);
 
  double dIssuePrice     = 1;
  double dParValue       = 110;
  double dRedemptionRate = 1;
  double dRecoveryRate   = 0.;  

   shared_ptr<BondTerms> 
      pBondTerms( new BondTerms(issueDate, dIssuePrice,
                         maturityDate, dParValue, dRedemptionRate,
                         dRecoveryRate)
                );

   Reset reset(pBondTerms,pResetConvSch);

}//ResetTest::LastResetDateDifferentFromBondMaturity()

void ResetTest::Dump()
{

  Date startDate(2001,Date::Jan,1);
  Date endDate(2002,Date::Jan,1);
  Date issueDate(2000,Date::Jan,1);
  Date maturityDate(2003,Date::Jan,1);

  ResetFlooredBy flooredBy = ResetFlooredBy_PrevailingConversionPrice;

  double dInitialConvPrice = 1.0;
  double dCurrentConvPrice = 1.0;

  shared_ptr<ResetConversionSchedule> 
    pResetConvSch( new ResetConversionSchedule
                   ( startDate, endDate, dInitialConvPrice, 
                    dCurrentConvPrice, flooredBy ) 
                 );

  double dFloorRate = 1.0;
  Date resetDate(2001, Date::Mar, 1);

  shared_ptr<ConversionPriceReset> 
    pConvPriceReset( new ConversionPriceReset(resetDate, dFloorRate) );

  pResetConvSch->AddConversionPriceReset(pConvPriceReset);

  double dIssuePrice     = 1;
  double dParValue       = 110;
  double dRedemptionRate = 1;
  double dRecoveryRate   = 0.;  

  shared_ptr<BondTerms> 
    pBondTerms( new BondTerms(issueDate, dIssuePrice,
                        maturityDate, dParValue, dRedemptionRate,
                        dRecoveryRate) 
              );

  Reset reset(pBondTerms, pResetConvSch);

  std::ostringstream oss;

  ExpectedXML expected(oss,
    "<?xml version=\"1.0\"?>"
    "<root>\n"
    "<reset>\n"
    "<bond_terms>\n"
    "<issue_date>2000-01-01</issue_date>\n"
    "<issue_price>1</issue_price>\n"
    "<maturity>2003-01-01</maturity>\n"
    "<nominal>110</nominal>\n"
    "<recovery_rate>0</recovery_rate>\n"
    "<redemption_price>1</redemption_price>\n"
    "</bond_terms>\n"
    "<new_share>0</new_share>\n"
    "<reset_conversion_schedule>\n"
    "<keep_accrued>0</keep_accrued>\n"
    "<forfeit_coupon>0</forfeit_coupon>\n"
    "<conversion_start>2001-01-01</conversion_start>\n"
    "<conversion_end>2002-01-01</conversion_end>\n"
    "<initial_conversion_price>1</initial_conversion_price>\n"
    "<current_conversion_price>1</current_conversion_price>\n"
    "<cash_value>0</cash_value>\n"
    "<floored_by>prevailing_conversion_price</floored_by>\n"
    "<conversion_price_resets>\n"
    "<conversion_price_reset>\n"
    "<date>2001-03-01</date>\n"
    "<cap_rate>1</cap_rate>\n"
    "<floor_rate>1</floor_rate>\n"
    "<multiplier>1</multiplier>\n"
    "</conversion_price_reset>\n"
    "</conversion_price_resets>\n"
    "</reset_conversion_schedule>\n"
    "</reset>\n"
    "</root>\n");


  ito33::XML::RootTag root("root",oss);

  reset.Dump(root);

}//ResetTest::Dump()
