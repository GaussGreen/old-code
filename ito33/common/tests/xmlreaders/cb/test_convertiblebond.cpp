/////////////////////////////////////////////////////////////////////////////
// Name:        tests/bondlike/test_convertiblebond.cpp
// Purpose:     test file for bond
// Author:      Zhang (converted to cppunit by David)
// Created:     24.06.04
// RCS-ID:      $Id: test_convertiblebond.cpp,v 1.19 2006/08/19 23:22:41 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------
#include <iostream>

#include "ito33/common.h"
#include "ito33/debug.h"
#include "ito33/exception.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/bondlike/bond.h"
#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/conversionschedule.h"

#include "ito33/finance/bondlike/bondliketerms.h"
#include "ito33/finance/bondlike/bondterms.h"

#include "ito33/list.h"
#include "ito33/vector.h"
#include "ito33/array.h"

#include "ito33/cppunit.h"

#include "ito33/arraycheckers.h"

#include "./testxmlreader_bondlike.h"
#include "./xmlreader_fortest.h"

#include "ito33/xml/read.h"

#include "ito33/xml/finance/bondlike/bondliketerms_all.h"
#include "ito33/xml/finance/bondlike/callschedule.h"
#include "ito33/xml/finance/bondlike/conversionschedule.h"

#include <xmlwrapp/init.h>
#include <xmlwrapp/document.h>
#include <xmlwrapp/tree_parser.h>

using namespace ito33;
using namespace ito33::XML;

using namespace ito33::finance;


void XMLReaderBondLikeTest::CheckBondLikeTerms()
{
  ito33::XML::TestingReader reader("xmlfiles/bondterm.xml");

  const xml::document& doc = reader.m_parser->get_document();
  const xml::node& node = doc.get_root_node();

  shared_ptr<BondLikeTerms> 
    pTerms = GetBondLikeTermsFromNode(*(node.find("bond_like_terms_test")));
  
  CPPUNIT_ASSERT_EQUAL( pTerms->GetIssueDate(), Date(2004, Date::Feb, 3) );

  CPPUNIT_ASSERT_DOUBLES_EQUAL( pTerms->GetIssuePrice(), 1.5, 1.e-10);

  CPPUNIT_ASSERT_EQUAL( pTerms->GetMaturityDate(), Date(2010, Date::Nov, 1) );
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( pTerms->GetRecoveryRate(), .2, 1.e-10);

  CPPUNIT_ASSERT( pTerms->GetCashDistribution() );
}


void XMLReaderBondLikeTest::CheckBondTerms()
{
  ito33::XML::TestingReader reader("xmlfiles/bondterm.xml");

  const xml::document& doc = reader.m_parser->get_document();
  const xml::node& node = doc.get_root_node();

  shared_ptr<BondTerms> pTerms = GetBondTermsFromNode
                                  (
                                    *(node.find("bond_terms_test"))
                                  );

  CPPUNIT_ASSERT_DOUBLES_EQUAL
    ( pTerms->GetNominal(), 58.6, 1.e-10);

  CPPUNIT_ASSERT_DOUBLES_EQUAL
    ( pTerms->GetYieldToMaturityOfAccretingBond(), .02, 1.e-10);

  CPPUNIT_ASSERT_EQUAL
    ( pTerms->GetMaturityDate(), Date(2010, Date::Nov, 1));
}

void XMLReaderBondLikeTest::CheckBondTermsWrongRedemption()
{
  ito33::XML::TestingReader reader("xmlfiles/bondterm.xml");

  const xml::document& doc = reader.m_parser->get_document();
  const xml::node& node = doc.get_root_node();

  shared_ptr<BondTerms>
    pTerms = GetBondTermsFromNode
                (
                  *(node.find("wrong_bond_terms_test_wrong_redemption"))
                );
}


void XMLReaderBondLikeTest::CheckCallSchedule()
{
  ito33::XML::TestingReader reader("xmlfiles/callschedule.xml");

  const xml::document& doc = reader.m_parser->get_document();
  const xml::node& node = doc.get_root_node();

  shared_ptr<CallSchedule> pCallSchedule;

  // we can restore a call schedule here
  CPPUNIT_ASSERT(
                  Restore
                    (
                      *(node.find("call_schedule_test")),
                      pCallSchedule
                    )
                );
  
  CPPUNIT_ASSERT( pCallSchedule->GetKeepAccrued() );

  CPPUNIT_ASSERT( !pCallSchedule->GetForfeitCoupon() );
  
  CPPUNIT_ASSERT_EQUAL( (int)pCallSchedule->GetNoticePeriod(), 20 );

  
  CPPUNIT_ASSERT_EQUAL( (int)pCallSchedule->GetTriggerPeriod(), 21 );
  
  CPPUNIT_ASSERT_EQUAL( (int)pCallSchedule->GetTriggerHistory(), 2 );
  
  CPPUNIT_ASSERT_EQUAL( pCallSchedule->GetMakeWholeType(), MakeWholeType_Coupon );
  
  CPPUNIT_ASSERT( !pCallSchedule->IsPVCouponMakeWhole() );

  shared_ptr<CallPeriod> pCallPeriod = pCallSchedule->GetAll().front();

  // check CallPeriod

  CPPUNIT_ASSERT_EQUAL( pCallPeriod->GetStartDate(), Date(2003, Date::Feb, 1));
  CPPUNIT_ASSERT_EQUAL( pCallPeriod->GetEndDate(), Date(2004, Date::May, 1));
  CPPUNIT_ASSERT_DOUBLES_EQUAL( pCallPeriod->GetStrike(), 1, 1.e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( pCallPeriod->GetTrigger(), 1.2, 1.e-10);
}

void XMLReaderBondLikeTest::CheckConversionSchedule()
{
  ito33::XML::TestingReader reader("xmlfiles/conversionschedule.xml");

  const xml::document& doc = reader.m_parser->get_document();
  const xml::node& node = doc.get_root_node();

  shared_ptr<ConversionSchedule> pConversionSchedule;

  // we can restore a call schedule here
  CPPUNIT_ASSERT(
                  Restore
                    (
                      *(node.find("conversion_schedule_test")),
                      pConversionSchedule
                    )
                );
  

  CPPUNIT_ASSERT( pConversionSchedule->GetKeepAccrued() );

  CPPUNIT_ASSERT( !pConversionSchedule->GetForfeitCoupon() );
  
  shared_ptr<ConversionPeriod> pConversionPeriod = pConversionSchedule->GetAll().front();

  // check ConversionPeriod

  CPPUNIT_ASSERT_EQUAL( pConversionPeriod->GetStartDate(), Date(2002, Date::Feb, 1));
  CPPUNIT_ASSERT_EQUAL( pConversionPeriod->GetEndDate(), Date(2004, Date::May, 1));
  CPPUNIT_ASSERT_DOUBLES_EQUAL( pConversionPeriod->GetRatio(), 1.5, 1.e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( pConversionPeriod->GetTrigger(), 1.2, 1.e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( pConversionPeriod->GetChangeRate(), 0.02, 1.e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( pConversionPeriod->GetExtremeTrigger(), 1.4, 1.e-10);
  CPPUNIT_ASSERT_EQUAL( pConversionPeriod->GetCoCoType(), CoCoType_CheckQuarterlyAndConvertAsOfCheckDate);
  
}
