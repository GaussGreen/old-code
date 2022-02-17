/////////////////////////////////////////////////////////////////////////////
// Name:        testresets.h
// Purpose:     acceptance tests for resets at the financial level
// Author:      Yann d'Halluin
// Created:     17/11/2004
// RCS-ID:      $Id: testresets.h,v 1.3 2006/05/20 16:17:50 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/exception.h"

// ----------------------------------------------------------------------------
// Conversion schedule reset tests
// ----------------------------------------------------------------------------

class ResetTest : public CppUnit::TestCase
{
public:
  ResetTest() { }


private:
 CPPUNIT_TEST_SUITE( ResetTest );
 
 CPPUNIT_TEST_EXCEPTION( ConversionScheduleExists,ito33::Exception );
 CPPUNIT_TEST_EXCEPTION( ContingentCallNotSupported,ito33::Exception );
 CPPUNIT_TEST_EXCEPTION( ConvSchStartDateAfterBondIssueDate,ito33::Exception );
 CPPUNIT_TEST_EXCEPTION( ConvSchEndDateBeforeBondMaturityDate,ito33::Exception );
 CPPUNIT_TEST_EXCEPTION( ConvSchContainsAtLeastOneResetDate,ito33::Exception );
 CPPUNIT_TEST_EXCEPTION( LastResetDateDifferentFromBondMaturity,ito33::Exception );

 CPPUNIT_TEST( Dump );


 CPPUNIT_TEST_SUITE_END();

 void ConversionScheduleExists();
 void ContingentCallNotSupported();
 void ConvSchStartDateAfterBondIssueDate();
 void ConvSchEndDateBeforeBondMaturityDate();
 void ConvSchContainsAtLeastOneResetDate();
 void LastResetDateDifferentFromBondMaturity();
 void Dump();

 
NO_COPY_CLASS(ResetTest);
};
