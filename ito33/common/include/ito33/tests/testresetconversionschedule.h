/////////////////////////////////////////////////////////////////////////////
// Name:        testresetconversionschedule.h
// Purpose:     acceptance tests for resets at the financial level
// Author:      Yann d'Halluin
// Created:     17/11/2004
// RCS-ID:      $Id: testresetconversionschedule.h,v 1.1 2004/11/18 18:55:19 yann Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/exception.h"

// ----------------------------------------------------------------------------
// Conversion schedule reset tests
// ----------------------------------------------------------------------------

class ConversionScheduleResetTest : public CppUnit::TestCase
{
public:
  ConversionScheduleResetTest() { }


private:
 CPPUNIT_TEST_SUITE( ConversionScheduleResetTest );
 
 CPPUNIT_TEST_EXCEPTION( StartDateBeforeEndDate,ito33::Exception );
 CPPUNIT_TEST_EXCEPTION( InitialRatioNegative,ito33::Exception );
 CPPUNIT_TEST_EXCEPTION( CurrentRatioNegative,ito33::Exception );
 CPPUNIT_TEST_EXCEPTION( ResetDateOutsideSchedule,ito33::Exception );

 CPPUNIT_TEST( Dump );


 CPPUNIT_TEST_SUITE_END();

 void StartDateBeforeEndDate();
 void CurrentRatioNegative();
 void InitialRatioNegative();
 void ResetDateOutsideSchedule();
 void Dump();

 
NO_COPY_CLASS(ConversionScheduleResetTest);
};
