/////////////////////////////////////////////////////////////////////////////
// Name:        testsharedependentconversion.h
// Purpose:     acceptance tests for shared dependent conversion
// Author:      ITO 33
// Created:     14/03/2005
// RCS-ID:      $Id: testsharedependentconversion.h,v 1.2 2005/04/21 17:34:18 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_COMMON_TESTS_SHAREDEPENDENTCONVERSION_H_
#define _ITO33_COMMON_TESTS_SHAREDEPENDENTCONVERSION_H_

#include "ito33/cppunit.h"
#include "ito33/exception.h"

// ----------------------------------------------------------------------------
// Shared dependent conversion
// ----------------------------------------------------------------------------

class ShareDependentConversionTest : public CppUnit::TestCase
{
public:
  ShareDependentConversionTest() { }

private:
 CPPUNIT_TEST_SUITE( ShareDependentConversionTest );
 
 CPPUNIT_TEST_EXCEPTION( EndDateBeforeStartDate, ito33::Exception);
 CPPUNIT_TEST_EXCEPTION( BaseRatioNegative, ito33::Exception);
 CPPUNIT_TEST_EXCEPTION( ResetDateBeforeStartDate, ito33::Exception);
 CPPUNIT_TEST_EXCEPTION( ResetDateAfterEndDate, ito33::Exception);
 CPPUNIT_TEST_EXCEPTION( IncrementalShareFactorNegative, ito33::Exception);
 CPPUNIT_TEST_EXCEPTION( CapRatioNegative, ito33::Exception);
 CPPUNIT_TEST_EXCEPTION( CapRatioLessThanBaseRatio, ito33::Exception);
 CPPUNIT_TEST_EXCEPTION( FixedStrikeNegative, ito33::Exception);
 CPPUNIT_TEST_EXCEPTION( TriggerRateNegative, ito33::Exception);
 CPPUNIT_TEST_EXCEPTION( ExtremeTriggerRateNegative, ito33::Exception);
 CPPUNIT_TEST_EXCEPTION( ChangeRatePositive, ito33::Exception);
 CPPUNIT_TEST_EXCEPTION( ChangeRateNegative, ito33::Exception);
 CPPUNIT_TEST( Dump );

 CPPUNIT_TEST_SUITE_END();

 
 void EndDateBeforeStartDate();
 void BaseRatioNegative();
 void ResetDateBeforeStartDate();
 void ResetDateAfterEndDate();
 void IncrementalShareFactorNegative();
 void CapRatioNegative();
 void CapRatioLessThanBaseRatio();
 void FixedStrikeNegative();
 void TriggerRateNegative();
 void ExtremeTriggerRateNegative();
 void ChangeRatePositive();
 void ChangeRateNegative();
 void Dump();
 
 
NO_COPY_CLASS(ShareDependentConversionTest);
};


#endif //_ITO33_COMMON_TESTS_SHAREDEPENDENTCONVERSION_H_
