/////////////////////////////////////////////////////////////////////////////
// Name:        testattachedwarrantconvertiblebond.h
// Purpose:     acceptance tests for attached warrant convertible bond
// Author:      ITO 33
// Created:     15/03/2005
// RCS-ID:      $Id: testattachedwarrantconvertiblebond.h,v 1.2 2005/04/21 17:34:18 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_COMMON_TESTS_ATTACHEDWARRANTCONVERTIBLEBOND_H_
#define _ITO33_COMMON_TESTS_ATTACHEDWARRANTCONVERTIBLEBOND_H_

#include "ito33/cppunit.h"
#include "ito33/exception.h"

// ----------------------------------------------------------------------------
// Shared dependent conversion
// ----------------------------------------------------------------------------

class AttachedWarrantConvertibleBondTest : public CppUnit::TestCase
{
public:
  AttachedWarrantConvertibleBondTest() { }

private:
 CPPUNIT_TEST_SUITE( AttachedWarrantConvertibleBondTest );
 
 CPPUNIT_TEST_EXCEPTION( ConversionStartDateBeforeBondStartDate, 
                         ito33::Exception);
 CPPUNIT_TEST_EXCEPTION( ConversionEndDateAfterBondEndDate, 
                         ito33::Exception);
 CPPUNIT_TEST_EXCEPTION( EmptySharedDependentConversion, ito33::Exception);
 CPPUNIT_TEST( Dump );

 CPPUNIT_TEST_SUITE_END();

 void ConversionStartDateBeforeBondStartDate();
 void ConversionEndDateAfterBondEndDate();
 void EmptySharedDependentConversion();
 void Dump();
 
 
NO_COPY_CLASS(AttachedWarrantConvertibleBondTest);
};


#endif //_ITO33_COMMON_TESTS_ATTACHEDWARRANTCONVERTIBLEBOND_H_
