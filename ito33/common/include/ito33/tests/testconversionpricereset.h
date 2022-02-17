/////////////////////////////////////////////////////////////////////////////
// Name:        testconversionpricereset.h
// Purpose:     acceptance tests for conversion price resets 
// Author:      Yann d'Halluin
// Created:     17/11/2004
// RCS-ID:      $Id: testconversionpricereset.h,v 1.1 2004/11/18 18:55:19 yann Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/exception.h"

// ----------------------------------------------------------------------------
// Conersion price reset tests
// ----------------------------------------------------------------------------

class ConversionPriceResetTest : public CppUnit::TestCase
{
public:
  ConversionPriceResetTest() { }


private:
 CPPUNIT_TEST_SUITE( ConversionPriceResetTest );
 
 CPPUNIT_TEST_EXCEPTION( CapRateTooLarge,ito33::Exception );
 CPPUNIT_TEST_EXCEPTION( CapRateTooSmall,ito33::Exception );
 CPPUNIT_TEST_EXCEPTION( FloorRateTooLarge,ito33::Exception );
 CPPUNIT_TEST_EXCEPTION( FloorRateTooSmall,ito33::Exception );
 CPPUNIT_TEST_EXCEPTION( MultiplierTooSmall,ito33::Exception );
 CPPUNIT_TEST_EXCEPTION( MultiplierTooLarge,ito33::Exception );
 CPPUNIT_TEST( Dump );
 CPPUNIT_TEST( CapRateDefaultToOne );
 CPPUNIT_TEST( MultiplierDefaultToOne );


 CPPUNIT_TEST_SUITE_END();

 void CapRateTooLarge();
 void CapRateTooSmall();
 void FloorRateTooLarge();
 void FloorRateTooSmall();
 void MultiplierTooLarge();
 void MultiplierTooSmall();
 void Dump();
 void CapRateDefaultToOne();
 void MultiplierDefaultToOne();


 
NO_COPY_CLASS(ConversionPriceResetTest);
};
