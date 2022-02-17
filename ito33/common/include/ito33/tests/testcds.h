/////////////////////////////////////////////////////////////////////////////
// Name:        testcds.h
// Purpose:     test for cds
// Author:      Yann d'Halluin
// Created:     22/11/2004
// RCS-ID:      $Id: testcds.h,v 1.2 2005/07/18 10:11:44 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/exception.h"

class CDSTest : public CppUnit::TestCase
{

public:
  CDSTest(){}
  

 private:
  CPPUNIT_TEST_SUITE( CDSTest );
    CPPUNIT_TEST_EXCEPTION( NegativeRecoveryRate, ito33::Exception );  
    CPPUNIT_TEST_EXCEPTION( RecoveryRateTooLarge, ito33::Exception );   
    CPPUNIT_TEST_EXCEPTION( NoCashFlowStream, ito33::Exception );   
    CPPUNIT_TEST( Dump );
    
    CPPUNIT_TEST( SlidingCDS );

  CPPUNIT_TEST_SUITE_END();

  void NegativeRecoveryRate();
  void RecoveryRateTooLarge();
  void NoCashFlowStream();
  void Dump();

  void SlidingCDS();

  NO_COPY_CLASS(CDSTest);
};
