/////////////////////////////////////////////////////////////////////////////
// Name:        tests/nonlinearsolver/testNLS.cpp
// Purpose:     unit tests for non linear equation 1D solvers
// Author:      Pedro Ferreira, Laurence Gozalo
// Created:     04.12.03
// RCS-ID:      $Id: testNLS.h,v 1.4 2006/03/01 14:42:05 yann Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_TEST_NLS_H_
#define _ITO33_TEST_NLS_H_

#include "ito33/common.h"
#include "ito33/exception.h"
#include "ito33/cppunit.h"


// Test class
class NLSTest : public CppUnit::TestCase {
public:
  NLSTest( ) {}

  void tearDown() {}

private:
  CPPUNIT_TEST_SUITE( NLSTest );
    CPPUNIT_TEST ( Sin );
    CPPUNIT_TEST ( Polynom );
    CPPUNIT_TEST ( IterPol);
    CPPUNIT_TEST ( RFSin );
    CPPUNIT_TEST ( RFPolynom );
    CPPUNIT_TEST ( RFPolynom3 );
    CPPUNIT_TEST ( RFIterPol);
    CPPUNIT_TEST ( RFLinearPolynom ); 
    CPPUNIT_TEST ( RFExp );
    CPPUNIT_TEST_EXCEPTION ( Tolerance, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION ( Div0 , ito33::Exception );
    CPPUNIT_TEST_EXCEPTION ( MaxIterations , ito33::Exception );

    CPPUNIT_TEST( RFFindMatch );
    CPPUNIT_TEST_EXCEPTION( RFFindMatchBad, ito33::Exception );
  CPPUNIT_TEST_SUITE_END();

  void Sin();
  void Polynom();
  void IterPol();
  void RFSin();
  void RFPolynom();
  void RFPolynom3();
  void RFIterPol();
  void RFFindMatch();
  void RFFindMatchBad();
  void RFLinearPolynom();
  void RFExp();
  void Tolerance();
  void Div0();
  void MaxIterations();

  NO_COPY_CLASS( NLSTest );
};

#endif
