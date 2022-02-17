/////////////////////////////////////////////////////////////////////////////
// Name:        common/include/tests/testbinarysearch.h
// Purpose:     Unit testing of binary search
// Author:      yann d'Halluin
// Created:     22/08/2004
// RCS-ID:      $Id: testbinarysearch.h,v 1.3 2004/10/05 09:13:39 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/common.h"
#include "ito33/exception.h"

#include "ito33/binarysearch.h"
#include "ito33/cppunit.h"


class BinarySearchTest : public CppUnit::TestCase
{
public:
  BinarySearchTest() { }

private:
  CPPUNIT_TEST_SUITE( BinarySearchTest );
    CPPUNIT_TEST( ElementTooLarge );
    CPPUNIT_TEST( ElementTooSmall );
    CPPUNIT_TEST( ExactLocation );
    CPPUNIT_TEST( HalfWayLocation );
  CPPUNIT_TEST_SUITE_END();

  void ElementTooLarge();
  void ElementTooSmall();
  void ExactLocation();
  void HalfWayLocation();

  NO_COPY_CLASS(BinarySearchTest);
};
