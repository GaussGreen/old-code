///////////////////////////////////////////////////////////////////////////////
// Name:        common/tests/binarysearch/testbinarysearch.cpp
// Purpose:     Unit testing of binary search
// Author:      yann d'Halluin
// Created:     22/08/2004
// RCS-ID:      $Id: testbinarysearch.cpp,v 1.4 2005/12/16 20:36:19 dave Exp $
// Copyright:   (c) 2004 Trilemma LLP
///////////////////////////////////////////////////////////////////////////////


#include "ito33/beforestd.h"
#include <cmath>
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/common.h"
#include "ito33/exception.h"
#include "ito33/cppunit.h"

#include "ito33/tests/testbinarysearch.h"

using namespace std;
using namespace ito33;

void BinarySearchTest::ElementTooLarge()
{
  size_t results = 0;
  
  size_t nSizeA = 4;
  double pdA[4]={0,1,2,3};

  double dx     = 5;

  results = BinSearch(pdA,nSizeA,dx);

  CPPUNIT_ASSERT( results == nSizeA-1 );

}

void BinarySearchTest::ElementTooSmall()
{

  size_t results = 0;
  double pdA[4]={10,20,30,300};
  size_t nSizeA = 4;
  double dx     = 1.5;

  results = BinSearch(pdA,nSizeA,dx);

  CPPUNIT_ASSERT( results == 1);

}


void BinarySearchTest::ExactLocation()
{

  size_t results = 0;
  double pdA[4]={10,20,30,300};
  size_t nSizeA = 4;
  double dx     = 20;

  results = BinSearch(pdA,nSizeA,dx);

  CPPUNIT_ASSERT( results == 2);

}

void BinarySearchTest::HalfWayLocation()
{
  size_t results = 0;
  double pdA[4]={10,20,30,300};
  size_t nSizeA = 4;
  double dx     = 10.1;

  results = BinSearch(pdA,nSizeA,dx);

  CPPUNIT_ASSERT( results == 1);

}
