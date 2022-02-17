/////////////////////////////////////////////////////////////////////////////
// Name:        test/dividends/main.cpp
// Purpose:     header file for unit test for Dividends class
// Author:      Vadim Zeitlin
// Created:     30.07.03
// RCS-ID:      $Id: testdividends.h,v 1.5 2006/01/03 17:14:14 zhang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_TEST_DIVIDENDS_H_
#define _ITO33_TEST_DIVIDENDS_H_

#include "ito33/common.h"
#include "ito33/exception.h"

#include "ito33/finance/dividends.h"

#include "ito33/cppunit.h"

// normally there should be no "using" in header but this is just a test...
using ito33::finance::Dividend;
using ito33::finance::Dividends;

class DivsTestCase : public CppUnit::TestCase
{
public:
  DivsTestCase() { }

  void tearDown() { m_divs.Clear(); }

private:

  CPPUNIT_TEST_SUITE( DivsTestCase );
    
    CPPUNIT_TEST_EXCEPTION( AddDuplicates, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( NegativeCash, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( LargeYield, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( LargePseudoYield, ito33::Exception );
    
    CPPUNIT_TEST( Add );
    CPPUNIT_TEST( AddReverse );
    
    CPPUNIT_TEST( SetCash );
    CPPUNIT_TEST( SetYield );
        
    CPPUNIT_TEST( Clear );
 
  CPPUNIT_TEST_SUITE_END();

  void AddDuplicates();
  void NegativeCash();
  void LargeYield();
  void LargePseudoYield();
  
  void Add();
  void AddReverse();
  
  void SetCash();
  void SetYield();

  void Clear();

  void DoCheckElement(const Dividend& div,
                      Dividend::Type type,
                      int date,
                      double value,
                      double pseudoYield = 0.);

  void DoCheckElements();

  Dividends m_divs;

  NO_COPY_CLASS(DivsTestCase);
};

#endif
