/////////////////////////////////////////////////////////////////////////////
// Purpose:     header file of XMLWrite test program
// Author:      Vadim Zeitlin
// Created:     22.03.04
// RCS-ID:      $Id: testxmlwrite.h,v 1.3 2004/10/05 09:13:39 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_TEST_XMLWRITE_H_
#define _ITO33_TEST_XMLWRITE_H_

#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/common.h"
#include "ito33/cppunit.h"


class XMLWriteTestCase : public CppUnit::TestCase
{
public:
  XMLWriteTestCase() { }

  virtual void tearDown() { m_oss.str().clear(); }

private:
  CPPUNIT_TEST_SUITE( XMLWriteTestCase );
    CPPUNIT_TEST( Empty );
    CPPUNIT_TEST( Close );
    CPPUNIT_TEST( Attr );
    CPPUNIT_TEST( AttrQuote );
    CPPUNIT_TEST( Value );
    CPPUNIT_TEST( ValueQuote );
    CPPUNIT_TEST( Comment );
    CPPUNIT_TEST( SubElements );
    CPPUNIT_TEST( SiblingElements );
    CPPUNIT_TEST( RandomExample );
  CPPUNIT_TEST_SUITE_END();

  void Empty();
  void Close();
  void Attr();
  void AttrQuote();
  void Value();
  void ValueQuote();
  void Comment();
  void SubElements();
  void SiblingElements();
  void RandomExample();

  std::ostringstream m_oss;
  std::string m_expected;

  NO_COPY_CLASS(XMLWriteTestCase);
};

#endif 
