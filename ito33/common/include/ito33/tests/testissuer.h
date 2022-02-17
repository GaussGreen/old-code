/////////////////////////////////////////////////////////////////////////////
// Name:        testissuer.h
// Purpose:     input validation for issuer
// Author:      Yann d'Halluin
// Created:     17/11/2004
// RCS-ID:      $Id: testissuer.h,v 1.2 2006/04/04 16:29:47 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/exception.h"

class IssuerTest : public CppUnit::TestCase
{

public:
  IssuerTest(){}
  

 private:
  CPPUNIT_TEST_SUITE( IssuerTest );
  CPPUNIT_TEST_SUITE_END();

  void NoMoneyMarket();

  NO_COPY_CLASS(IssuerTest);
};
