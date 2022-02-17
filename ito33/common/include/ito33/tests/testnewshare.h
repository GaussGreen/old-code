/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testnewshare.h
// Purpose:     header file for new share test
// Author:      Nabil
// Created:     2005/04/11
// RCS-ID:      $Id: testnewshare.h,v 1.5 2006/08/19 22:28:08 wang Exp $
// Copyright:   (c) 2005- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_TEST_NEWSHARE_H_
#define _ITO33_TEST_NEWSHARE_H_

#include "ito33/cppunit.h"
#include "ito33/common.h"
#include "ito33/exception.h"

#include "ito33/finance/sessiondata.h"

class NewShareTest : public CppUnit::TestCase 
{

public:

  NewShareTest() {}

  void tearDown() {}

private:

  CPPUNIT_TEST_SUITE( NewShareTest );
  
    CPPUNIT_TEST_EXCEPTION ( NewShareAndExchangeable, ito33::Exception);
    CPPUNIT_TEST_EXCEPTION ( ExchangeableAndNewShare, ito33::Exception);

  CPPUNIT_TEST_SUITE_END();

  void NewShareAndExchangeable();
  void ExchangeableAndNewShare();
  
  // common functions
  ito33::shared_ptr<ito33::finance::SessionData> InitSessionData();

  NO_COPY_CLASS( NewShareTest );

}; // Class NewShareTest

#endif // _ITO33_TEST_NEWSHARE_H_
