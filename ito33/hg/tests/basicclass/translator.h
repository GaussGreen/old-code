/////////////////////////////////////////////////////////////////////////////
// Name:        hg/tests/basicclass/translator.h
// Purpose:     testing Translator class
// Created:     2005/04/26
// RCS-ID:      $Id: translator.h,v 1.1 2005/04/27 13:06:35 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <math.h>

#include "ito33/cppunit.h"
#include "ito33/exception.h"

namespace ito33
{

namespace hg
{
  class Translator;
}

}

// ----------------------------------------------------------------------------
// TranslatorTest tests
// ----------------------------------------------------------------------------

class TranslatorTest : public CppUnit::TestCase
{
public:
  TranslatorTest() { }

private:
  CPPUNIT_TEST_SUITE( TranslatorTest );
    CPPUNIT_TEST( ExpecetedUnderlyingProcess );
    CPPUNIT_TEST( ExpectedParameters);
  CPPUNIT_TEST_SUITE_END();

  // Test the individual functions
  void ExpecetedUnderlyingProcess(); 
  
  /**
     Test translation of a given underlying process against parameters obtained
     manually
   */
  void ExpectedParameters();

  NO_COPY_CLASS(TranslatorTest);
};
