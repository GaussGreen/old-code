/////////////////////////////////////////////////////////////////////////////
// Name:        hg/tests/basicclass/theoreticalmodel.h
// Purpose:     testing TheoreticalModel class
// Created:     2005/04/28
// RCS-ID:      $Id: theoreticalmodel.h,v 1.2 2006/07/06 16:04:28 wang Exp $
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

class TheoreticalModelTest : public CppUnit::TestCase
{
public:
  TheoreticalModelTest() { }

private:
  CPPUNIT_TEST_SUITE( TheoreticalModelTest );
    CPPUNIT_TEST( ZeroSharpeRatio );
    CPPUNIT_TEST( AnalyticCheck );
    CPPUNIT_TEST( CumulativeDefaultProba );
  CPPUNIT_TEST_SUITE_END();

  /**
     when the sharpe ratio is zero, the computed underlying process 
     should be the same
   */
  void ZeroSharpeRatio(); 
  
  /**
     Checks the solution against the original equation
   */
  void AnalyticCheck();

  /** 
     Checks the cumulative proba against analytical solution.
   */
  void CumulativeDefaultProba();
    
  NO_COPY_CLASS(TheoreticalModelTest);
};
