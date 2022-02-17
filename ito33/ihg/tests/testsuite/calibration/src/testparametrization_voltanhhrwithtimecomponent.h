/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testsuite/calibration/src/testparametrization_voltanhhrwithtimecomponent.h
// Purpose:     header file for testing parametrizationVolTanhHRWithTimeComponent
// Author:      ITO33
// Created:     2006/01/27
// RCS-ID:      $Id: testparametrization_voltanhhrwithtimecomponent.h,v 1.3 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"
#include "ito33/sharedptr.h"

namespace ito33
{
namespace finance
{
  class Derivatives;
  class SessionData;
}
namespace ihg
{
  class ParametrizationVolTanhHRWithTimeComponent;
}
}



class Parametrization_VolTanhHRWithTimeComponentTest : public CppUnit::TestCase
{
public:
  Parametrization_VolTanhHRWithTimeComponentTest() { }

  virtual void tearDown() { }

private:
  CPPUNIT_TEST_SUITE( Parametrization_VolTanhHRWithTimeComponentTest );
    CPPUNIT_TEST( TestAsSpecified );
    CPPUNIT_TEST( TestGeneratedValues );
  CPPUNIT_TEST_SUITE_END();

  // Use the parametrization as specified in the xml file
  void TestAsSpecified();

  // Generate some time only hazard rates. Should fit exactly
  void TestGeneratedValues();
  

  /// Setup the member variables by reading global xml file
  void Setup();

  /// The CDS termstructure to fit
  ito33::shared_ptr<ito33::finance::TermStructureCDS> m_ptsCDS;

  /// The first option to fit
  ito33::shared_ptr<ito33::finance::Option> m_pOption1;

  /// The 2nd option to fit
  ito33::shared_ptr<ito33::finance::Option> m_pOption2;

  /// The session used for calibration and pricing
  ito33::shared_ptr<ito33::finance::SessionData> m_pSessionData;

  /// The parametrization object
  ito33::shared_ptr<ito33::ihg::ParametrizationVolTanhHRWithTimeComponent> 
    m_pParametrization;

  NO_COPY_CLASS(Parametrization_VolTanhHRWithTimeComponentTest);
};

