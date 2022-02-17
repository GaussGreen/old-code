/////////////////////////////////////////////////////////////////////////////
// Name:        hg/tests/hedging.h
// Purpose:     header file for hedging tests
// Created:     2005/10/01
// RCS-ID:      $Id: hedging.h,v 1.4 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

// ----------------------------------------------------------------------------
// Tests for hedging (hero is done separately)
// ----------------------------------------------------------------------------

namespace ito33
{

  namespace finance
  {
    class Derivative;
    class Derivatives;
  }

namespace hg
{

  class TheoreticalModel;

namespace HedgingTest
{


class HedgingTest : public CppUnit::TestCase
{
public:
  
  HedgingTest() { }


private:

  CPPUNIT_TEST_SUITE( HedgingTest );
    CPPUNIT_TEST( BSDeltaTest );
    CPPUNIT_TEST( HedgeWithTargetTest );
    CPPUNIT_TEST( GenericTests );
  CPPUNIT_TEST_SUITE_END();
   
  // One regime, no jumps, should give hedge = delta as in Black-Scholes
  void BSDeltaTest();

  // If the target contract is one of the hedge contracts, ratio should be 1
  void HedgeWithTargetTest();
       
  // Just run the hedging code on a list of xml files
  void GenericTests();


  // Read the specified xml file
  void Setup(const std::string& sFileName);

  // All of these are read from the xml file
  shared_ptr<hg::TheoreticalModel> m_pModel;
  shared_ptr<finance::Derivative> m_pTarget;
  shared_ptr<finance::Derivatives> m_pHedgeInstruments;

  NO_COPY_CLASS(HedgingTest);
};


} // namespace HedgingTest

} // namespace hg

} // namespace ito33
