/////////////////////////////////////////////////////////////////////////////
// Name:        hg/tests/hero.h
// Purpose:     header file for HERO tests
// Created:     2005/10/01
// RCS-ID:      $Id: hero.h,v 1.4 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

// ----------------------------------------------------------------------------
// Tests for HERO (hedging is done separately)
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

namespace HeroTest
{


class HeroTest : public CppUnit::TestCase
{
public:
  
  HeroTest() { }


private:

  CPPUNIT_TEST_SUITE( HeroTest );
    CPPUNIT_TEST( ZeroHeroTest );
    CPPUNIT_TEST( AddHedgingContractsTest );    
    CPPUNIT_TEST( CompareSVDandNAGTest );
    CPPUNIT_TEST( GenericTests );
    //CPPUNIT_TEST_EXCEPTION( MissingDataTest, ito33::Exception );
  CPPUNIT_TEST_SUITE_END();
   
  // One regime, no jumps, should give hedge = delta as in Black-Scholes
  void ZeroHeroTest();
      
  // Add hedging contracts.  Hero should decrease
  void AddHedgingContractsTest();

  // Add hedging contracts.  Hero should decrease
  void CompareSVDandNAGTest();

  // Just run the hedging code on a list of xml files
  void GenericTests();

  // Hedging contract has insufficient surface data
  // Add this test when hg supports path-dependent contracts
  //void MissingDataTest();


  // _________ Helper functions and member data ___________
  void Setup(const std::string& sFileName);

  shared_ptr<hg::TheoreticalModel> m_pModel;
  shared_ptr<finance::Derivative> m_pTarget;
  shared_ptr<finance::Derivatives> m_pHedgeInstruments;

  NO_COPY_CLASS(HeroTest);
};


} // namespace HeroTest

} // namespace hg

} // namespace ito33
