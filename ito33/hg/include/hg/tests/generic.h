/////////////////////////////////////////////////////////////////////////////
// Name:        hg/tests/generic.h
// Purpose:     header file for generic tests on HG model
// Created:     2005/04/19
// RCS-ID:      $Id: generic.h,v 1.6 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

// ----------------------------------------------------------------------------
// Generic tests with HG Model
// ----------------------------------------------------------------------------

namespace ito33
{
  
namespace finance
{
  class Derivative; 
}

namespace hg
{

namespace GenericTest
{


void Setup(shared_ptr<finance::Derivative> pDerivative);

class GenericTest : public CppUnit::TestCase
{
public:
  
  GenericTest() { }

private:

  CPPUNIT_TEST_SUITE( GenericTest );
    CPPUNIT_TEST( IHG );
    CPPUNIT_TEST( NoInterRegimeJump );
    CPPUNIT_TEST( IdenticalRegimes );
    CPPUNIT_TEST( Symmetry );
    CPPUNIT_TEST( TransitionProbability );

  CPPUNIT_TEST_SUITE_END();
   
  // Test the individual functions
  void IHG();
  void NoInterRegimeJump();
  void IdenticalRegimes();
  void Symmetry();

  void TransitionProbability();

  NO_COPY_CLASS(GenericTest);
};


} // namespace GenericTest

} // namespace hg

} // namespace ito33
