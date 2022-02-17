/////////////////////////////////////////////////////////////////////////////
// Name:        hg/tests/onetouch.h
// Purpose:     tests on homogeneity applied to one touch
// Created:     2006/02/20
// RCS-ID:      $Id: onetouch.h,v 1.3 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

// ----------------------------------------------------------------------------
// Homogeneity tests with HG Model
// ----------------------------------------------------------------------------

namespace ito33
{
  namespace finance
  {
    class OneTouch;
  }

namespace hg
{

namespace OneTouchTest
{


void Setup(shared_ptr<TheoreticalModel> pModel,
           shared_ptr<finance::OneTouch> pOneTouch);

class OneTouchTest : public CppUnit::TestCase
{
public:
  
  OneTouchTest() { }

private:

  CPPUNIT_TEST_SUITE( OneTouchTest );
    CPPUNIT_TEST( Homogeneity );
    CPPUNIT_TEST( HomogeneityImplementation );
  CPPUNIT_TEST_SUITE_END();
  
  // Test the individual functions
  void Homogeneity();

  void HomogeneityImplementation();
  
  bool HomogeneityCanApply();

  NO_COPY_CLASS(OneTouchTest);
};


} // namespace OneTouchTest

} // namespace hg

} // namespace ito33
