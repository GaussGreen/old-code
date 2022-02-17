/////////////////////////////////////////////////////////////////////////////
// Name:        hg/tests/varianceswap.h
// Purpose:     tests on variance swaps
// Created:     2006/03/090
// RCS-ID:      $Id: varianceswap.h,v 1.6 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

// ----------------------------------------------------------------------------
// Variance swap tests
// ----------------------------------------------------------------------------

namespace ito33
{
  namespace finance
  {
    class VarianceSwap;
  }

namespace hg
{

namespace VarianceSwapTest
{


void Setup(shared_ptr<TheoreticalModel> pModel,
           shared_ptr<finance::VarianceSwap> pVarianceSwap);


class VarianceSwapTest : public CppUnit::TestCase
{
public:
  
  VarianceSwapTest() { }

private:

  CPPUNIT_TEST_SUITE( VarianceSwapTest );
    CPPUNIT_TEST( BasicVarianceSwap );
    CPPUNIT_TEST( BasicVolatilitySwap );
    CPPUNIT_TEST( ChangeParams );
    CPPUNIT_TEST( SetCurrentValues );
    CPPUNIT_TEST( Compare2Dand3D );
    CPPUNIT_TEST( ClosedForm );
  CPPUNIT_TEST_SUITE_END();
  
  // Basic variance swap
  void BasicVarianceSwap();

  // Basic volatility swap
  void BasicVolatilitySwap();

  // Test effect of changing var swap constructor params
  void ChangeParams();

  // Test setting the volatility when valuation is after the sampling start
  void SetCurrentValues();

  // Compare full pricing (3D) and similarity reduction (2D)
  void Compare2Dand3D();

  // compare with closed form
  void ClosedForm();

  // Relative error comparison
  void RelErrorCheck(double dVal1, double dVal2, double dPrecision = 1.e-3);

   NO_COPY_CLASS(VarianceSwapTest);
};


} // namespace VarianceSwapTest

} // namespace hg

} // namespace ito33
