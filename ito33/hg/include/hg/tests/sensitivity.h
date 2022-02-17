/////////////////////////////////////////////////////////////////////////////
// Name:        hg/tests/sensitivity.h
// Purpose:     header file for Sensitivity tests
// Created:     2005/05/11
// RCS-ID:      $Id: sensitivity.h,v 1.4 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

// ----------------------------------------------------------------------------
// Generic tests on sensitivities
// ----------------------------------------------------------------------------

namespace ito33
{

  namespace finance
  {
    class Derivative;
  }

namespace hg
{

  class TheoreticalModel;

namespace SensitivityTest
{

// Set global (for sensitivity tests) flags
void Setup(shared_ptr<finance::Derivative> pDerivative,
           shared_ptr<TheoreticalModel> pModel);

// compare pdSensFull to pdSensComputed for all true entries in pbSensFlags.
// pdSensFull is assumed to contain all sensitivities. 
// pdSensComputed can be smaller (equal to count of trues in pbSensFlags).
void CompareSensitivities(const std::vector<double>& pdSensFull,
                          const std::vector<double>& pdSensComputed,
                          const std::vector<bool>& pbSensFlags);


class SensitivityTest : public CppUnit::TestCase
{
public:
  
  SensitivityTest() { }

private:

  CPPUNIT_TEST_SUITE( SensitivityTest );
    CPPUNIT_TEST( PDEAgainstFD );
    CPPUNIT_TEST( AdjointAgainstFD );
  CPPUNIT_TEST_SUITE_END();
   
  // Compare PDE sensitivities versus FD sensitivities
  void PDEAgainstFD();

  // Compare adjoint sensitivities versus FD sensitivities
  void AdjointAgainstFD();

  NO_COPY_CLASS(SensitivityTest);
};


class PartialSensitivityTest : public CppUnit::TestCase
{
public:
  
  PartialSensitivityTest() { }

private:

  CPPUNIT_TEST_SUITE( PartialSensitivityTest );
    CPPUNIT_TEST( FreezeNone );
    CPPUNIT_TEST( FreezeFirst );
    CPPUNIT_TEST( FreezeLast );
    CPPUNIT_TEST( FreezeMultiple );
    CPPUNIT_TEST( FreezeAll );
  CPPUNIT_TEST_SUITE_END();
   
  // Compute all sensitivities by setting all flags to true
  void FreezeNone();

  // Compute all but the first sensitivity
  void FreezeFirst();

  // Compute all but the last sensitivity
  void FreezeLast();

  // Only compute some of the sensitivities (freeze 2 or more)
  void FreezeMultiple();

  // Freeze all of the sensitivities
  void FreezeAll();

  // Helper function to run the tests
  void RunTest(const std::vector<bool>& pbSensFlags);

  NO_COPY_CLASS(PartialSensitivityTest);
};


} // namespace SensitivityTest

} // namespace hg

} // namespace ito33
