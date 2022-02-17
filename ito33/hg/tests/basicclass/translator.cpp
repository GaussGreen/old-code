/////////////////////////////////////////////////////////////////////////////
// Name:        hg/tests/basicclass/translator.cpp
// Purpose:     testing Translator
// Created:     2005/04/27
// RCS-ID:      $Id: translator.cpp,v 1.3 2006/08/19 23:47:17 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "hg/translator.h"

#include "translator.h"

using namespace ito33;
using namespace ito33::hg;

// ----------------------------------------------------------------------------
// TranslatorTest tests
// ----------------------------------------------------------------------------

void TranslatorTest::ExpecetedUnderlyingProcess()
{
  // Create an underlying process
  size_t nNbRegimes = 2;

  std::vector<double> pdVols(nNbRegimes);
  std::vector<double> pdIntensities(nNbRegimes);

  UnderlyingProcess underlyingProcess(nNbRegimes, pdVols, pdIntensities);

  // Just define the structure of the underlying process.
  std::vector<double> pdIntensity(1);
  std::vector<double> pdAmplitude(1);

  underlyingProcess.SetJumps(0, 0, pdIntensity, pdAmplitude);

  underlyingProcess.SetJumps(0, 1, pdIntensity, pdAmplitude);

  underlyingProcess.SetJumps(1, 0, pdIntensity, pdAmplitude);

  underlyingProcess.SetJumps(1, 1, pdIntensity, pdAmplitude);

  Translator trans(underlyingProcess);

  // Get the parameters corresponding to the underlying process
  std::vector<double> parameters(trans.GetNbParameters());

  for (size_t nIdx = 0; nIdx < parameters.size(); nIdx++)
    parameters[nIdx] = 0.01 * (nIdx + 1);

  // Define the expected underlying process
  pdVols[0] = 0.01;
  pdVols[1] = 0.02;
  
  pdIntensities[0] = 0.03;
  pdIntensities[1] = 0.04;

  UnderlyingProcess expected(nNbRegimes, pdVols, pdIntensities);

  pdIntensity[0] = 0.05;
  pdAmplitude[0] = 0.06;
  expected.SetJumps(0, 0, pdIntensity, pdAmplitude);

  pdIntensity[0] = 0.07;
  pdAmplitude[0] = 0.08;
  expected.SetJumps(0, 1, pdIntensity, pdAmplitude);

  pdIntensity[0] = 0.09;
  pdAmplitude[0] = 0.10;
  expected.SetJumps(1, 0, pdIntensity, pdAmplitude);

  pdIntensity[0] = 0.11;
  pdAmplitude[0] = 0.12;
  
  expected.SetJumps(1, 1, pdIntensity, pdAmplitude);
 
  shared_ptr<UnderlyingProcess> pUnderlyingProcess = trans(parameters);

  // compare the two underlying process
  CPPUNIT_ASSERT(pUnderlyingProcess->GetNbRegimes() == nNbRegimes);

  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
  {
    CPPUNIT_ASSERT(    pUnderlyingProcess->GetVolatilities() 
                    == expected.GetVolatilities() );

    CPPUNIT_ASSERT(    pUnderlyingProcess->GetJumpsToDefault() 
                    == expected.GetJumpsToDefault());

    for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
    {
      const Jumps& jumps1 = expected.GetJumps(nIdxR1, nIdxR2);
      const Jumps& jumps2 = pUnderlyingProcess->GetJumps(nIdxR1, nIdxR2);

      CPPUNIT_ASSERT(jumps1.size() == jumps2.size());

      Jumps::const_iterator pJump1, pJump2;
      for (pJump1 = jumps1.begin(), pJump2 = jumps2.begin();
           pJump1 != jumps1.end() && pJump2 != jumps2.end();
           ++pJump1, ++pJump2)
      {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(pJump1->GetIntensity(), pJump2->GetIntensity(), 1.e-8);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(pJump1->GetAmplitude(), pJump2->GetAmplitude(), 1.e-8);
      }
    }
  }
}

void TranslatorTest::ExpectedParameters()
{
  // Create an underlying process
  size_t nNbRegimes = 2;

  std::vector<double> pdVols(nNbRegimes);
  std::vector<double> pdIntensities(nNbRegimes);

  pdVols[0] = 0.01;
  pdVols[1] = 0.02;
  
  pdIntensities[0] = 0.03;
  pdIntensities[1] = 0.04;

  UnderlyingProcess underlyingProcess(nNbRegimes, pdVols, pdIntensities);

  std::vector<double> pdIntensity(1);
  std::vector<double> pdAmplitude(1);

  pdIntensity[0] = 0.05;
  pdAmplitude[0] = 0.06;
  underlyingProcess.SetJumps(0, 0, pdIntensity, pdAmplitude);

  pdIntensity[0] = 0.07;
  pdAmplitude[0] = 0.08;
  underlyingProcess.SetJumps(0, 1, pdIntensity, pdAmplitude);

  pdIntensity[0] = 0.09;
  pdAmplitude[0] = 0.10;
  underlyingProcess.SetJumps(1, 0, pdIntensity, pdAmplitude);

  pdIntensity[0] = 0.11;
  pdAmplitude[0] = 0.12;
  
  underlyingProcess.SetJumps(1, 1, pdIntensity, pdAmplitude);

  Translator trans(underlyingProcess);

  // Get the parameters corresponding to the underlying process
  std::vector<double> parameters(trans.GetNbParameters());

  trans.GetParameters(&parameters[0]);

  for (size_t nIdx = 0; nIdx < parameters.size(); nIdx++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(parameters[nIdx], 0.01 * (nIdx + 1), 1.e-8);
}
