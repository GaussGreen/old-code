/////////////////////////////////////////////////////////////////////////////
// Name:        hg/tests/basicclass/theoreticalmodel.cpp
// Purpose:     testing TheoreticalModel
// Created:     2005/04/28
// RCS-ID:      $Id: theoreticalmodel.cpp,v 1.5 2006/08/19 23:47:17 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @todo Add more tests for regime switching proba:symmetrical jumps should
          give identical probabilities, jump to default and no regime jumps
          should be same as no jump to default and one regime jump, etc
 */

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/numeric/densematrix.h"

#include "ito33/hg/theoreticalmodel.h"

#include "hg/translator.h"

#include "theoreticalmodel.h"

using namespace ito33;
using namespace ito33::hg;

// ----------------------------------------------------------------------------
// TranslatorTest tests
// ----------------------------------------------------------------------------

void TheoreticalModelTest::ZeroSharpeRatio()
{
  // Create an underlying process
  size_t nNbRegimes = 2;

  std::vector<double> pdVols(nNbRegimes);
  std::vector<double> pdDefaultIntensities(nNbRegimes);

  // Define a underlying process
  std::vector<double> pdIntensity(1);
  std::vector<double> pdAmplitude(1);

  // Define the  underlying process
  pdVols[0] = 0.01;
  pdVols[1] = 0.02;
  
  pdDefaultIntensities[0] = 0.03;
  pdDefaultIntensities[1] = 0.04;

  shared_ptr<UnderlyingProcess>
    pUnderlyingProcess( new UnderlyingProcess
                            (nNbRegimes, pdVols, pdDefaultIntensities) );

  pdIntensity[0] = 0.05;
  pdAmplitude[0] = 0.06;
  pUnderlyingProcess->SetJumps(0, 0, pdIntensity, pdAmplitude);

  pdIntensity[0] = 0.07;
  pdAmplitude[0] = 0.08;
  pUnderlyingProcess->SetJumps(0, 1, pdIntensity, pdAmplitude);

  pdIntensity[0] = 0.09;
  pdAmplitude[0] = 0.10;
  pUnderlyingProcess->SetJumps(1, 0, pdIntensity, pdAmplitude);

  pdIntensity[0] = 0.11;
  pdAmplitude[0] = 0.12;

  pUnderlyingProcess->SetJumps(1, 1, pdIntensity, pdAmplitude);

  TheoreticalModel model(pUnderlyingProcess);

  model.SetSharpeRatio(0);

  shared_ptr<UnderlyingProcess> pRealUP = model.GetRealUnderlyingProcess();

  // A hack, but should be ok, just compare the translated parameters
  Translator trans(*pUnderlyingProcess);

  // Get the parameters corresponding to the underlying process
  std::vector<double> parameters(trans.GetNbParameters());
  trans.GetParameters(&parameters[0]);
 
  Translator trans1(*pRealUP);

  // Get the parameters corresponding to the underlying process
  std::vector<double> parameters1(trans1.GetNbParameters());
  trans1.GetParameters(&parameters1[0]);

  CPPUNIT_ASSERT( parameters.size() == parameters1.size() );

  for (size_t nIdx = 0; nIdx < parameters.size(); nIdx++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(parameters[nIdx], parameters1[nIdx], 1.e-12);
}


void TheoreticalModelTest::AnalyticCheck()
{
  // Create an underlying process
  size_t nNbRegimes = 2;

  std::vector<double> pdVols(nNbRegimes);
  std::vector<double> pdDefaultIntensities(nNbRegimes);

  // Define a underlying process
  std::vector<double> pdIntensity(1);
  std::vector<double> pdAmplitude(1);

  // Define the  underlying process
  pdVols[0] = 0.01;
  pdVols[1] = 0.02;
  
  pdDefaultIntensities[0] = 0.03;
  pdDefaultIntensities[1] = 0.04;

  shared_ptr<UnderlyingProcess>
    pUnderlyingProcess( new UnderlyingProcess
                            (nNbRegimes, pdVols, pdDefaultIntensities) );

  pdIntensity[0] = 0.05;
  pdAmplitude[0] = 0.06;
  pUnderlyingProcess->SetJumps(0, 0, pdIntensity, pdAmplitude);

  pdIntensity[0] = 0.07;
  pdAmplitude[0] = 0.08;
  pUnderlyingProcess->SetJumps(0, 1, pdIntensity, pdAmplitude);

  pdIntensity[0] = 0.09;
  pdAmplitude[0] = 0.10;
  pUnderlyingProcess->SetJumps(1, 0, pdIntensity, pdAmplitude);

  pdIntensity[0] = 0.11;
  pdAmplitude[0] = 0.12;

  pUnderlyingProcess->SetJumps(1, 1, pdIntensity, pdAmplitude);

  TheoreticalModel model(pUnderlyingProcess);

  double dSR = 0.5;

  model.SetSharpeRatio(dSR);

  shared_ptr<UnderlyingProcess> pRealUP = model.GetRealUnderlyingProcess();

  // Check against the original equations
  
  for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
  {
    // compute first the total vol of the real underlying process
    double dVol = pRealUP->GetVolatilities()[nIdxR1];
    double dTotalVol = dVol * dVol;
    for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
    {
      const Jumps& jumps = pRealUP->GetJumps(nIdxR1, nIdxR2);
      Jumps::const_iterator pJump;
      for (pJump = jumps.begin(); pJump != jumps.end(); ++pJump)
        dTotalVol += pJump->GetIntensity() * pJump->GetAmplitude()
                   * pJump->GetAmplitude();
    }

    dTotalVol += pRealUP->GetJumpsToDefault()[nIdxR1];

    dTotalVol = sqrt(dTotalVol);

    // check each jumps
    for (size_t nIdxR2 = 0; nIdxR2 < nNbRegimes; nIdxR2++)
    {      
      const Jumps& jumps1 = pRealUP->GetJumps(nIdxR1, nIdxR2);
      Jumps::const_iterator pJump1;
      const Jumps& jumps2 = pUnderlyingProcess->GetJumps(nIdxR1, nIdxR2);
      Jumps::const_iterator pJump2;

      CPPUNIT_ASSERT( jumps1.size() == jumps2.size() );

      for (pJump1 = jumps1.begin(), pJump2 = jumps2.begin(); 
           pJump1 != jumps1.end() && pJump2 != jumps2.end();
           ++pJump1, ++pJump2)
      {
        double dIntensity = pJump1->GetIntensity();
        dIntensity *= (1. - dSR * pJump1->GetAmplitude() / dTotalVol);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(dIntensity, pJump2->GetIntensity(), 1.e-6);
      }
    }

    // check the default intensity
    double dDefaultIntensity = pRealUP->GetJumpsToDefault()[nIdxR1];
    dDefaultIntensity *= (1. - dSR * (- 1.) / dTotalVol);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(dDefaultIntensity, pdDefaultIntensities[nIdxR1], 1.e-6);
  }
}

void TheoreticalModelTest::CumulativeDefaultProba()
{
  // Create an underlying process
  size_t nNbRegimes = 2;

  std::vector<double> pdVols(nNbRegimes, 0.2);
  std::vector<double> pdDefaultIntensities(nNbRegimes);
  for (size_t nIdxR = 0; nIdxR < nNbRegimes; nIdxR++)
    pdDefaultIntensities[nIdxR] = 0.02 * (nIdxR + 1);

  UnderlyingProcess 
    underlyingProcess(nNbRegimes, pdVols, pdDefaultIntensities);

  double dT = 1;
  numeric::DenseMatrix 
    matrix(underlyingProcess.ComputeRegimeTransitionProba(dT));

  for (size_t nIdxR = 0; nIdxR < nNbRegimes; nIdxR++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1. - exp(- dT * pdDefaultIntensities[nIdxR]),
                                 matrix[nIdxR][nNbRegimes],
                                 1.e-4);
}
