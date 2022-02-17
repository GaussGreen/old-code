/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/tests/sensitivity.cpp
// Purpose:     Implementation of sensitivity tests
// Created:     2005/05/11
// RCS-ID:      $Id: sensitivity.cpp,v 1.15 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/beforestd.h"
#include <iostream>
#include <vector>
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/finance/modeloutput.h"

#include "ito33/hg/theoreticalmodel.h"

#include "hg/computesensitivity.h"
#include "hg/numoutput.h"

#include "hg/tests/sensitivity.h"

namespace ito33
{

namespace hg
{

namespace SensitivityTest
{


// Global object to work around CppUnit limitation
shared_ptr<finance::Derivative> m_pDerivative;
shared_ptr<TheoreticalModel> m_pModel;

void Setup(shared_ptr<finance::Derivative> pDerivative,
           shared_ptr<TheoreticalModel> pModel)
{
  m_pDerivative = pDerivative;
  m_pModel = pModel;
}

// Called by both SensitivityTest and PartialSensitivityTest
// pdSensFull is assumed to contain all sensitivities. pdSensComputed
// can be smaller (equal to count of trues in pbSensFlags)
void CompareSensitivities(const std::vector<double>& pdSensFull,
                          const std::vector<double>& pdSensComputed,
                          const std::vector<bool>& pbSensFlags)
{

  size_t nIdxComputed = 0;
  for (size_t nIdx = 0; nIdx < pbSensFlags.size(); nIdx++)
  {  
    // only compare if computed
    if (pbSensFlags[nIdx])
    {
      const double dScale = 1. / ( 1. + fabs(pdSensFull[nIdx]) );

      double dError = fabs(pdSensComputed[nIdxComputed] - pdSensFull[nIdx]) 
                    * dScale;

      // if the absolute value is small, loosen the tolerance
      double dTol = 5.e-4;
      if ( fabs(pdSensFull[nIdx]) < 1.0 )
        dTol = 1.e-3;

      CPPUNIT_ASSERT_DOUBLES_EQUAL( 0, dError, dTol );

      nIdxComputed++;
    } // if this sensitivity was computed
  } // loop over all sensitivities

}


void SensitivityTest::PDEAgainstFD()
{
  // Create new model, and force the same mesh
  TheoreticalModel model(*m_pModel);
  model.SetUnderlyingProcessForMesh( m_pModel->GetUnderlyingProcess() );

  // use PDEs to compute sensitivity
  shared_ptr<finance::ComputationalFlags> 
    pFlags(new finance::ComputationalFlags);
  pFlags->ActivateAllSensitivities(true);
  pFlags->SetSensitivityMethod(1);

  m_pDerivative->SetComputationalFlags(pFlags);

  shared_ptr<finance::ModelOutput> pMO = model.Compute(*m_pDerivative);

  shared_ptr<hg::NumOutput>
    pNumOutput( dynamic_pointer_cast<hg::NumOutput>( pMO->GetNumOutput() ) );

  if ( !pNumOutput->HasSensitivities() )
    return;

  std::vector<double> pdSensitivities( pNumOutput->GetSensitivities() );
 
  std::vector<double> 
    sensitivitiesFD( ComputeSensitivity(model, *m_pDerivative) );

  CPPUNIT_ASSERT( pdSensitivities.size() == sensitivitiesFD.size() );

  // Compare all sensitivity values
  std::vector<bool> pbSensFlags(pdSensitivities.size(), true);

  CompareSensitivities(sensitivitiesFD, pdSensitivities, pbSensFlags);

}


void SensitivityTest::AdjointAgainstFD()
{
  // Create new model, and force the same mesh
  TheoreticalModel model(*m_pModel);
  model.SetUnderlyingProcessForMesh( m_pModel->GetUnderlyingProcess() );

  // use adjoint to compute sensitivity
  shared_ptr<finance::ComputationalFlags> 
    pFlags(new finance::ComputationalFlags);
  
  pFlags->ActivateAllSensitivities(true);
  pFlags->SetSensitivityMethod(2);

  m_pDerivative->SetComputationalFlags(pFlags);

  shared_ptr<finance::ModelOutput> pMO = model.Compute(*m_pDerivative);

  shared_ptr<hg::NumOutput>
    pNumOutput( dynamic_pointer_cast<hg::NumOutput>( pMO->GetNumOutput() ) );

  if ( !pNumOutput->HasSensitivities() )
    return;

  std::vector<double> pdSensitivities( pNumOutput->GetSensitivities() );

  std::vector<double> 
    sensitivitiesFD( ComputeSensitivity(model, *m_pDerivative) );

  CPPUNIT_ASSERT( pdSensitivities.size() == sensitivitiesFD.size() );
  
  // Compare all sensitivity values
  std::vector<bool> pbSensFlags(pdSensitivities.size(), true);

  CompareSensitivities(sensitivitiesFD, pdSensitivities, pbSensFlags);
  
}


void PartialSensitivityTest::RunTest(const std::vector<bool>& pbSensFlags)                             
{

  // Create new model, and force the same mesh
  TheoreticalModel model(*m_pModel);
  shared_ptr<UnderlyingProcess> pProcess = m_pModel->GetUnderlyingProcess();
  model.SetUnderlyingProcessForMesh( pProcess );
  
  shared_ptr<finance::ComputationalFlags> 
    pFlags(new finance::ComputationalFlags);
  
  m_pDerivative->SetComputationalFlags(pFlags);

  // Compute all of the sensitivities by finite difference
  std::vector<double> 
    sensitivitiesFD( ComputeSensitivity(model, *m_pDerivative) );

  // Set the partial flags for this test, and count how many entries
  // are frozen so we can check computed array sizes
  pFlags->SetSensitivityFlags(pbSensFlags);
  size_t nNbFrozen = 0;
  for (size_t nIdx = 0; nIdx < pbSensFlags.size(); nIdx++)
  {
    if (pbSensFlags[nIdx] == false)
      nNbFrozen++;
  }
  
  // Test with both adjoint and PDE.  For contract types that do not support
  // one or both of these methods, the tests will be duplicated. 
  for (size_t nIdxMethod = 1; nIdxMethod <= 2; nIdxMethod++)
  {
    // Set method
    pFlags->SetSensitivityMethod(nIdxMethod);    
    m_pDerivative->SetComputationalFlags(pFlags);

    // Actually compute
    shared_ptr<finance::ModelOutput> pMO = model.Compute(*m_pDerivative);

    shared_ptr<NumOutput>
      pNumOutput( static_pointer_cast<NumOutput>(pMO->GetNumOutput()) );

    // Get the computed sensitivities and compare against FD
    std::vector<double> pdSensitivities( pNumOutput->GetSensitivities() );
 
    CPPUNIT_ASSERT( pdSensitivities.size() == 
                    sensitivitiesFD.size() - nNbFrozen );
  
    CompareSensitivities(sensitivitiesFD, pdSensitivities, pbSensFlags);
  
  } // loop over sensitivity methods

}

void PartialSensitivityTest::FreezeNone()
{

  // Compute all sensitivities via the flags
  shared_ptr<UnderlyingProcess> pProcess = m_pModel->GetUnderlyingProcess();
  std::vector<bool> pbSensFlags(pProcess->GetNbParameters(), true);
  RunTest(pbSensFlags);

}

void PartialSensitivityTest::FreezeFirst()
{

  // Compute all but the first sensitivity
  shared_ptr<UnderlyingProcess> pProcess = m_pModel->GetUnderlyingProcess();
  std::vector<bool> pbSensFlags(pProcess->GetNbParameters(), true);
  pbSensFlags[0] = false;
  RunTest(pbSensFlags);

}

void PartialSensitivityTest::FreezeLast()
{

  // Compute all but last sensitivity
  shared_ptr<UnderlyingProcess> pProcess = m_pModel->GetUnderlyingProcess();
  size_t nNbParams = pProcess->GetNbParameters();
  std::vector<bool> pbSensFlags(nNbParams, true);
  pbSensFlags[nNbParams - 1] = false;
  RunTest(pbSensFlags);

}


void PartialSensitivityTest::FreezeMultiple()
{

  // Freeze some pseudo-random sensitivities
  shared_ptr<UnderlyingProcess> pProcess = m_pModel->GetUnderlyingProcess();
  size_t nNbParams = pProcess->GetNbParameters();
  std::vector<bool> pbSensFlags(nNbParams, true);

  if (nNbParams > 8)
  {
    pbSensFlags[8] = false;
    pbSensFlags[7] = false;
    pbSensFlags[3] = false;
  }
  else if (nNbParams > 4)
  {
    pbSensFlags[1] = false;
    pbSensFlags[2] = false;
    pbSensFlags[4] = false;
  }
  else
  {
    pbSensFlags[0] = false;
  }

  RunTest(pbSensFlags);

}


void PartialSensitivityTest::FreezeAll()
{

  // Compute none of the sensitivities
  shared_ptr<UnderlyingProcess> pProcess = m_pModel->GetUnderlyingProcess();
  size_t nNbParams = pProcess->GetNbParameters();
  std::vector<bool> pbSensFlags(nNbParams, false);
  RunTest(pbSensFlags);

}


} // namespace SensitivityTest

} // namespace ito33

} // namespace hg
