/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/tests/generic.cpp
// Purpose:     Implementation of generic tests on HG model
// Created:     2005/04/19
// RCS-ID:      $Id: generic.cpp,v 1.17 2006/08/19 23:44:27 wang Exp $
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

#include "ito33/date.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/predicatedouble.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardrateflat.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/transitionprobabilityoutput.h"

#include "hg/backwardnumoutput.h"

#include "hg/tests/generic.h"

namespace ito33
{

namespace hg
{

namespace GenericTest
{


// Global object to work around CppUnit limitation
shared_ptr<finance::Derivative> m_pDerivative;
void Setup(shared_ptr<finance::Derivative> pDerivative)
{
  m_pDerivative = pDerivative;
}

const double dVol = 0.5;
const double dIntensity = 0.2;

shared_ptr<ihg::TheoreticalModel> 
CreateIHGModelFromHGModel(const shared_ptr<TheoreticalModel>& pHG)
{
  ihg::TheoreticalModel* pIHG = new ihg::TheoreticalModel;
  
  shared_ptr<UnderlyingProcess> pUP = pHG->GetUnderlyingProcess();

  pIHG->SetVolatility( shared_ptr<ihg::VolatilityFlat>
                       ( new ihg::VolatilityFlat(pUP->GetVolatilities()[0]) ) );
  
  pIHG->SetHazardRate( shared_ptr<ihg::HazardRateFlat>
                       ( new ihg::HazardRateFlat(pUP->GetJumpsToDefault()[0]) ) );

  return shared_ptr<ihg::TheoreticalModel>(pIHG);
}

shared_ptr<TheoreticalModel> CreateIHGCompatibleModel()
{
  size_t nNbRegimes = 1;

  std::vector<double> pdVols(nNbRegimes, dVol);

  std::vector<double> pdIntensities(nNbRegimes, dIntensity);
  
  shared_ptr<UnderlyingProcess> 
    pUP( new UnderlyingProcess(nNbRegimes, pdVols, pdIntensities) );

  return shared_ptr<TheoreticalModel>( new TheoreticalModel(pUP) );
}

shared_ptr<TheoreticalModel> CreateIdenticalRegimes()
{
  size_t nNbRegimes = 2;

  std::vector<double> pdVols(nNbRegimes, dVol);

  std::vector<double> pdIntensities(nNbRegimes, dIntensity);
  
  shared_ptr<UnderlyingProcess> 
    pUP( new UnderlyingProcess(nNbRegimes, pdVols, pdIntensities) );

  return shared_ptr<TheoreticalModel>( new TheoreticalModel(pUP) );
}

shared_ptr<TheoreticalModel> CreateIndependentRegimes()
{
  size_t nNbRegimes = 2;

  double dVol1 = 0.2;
  std::vector<double> pdVols(nNbRegimes);
  pdVols[0] = dVol;
  pdVols[1] = dVol1;

  double dIntensity1 = 0.15;
  std::vector<double> pdIntensities(nNbRegimes);
  pdIntensities[0] = dIntensity;
  pdIntensities[1] = dIntensity1;

  shared_ptr<UnderlyingProcess> 
    pUP( new UnderlyingProcess(nNbRegimes, pdVols, pdIntensities) );

  return shared_ptr<TheoreticalModel>( new TheoreticalModel(pUP) );
}

shared_ptr<TheoreticalModel> CreateSymmetricalRegimes()
{
  shared_ptr<TheoreticalModel> pModel( CreateIdenticalRegimes() );

  shared_ptr<UnderlyingProcess> pUP = pModel->GetUnderlyingProcess();

  Jumps jumps;
  jumps.push_back(Jump(0.2, -0.2));

  pUP->SetJumps(0, 1, jumps);
  pUP->SetJumps(1, 0, jumps);

  return pModel;
}

void CompareOutputs(const finance::ModelOutput& mo1, 
                    const finance::ModelOutput& mo2)
{
  // For the moment, just check with price. 
  double dTmp = fabs( mo1.GetPrice() );
  double dDenom = ( dTmp > 1.0 ) ? dTmp : 1.0;
  double dError = fabs( mo1.GetPrice() - mo2.GetPrice() ) / dDenom;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, dError, 1.e-3);

  //CPPUNIT_ASSERT_DOUBLES_EQUAL(mo1.GetPrice(), mo2.GetPrice(), 5 * 1.e-3);

  // there is data at analysis date. Compare them
  if ( mo1.HasSpotAtAnalysisDate() )
  {
    std::vector<double> spots1 = mo1.GetSpotsAtAnalysisDate();
    std::vector<double> prices1 = mo1.GetPricesAtAnalysisDate();

    std::vector<double> spots2 = mo2.GetSpotsAtAnalysisDate();
    std::vector<double> prices2 = mo2.GetPricesAtAnalysisDate();

    std::vector<double> pdPricesTmp1;
    pdPricesTmp1.resize(spots1.size());

    numeric::Interpolate( &spots2[0], &prices2[0], spots2.size(),
                          &spots1[0], &pdPricesTmp1[0], spots1.size() );

    double dTol = 1.e-2;
    for (size_t nIdx = 0; nIdx < spots1.size(); nIdx++)
    {
      double dTmp = fabs(prices1[nIdx]);
      double dDenom = ( dTmp > 1.0 ) ? dTmp : 1.0;
      double dError = fabs( prices1[nIdx] - pdPricesTmp1[nIdx] )/ dDenom;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, dError, dTol);

      // accuracy decreases when far from the spot
      if (nIdx == 3 * spots1.size() / 4)
        dTol = 2.e-2;
    }

    std::vector<double> pdPricesTmp2;
    pdPricesTmp2.resize(spots2.size());

    numeric::Interpolate( &spots1[0], &prices1[0], spots1.size(),
                          &spots2[0], &pdPricesTmp2[0], spots2.size() );

    dTol = 1.e-2;
    for (size_t nIdx = 0; nIdx < spots2.size(); nIdx++)
    {
      double dTmp = fabs(prices2[nIdx]);
      double dDenom = ( dTmp > 1.0 ) ? dTmp : 1.0;
      double dError = fabs( prices2[nIdx] - pdPricesTmp2[nIdx] )/ dDenom;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, dError, dTol);

      // accuracy decreases when far from the spot
      if (nIdx == 3 * spots1.size() / 4)
        dTol = 2.e-2;
    }
  }

}


void CompareRegimesInOutput(const finance::ModelOutput& /* mo1 */)
{
  // Todo. Add implementation
}

void GenericTest::IHG()
{
  // Create a one regime model
  shared_ptr<TheoreticalModel> pModel1( CreateIHGCompatibleModel() );

  shared_ptr<finance::ComputationalFlags> 
    pFlagsHG(new finance::ComputationalFlags);

  pFlagsHG->SetAnalysisDate
      (m_pDerivative->GetSessionData()->GetValuationDate());

  m_pDerivative->SetComputationalFlags(pFlagsHG);

  // Compute with HG
  shared_ptr<finance::ModelOutput> pMO1( pModel1->Compute(*m_pDerivative) );

  // Create a IHG model
  shared_ptr<ihg::TheoreticalModel> 
    pModel2( CreateIHGModelFromHGModel(pModel1) );

  shared_ptr<finance::ComputationalFlags>
    pFlagsIHG(new finance::ComputationalFlags);

  pFlagsIHG->SetAnalysisDate
      (m_pDerivative->GetSessionData()->GetValuationDate());

  m_pDerivative->SetComputationalFlags(pFlagsIHG);

  // Compute with IHG
  shared_ptr<finance::ModelOutput> pMO2( pModel2->Compute(*m_pDerivative) );

  CompareOutputs(*pMO1, *pMO2);
}

void GenericTest::NoInterRegimeJump()
{
  // Create a one regime model
  shared_ptr<TheoreticalModel> pModel1( CreateIndependentRegimes() );

  // Compute with HG
  shared_ptr<finance::ModelOutput> pMO1( pModel1->Compute(*m_pDerivative) );
 
  // Create one regime model from the independent regimes model
  const shared_ptr<UnderlyingProcess>& pUP1 = pModel1->GetUnderlyingProcess();

  size_t nNbRegimes = 1;
  size_t nIdxR = 0;
  std::vector<double> pdVols2(1, pUP1->GetVolatilities()[nIdxR]);
  std::vector<double> pdIntensities2(1, pUP1->GetJumpsToDefault()[nIdxR]);
 
  shared_ptr<UnderlyingProcess> 
    pUP2( new UnderlyingProcess(nNbRegimes, pdVols2, pdIntensities2) );

  pUP2->SetJumps( nIdxR, nIdxR, pUP1->GetJumps(nIdxR, nIdxR) );
  
  shared_ptr<TheoreticalModel> pModel2( new TheoreticalModel(pUP2) );

  shared_ptr<finance::ModelOutput> pMO2( pModel2->Compute(*m_pDerivative) );

  CompareOutputs(*pMO1, *pMO2);
}

void GenericTest::IdenticalRegimes()
{
  // Create a one regime model
  shared_ptr<TheoreticalModel> pModel1( CreateIdenticalRegimes() );

  // Compute with HG
  shared_ptr<finance::ModelOutput> pMO1( pModel1->Compute(*m_pDerivative) );

  CompareRegimesInOutput(*pMO1);
}

void GenericTest::Symmetry()
{
  // Create a one regime model
  shared_ptr<TheoreticalModel> pModel1( CreateSymmetricalRegimes() );

  // Compute with HG
  shared_ptr<finance::ModelOutput> pMO1( pModel1->Compute(*m_pDerivative) );

  CompareRegimesInOutput(*pMO1);
}

static void 
CheckSum(const shared_ptr<TransitionProbabilityOutput>& pOutput, bool bStraight)
{
  std::vector< shared_ptr<finance::ModelOutput> > 
    probas( pOutput->GetProbabilities() );

  const std::vector<double>& returns(pOutput->GetSpaceMesh());
  size_t nNbS = returns.size();

  size_t nNbRegimes = probas.size();
  for (size_t nIdxR = 0; nIdxR < nNbRegimes; nIdxR++)
  {
    shared_ptr<finance::ModelOutput> probaR(probas[nIdxR]);
    shared_ptr<BackwardNumOutput> 
      pNumOutput( static_pointer_cast<BackwardNumOutput>
                  ( probaR->GetNumOutput() ) );
    
    double dSum1 = 0.;
    double dSum2 = 0.;
    for (size_t nIdxR1 = 0; nIdxR1 < nNbRegimes; nIdxR1++)
    {
      const double* pdTmp = &pNumOutput->GetFinalPrices()[0] + nNbS * nIdxR1;

      for (size_t nIdxS = 0; nIdxS < returns.size(); nIdxS++)
      {
        // probability can't be negative
        CPPUNIT_ASSERT(numeric::IsEqualOrGreater(pdTmp[nIdxS], 0));

        if ( bStraight )
          dSum1 += pdTmp[nIdxS];

        dSum2 += pdTmp[nIdxS] * returns[nIdxS];
      } 
    }
    
    if ( bStraight )
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1, dSum1, 1.e-6);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1, dSum2, 1.e-6);
  }
}

void GenericTest::TransitionProbability()
{
  shared_ptr<TheoreticalModel> pModel;
  shared_ptr<TransitionProbabilityOutput> pOutput;

  pModel = CreateSymmetricalRegimes();

  pOutput = pModel->ComputeTransitionProbability(1);

  CheckSum(pOutput, false);
  
  size_t nNbRegimes = 1;

  std::vector<double> pdVols;
  pdVols.resize(nNbRegimes, 0.2);

  std::vector<double> pdDefaultIntensities;
  pdDefaultIntensities.resize(nNbRegimes, 0.);

  shared_ptr<hg::UnderlyingProcess>
    pUnderlyingProcess( new hg::UnderlyingProcess
                            (nNbRegimes, pdVols, pdDefaultIntensities) );

  Jumps jumps;
    
  jumps.push_back(Jump(0.1, 0.2));
  pUnderlyingProcess->SetJumps(0, 0, jumps);

  pModel = make_ptr( new TheoreticalModel(pUnderlyingProcess) );
  
  double dDeltaT = 1./ 252;

  pOutput = pModel->ComputeTransitionProbability(dDeltaT);

  CheckSum(pOutput, true);

  pModel = CreateIHGCompatibleModel();

  std::vector<double> pdIntensities(1, 0);

  pModel->GetUnderlyingProcess()->SetJumpsToDefault(pdIntensities);

  pOutput = pModel->ComputeTransitionProbability(1);

  CheckSum(pOutput, true);
}


} // namespace GenericTest

} // namespace ito33

} // namespace hg
