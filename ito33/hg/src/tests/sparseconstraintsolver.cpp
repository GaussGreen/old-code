/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/tests/sparseconstraintsolver.cpp
// Purpose:     Implementation of sparse constraint solver tests
// Created:     2005/05/11
// RCS-ID:      $Id: sparseconstraintsolver.cpp,v 1.7 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include <iostream>

#include "ito33/finance/computationalflags.h"
#include "ito33/hg/theoreticalmodel.h"

#include "hg/optionstepper.h"

#include "hg/tests/sparseconstraintsolver.h"

namespace ito33
{

namespace hg
{

namespace SparseConstraintSolverTest
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


void SparseConstraintSolverTest::FrozenAgainstPenalty()
{
  TheoreticalModel model(*m_pModel);
  shared_ptr<finance::ComputationalFlags> 
    pFlags(new finance::ComputationalFlags);

  // Use frozen method
  pFlags->SetSolverType(1);

  m_pDerivative->SetComputationalFlags(pFlags);

  shared_ptr<finance::ModelOutput> pMO = model.Compute(*m_pDerivative);

  double dFrozenPrice = pMO->GetPrice();

  // Use penalty method
  pFlags->SetSolverType(0);

  m_pDerivative->SetComputationalFlags(pFlags);

  pMO = model.Compute(*m_pDerivative);

  double dPenaltyPrice = pMO->GetPrice();

  // Compare
  double dError = (dFrozenPrice - dPenaltyPrice) / dFrozenPrice;
  
  //std::cout << "frozen = " << dFrozenPrice
  //          << ", penalty = " << dPenaltyPrice
  //          << ", error = " << dError
  //          << std::endl;

  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0, dError, 1.e-4 );

}


} // namespace SparseConstraintcSolverTest

} // namespace ito33

} // namespace hg
