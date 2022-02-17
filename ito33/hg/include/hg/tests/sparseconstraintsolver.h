/////////////////////////////////////////////////////////////////////////////
// Name:        hg/tests/sparseconstraintsolver.h
// Purpose:     header file for sparse constraint solvers
// Created:     2005/06/15
// RCS-ID:      $Id: sparseconstraintsolver.h,v 1.2 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

// ----------------------------------------------------------------------------
// Generic tests on sparse consraint solvers (frozen, penalty)
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

namespace SparseConstraintSolverTest
{


void Setup(shared_ptr<finance::Derivative> pDerivative,
           shared_ptr<TheoreticalModel> pModel);


class SparseConstraintSolverTest : public CppUnit::TestCase
{
public:
  
  SparseConstraintSolverTest() { }

private:

  CPPUNIT_TEST_SUITE( SparseConstraintSolverTest );
    CPPUNIT_TEST( FrozenAgainstPenalty );
  CPPUNIT_TEST_SUITE_END();
   
  // Compare frozen method against the penalty method
  void FrozenAgainstPenalty();
       

  NO_COPY_CLASS(SparseConstraintSolverTest);
};


} // namespace SparseConstraintSolverTest

} // namespace hg

} // namespace ito33
