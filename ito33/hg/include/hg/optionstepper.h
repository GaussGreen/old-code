/////////////////////////////////////////////////////////////////////////////
// Name:        hg/optionstepper.h
// Purpose:     HG stepper class for regular options
// Created:     2005/01/13
// RCS-ID:      $Id: optionstepper.h,v 1.11 2006/08/04 16:09:13 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/optionstepper.h
   @brief HG stepper class for regular options
 */

#ifndef _HG_OPTIONSTEPPER_H_
#define _HG_OPTIONSTEPPER_H_

#include "ito33/autoptr.h"

#include "ito33/numeric/tridiagonalsolver.h"
#include "ito33/numeric/trisparseconstraintsolver.h"
#include "ito33/numeric/trisparselinearsolver.h"
#include "ito33/numeric/trisparselinearsolver_fixedpoint.h"
#include "ito33/numeric/trisparselinearsolver_gmres.h"
#include "ito33/numeric/trisparseconstraintsolver_penalty.h"
#include "ito33/numeric/trisparseconstraintsolver_frozen.h"

#include "hg/optioninstdata.h"
#include "hg/stepper_fix.h"
#include "hg/stepper_fe_fix.h"

namespace ito33
{

namespace hg
{


/// HG stepper class for regular option
template <class OptionStepperBase>
class OptionStepper : public OptionStepperBase
{

public:

  OptionStepper(OptionInstData& instdata,
                const finance::ComputationalFlags& flags) 
              : OptionStepperBase(instdata, flags), m_instdata(instdata)
  {
  }

  virtual ~OptionStepper() { }

  /**
      Initialization of the class data

      Most data is initialized by the base class Init() function, but the
      computational grid must still be defined (and delta, gamma).
   */
  void Init()
  {
    OptionStepperBase::Init();

    using namespace numeric;

    if (m_instdata.m_pConstraints == 0)
    {
      m_pLinearSolver = AutoPtr<TriSparseLinearSolver>
        //   ( new TriSparseLinearSolverFixedPoint(m_nNbX) );  
        ( new TriSparseLinearSolverGmres(m_nNbX) );
    }
    else
    {
      if ( m_iSolverType ) 
        m_pConstraintSolver = AutoPtr<TriSparseConstraintSolver>
          ( new TriSparseConstraintSolverFrozen(m_nNbRegimes, m_nNbX) );
      else
        m_pConstraintSolver = AutoPtr<numeric::TriSparseConstraintSolver>
          ( new TriSparseConstraintSolverPenalty(m_nNbRegimes, m_nNbX) ); 
    }
  }

  /**
      Makes a single timestep.
      
      Called by engine to make a single (possibly iterative) timestep.
      This function is pure virtual in the base classes.
   */
  void Run()
  {
    // Build the main system for price
    OptionStepperBase::BuildMainSystem();

    // Check if constraints are active, and choose the solver accordingly
    if (m_instdata.m_pConstraints == 0) 
      m_pLinearSolver->Solve(m_tridiagonalMatrix, 
                            m_sparseMatrix,
                            m_pdRHS.Get(),
                            m_instdata.m_pdPrices.Get());
    else
      m_pConstraintSolver->Solve(m_tridiagonalMatrix,
                                m_sparseMatrix,
                                m_pdRHS.Get(),
                                m_instdata.m_pConstraints,
                                m_instdata.m_pdPrices.Get(),
                                m_instdata.m_piFrozenFlags.Get());

    // Calculate sensitivities if required. Note that gamma is needed.
    // Pricing data comes from the solution computed above
    size_t nNbSensitivity = m_instdata.m_pSensitivityData.size();
    if ( nNbSensitivity 
      && m_instdata.m_sensitivityMethod == SensitivityMethod_PDE)
    {
      OptionStepperBase::ComputeHelperArrays();

      for (size_t nIdxD = 0; nIdxD < nNbSensitivity; nIdxD++)
      {
        OptionStepperBase::MakeSensitivityCoefficients
                           ( m_instdata.m_pSensitivityData[nIdxD] );

        OptionStepperBase::BuildRHS
            ( m_instdata.m_ppdOldSensitivities[nIdxD].Get(), 
              m_instdata.m_ppdOldOldSensitivities[nIdxD].Get() );

        if (m_instdata.m_pConstraints == 0)
          m_pLinearSolver->Solve(
                          m_tridiagonalMatrix,
                          m_sparseMatrix,
                          m_pdRHS.Get(),
                          m_instdata.m_ppdSensitivities[nIdxD].Get());
        else
          m_pConstraintSolver->SolveGreek(
                          m_tridiagonalMatrix,
                          m_sparseMatrix,
                          m_pdRHS.Get(),
                          m_instdata.m_piFrozenFlags.Get(),
                          m_instdata.m_ppdSensitivities[nIdxD].Get());

      } // loop over the model params
    } // if computing sensitivities

    if ( m_instdata.m_bDualSystemRequired )
      OptionStepperBase::SetupDualSystemData();

    if ( m_instdata.m_sensitivityMethod == SensitivityMethod_Adjoint )
    {
      OptionStepperBase::SetupSensivityByAdjointData();

      // constant used to modify the matrix according to the constraint flags
      const double dLarge = 1.e8;

      if ( m_instdata.m_pConstraints )
      {
        double* pdB = m_instdata.m_aData.m_pMatrix->GetB();

        for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
          if ( m_instdata.m_piFrozenFlags[nIdx] )
            pdB[nIdx] += dLarge;
      }

      if ( m_instdata.m_pConstraints )
      {
        Array<int> piConstraintFlags(m_nNbX);
        for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
          piConstraintFlags[nIdx] = m_instdata.m_piFrozenFlags[nIdx];
        
        m_instdata.m_aData.m_piConstraintFlags = piConstraintFlags;
      }
    }
  }

protected:
  
  AutoPtr<numeric::TriSparseConstraintSolver> m_pConstraintSolver;

  /// The type info of instdata is needed (e.g. for Greeks)
  OptionInstData& m_instdata;

  using OptionStepperBase::m_nNbRegimes;
  using OptionStepperBase::m_nNbX;
  using OptionStepperBase::m_pdRHS;
  using OptionStepperBase::m_tridiagonalMatrix;
  using OptionStepperBase::m_sparseMatrix;
  using OptionStepperBase::m_pLinearSolver;
  using OptionStepperBase::m_iSolverType;

private:

  NO_COPY_TEMPLATE_CLASS(OptionStepper, OptionStepperBase);

}; // class OptionStepper;


} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_OPTIONSTEPPER_H_
