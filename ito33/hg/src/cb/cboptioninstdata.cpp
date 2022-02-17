/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/cb/cboptioninstdata.cpp
// Purpose:     implementations for cb option instdata class in HG
// Created:     2006/01/19
// RCS-ID:      $Id: cboptioninstdata.cpp,v 1.8 2006/05/28 15:36:52 zhang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#if ( defined(_MSC_VER) && _MSC_VER >= 1400 )
// still want to use fill_n but don't know how to avoid the compile warning
// warning C4996: 'std::_Fill_n' was declared deprecated
// so just disable the warning
#pragma warning(disable:4996)
#endif

#include "ito33/array.h"

#include "ito33/numeric/interpolation.h"

#include "ito33/pricing/event.h"
#include "ito33/pricing/cbevent.h"
#include "ito33/pricing/cbconstraints.h"
#include "ito33/pricing/cbmeshmanager.h"
#include "ito33/pricing/cboptionparams.h"

#include "ito33/finance/computationalflags.h"

#include "hg/cboptioninstdata.h"
#include "hg/cbinstdata.h"

namespace ito33
{

using namespace numeric;

namespace hg
{

CBOptionInstData::CBOptionInstData(pricing::CBOptionParams& cboptionparams, 
      Model& model, pricing::CBMeshManager& cbmeshes)
    : InstDataWithConstraints(cboptionparams, model, cbmeshes),
      m_cboptionparams(cboptionparams), 
      m_cbmeshes(cbmeshes),
      m_cbinstdata(cboptionparams, model, cbmeshes)
{
  // No recovery value for a convertible bond option in case of default
  m_dRecoveryValue = 0; 
}

void CBOptionInstData::Alloc(size_t nNbX)
{
  // Allocate the cb arrays
  m_cbinstdata.Alloc(nNbX);

  // Allocate the cb option arrays
  InstDataWithConstraints::Alloc(nNbX);

  m_pdConstriantsTmp = Array<double>(nNbX);
}

void CBOptionInstData::ApplyEvent(const pricing::Event* pEvent)
{
  // Apply event for the cb
  m_cbinstdata.ApplyEvent( pEvent );

  if ( IsCBOptionWindow() && pEvent->GetType() == pricing::ET_Dividend )
    for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
      pEvent->ApplyToPrice(m_pdS, m_pdPrices.Get() + nIdxR * m_nNbS, m_nNbS);
}

void CBOptionInstData::Swap()
{
  // we don't need to swap after the maturity of the ASW, 
  // as the price is just 0
  if ( IsCBOptionWindow() )
    InstDataWithConstraints::Swap();
}

void CBOptionInstData::SetupFlags(const finance::ComputationalFlags& flags)
{
  InstDataWithConstraints::SetupFlags(flags);

  m_sensitivityMethod = SensitivityMethod_None;

  // m_sensitivityMethod = SensitivityMethod_PDE;

  m_cbinstdata.SetupFlags(flags);
}

void CBOptionInstData::InterpWithPassageOfSpaceMesh
                 (const double *pdOldS, size_t nNbOldS, 
                  const double *pdS, size_t nNbS)
{
  // Treatment for the cb
  m_cbinstdata.InterpWithPassageOfSpaceMesh(pdOldS, nNbOldS, pdS, nNbS);

  // We have to update cb option prices only if it is needed. Then, only if we 
  // are in the life window of the cb option. Then, we have to test with the 
  // function IsCBOptionWindow(). 
  if ( IsCBOptionWindow() )
  {
    // Interpolation of price vector on the new grid. 
    for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
      Interpolate(pdOldS, m_pdOldPrices.Get() + nIdxR * nNbOldS,
                  nNbOldS, 
                  pdS, m_pdPrices.Get() + nIdxR * nNbS, nNbS, 
                  ExtrapolationMode_Linear, ExtrapolationMode_Linear);
  }
}

void CBOptionInstData::Init()
{
  // Init for the cb
  m_cbinstdata.Init();

  m_nNbSpotsMax = m_cbmeshes.GetNbSMax();
 
  Alloc(m_nNbSpotsMax * m_nNbRegimes);
  
  // Some initializations.
  m_pdS    = m_cbmeshes.GetS();
  m_pdLogS = m_cbmeshes.GetLogS();
}

void CBOptionInstData::SetInitialValue()
{
  // Initial values for cb
  m_cbinstdata.SetInitialValue();

  m_nNbS = m_cbmeshes.GetNbS();
  m_nNbX = m_nNbS * m_nNbRegimes;

  std::fill_n(m_pdPrices.Get(), m_nNbX, 0.);
  std::fill_n(m_pdOldPrices.Get(), m_nNbX, 0.);
  std::fill_n(m_pdOldOldPrices.Get(), m_nNbX, 0.);
  
  std::fill_n(m_piFrozenFlags.Get(), m_nNbX, 0);

  m_aData.m_bIsValid = false;

  DoEvents();
}

double CBOptionInstData::GetSpeedCorrection() const
{
  return m_cbmeshes.GetSpeedCorrection();
}

void CBOptionInstData::UpdateBeforeStep()
{  
  if ( m_cbmeshes.IsEndOfGrid() )
  {
    /*
      Get initial condition for the new subgrid from the followed subgrid.
      
      IMPORTANT: Note that the UpdateAtEndOfGrid() for the underlying CB
                 is called throughout the UpdateAtEndOfGrid() function
                 of the cb option. Then, all the updates at an end of grid
                 are made throughout this function.
      
      Time index is not changed, and is not the one at which the system will 
      be solved at this step.
    */
    UpdateAtEndOfGrid(); 

    /*
      Really go ahead, and change the time index, so constraints need to be 
      updated again.
    */
    m_cbmeshes.TryGoAhead();
  }

  /*
    Update for the cb.
    Note the the update at the end of a grid is made throughout the 
    function UpdateAtEndOfGrid() of the cb option.
  */
  m_cbinstdata.UpdateBeforeStep();

  InstDataWithConstraints::UpdateBeforeStep();
}

void CBOptionInstData::ApplyConstraintsAfterEventsOrConstraintsUpdate()
{
  // treatment for the underlying cb
  m_cbinstdata.ApplyConstraintsAfterEventsOrConstraintsUpdate();
  
  if ( IsCBOptionWindow() )
    UpdateCBOptionResults();
}

void CBOptionInstData::UpdateAtEndOfGrid()
{ 
  // Update at end of grid for the underlying cb
  m_cbinstdata.UpdateAtEndOfGrid();

  Swap();

  /*
    Maturity for the moment is considered as an end of grid, although 
    it's not necessary to do the interpolation, it's rather an exception.
    "Real" end of grid can have different space mesh at t+ and t-

    The interpolated values(after having applied the possible constraints)
    will be used as the initial condition of the previous subgrid 
    (we go backward on time).
  */
  InterpWithPassageOfSpaceMesh( m_cbmeshes.GetOldS(), m_cbmeshes.GetOldNbS(), 
                                m_cbmeshes.GetS(), m_cbmeshes.GetNbS() );
  
  m_nNbS = m_cbmeshes.GetNbS();

  if ( IsCBOptionWindow() )
    UpdateCBOptionResults();
}

void CBOptionInstData::UpdateCBOptionConstraint()
{  
  ASSERT( IsCBOptionWindow() );

  double dStrike = m_cboptionparams.GetCBOptionStrike();

  for (size_t nIdx = 0; nIdx < m_nNbX; ++nIdx)
  {
    m_pdConstriantsTmp[nIdx] = m_cbinstdata.m_pdPrices[nIdx] - dStrike;

    if ( m_pdConstriantsTmp[nIdx] < 0. )
      m_pdConstriantsTmp[nIdx] = 0.;
  }

  m_CBOptionConstraint.Update(m_pdConstriantsTmp.Get(), m_nNbX);
}

void CBOptionInstData::UpdateCBOptionResults()
{
  // Update cb option constraints.
  UpdateCBOptionConstraint();

  // Apply Constraints to the price.
  m_CBOptionConstraint.Apply(m_pdPrices.Get(), m_piFrozenFlags.Get(), m_nNbX);
}

} // namespace hg

} // namespace ito33
