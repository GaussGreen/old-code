/////////////////////////////////////////////////////////////////////////////
// Name:        cb/cboptioninstdata.cpp
// Purpose:     implementations for cb option instdata class in the ihg model
// Author:      Nabil
// Created:     2004/10/14
// RCS-ID:      $Id: cboptioninstdata.cpp,v 1.8 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004-2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#if ( defined(_MSC_VER) && _MSC_VER >= 1400 )
// still want to use fill_n but don't know how to avoid the compile warning
// warning C4996: 'std::_Fill_n' was declared deprecated
// so just disable the warning
#pragma warning(disable:4996)
#endif

#include "ito33/sharedptr.h"
#include "ito33/array.h"

#include "ito33/finance/computationalflags.h"

#include "ito33/numeric/interpolation.h"

#include "ito33/pricing/event.h"
#include "ito33/pricing/cbevent.h"
#include "ito33/pricing/cbconstraints.h"
#include "ito33/pricing/cbmeshmanager.h"
#include "ito33/pricing/cboptionparams.h"


#include "ihg/cboptioninstdata.h"
#include "ihg/cbinstdata.h"
#include "ihg/model.h"

namespace ito33
{

using namespace numeric;
using namespace pricing;

namespace ihg
{

CBOptionInstData::CBOptionInstData( pricing::CBOptionParams& cboptionparams, 
                                    Model& model, 
                                    pricing::CBMeshManager& cbmeshes)
                                  : InstDataMulti(cboptionparams, model, 
                                                  cbmeshes),
                                    m_cboptionparams(cboptionparams), 
                                    m_cbmeshes(cbmeshes),
                                    m_cbinstdata(CBInstData(cboptionparams, model, cbmeshes))
{ 
  // No recovery value for a convertible bond option in case of default
  m_dRecoveryValue = 0; 
}

void CBOptionInstData::Alloc(size_t nNbS)
{
  // Allocate the cb arrays
  m_cbinstdata.Alloc(nNbS);

  // Allocate the cb option arrays
  InstDataMulti::Alloc(nNbS);

  m_pdTmp1 = Array<double>(nNbS);
}

void CBOptionInstData::ApplyEvent(const pricing::Event *pEvent)
{
  // Apply event for the cb
  m_cbinstdata.ApplyEvent( pEvent );

  if( IsCBOptionWindow() && pEvent->GetType() == ET_Dividend )
  {
    pEvent->ApplyToPrice( m_pdS, m_pdPrices.Get(), m_nNbS );

    if( m_bComputeVega )
      pEvent->ApplyToGreek( m_pdS, m_pdVegas.Get(), m_nNbS );
  }
}

void CBOptionInstData::Swap()
{
  // we don't need to swap after the maturity of the ASW, 
  // as the price is just 0
  if( IsCBOptionWindow() )
    InstDataMulti::Swap();
}

void CBOptionInstData::InterpWithPassageOfSpaceMesh
                 (const double *pdOldSpots, size_t nNbOldSpots, 
                  const double *pdSpots, size_t nNbSpots)
{
  // Treatment for the cb
  m_cbinstdata.InterpWithPassageOfSpaceMesh
                ( pdOldSpots, nNbOldSpots, pdSpots, nNbSpots );

  // We have to update cb option prices only if it is needed. Then, only if we 
  // are in the life window of the cb option. Then, we have to test with the 
  // function IsCBOptionWindow(). 
  if( IsCBOptionWindow() )
    InstDataMulti::InterpWithPassageOfSpaceMesh
                  ( pdOldSpots, nNbOldSpots, pdSpots, nNbSpots );
}

void CBOptionInstData::Init()
{
  // Init for the cb
  m_cbinstdata.Init();

  m_nNbSpotsMax   = m_cbmeshes.GetNbSMax();
 
  Alloc(m_nNbSpotsMax);
  
  // Some initializations.
  m_pdS    = m_cbmeshes.GetS();
  m_pdLogS = m_cbmeshes.GetLogS();  

  bool bTmp;
  m_bDerivativeHasOwnHR = m_cboptionparams.GetCBLike().IsExchangeable(bTmp);
}

void CBOptionInstData::SetInitialValue()
{
  // Initial values for cb
  m_cbinstdata.SetInitialValue();

  m_nNbS = m_cbmeshes.GetNbS();

  std::fill_n(m_pdPrices.Get(), m_nNbSpotsMax, 0.);

  std::fill_n(m_piFrozenFlags.Get(), m_nNbSpotsMax, 0);

  if( m_bComputeVega )
    std::fill_n(m_pdVegas.Get(), m_nNbSpotsMax, 0.);
  
  if ( m_bComputeFugit )
    std::fill_n(m_pdFugits.Get(), m_nNbS, 0.);

  DoEvents();
}

double CBOptionInstData::GetSpeedCorrection() const
{
  return m_cbmeshes.GetSpeedCorrection();
}

void CBOptionInstData::UpdateBeforeStep()
{  
  if( m_cbmeshes.IsEndOfGrid() )
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

  InstDataMulti::UpdateBeforeStep();

  if(m_bDerivativeHasOwnHR)
    m_dHROfDerivative = m_model.GetHROfDerivative();
}

void CBOptionInstData::ApplyConstraintsAfterEventsOrConstraintsUpdate()
{
  // treatment for the underlying cb
  m_cbinstdata.ApplyConstraintsAfterEventsOrConstraintsUpdate();
  
  if( IsCBOptionWindow() )
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
  InterpWithPassageOfSpaceMesh( m_meshes.GetOldS(), m_meshes.GetOldNbS(), 
                                m_meshes.GetS(), m_meshes.GetNbS() );
  
  m_nNbS = m_meshes.GetNbS();

  if( IsCBOptionWindow() )
    UpdateCBOptionResults();
}

void CBOptionInstData::UpdateCBOptionConstraint()
{  
  ASSERT( IsCBOptionWindow() );

  double 
    dStrike = m_cboptionparams.GetCBOptionStrike();

  for(size_t nIdxS = 0; nIdxS < m_nNbS; ++nIdxS)
  {
    m_pdTmp1[nIdxS] = m_cbinstdata.m_pdPrices[nIdxS] - dStrike;

    if( m_pdTmp1[nIdxS] < 0. )
      m_pdTmp1[nIdxS] = 0.;
  }
  m_CBOptionConstraint.Update(m_pdTmp1.Get(), m_nNbS);
}

void CBOptionInstData::UpdateCBOptionResults()
{
  // Update cb option constraints.
  UpdateCBOptionConstraint();

  // Apply Constraints to the price.
  m_CBOptionConstraint.Apply(m_pdPrices.Get(), 
                              m_piFrozenFlags.Get(), 
                              m_nNbS);
          
  // Apply Constraints to the Vega.   
  if( m_bComputeVega )
    ApplyConstraintsToCBOptionGreek(m_pdVegas.Get(),
                                    m_cbinstdata.m_pdVegas.Get(),
                                    m_piFrozenFlags.Get(), 
                                    m_nNbS);
}

void CBOptionInstData::ApplyConstraintsToCBOptionGreek(double* pdGreeks,
                                                 const double* pdCBGreeks, 
                                                 const int* piFlagConstraints, 
                                                 size_t nNbValues)
{
  size_t
    nIdxSpot;

  for(nIdxSpot = 0; nIdxSpot < nNbValues; nIdxSpot++)
  {
    if( piFlagConstraints[nIdxSpot] )//cb_option_price[i] = cb_price[i] - K(t)
       pdGreeks[nIdxSpot] = pdCBGreeks[nIdxSpot];
  }
}

} // namespace ihg

} // namespace ito33
