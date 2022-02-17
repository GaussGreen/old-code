/////////////////////////////////////////////////////////////////////////////
// Name:        instdatamulti.cpp
// Purpose:     instdatamulti class for ihg project
// Author:      Nabil
// Created:     2004/03/29
// RCS-ID:      $Id: instdatamulti.cpp,v 1.12 2004/10/04 18:04:07 pedro Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/numeric/interpolation.h"

#include "ito33/pricing/constraints.h"

#include "ihg/instdatamulti.h"

namespace ito33
{

using namespace numeric;

namespace ihg
{
  
void InstDataMulti::InterpWithPassageOfSpaceMesh(const double *pdOldSpots, 
                        size_t nNbOldSpots, const double *pdSpots, 
                        size_t nNbSpots)
{
  
  // Interpolation of price vector on the new grid. 
  Interpolate(pdOldSpots, m_pdOldPrices.Get(), nNbOldSpots, 
              pdSpots, m_pdPrices.Get(), nNbSpots, 
              ExtrapolationMode_Linear, ExtrapolationMode_Linear);
  
  // Interpolation of vega vector on the new grid.
  if( m_bComputeVega )
    Interpolate(pdOldSpots, m_pdOldVegas.Get(), nNbOldSpots, 
                pdSpots, m_pdVegas.Get(), nNbSpots, 
                ExtrapolationMode_Linear, ExtrapolationMode_Linear);
  
  // Interpolation of fugit vector on the new grid.
  if( m_bComputeFugit )
    Interpolate(pdOldSpots, m_pdOldFugits.Get(), nNbOldSpots, 
                pdSpots, m_pdFugits.Get(), nNbSpots, 
                ExtrapolationMode_Linear, ExtrapolationMode_Linear);
}

void InstDataMulti::ApplyConstraintsToAll()
{
  //Apply Constraints to the price.
  m_pConstraints->Apply(m_pdPrices.Get(), m_piFrozenFlags.Get(), m_nNbS);
          
  //Apply Constraints to the Vega.   
  if( m_bComputeVega )
    pricing::ApplyConstraintsToGreek(m_pdVegas.Get(), m_piFrozenFlags.Get(), 
      m_nNbS);
   
  //Apply Constraints to the fugit.   
  if( m_bComputeFugit )
    pricing::ApplyConstraintsToGreek(m_pdFugits.Get(), m_piFrozenFlags.Get(), 
      m_nNbS);
}

void InstDataMulti::UpdateAtEndOfGrid()
{ 
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
  
  /*
    Update the constraints at t-, it can be different from constraints at t+
    because of keep accrued or forfeit coupon flags and/or change on index 
    on continuous call or conversion.

    The constraints will not be used for system solving. the constraints used
    to solve the non linear system is done in InstData::UpdateBeforeStep.
  */
  UpdateConstraints();
  
  if( m_pConstraints )
    ApplyConstraintsToAll();
}


} // namespace ihg

} // namespace ito33
