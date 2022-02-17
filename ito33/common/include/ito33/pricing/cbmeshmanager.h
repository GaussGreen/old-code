/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/cbmeshmanager.h
// Purpose:     cb multi-grid mesh manager for backward PDE problems
// Author:      Nabil
// Created:     2004/03/17
// RCS-ID:      $Id: cbmeshmanager.h,v 1.34 2006/03/23 09:27:14 yann Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/cbmeshmanager.h
    @brief cb multi-grid mesh manager class

    Implementation of the meshes (space and time) manager for cb.
    
    The class contains not only the meshes but also all data related to
    this meshes.
 */

#ifndef _ITO33_PRICING_CBMESHMANAGER_H_
#define _ITO33_PRICING_CBMESHMANAGER_H_

#include "ito33/debug.h"

#include "ito33/numeric/mesh/roots.h"

#include "ito33/pricing/backwardmeshmanager_multi.h"
#include "ito33/pricing/cblikeparams.h"

namespace ito33
{

namespace pricing
{

/// Declaration of the CBMeshManager class
class CBMeshManager: public BackwardMeshManagerMulti
{
public:

  CBMeshManager(CBLikeParams &params, Model &model)
              : BackwardMeshManagerMulti(params, model), m_cbParams(params) 
  {
  }

  /// Dummy virtual dctor
  virtual ~CBMeshManager() { };

  // SetupMe already implemented in the base class BackwardMeshManagerMulti
   
  /**
     Modify the rates so that constant can be an exact solution for 
     the discretized PDE.
   */
  virtual void SetupRates();

  /** 
     Gets the correction for the convection term due to the change of
     variable.

     @return The value of this correction.
   */
  double GetSpeedCorrection()
  {
    if( IsChangeOfVar() )    
      return ( m_pGrids[m_nIdxGrid].m_pdCorrectionS[m_nIdxTimeGrid + 1] /
        m_pGrids[m_nIdxGrid].m_pdCorrectionS[m_nIdxTimeGrid] - 1. ) /
        GetTimeStep();
    else
      return 0.;
  }  

  /// Update the foreign exchange rate
  void UpdateFXRate()
  {
    m_cbParams.SetFXRate( m_pdFXRates[m_nIdx] );
  }


protected:
  
  /// Construct the space mesh 
  virtual void ConstructSpaceMesh();

  /// Construct adaptative grid
  void AdaptativeGrid();

  /**
     Computes fixed boundary for each time point if any and determines
     definitifly the begin time and end time for all grids.

     called by AdaptativeGrid().
   */
  void ComputeRootsAndUpdateGrids();

  /**
     Gets special points of this grid. 

     called by AdaptativeGrid().

     @param grid working grid
     @param pdEventLogSpot (output) array of special points
     @param nNbEventSpot (output) number of special points
   */
  void GetSpecialSpacePointsForGrid
       (
         numeric::mesh::GridWithCV& grid,
         double* pdEventLogSpot,
         size_t& nNbEventSpot
       );

  /**
     builds space meshes for given grid

     called by AdaptativeGrid().

     @param grid given grid
     @param dDelta given standard mesh size
     @param pdEventLogSpot array of special points to be put in the mesh
     @param nNbEventSpot number of special points
     @param dLogSpotInitial log(InitialSpot)
   */
  void MakeSpaceMeshForGrid
       (
         numeric::mesh::GridWithCV& grid, double dDelta, 
         double* pdEventLogSpot, size_t nNbEventSpot,
         double dLogSpotInitial
       );

  /**
     Calculates some useful mesh data

     called by AdaptativeGrid()

     @param dLogSpotMin left boundary of main mesh domain
     @param dLogSpotMax right boundary of main mesh domain
     @param dLogSpotInitial log(InitialSpot)
     @param dDelta allowed largest mesh size in main mesh domain
   */
  void ComputeMeshHelperData
       ( 
         double& dLogSpotMin,double& dLogSpotMax,double& dLogSpotInitial,
         double& dDelta
       );

  /// Compute Roots.
  void ComputeRoots(double dTime, size_t &nNbRoots, numeric::mesh::Root *pRoots,
                    bool bPlus = true)
  {
    m_cbParams.GetConversions()->ComputeRoots(dTime, nNbRoots, pRoots, bPlus);
  }

  /// the foreign exchange rates on the time mesh.
  Array<double> m_pdFXRates;
  
  /// the params for cb
  CBLikeParams& m_cbParams;


private:

  NO_COPY_CLASS(CBMeshManager);

}; // class CBMeshManager


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_CBMESHMANAGER_H_

