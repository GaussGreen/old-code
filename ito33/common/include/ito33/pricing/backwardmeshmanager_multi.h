/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/backwardmeshmanager_multi.h
// Purpose:     base multi-grid mesh manager for backward PDE problems
// Author:      Nabil
// Created:     2004/03/16
// RCS-ID:      $Id: backwardmeshmanager_multi.h,v 1.17 2006/01/03 17:18:48 zhang Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_PRICING_BACKWARDMESHMANAGERMULTI_H_
#define _ITO33_PRICING_BACKWARDMESHMANAGERMULTI_H_

#include "ito33/vector.h"
#include "ito33/array.h"

#include "ito33/pricing/backwardmeshmanager_withspacemesh.h"

#include "ito33/numeric/mesh/gridwithcv.h"


namespace ito33
{

namespace pricing
{

/// Backward time mesh manager
class BackwardMeshManagerMulti : public BackwardMeshManagerWithSpaceMesh
{

public:
  
  /**
     Ctor, call the ctor of the class BackwardMeshManager

     @param params The underlying contract paramaters
     @param model The model used for the pricing
   */  
  BackwardMeshManagerMulti(Params &params, Model &model) 
                         : BackwardMeshManagerWithSpaceMesh(params, model) { }

  /// dummy virtual dtor
  virtual ~BackwardMeshManagerMulti() { }

  /// Gets the m_nNbSpotsMax member
  size_t GetNbSMax() const { return m_nNbSMax; }

  /**
     @name Functions for grid transition.
      
     @attention{They must be used only in a transition grid time}
   */
  //@{

  /// Gets the m_nOldNbSpots member
  size_t GetOldNbS() const { return m_nOldNbS; }
  
  /// Gets a pointer to the old space mesh
  const double *GetOldS() const 
  {
    return m_pdOldS.Get();
  }

  /// Gets a pointer to the old space mesh
  const double *GetOldLogS() const 
  {
    return m_pdOldLogS.Get();
  }

  //@}

  /** 
     Gets the m_bChangeOfVar member of the current grid.
   
     @return a bool : the m_bChangeOfVar member of the current grid.
   */
  bool IsChangeOfVar() const { return m_pGrids[m_nIdxGrid].m_bChangeOfVar; }

  /** 
     Gets the m_bIsEndOfGrid member.
   
     @return a bool : The m_bIsEndOfGrid member.
   */
  bool IsEndOfGrid() const { return m_bIsEndOfGrid; }

  /// Set up the meshes and the data 
  virtual void SetupMe();

  /// set the initial state of the mesh manager
  virtual void SetInitialState() 
  { 
    BackwardMeshManager::SetInitialState();

    m_nIdxGrid = m_nNbGrids - 1;
    m_nIdxTimeGrid = m_pGrids[m_nIdxGrid].m_nNbTimes - 1;
    m_nNbS = m_pGrids[m_nIdxGrid].m_nNbS;
    UpdateSpaceMeshes();

    m_bIsEndOfGrid = false;
    m_bIsInitialTime = true;
  }


protected:
 
  /// Construct the space meshes
  virtual void ConstructSpaceMesh() = 0;

  /// Go ahead in time
  virtual void GoAhead();

  /// Update space meshes
  void UpdateSpaceMeshes();
 
  /// maximum size of all the space meshes
  size_t m_nNbSMax;

  /// the number of subgrids
  size_t m_nNbGrids;
  
  /// the subgrid
  std::vector<numeric::mesh::GridWithCV> m_pGrids;  

  /**
     @name state variable.
     
     Pay attention to update them correctly together with m_nNbS
   */
  //@{

  /// check if it's at end of a subgrid (can be different at t or t-)
  bool m_bIsEndOfGrid;

  /// check if the current time is initial time (need to think more on it)
  bool m_bIsInitialTime;

  /// index of subgrid for current time (can be different at t or t-)
  size_t m_nIdxGrid;

  /// local time index inside of subgrid (can be different at t or t-)
  size_t m_nIdxTimeGrid;

  /// size of the old space mesh, updated during grid transition.
  size_t m_nOldNbS;
  
  //@}  
  
  /// storage for old space mesh, shouldn't be changed once initialized.
  Array<double> m_pdOldS;

  /// storage for old log space mesh, shouldn't be changed once initialized.
  Array<double> m_pdOldLogS;


private:

  NO_COPY_CLASS(BackwardMeshManagerMulti);

}; // class BackwardMeshManagerMulti


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_BACKWARDMESHMANAGERMULTI_H_

