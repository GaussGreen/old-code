/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/instdatamulti.h
// Purpose:     instdatamulti class for ihg project
// Author:      Nabil
// Created:     2004/03/26
// RCS-ID:      $Id: instdatamulti.h,v 1.11 2005/12/30 11:51:14 nabil Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/instdatamulti.h
    @brief ihg instdatamulti class
 */

#ifndef _IHG_INSTDATAMULTI_H_
#define _IHG_INSTDATAMULTI_H_

#include "ito33/pricing/backwardmeshmanager_multi.h"

#include "ihg/instdatawithconstraints.h"

namespace ito33
{

namespace ihg
{

/**
   Base instadata class for cb-like instruments
 */
class InstDataMulti : public InstDataWithConstraints
{

public:

  InstDataMulti(pricing::Params &params, 
                Model &model, 
                pricing::BackwardMeshManagerMulti &meshes)
    : InstDataWithConstraints(params, model, meshes), m_meshes(meshes)
  {  }

  virtual ~InstDataMulti() { }

  /// @name Functions required by the Engine
  //@{
   
  // Init() : No implementation at this level, so no need to redeclare

  // SetInitialValue(): No implementation at this level, so no need to redeclare

  //@}

  /**
    Gets the maximum mesh size in multi-grid case.
    
    */
  size_t GetNbSpotsMax() const { return m_nNbSpotsMax; }
    
  /**
    This function updates the data at the end of the grid 
    (except for the last one)

    @attention{This function updates in particular the member m_nNbS because
              the size of the space mesh can change only if the grid change:
              i.e For a given grid, the size of the space mesh is unchanged}
    
    */
  virtual void UpdateAtEndOfGrid();
protected:
 

  /**
    This function updates the data depending on space mesh at a change of grid
    time.

    @attention{ The function Swap() must be called just before }
    
    @param pdOldSpots Pointer to the space mesh of the old grid
    @param nNbOldSpots size of the space mesh of the old grid
    @param pdSpots Pointer to the space mesh of the new grid
    @param nNbSpots size of the space mesh of the new grid

   */
  virtual void InterpWithPassageOfSpaceMesh(const double *pdOldSpots, 
                  size_t nNbOldSpots, const double *pdSpots, size_t nNbSpots);

   /**
     Updates constraints.
  
     Query params to update the constraints
   */
  virtual void UpdateConstraints() = 0;

  /**
    Apply the constraints to all arrays depends on constraints.
    */
  virtual void ApplyConstraintsToAll();

  /// Maximum mesh size in multi-grid case
  size_t m_nNbSpotsMax;

  /// The backward mesh manager for multi-grids case.
  pricing::BackwardMeshManagerMulti &m_meshes;

private:

  NO_COPY_CLASS(InstDataMulti);

}; // class InstDataMulti

} // namespace ihg

} // namespace ito33 

#endif // #ifndef _IHG_INSTDATAMULTI_H_

