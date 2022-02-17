/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/backwardmeshmanager_withspacemesh.h
// Purpose:     mesh manager for backward PDE problems with space mesh
// Author:      Wang
// Created:     2004/09/02
// RCS-ID:      $Id: backwardmeshmanager_withspacemesh.h,v 1.4 2006/05/05 14:36:48 dave Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/pricing/backwardmeshmanager_withspacemesh.h
   @brief mesh manager for backward problem having space mesh
 */

#ifndef _ITO33_PRICING_BACKWARDMESHMANAGER_WITHSPACEMESH_H_
#define _ITO33_PRICING_BACKWARDMESHMANAGER_WITHSPACEMESH_H_

#include "ito33/pricing/backwardmeshmanager.h"

namespace ito33
{

namespace pricing
{


/// Backward time mesh manager with space mesh
class BackwardMeshManagerWithSpaceMesh : public BackwardMeshManager                                                          
{

public:
  
  /**
     Ctor, call the ctor of the base class MeshManagerFix

     @param params The underlying contract paramaters
     @param model The model used for the pricing
   */  
  BackwardMeshManagerWithSpaceMesh(Params& params, Model& model) 
                                 : BackwardMeshManager(params, model)
  { }

  /// dummy virtual dtor
  virtual ~BackwardMeshManagerWithSpaceMesh() { }

  /// Set up the meshes and the data 
  virtual void SetupMe();

  /// Get the current number of space points
  size_t GetNbS() const { return m_nNbS; }
  
  /// Get a pointer to the space mesh storage
  const double* GetS() const { return m_pdS.Get(); }

  /// Get a pointer to the space mesh storage
  const double* GetLogS() const { return m_pdLogS.Get(); }
  
  /// Check for end of subgrid (default to false, implement in derived classes)
  virtual bool IsEndOfGrid() { return false; }

protected:

  /// construct the space mesh (implement in derived classes)
  virtual void ConstructSpaceMesh() = 0;

  /// current number of space points.
  size_t m_nNbS;

  /// storage for current space mesh, shouldn't be changed once initialized.
  Array<double> m_pdS;

  /// storage for current log space mesh, shouldn't be changed once initialized.
  Array<double> m_pdLogS;


private:

  NO_COPY_CLASS(BackwardMeshManagerWithSpaceMesh);

}; // class BackwardMeshManagerWithSpaceMesh


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_BACKWARDMESHMANAGER_WITHSPACEMESH_H_

