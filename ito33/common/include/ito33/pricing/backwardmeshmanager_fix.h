/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/backwardmeshmanager_fix.h
// Purpose:     mesh manager for backward PDE problems with fixed space mesh
// Author:      Wang
// Created:     2004/02/11
// RCS-ID:      $Id: backwardmeshmanager_fix.h,v 1.7 2004/10/05 09:13:38 pedro Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/pricing/backwardmeshmanager_fix.h
   @brief mesh manager for fixed space mesh backward problem

   @todo no more needed?
 */

#ifndef _ITO33_PRICING_BACKWARDMESHMANAGER_FIX_H_
#define _ITO33_PRICING_BACKWARDMESHMANAGER_FIX_H_

#include "ito33/pricing/backwardmeshmanager_withspacemesh.h"

namespace ito33
{

namespace pricing
{


/// Backward time mesh manager
class BackwardMeshManagerFix : public BackwardMeshManagerWithSpaceMesh
{

public:
  
  /**
     Ctor, call the ctor of the base class MeshManagerFix

     @param params The underlying contract paramaters
     @param model The model used for the pricing
   */  
  BackwardMeshManagerFix(Params &params, Model &model) 
                       : BackwardMeshManagerWithSpaceMesh(params, model)
  { }

  /// dummy virtual dtor
  virtual ~BackwardMeshManagerFix() { }


  /// Set up the meshes and the data 
  virtual void SetupMe();


private:

  NO_COPY_CLASS(BackwardMeshManagerFix);

}; // class BackwardMeshManagerFix


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_BACKWARDMESHMANAGER_FIX_H_

