/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/forwardmeshmanager_fix.h
// Purpose:     mesh manager for forward PDE problems with fixed space mesh
// Author:      Wang
// Created:     2004/02/11
// RCS-ID:      $Id: forwardmeshmanager_fix.h,v 1.4 2004/10/05 09:13:39 pedro Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/pricing/forwardmeshmanager_fix.h
   @brief mesh manager for fixed space mesh forward problem
 */

#ifndef _ITO33_PRICING_FORWARDMESHMANAGER_FIX_H_
#define _ITO33_PRICING_FORWARDMESHMANAGER_FIX_H_

#include "ito33/pricing/forwardmeshmanager.h"
#include "ito33/pricing/fixedspacemeshmanager.h"

namespace ito33
{

namespace pricing
{


/// Backward time mesh manager
class ForwardMeshManagerFix : public ForwardMeshManager, 
                              public FixedSpaceMeshManager
                             
{

public:
  
  /**
     Ctor, call the ctor of the base class MeshManagerFix

     @param params The underlying contract paramaters
     @param model The model used for the pricing
   */  
  ForwardMeshManagerFix(Params &params, Model &model) 
                     : ForwardMeshManager(params, model)
  {
  }

  /// dummy virtual dtor
  virtual ~ForwardMeshManagerFix() { }

  /// Set up the meshes and the data 
  virtual void SetupMe();

  
private:

  NO_COPY_CLASS(ForwardMeshManagerFix);

}; // class ForwardMeshManagerFix


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_FORWARDMESHMANAGER_FIX_H_

