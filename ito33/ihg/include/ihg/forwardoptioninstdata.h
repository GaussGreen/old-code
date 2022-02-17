/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/forwardoptioninstdata.h
// Purpose:     optioin instdata class for forward PDE
// Author:      Wang
// Created:     2004/03/10
// RCS-ID:      $Id: forwardoptioninstdata.h,v 1.5 2004/10/04 18:04:04 pedro Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/forwardoptioninstdata.h
    @brief option instdata class for forward PDE

    Implementation of the Option InstData class for forward pricing.
 */

#ifndef _IHG_FORWARDOPTIONINSTDATA_H_
#define _IHG_FORWARDOPTIONINSTDATA_H_

#include "ito33/pricing/forwardoptionparams.h"
#include "ito33/pricing/forwardoptionmeshmanager.h"

#include "ihg/model.h"
#include "ihg/instdata.h"

namespace ito33
{

namespace ihg
{

/**
   Option instdata class for forward equation
 */
class ForwardOptionInstData : public InstData
{
public:

  ForwardOptionInstData(pricing::ForwardOptionParams& params, 
	                      Model& model,
					              pricing::ForwardOptionMeshManager& meshes)
    : InstData(params, model, meshes), 
  	  m_forwardOptionParams(params), m_forwardOptionMeshes(meshes)
  {
    m_bIsHRTimeOnly = false;
  }

  // Default dtor is ok

  /** 
     Allocates memory when the maximum size of the space mesh is known.
     It does necessary initialization work as well.

     m_meshes.SetupMe() must be called before calling this function.   
   */
  virtual void Init();


  /**
     Update variables at the beginning of each time step
    
     Required by Engine class. At each time step, this function 
     must update the variables required by the stepper class. The
     stepper class will update the solution values.
   */
  void UpdateBeforeStep();

  /**
     Set the initial values of myself from the mesh manager.
    
     At the start period, the mesh manager tell me what are the start 
     conditions for the Data part. It also works on the InstParam part.
   */
  virtual void SetInitialValue();

  /// Is the hazard rate time only? If so, use linear solver
  bool m_bIsHRTimeOnly;


protected:

  void Alloc(size_t nNbS);

  /// The params related to option
  pricing::ForwardOptionParams& m_forwardOptionParams;
  
  /// the option mesh manager
  pricing::ForwardOptionMeshManager& m_forwardOptionMeshes;


private:

  NO_COPY_CLASS(ForwardOptionInstData);

}; // class ForwardOptionInstData


} // namespace ihg

} // namespace ito33 

#endif // #ifndef _IHG_FORWARDOPTIONINSTDATA_H_

