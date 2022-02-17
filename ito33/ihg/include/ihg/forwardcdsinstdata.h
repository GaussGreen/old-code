/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/forwardcdsinstdata.h
// Purpose:     cds instdata class for forward PDE
// Author:      David
// Created:     2004/03/29
// RCS-ID:      $Id: forwardcdsinstdata.h,v 1.8 2006/04/21 09:26:05 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/forwardcdsinstdata.h
    @brief cds instdata class for forward PDE

    Implementation of the CDS InstData class for forward pricing.
 */

#ifndef _IHG_FORWARDCDSINSTDATA_H_
#define _IHG_FORWARDCDSINSTDATA_H_

#include "ito33/array.h"

#include "ihg/model.h"
#include "ihg/instdata.h"

namespace ito33
{

namespace pricing
{
  class ForwardCDSParams;
  class ForwardCDSMeshManager;
}

namespace ihg
{


class ForwardCDSInstData : public InstData
{
public:

ForwardCDSInstData(pricing::ForwardCDSParams &params, 
                   Model &model,
                   pricing::ForwardCDSMeshManager &meshes);

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
    
     At the start period, the mesh manager tell me what are the start conditions
     for the Data part. It also works on the InstParam part.
   */
  virtual void SetInitialValue();

  /// The CDS price at the current timestep (ie the desired solution)
  double m_dCDSPrice;

  /// The recovery part of the current CDS price
  double m_dRecoveryTerm;

  /// The spread part of the current CDS price
  double m_dSpreadTerm;

  /// The accrued part of the current CDS price
  double m_dAccruedTerm;

  /// The recovery term (X_d) at previous timestep
  double m_dRecoveryTermOld;

  /// The sum of the spread term up to the last spread date
  double m_dSpreadTermOld;

  /// The accrued spread term (Y_d) at previous timestep
  double m_dAccruedTermOld;

  /// The rhs of the recovery ODE (X_d) at previous timestep
  double m_dRecoveryRHSOld;

  /// The rhs of the accrued spread ODE (Y_d) at previous timestep
  double m_dAccruedRHSOld;

  /// The current timestepping index
  size_t m_nCurrentIndex;

  /// The recovery amount
  double m_dRecovery;

  /// The spread amount
  double m_dSpread;

  /// The fraction of the spread owned (accrued part, as fraction of spread)
  double m_dAccruedFraction;
  
  /// Is the hazard rate time only? If so, use linear solver
  bool m_bIsHRTimeOnly;


protected:

  /**
    Allocate memory

    @param nNbS The grid size
  */
  void Alloc(size_t nNbS);

  /// The params related to cds
  pricing::ForwardCDSParams &m_forwardCDSParams;
  
  /// The forward cds mesh manager
  pricing::ForwardCDSMeshManager &m_forwardCDSMeshes;

private:

  NO_COPY_CLASS(ForwardCDSInstData);

}; // class ForwardCDSInstData


} // namespace ihg

} // namespace ito33 

#endif // #ifndef _IHG_FORWARDCDSINSTDATA_H_
