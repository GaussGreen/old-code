/////////////////////////////////////////////////////////////////////////////
// Name:        hg/backwardinstdata.h
// Purpose:     HG instdata class for backward solving
// Created:     2005/01/13
// RCS-ID:      $Id: backwardinstdata.h,v 1.15 2006/03/31 17:43:59 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/backwardinstdata.h
   @brief HG instdata class for backward solving
 */

#ifndef _HG_BACKWARDINSTDATA_H_
#define _HG_BACKWARDINSTDATA_H_

#include "ito33/vector.h"
#include "ito33/array.h"
#include "ito33/dlldecl.h"

#include "hg/instdata.h"
#include "hg/sensitivitydata.h"
#include "hg/sensitivitymethod.h"
#include "hg/sensitivitybyadjointdata.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL ComputationalFlags;
}
  
namespace numeric 
{   
  class TridiagonalMatrix;   
  class MorseMatrix;
}

namespace hg
{
  

/// HG instdata class for backward solving
class BackwardInstData : public InstData
{

public:

  BackwardInstData(pricing::Params& params, 
                   Model& model, 
                   pricing::MeshManager& meshes)
                 : InstData(params, model, meshes)
  {
    // Set all HG model dependant compute flags as false
    m_bComputeFugit = false;

    // Initial recovery value may be referenced by numoutput
    m_dRecoveryValue = 0.0;

    // by default, don't store dual system
    m_bDualSystemRequired = false;

    // the flags should be set manually and explicitly later.
  }

  virtual ~BackwardInstData() { }

  // Functions required by the Engine

  // Init() : No implementation at this level, so no need to redeclare

  // UpdateBeforeStep(): Same implementation as in base, no need to reimplement

  // DoEvents(): Same implementation as in base, no need to reimplement

  // SetInitialValue(): No implementation at this level, so no need to redeclare

  /**
     Sets the relevant flags (greeks, sensitivities, etc) for this class.

     Derived classes are expected to overload this function to setup
     contract specific flags in the derived classes.

     @param flags the computational flags (usually set by the user)
   */
  virtual void SetupFlags(const finance::ComputationalFlags& flags);
  
  /// A dummy function just to work around for BackwardNumOutput
  virtual const int* GetConstraintFlags() const { return 0; }

  /**
     Gets initial spot.

     @return initial spot in params
   */
  double GetInitialSpot() const;

  /**
     Gets space mesh information for numOutput.

     @param nNbS (output) number of space points
     @return space points
   */
  const double* GetSpaceMesh(size_t& nNbS)
  {
    nNbS = m_nNbS;
    return m_pdS;
  }

  /**
     Gets the recovery value at initial spot.

     This function is reimplemented as the recovery value function
     is not constant for exchangeable.

     @return the recovery value at initial spot.
   */
  virtual double GetValueAfterDefaultAtInitialSpot() 
  {
    return m_dRecoveryValue; 
  }

  /**
     Sets boundary conditions for sensitivity pdes.

     Primarily used to reset Dirichlet nodes. Assume Dirichlet values for
     doesn't depend on the parameters, the Dirichlet boundary condition of 
     sensitivity PDE should be zero.
   */
  void SetSensitivityBoundary();


public:

  /// Default value at current time step
  double m_dRecoveryValue;

  /// @name sensitivity related data
  //@{

  /// Flags indicate if we compute sensitivities
  std::vector<bool> m_pbComputeSensitivities;

  /// Details of each model param
  std::vector< SensitivityData > m_pSensitivityData;

  /// Sensitivities at current step
  std::vector< Array<double> > m_ppdSensitivities;

  /// Sensitivities at last time step
  std::vector< Array<double> > m_ppdOldSensitivities;

  /// Sensitivities at the time step before last (for multistep methods)
  std::vector< Array<double> > m_ppdOldOldSensitivities;

  /// The sparse matrix pointer for sensitivity by adjoint method
  numeric::MorseMatrix* m_pSparseMatrix;

  /// Mass matrix for sensitivity by adjoint method using FE discretization
  numeric::TridiagonalMatrix* m_pMassMatrix;

  /// Data required by the adjoint method at the current timestep
  SensitivityByAdjointData m_aData;

  /// Method used to compute sensitivities
  SensitivityMethod m_sensitivityMethod;

  /**
      If the dual system is required: sensitivity by adjoint or 
      forward/(backward) by backward(forward) transposed.
   */
  bool m_bDualSystemRequired;

  //@}
 
  /// @name fugit related data
  //@{
   
  /// Flag indicates if we compute fugit
  bool m_bComputeFugit;

  /// Fugit at current time step (to be calculated).
  Array<double> m_pdFugits;

  /// Fugit at last time step.
  Array<double> m_pdOldFugits;

  /// Fugit at the time step before last (for multistep methods).
  Array<double> m_pdOldOldFugits;
  
  //@}


protected:

  /**
     Allocate the memory for the price, fugit arrays

     @param nNbS the size of the arrays to be allocated
   */
  void Alloc(size_t nNbS);

  /**
     Apply an event to prices and Greeks
     
     @param pEvent A pointer to a to be applied event
   */
  void ApplyEvent(const pricing::Event* pEvent);

  /** 
     Swap the price and Greek (if any) arrays
   */
   virtual void Swap();

  /**
     Sets initial value for sensitivities.
   */
  void SetInitialSensitivityValue();


private:

  NO_COPY_CLASS(BackwardInstData);

}; // class BackwardInstData


} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_BACKWARDINSTDATA_H_
