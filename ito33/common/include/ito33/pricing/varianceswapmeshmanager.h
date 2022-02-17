/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/varianceswapmeshmanager.h
// Purpose:     variance swap mesh manager class
// Created:     2006/02/21
// RCS-ID:      $Id: varianceswapmeshmanager.h,v 1.8 2006/08/19 22:01:27 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/varianceswapmeshmanager.h
    @brief variance swap mesh manager class

    Implementation of the meshes (space and time) manager for variance swaps.
    
    The class contains not only the meshes but also all data related to
    the meshes.
 */

#ifndef _ITO33_PRICING_VARIANCESWAPMESHMANAGER_H_
#define _ITO33_PRICING_VARIANCESWAPMESHMANAGER_H_

#include "ito33/array.h"

#include "ito33/pricing/backwardmeshmanager_fix.h"

namespace ito33
{

namespace pricing
{

  class VarianceSwapParams;

/// Class for managing the space and time meshes for variance swaps.
class VarianceSwapMeshManager: public BackwardMeshManagerFix
{
public:

  /**
     Constructor.

     An external space mesh can be specified so that the meshes in different
     paths of a path dependent pricing can be coordinated.  If the external 
     space mesh is empty, a new space mesh is constructed when SetupMe() is 
     called.

     @param params The variance swap params related to this mesh
     @param model The model used for pricing
     @param pdSpaceMesh An external space mesh to use
   */
  VarianceSwapMeshManager(VarianceSwapParams& params, 
                          Model& model,
                          const std::vector<double>& pdSpaceMesh);

  /// Dummy virtual dctor
  virtual ~VarianceSwapMeshManager() { }

  /// Set as the fixed leg of variance swap for computing recovery values
  void SetIsConditionalFixed(bool bIsConditionalFixed)
  {
    m_bIsConditionalFixed = bIsConditionalFixed;
  }

  /// Get the current recovery value
  double GetRecoveryValue() const { return m_pdRecoveryValues[m_nIdx]; }

  /// Setup this variance swap meshmanager
  virtual void SetupMe();
  
  /// Pre-computes recovery values
  virtual void ComputeRecoveryValues();

  /// Return the number of times requested by params, adjusted for var swaps.
  size_t GetNbRequestedTimes() const;

protected:
  
  /// Constructs a non-uniform space mesh 
  virtual void ConstructSpaceMesh();

  /// Constructs a uniform space mesh 
  virtual void ConstructUniformSpaceMesh();

  /// Pre-computed recovery values
  Array<double> m_pdRecoveryValues;

  /// The params for the variance swap
  VarianceSwapParams& m_varianceSwapParams;

  /// Get the number of observations up to (but not including) dTime
  size_t GetObservationNumber(double dTime, bool &bIsTimeOnEvent);

  /// Is this the fixed leg of a conditional swap
  bool m_bIsConditionalFixed;

private:

  NO_COPY_CLASS(VarianceSwapMeshManager);

}; // class VarianceSwapMeshManager


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_VARIANCESWAPMESHMANAGER_H_
