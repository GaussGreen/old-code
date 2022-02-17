/////////////////////////////////////////////////////////////////////////////
// Name:        hr/heromeshmanager.h
// Purpose:     HERO mesh namager class
// Created:     2005/09/26
// RCS-ID:      $Id: heromeshmanager.h,v 1.2 2005/12/01 18:43:00 dave Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/heromeshmanager.h
    @brief HERO mesh manager class

    Implementation of the meshes (space and time) manager for HERO.
    
    The class contains not only the meshes but also all data related to
    this meshes.
 */

#ifndef _HG_HEROMESHMANAGER_H_
#define _HG_HEROMESHMANAGER_H_

#include "ito33/pricing/backwardmeshmanager_fix.h"


namespace ito33
{

namespace hg
{

  class HeroParams;
  class Model;

/**
  @brief Class for managing the space and time meshes for HERO.

  Based on the backward, fixed mesh manager class. Notable behaviour 
  for HERO include 
  - function to compute the 'h' values needed for the HERO PDE
  
*/

class HeroMeshManager: public pricing::BackwardMeshManagerFix
{
public:

  /// Constructor
  HeroMeshManager(HeroParams& params, Model& model);

  /// Dummy virtual dctor
  virtual ~HeroMeshManager() { };

  /// Get the current recovery value
  double GetRecoveryValue() const { return m_pdRecoveryValues[m_nIdx]; }

  /// Get the current 'h' value
  double GetHValue() const { return m_pdHValues[m_nIdx]; }

  /// Setup this meshmanager
  virtual void SetupMe();
  
  /// Pre-computes recovery values
  virtual void ComputeRecoveryValues();

  /// Pre-computes the 'h' values
  virtual void ComputeHValues();
  
  
protected:
  
  /// Construct a non-uniform space mesh 
  virtual void ConstructSpaceMesh();

  /// Construct a uniform space mesh 
  virtual void ConstructUniformSpaceMesh();

  /// Pre-computed recovery values
  Array<double> m_pdRecoveryValues;

  /// Pre-computed 'h' values
  Array<double> m_pdHValues;

  /// The params for HERO
  HeroParams& m_heroParams;

  /// The HG model
  Model& m_HGmodel;

private:

  NO_COPY_CLASS(HeroMeshManager);

}; // class HeroMeshManager


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_HEROMESHMANAGER_H_
