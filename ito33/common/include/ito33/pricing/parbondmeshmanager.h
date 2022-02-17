/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/parbondmeshmanager.h
// Purpose:     parbond mesh namager class (backward)
// Author:      ZHANG
// Created:     2005/05/20
// RCS-ID:      $Id: parbondmeshmanager.h,v 1.1 2005/06/08 15:41:22 zhang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/pricing/parbondmeshmanager.h
   @brief parbond mesh manager class (backward)

   Implementation of the meshes (space and time) manager for parbond.
 */

#ifndef _ITO33_PRICING_ParBondMESHMANAGER_H_
#define _ITO33_PRICING_ParBondMESHMANAGER_H_

#include "ito33/common.h"

#include "ito33/pricing/backwardmeshmanager_fix.h"

namespace ito33
{

namespace pricing
{

  class ParBondParams;

/**
   Class for managing the space and time meshes for ParBond contracts
 */
class ParBondMeshManager: public BackwardMeshManagerFix
{
public:

  ParBondMeshManager(ParBondParams &params, Model &model);

  // default dctor is ok

  /// Get the current recovery value
  double GetRecoveryValue() const { return m_pdRecoveryValues[m_nIdx]; }

  /// Setup this parbond meshmanager
  virtual void SetupMe();

  /// Time mesh construction only for time only pricer
  void SetupMeTimeOnly();


protected:

  /// Compute the recovery values
  void ComputeRecoveryValues();
  
  /// Construct the space mesh 
  virtual void ConstructSpaceMesh();

  /// Pre-computed recovery values
  Array<double> m_pdRecoveryValues;

  /// The params for option
  ParBondParams& m_parbondParams;


private:

  NO_COPY_CLASS(ParBondMeshManager);

}; // class ParBondMeshManager


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_ParBondMESHMANAGER_H_

