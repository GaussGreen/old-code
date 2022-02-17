/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/cdsmeshmanager.h
// Purpose:     cds mesh namager class (backward)
// Author:      Wang
// Created:     2004/03/02
// RCS-ID:      $Id: cdsmeshmanager.h,v 1.8 2005/02/17 10:38:38 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/pricing/cdsmeshmanager.h
   @brief cds mesh manager class (backward)

   Implementation of the meshes (space and time) manager for cds.
 */

#ifndef _ITO33_PRICING_CDSMESHMANAGER_H_
#define _ITO33_PRICING_CDSMESHMANAGER_H_

#include "ito33/common.h"

#include "ito33/pricing/backwardmeshmanager_fix.h"

namespace ito33
{

namespace pricing
{

  class CDSParams;

/**
   Class for managing the space and time meshes for CDS contracts
 */
class CDSMeshManager: public BackwardMeshManagerFix
{
public:

  CDSMeshManager(CDSParams &params, Model &model);

  // default dctor is ok

  /// Get the current recovery value
  double GetRecoveryValue() const { return m_pdRecoveryValues[m_nIdx]; }

  /// Setup this cds meshmanager
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
  CDSParams& m_cdsParams;


private:

  NO_COPY_CLASS(CDSMeshManager);

}; // class CDSMeshManager


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_CDSMESHMANAGER_H_

