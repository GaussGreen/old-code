/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/edsmeshmanager.h
// Purpose:     EDS mesh namager class (backward)
// Created:     2005/01/26
// RCS-ID:      $Id: edsmeshmanager.h,v 1.1 2005/01/26 18:24:14 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/pricing/edsmeshmanager.h
   @brief Implementation of the meshes (space and time) manager for EDS.
 */

#ifndef _ITO33_PRICING_EDSMESHMANAGER_H_
#define _ITO33_PRICING_EDSMESHMANAGER_H_

#include "ito33/common.h"

#include "ito33/pricing/backwardmeshmanager_fix.h"

namespace ito33
{

namespace pricing
{

  class EDSParams;

/**
   Class for managing the space and time meshes for CDS contracts
 */
class EDSMeshManager: public BackwardMeshManagerFix
{
public:

  EDSMeshManager(EDSParams& params, Model& model);

  /// Get the current recovery value
  double GetRecoveryValue() const { return m_pdRecoveryValues[m_nIdx]; }

  /// Setup this EDS meshmanager
  virtual void SetupMe();

  // default dctor is ok


protected:
  
  /// Construct the space mesh 
  virtual void ConstructSpaceMesh();

  /// Pre-computed recovery values
  Array<double> m_pdRecoveryValues;

  /// The params for option
  EDSParams& m_edsParams;


private:

  NO_COPY_CLASS(EDSMeshManager);

}; // class EDSMeshManager


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_EDSMESHMANAGER_H_

