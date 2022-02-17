/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/asianoptionmeshmanager.h
// Purpose:     Asian option mesh manager class
// Author:      Ito33 Canada
// Created:     2004/06/02
// RCS-ID:      $Id: asianoptionmeshmanager.h,v 1.4 2006/08/19 22:01:27 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/asianoptionmeshmanager.h
    @brief Asian option mesh manager class

    Implementation of the meshes (space and time) manager for asian options.
 */

#ifndef _ITO33_PRICING_ASIANOPTIONMESHMANAGER_H_
#define _ITO33_PRICING_ASIANOPTIONMESHMANAGER_H_

#include "ito33/pricing/optionmeshmanager.h"

namespace ito33
{

namespace pricing
{

  class AsianOptionParams;

/**
   Class for managing the space and time meshes for asian options
 */
class AsianOptionMeshManager: public OptionMeshManager
{

public:

  AsianOptionMeshManager(AsianOptionParams& params, Model& model);

  /// Pre-computes recovery values
  void ComputeRecoveryValues();
  
  /// Dummy virtual dctor
  virtual ~AsianOptionMeshManager() { }

protected:

  /// The params for asian option
  AsianOptionParams& m_asianOptionParams;

  /// Gets the number of observations up to (but not including) dTime
  size_t GetObservationNumber(double dTime);

private:

  NO_COPY_CLASS(AsianOptionMeshManager);

}; // class AsianOptionMeshManager


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_ASIANOPTIONMESHMANAGER_H_
