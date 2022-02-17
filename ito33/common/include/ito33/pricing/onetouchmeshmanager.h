/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/onetouchmeshmanager.h
// Purpose:     OneTouch mesh namager class (backward)
// Created:     2005/07/04
// RCS-ID:      $Id: onetouchmeshmanager.h,v 1.1 2005/07/05 10:03:23 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/pricing/onetouchmeshmanager.h
   @brief Implementation of the meshes (space and time) manager for OneTouch.
 */

#ifndef _ITO33_PRICING_ONETOUCHMESHMANAGER_H_
#define _ITO33_PRICING_ONETOUCHMESHMANAGER_H_

#include "ito33/pricing/backwardmeshmanager_fix.h"

namespace ito33
{

namespace pricing
{

  class OneTouchParams;


/**
   Class for managing the space and time meshes for OneTouch contract
 */
class OneTouchMeshManager: public BackwardMeshManagerFix
{
public:

  OneTouchMeshManager(OneTouchParams& params, Model& model);

  /// Get the current recovery value
  double GetRecoveryValue() const { return m_pdRecoveryValues[m_nIdx]; }

  /// Get the current boundary value
  double GetBoundaryValue() const
  {
    return m_pdBoundaryValues[m_nIdx]; 
  }

  /// Setup this OneTouch meshmanager
  virtual void SetupMe();

  // default dctor is ok


protected:
  
  /// Construct the space mesh 
  virtual void ConstructSpaceMesh();

  /// Pre-computed recovery values
  Array<double> m_pdRecoveryValues;

  /// Pre-computed boundary values
  Array<double> m_pdBoundaryValues;

  /// The params for one touch
  OneTouchParams& m_oneTouchParams;


private:

  NO_COPY_CLASS(OneTouchMeshManager);

}; // class OneTouchMeshManager


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_ONETOUCHMESHMANAGER_H_
