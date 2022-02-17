/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/forwardoptionmeshmanager.h
// Purpose:     option mesh mamager class for forward PDE
// Author:      Wang
// Created:     2004/03/04
// RCS-ID:      $Id: forwardoptionmeshmanager.h,v 1.5 2004/10/05 09:13:39 pedro Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/forwardoptionmeshmanager.h
    @brief option mesh manager class for forward PDE
 */

#ifndef _ITO33_PRICING_FORWARDOPTIONMESHMANAGER_H_
#define _ITO33_PRICING_FORWARDOPTIONMESHMANAGER_H_

#include "ito33/pricing/forwardoptionparams.h"
#include "ito33/pricing/forwardmeshmanager_fix.h"

namespace ito33
{

namespace pricing
{

/**
   Class for managing the space and time meshes for forward options.

   Based on the forward, fixed mesh classes.
 */
class ForwardOptionMeshManager: public ForwardMeshManagerFix
{
public:

  ForwardOptionMeshManager(ForwardOptionParams &params, Model &model)
                         : ForwardMeshManagerFix(params, model),
                           m_forwardOptionParams(params)  { }

  /// get the current recovery value
  double GetRecoveryValue() const { return m_pdRecoveryValues[m_nIdx]; }

  /// setup this option meshmanager
  virtual void SetupMe();

  // Default dtor is ok


protected:
  
  /// Construct the space mesh 
  virtual void ConstructSpaceMesh();

  /// Construct a uniform space mesh 
  virtual void ConstructUniformSpaceMesh();

  /// pre-computed recovery values
  Array<double> m_pdRecoveryValues;

  /// the params for option (forward PDE)
  ForwardOptionParams &m_forwardOptionParams;


private:

  NO_COPY_CLASS(ForwardOptionMeshManager);

}; // class ForwardOptionMeshManager


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_FORWARDOPTIONMESHMANAGER_H_

