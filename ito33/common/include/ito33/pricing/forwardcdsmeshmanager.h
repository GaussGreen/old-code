/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/forwardcdsmeshmanager.h
// Purpose:     cds mesh mamager class for forward PDE pricing
// Author:      David
// Created:     2004/03/31
// RCS-ID:      $Id: forwardcdsmeshmanager.h,v 1.8 2006/08/19 22:01:27 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/forwardcdsmeshmanager.h
    @brief cds mesh manager class for forward PDE pricing
 */

#ifndef _ITO33_PRICING_FORWARDCDSMESHMANAGER_H_
#define _ITO33_PRICING_FORWARDCDSMESHMANAGER_H_

#include "ito33/pricing/forwardmeshmanager_fix.h"
#include "ito33/pricing/forwardcdsparams.h"

namespace ito33
{

namespace pricing
{

/**
    Class for managing the space and time meshes for forward cds pricing.

    Based on the forward, fixed mesh classes.
 */
class ForwardCDSMeshManager: public ForwardMeshManagerFix
{
public:

  ForwardCDSMeshManager(ForwardCDSParams& params, Model& model)
                      : ForwardMeshManagerFix(params, model),
                        m_forwardCDSParams(params)  { }

  // Default dtor is ok

  /// Get the current accrued spread fraction
  double GetAccruedFraction() const { return m_pdAccruedFraction[m_nIdx]; }

  /// Get the current spread
  double GetSpread() const { return m_pdSpreads[m_nIdx]; }

  /// Setup this forward cds meshmanager
  virtual void SetupMe();

  /// Get the current timestepping index
  size_t GetCurrentIndex() { return m_nIdx; }

protected:
  
  /// Construct the space mesh 
  virtual void ConstructSpaceMesh();

  /// Construct a uniform space mesh 
  virtual void ConstructUniformSpaceMesh();

  /// The spread payment of the period in which this timestep occurs
  Array<double> m_pdSpreads;

  /// pre-computed accrued fractions at each timestep
  Array<double> m_pdAccruedFraction;

  /// the params for the cds contracts (forward PDE)
  ForwardCDSParams &m_forwardCDSParams;


private:

  NO_COPY_CLASS(ForwardCDSMeshManager);

}; // class ForwardCDSMeshManager


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_FORWARDCDSMESHMANAGER_H_

