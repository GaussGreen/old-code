/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/optionmeshmanager.h
// Purpose:     option mesh namager class
// Created:     2004/02/11
// RCS-ID:      $Id: optionmeshmanager.h,v 1.14 2006/08/19 22:01:27 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/optionmeshmanager.h
    @brief option mesh manager class

    Implementation of the meshes (space and time) manager for options.
    
    The class contains not only the meshes but also all data related to
    this meshes.
 */

#ifndef _ITO33_PRICING_OPTIONMESHMANAGER_H_
#define _ITO33_PRICING_OPTIONMESHMANAGER_H_

#include "ito33/pricing/backwardmeshmanager_fix.h"

namespace ito33
{

namespace pricing
{

  class OptionParams;

/**
    Class for managing the space and time meshes for regular options.

    Based on the backward, fixed mesh manager class. Notable functions 
    include 
    - constructing the space mesh (uniform or non-uniform)
    - getting the recovery value in case of default
    - general SetupMe() function
 */
class OptionMeshManager: public BackwardMeshManagerFix
{
public:

  /// Constructor
  OptionMeshManager(OptionParams& params, Model& model);

  /// Get the current recovery value
  double GetRecoveryValue() const { return m_pdRecoveryValues[m_nIdx]; }

  /// Setup this option meshmanager
  virtual void SetupMe();
  
  /// Pre-computes recovery values
  virtual void ComputeRecoveryValues();
  
  /// Dummy virtual dctor
  virtual ~OptionMeshManager() { }


protected:
  
  /// Construct a non-uniform space mesh 
  virtual void ConstructSpaceMesh();

  /// Construct a uniform space mesh 
  virtual void ConstructUniformSpaceMesh();

  /// Pre-computed recovery values
  Array<double> m_pdRecoveryValues;

  /// The params for option
  OptionParams& m_optionParams;


private:

  NO_COPY_CLASS(OptionMeshManager);

}; // class OptionMeshManager


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_OPTIONMESHMANAGER_H_

