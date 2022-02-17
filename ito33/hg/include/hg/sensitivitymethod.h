/////////////////////////////////////////////////////////////////////////////
// Name:        hg/sensitivitytype.h
// Purpose:     enum for the method used to compute sensitivities
// Created:     2005/06/09
// RCS-ID:      $Id: sensitivitymethod.h,v 1.1 2005/06/10 16:13:07 dave Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/sensitivitymethod.h
   @brief enum for the method used to compute sensitivities
 */

#ifndef _HG_SENSITIVITYMETHOD_H_
#define _HG_SENSITIVITYMETHOD_H_


namespace ito33
{

namespace hg
{


/**
   The method used to compute sensitvities.
 */
enum SensitivityMethod
{
  // Sensitivities are not computed
  SensitivityMethod_None,

  // Compute via PDEs
  SensitivityMethod_PDE,

  // Compute via the adjoint method
  SensitivityMethod_Adjoint
};


} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_SENSITIVITYMETHOD_H_
