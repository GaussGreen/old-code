/////////////////////////////////////////////////////////////////////////////
// Name:        hg/sensitivitytype.h
// Purpose:     enum for the sensitivity on HG model parameter
// Created:     2005/05/27
// RCS-ID:      $Id: sensitivitytype.h,v 1.1 2005/05/27 07:37:02 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/sensitivitytype.h
   @brief enum for the sensitivity on HG model parameter
 */

#ifndef _HG_SENSITIVITYTYPE_H_
#define _HG_SENSITIVITYTYPE_H_


namespace ito33
{

namespace hg
{


/**
   The type of sensitvity.

   Makes it easy to differentiate between the different model parameters.
 */
enum SensitivityType
{
  SensitivityType_Volatility,

  SensitivityType_DefaultIntensity,

  SensitivityType_JumpIntensity,

  SensitivityType_JumpAmplitude,

  SensitivityType_Undefined
};


} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_SENSITIVITYTYPE_H_
