/////////////////////////////////////////////////////////////////////////////
// Name:        hg/sensitivitydata.h
// Purpose:     class storing sensitivity data for the HG model parameters
// Created:     2005/04/19
// RCS-ID:      $Id: sensitivitydata.h,v 1.3 2005/05/27 07:37:52 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/sensitivitydata.h
   @brief class storing sensitivity data for the HG model parameters
 */

#ifndef _HG_SENSITIVITYDATA_H_
#define _HG_SENSITIVITYDATA_H_

#include "ito33/autoptr.h"
#include "ito33/constants.h"

#include "hg/sensitivitytype.h"

namespace ito33
{

namespace numeric
{
  class MorseMatrix;
}

namespace hg
{


/**
   Data for each model parameter.

   Used when computing the sensitivities, since each sensitivity type has
   a different form of PDE.  Extra data for parameter can be stored, such
   as jump amplitudes and intensity.  Sometimes this extra data is needed
   when constructing the appropriate PDE.  Sometimes this data is not
   needed, and doesn't even make sense for the type of parameter.  It may
   be better to split this into different classes, or a class hierarchy, but
   with the limited number of parameter types, storing everything in one 
   struct is just easier. 
 */
struct SensitivityData
{
  SensitivityData() :
    m_sensitivityType(SensitivityType_Undefined),
    m_nRegime1(INVALIDINDEX),
    m_nRegime2(INVALIDINDEX),
    m_dAmplitude(0.0),
    m_dIntensity(0.0)
    {}

  /// The type of model parameter
  SensitivityType m_sensitivityType;

  /// The regime to which the parameter applies
  size_t m_nRegime1;

  /// If a jump, the regime it jumps to
  size_t m_nRegime2;

  /// If a jump, the jump amplitude
  double m_dAmplitude;

  /// If a jump, the jump intensity
  double m_dIntensity;

  /// The amplitude and intensity PDEs need a sparse matrix
  mutable AutoPtr<numeric::MorseMatrix> m_pSparseMatrix;

}; // struct sensitivitydata


} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_SENSITIVITYDATA_H_
