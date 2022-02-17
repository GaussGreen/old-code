/////////////////////////////////////////////////////////////////////////////
// Name:        hg/sensitivitybyadjointdata.h
// Purpose:     class storing sensitivity data for the HG model parameters
// Created:     2005/05/27
// RCS-ID:      $Id: sensitivitybyadjointdata.h,v 1.7 2006/01/24 10:49:29 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/sensitivitybyadjointdata.h
   @brief class storing sensitivity data using adjoint method for the HG model
          parameters
 */

#ifndef _HG_SENSITIVITYBYADJOINTDATA_H_
#define _HG_SENSITIVITYBYADJOINTDATA_H_

#include "ito33/autoptr.h"
#include "ito33/array.h"

#include "hg/sensitivitytype.h"

namespace ito33
{

namespace numeric
{
  class TridiagonalMatrix;
  class MorseMatrix;
  class InterpolationMatrix;
}

namespace hg
{


/**
   datas of the non trivial systems.

   Note that we should try the best to avoid this class bloated
 */
struct SensitivityByAdjointData
{
  SensitivityByAdjointData() : m_bIsValid(false) { }

  /// Boolean helps Numoutput to take or not the current data
  bool m_bIsValid;

  /// The coefficient om the old prices/matrix
  double m_dOldTimeWeight;
  
  /// The coefficient on the old old prices/matrix
  double m_dOldOldTimeWeight;
  
  /**
     The recovery value used to do \frac{\partial B}{\partial p} 
     where p is one of the default intensity.
   */
  double m_dRecoveryValue;

  /// The tridiagonal matrix at each time step
  AutoPtr<numeric::TridiagonalMatrix> m_pMatrix;

  /// Interpolation due to event(mainly dividend) if any
  AutoPtr<numeric::InterpolationMatrix> m_pInterpMatrix;

  /**
     Interpolation on the current observation points if any.
     It can be merged with the residuals so that only the product of 
     transpose matrix vector is stored(require more memory in general).
   */
  AutoPtr<numeric::InterpolationMatrix> m_pRMatrix;

  /// The (weighted) residuals at the current observation points if any
  Array<double> m_pdResiduals;

  /**
     The discretization method used by finite difference at each point,
     not used by finite element implementation.
   */
  Array<int> m_piFDFlags;

  /**
     The constraint flags if any, current implementation doesn't use it,
     might be removed once the implementation is well tested.
   */
  Array<int> m_piConstraintFlags;

  /// The solution of the system at current time
  Array<double> m_pdPrices;

}; // struct SensitivityByAdjointData


} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_SENSITIVITYBYADJOINTDATA_H_
