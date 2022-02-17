/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/hg/hedgeratiodata.h
// Purpose:     hg output structure for individual hedge contracts
// Created:     2005/12/01
// RCS-ID:      $Id: hedgeratiodata.h,v 1.5 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/hg/hedgeratiodata.h

    Hedging output class for a single hedge contract in the 
    homogeneous (HG) model.
 */

#ifndef _ITO33_HG_HEDGERATIODATA_H_
#define _ITO33_HG_HEDGERATIODATA_H_

#include "ito33/sharedptr.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/hg/common.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL Derivative;
}

namespace hg
{

/**
    Container for a hedge ratio with associated price output and derivative.

    @nocreate
 */
class ITO33_HG_DLLDECL HedgeRatioData
{
public:

  /// Default constructor.
  HedgeRatioData() { }

  /**
      @internal
      @brief The hedge ratio of the associated derivative.

      @param dRatio the hedge ratio
      @noexport
   */
  void SetRatio(double dRatio)
  {
    m_dRatio = dRatio;
  }

  /**
      @internal
      @brief The derivative used for hedging.

      @param pDerivative the hedging derivative
      @noexport
   */
  void SetDerivative(shared_ptr<finance::Derivative> pDerivative)
  {
    m_pDerivative = pDerivative;
  }
  
  /**
      @internal
      @brief The pricing model output of the associated derivative.

      @param pModelOutput the pricing model output
      @noexport
   */
  void SetModelOutput(shared_ptr<finance::ModelOutput> pModelOutput)
  {
    m_pModelOutput = pModelOutput;
  }

  /**
      The hedge ratio of the associated derivative.

      @return The hedge ratio
   */
  double GetRatio() const
  {
    return m_dRatio;
  }

  /**
      The derivative used for hedging.

      @return The hedging derivative
   */
  shared_ptr<finance::Derivative> GetDerivative() const
  {
    return m_pDerivative;
  }
  
  /**
      The pricing model output of the associated derivative.

      @return The pricing model output
   */
  shared_ptr<finance::ModelOutput> GetModelOutput() const
  {
    return m_pModelOutput;
  }


protected:

  /// The hedge ratio
  double m_dRatio;

  /// The hedging derivative
  shared_ptr<finance::Derivative> m_pDerivative;

  /// The pricing model output of the associated derivative
  shared_ptr<finance::ModelOutput> m_pModelOutput;

}; // class HedgeRatioData

} // namespace hg

} // namespace ito33

#endif // #ifndef _ITO33_HG_HEDGERATIODATA_H_
