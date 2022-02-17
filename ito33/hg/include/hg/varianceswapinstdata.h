/////////////////////////////////////////////////////////////////////////////
// Name:        hg/varianceswapinstdata.h
// Purpose:     HG variance swap instdata class
// Created:     2006/03/05
// RCS-ID:      $Id: varianceswapinstdata.h,v 1.2 2006/04/10 12:02:12 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/varianceswapinstdata.h
   @brief HG variance swap instdata class
 */

#ifndef _HG_VARIANCESWAPINSTDATA_H_
#define _HG_VARIANCESWAPINSTDATA_H_

#include "hg/backwardinstdata.h"

namespace ito33
{

namespace pricing
{
  class VarianceSwapParams;
  class VarianceSwapMeshManager;
}

namespace hg
{


class VarianceSwapInstData : public BackwardInstData
{

public:

  VarianceSwapInstData(pricing::VarianceSwapParams& params,
                       Model& model,
                       pricing::VarianceSwapMeshManager& meshes);

  virtual ~VarianceSwapInstData() { }

  virtual void Init();

  void UpdateBeforeStep();

  virtual void SetInitialValue();

protected:

  /// The params related to variance swap
  pricing::VarianceSwapParams& m_varianceSwapParams;
  
  /// The variance swap mesh manager
  pricing::VarianceSwapMeshManager& m_varianceSwapMeshes;

private:

  NO_COPY_CLASS(VarianceSwapInstData);

}; // class VarianceSwapInstData


} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_VARIANCESWAPINSTDATA_H_

