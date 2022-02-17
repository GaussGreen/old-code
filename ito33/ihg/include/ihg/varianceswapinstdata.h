/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/varianceswapinstdata.h
// Purpose:     variance swap instdata class
// Created:     2006/02/21
// RCS-ID:      $Id: varianceswapinstdata.h,v 1.2 2006/04/10 11:38:55 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/varianceswapinstdata.h
    @brief Variance swap instdata class
 */

#ifndef _IHG_VARIANCESWAPINSTDATA_H_
#define _IHG_VARIANCESWAPINSTDATA_H_

#include "ihg/backwardinstdata.h"

namespace ito33
{

namespace pricing
{
  class VarianceSwapParams;
  class VarianceSwapMeshManager;
}

namespace ihg
{


class VarianceSwapInstData : public BackwardInstData
{

public:

  VarianceSwapInstData(pricing::VarianceSwapParams& params,
                       Model& model,
                       pricing::VarianceSwapMeshManager& meshes);

  virtual ~VarianceSwapInstData() { }


  // base class virtual functions
  virtual void Init();

  void UpdateBeforeStep();

  virtual void SetInitialValue();


protected:

  /// The variance swap params
  pricing::VarianceSwapParams& m_varianceSwapParams;
  
  /// The variance swap mesh manager
  pricing::VarianceSwapMeshManager& m_varianceSwapMeshes;


private:

  NO_COPY_CLASS(VarianceSwapInstData);

}; // class VarianceSwapInstData


} // namespace ihg

} // namespace ito33 

#endif // #ifndef _IHG_VARIANCESWAPINSTDATA_H_
