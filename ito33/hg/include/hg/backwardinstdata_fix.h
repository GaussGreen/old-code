/////////////////////////////////////////////////////////////////////////////
// Name:        hg/backwardinstdata_fix.h
// Purpose:     instdata class for backward problem with fixed space mesh
// Created:     2005/01/31
// RCS-ID:      $Id: backwardinstdata_fix.h,v 1.1 2005/01/31 17:00:10 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _HG_BACKWARDINSTDATA_FIX_H_
#define _HG_BACKWARDINSTDATA_FIX_H_

#include "hg/backwardinstdata.h"

namespace ito33
{

namespace pricing
{
  class Params;
  class BackwardMeshManagerFix;
}

namespace hg
{

/// Inst data class for backward pricing with fixed space mesh in the HG model
class BackwardInstDataFix : public BackwardInstData
{
public:

  BackwardInstDataFix(pricing::Params& params,
                      Model& model,
                      pricing::BackwardMeshManagerFix& meshes);

  // Default dtor is ok
  
  void Init();

  // DoEvents is the same as in the base class BackwardInstData

  /// Default value at current time step
  double m_dRecoveryValue;


protected:

  pricing::BackwardMeshManagerFix& m_fixedMeshes;


private:

  NO_COPY_CLASS(BackwardInstDataFix);

}; // class BackwardInstDataFix


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_BACKWARDINSTDATA_FIX_H_


