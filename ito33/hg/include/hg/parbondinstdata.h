/////////////////////////////////////////////////////////////////////////////
// Name:        hg/parbondinstdata.h
// Purpose:     ParBond instdata class
// Created:     2005/06/09
// RCS-ID:      $Id: parbondinstdata.h,v 1.1 2005/06/09 15:36:37 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/parbondinstdata.h
 */

#ifndef _HG_PARBONDINSTDATA_H_
#define _HG_PARBONDINSTDATA_H_

#include "hg/instdatatimeonly.h"

namespace ito33
{

namespace pricing
{
  class ParBondParams;
  class ParBondMeshManager;
}

namespace hg
{


/// Inst data class for ParBond in the HG model
class ParBondInstData : public InstDataTimeOnly
{
public:

  ParBondInstData(pricing::ParBondParams& params,
                  Model& model,
                  pricing::ParBondMeshManager& meshes);

  // Default dtor is ok

  virtual void UpdateRecoveryValue(); 


protected:

  pricing::ParBondParams& m_parBondParams;

  pricing::ParBondMeshManager& m_parBondMeshes;


private:
  
  NO_COPY_CLASS(ParBondInstData);

}; // class ParBondInstData


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_PARBONDINSTDATA_H_
