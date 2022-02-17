/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/parbondinstdata.h
// Purpose:     parbond instdata class
// Author:      Nabil, ZHANG
// Created:     2005/05/20
// RCS-ID:      $Id: parbondinstdata.h,v 1.1 2005/06/08 15:53:35 zhang Exp $
// Copyright:   (c) 2003-2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _IHG_PARBONDINSTDATA_H_
#define _IHG_PARBONDINSTDATA_H_

#include "ihg/backwardinstdata.h"

namespace ito33
{

namespace pricing
{
  class ParBondParams;
  class ParBondMeshManager;
}

namespace ihg
{

/// Inst data class for parbond in the ihg model
class ParBondInstData : public BackwardInstData
{
public:

  ParBondInstData(pricing::ParBondParams& params,
              Model& model,
              pricing::ParBondMeshManager &meshes);

  // Default dtor is ok
  
  void Init();

  void UpdateBeforeStep();

  // DoEvents is the same as in the base class BackwardInstData

  void SetInitialValue();


protected:
  
  pricing::ParBondParams& m_parbondParams;

  pricing::ParBondMeshManager& m_parbondMeshes;


private:

  NO_COPY_CLASS(ParBondInstData);

}; // class ParBondInstData


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_PARBONDINSTDATA_H_


