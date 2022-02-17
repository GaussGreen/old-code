/////////////////////////////////////////////////////////////////////////////
// Name:        hg/onetouchinstdata.h
// Purpose:     OneTouch instdata class
// Created:     2005/07/04
// RCS-ID:      $Id: onetouchinstdata.h,v 1.1 2005/07/05 10:18:20 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _HG_ONETOUCHINSTDATA_H_
#define _HG_ONETOUCHINSTDATA_H_

#include "hg/backwardinstdata.h"

namespace ito33
{

namespace pricing
{
  class OneTouchParams;
  class OneTouchMeshManager;
}

namespace hg
{


/// Inst data class for EDS in the HG model
class OneTouchInstData : public BackwardInstData
{
public:

  OneTouchInstData(pricing::OneTouchParams& params,
                   Model& model,
                   pricing::OneTouchMeshManager& meshes);

  // Default dtor is ok

  void SetInitialValue();
  
  void UpdateBeforeStep();

  void Init();


protected:
  
  pricing::OneTouchParams& m_oneTouchParams;

  pricing::OneTouchMeshManager& m_oneTouchMeshes;


private:

  NO_COPY_CLASS(OneTouchInstData);

}; // class OneTouchInstData


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_ONETOUCHINSTDATA_H_
