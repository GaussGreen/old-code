/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/onetouchinstdata.h
// Purpose:     OneTouch instdata class
// Created:     2005/01/26
// RCS-ID:      $Id: onetouchinstdata.h,v 1.1 2006/08/10 23:10:42 wang Exp $
// Copyright:   (c) 2006  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _IHG_ONETOUCHINSTDATA_H_
#define _IHG_ONETOUCHINSTDATA_H_

#include "ihg/backwardinstdata.h"

namespace ito33
{

namespace pricing
{
  class OneTouchParams;
  class OneTouchMeshManager;
}

namespace ihg
{

/// Inst data class for OneTouch in the IHG model
class OneTouchInstData : public BackwardInstData
{
public:

  OneTouchInstData(pricing::OneTouchParams& params,
                   Model& model,
                   pricing::OneTouchMeshManager& meshes);

  // Default dtor is ok
  
  void Init();

  void SetInitialValue();
  
  void UpdateBeforeStep();

  // DoEvents is the same as in the base class BackwardInstData

protected:
  
  pricing::OneTouchParams& m_oneTouchParams;

  pricing::OneTouchMeshManager& m_oneTouchMeshes;


private:

  NO_COPY_CLASS(OneTouchInstData);

}; // class OneTouchInstData


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_ONETOUCHINSTDATA_H_
