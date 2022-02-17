/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/edsinstdata.h
// Purpose:     EDS instdata class
// Created:     2005/01/26
// RCS-ID:      $Id: edsinstdata.h,v 1.2 2005/02/04 15:11:07 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _IHG_EDSINSTDATA_H_
#define _IHG_EDSINSTDATA_H_

#include "ihg/backwardinstdata.h"

namespace ito33
{

namespace pricing
{
  class EDSParams;
  class EDSMeshManager;
}

namespace ihg
{

/// Inst data class for EDS in the ihg model
class EDSInstData : public BackwardInstData
{
public:

  EDSInstData(pricing::EDSParams& params,
              Model& model,
              pricing::EDSMeshManager& meshes);

  // Default dtor is ok
  
  void Init();

  void UpdateBeforeStep();

  // DoEvents is the same as in the base class BackwardInstData

  void SetInitialValue();


protected:
  
  pricing::EDSParams& m_edsParams;

  pricing::EDSMeshManager& m_edsMeshes;


private:

  NO_COPY_CLASS(EDSInstData);

}; // class EDSInstData


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_EDSINSTDATA_H_


