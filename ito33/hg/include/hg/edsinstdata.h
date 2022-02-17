/////////////////////////////////////////////////////////////////////////////
// Name:        hg/edsinstdata.h
// Purpose:     EDS instdata class
// Created:     2005/01/31
// RCS-ID:      $Id: edsinstdata.h,v 1.3 2005/04/21 15:21:15 dave Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _HG_EDSINSTDATA_H_
#define _HG_EDSINSTDATA_H_

#include "hg/backwardinstdata.h"

namespace ito33
{

namespace pricing
{
  class EDSParams;
  class EDSMeshManager;
}

namespace hg
{

/// Inst data class for EDS in the HG model
class EDSInstData : public BackwardInstData
{
public:

  EDSInstData(pricing::EDSParams& params,
              Model& model,
              pricing::EDSMeshManager& meshes);

  // Default dtor is ok

  void SetInitialValue();
  
  void UpdateBeforeStep();

  void Init();

protected:
  
  pricing::EDSParams& m_edsParams;

  pricing::EDSMeshManager& m_edsMeshes;


private:

  NO_COPY_CLASS(EDSInstData);

}; // class EDSInstData


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_EDSINSTDATA_H_


