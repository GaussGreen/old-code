/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/cdsinstdata.h
// Purpose:     cds instdata class
// Author:      Nabil, Wang
// Created:     2003/10/29
// RCS-ID:      $Id: cdsinstdata.h,v 1.9 2005/02/01 12:16:24 nabil Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _IHG_CDSINSTDATA_H_
#define _IHG_CDSINSTDATA_H_

#include "ihg/backwardinstdata.h"

namespace ito33
{

namespace pricing
{
  class CDSParams;
  class CDSMeshManager;
}

namespace ihg
{

/// Inst data class for cds in the ihg model
class CDSInstData : public BackwardInstData
{
public:

  CDSInstData(pricing::CDSParams& params,
              Model& model,
              pricing::CDSMeshManager &meshes);

  // Default dtor is ok
  
  void Init();

  void UpdateBeforeStep();

  // DoEvents is the same as in the base class BackwardInstData

  void SetInitialValue();


protected:
  
  pricing::CDSParams& m_cdsParams;

  pricing::CDSMeshManager& m_cdsMeshes;


private:

  NO_COPY_CLASS(CDSInstData);

}; // class CDSInstData


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_CDSINSTDATA_H_


