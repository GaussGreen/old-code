/////////////////////////////////////////////////////////////////////////////
// Name:        hg/cdsinstdata.h
// Purpose:     CDS instdata class
// Created:     2005/02/16
// RCS-ID:      $Id: cdsinstdata.h,v 1.8 2005/06/09 14:14:24 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/cdsinstdata.h
 */

#ifndef _HG_CDSINSTDATA_H_
#define _HG_CDSINSTDATA_H_

#include "hg/instdatatimeonly.h"

namespace ito33
{

namespace pricing
{
  class CDSParams;
  class CDSMeshManager;
}

namespace hg
{


/// Inst data class for CDS in the HG model
class CDSInstData : public InstDataTimeOnly
{
public:

  CDSInstData(pricing::CDSParams& params,
              Model& model,
              pricing::CDSMeshManager& meshes);

  // Default dtor is ok

  virtual void UpdateRecoveryValue(); 


protected:

  pricing::CDSParams& m_cdsParams;

  pricing::CDSMeshManager& m_cdsMeshes;


private:
  
  NO_COPY_CLASS(CDSInstData);

}; // class CDSInstData


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_CDSINSTDATA_H_
