/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/cdsparams.h
// Purpose:     cds params class
// Author:      Wang
// Created:     2004/03/02
// RCS-ID:      $Id: cdsparams.h,v 1.7 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/cdsparams.h
    @brief cds params class

    Implementation of the params class for cds.
 */

#ifndef _ITO33_PRICING_CDSPARAMS_H_
#define _ITO33_PRICING_CDSPARAMS_H_

#include "ito33/autoptr.h"
#include "ito33/list.h"

#include "ito33/pricing/params.h"
#include "ito33/pricing/cds.h"

namespace ito33 
{
 
namespace finance
{
  class ITO33_DLLDECL SessionData;
}

namespace pricing
{


/**
   Option params class

   Pricing parameters for European and American options
 */
class CDSParams : public Params
{

public: 

  CDSParams(CDS& cds) : Params(cds), m_cds(cds)                                 
  { 
  }

  CDSParams(CDS& cds,
            const finance::SessionData& sessionData,
            const shared_ptr<numeric::NumParams>& pNumParams,
            const shared_ptr<numeric::MeshParams>& pMeshParams)
    : Params(cds, sessionData, pNumParams, pMeshParams),
      m_cds(cds)
  {
  }

  // Default dtor is ok

  virtual void Init();

  /// Get a reference to the cds contract
  CDS& GetCDS() const { return m_cds; }


protected:

  CDS& m_cds;


private:

  NO_COPY_CLASS(CDSParams);

}; // class CDSParams


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_CDSPARAMS_H_
