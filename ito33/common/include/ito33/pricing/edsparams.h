/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/edsparams.h
// Purpose:     EDS params class
// Created:     2005/01/26
// RCS-ID:      $Id: edsparams.h,v 1.3 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/edsparams.h
    @brief Implementation of the params class for EDS.
 */

#ifndef _ITO33_PRICING_EDSPARAMS_H_
#define _ITO33_PRICING_EDSPARAMS_H_

#include "ito33/pricing/params.h"
#include "ito33/pricing/eds.h"

namespace ito33 
{
 
namespace pricing
{


/**
   Option params class

   Pricing parameters for European and American options
 */
class EDSParams : public Params
{

public: 

  EDSParams(EDS& eds) : Params(eds), m_eds(eds)                                 
  { 
  }

  EDSParams(EDS& eds,
            const finance::SessionData& sessionData,
            const shared_ptr<numeric::NumParams>& pNumParams,
            const shared_ptr<numeric::MeshParams>& pMeshParams)
    : Params(eds, sessionData, pNumParams, pMeshParams),
      m_eds(eds)
  {
  }

  // Default dtor is ok

  void Init();

  /// Get a reference to the EDS contract
  EDS& GetEDS() const { return m_eds; }


private:
  
  EDS& m_eds; 

  NO_COPY_CLASS(EDSParams);

}; // class EDSParams


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_EDSPARAMS_H_
