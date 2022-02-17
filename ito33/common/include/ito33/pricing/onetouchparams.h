/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/onetouchparams.h
// Purpose:     OneTouch params class
// Created:     2005/07/04
// RCS-ID:      $Id: onetouchparams.h,v 1.3 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/onetouchparams.h
    @brief Implementation of the params class for OneTouch.
 */

#ifndef _ITO33_PRICING_ONETOUCHPARAMS_H_
#define _ITO33_PRICING_ONETOUCHPARAMS_H_

#include "ito33/pricing/params.h"
#include "ito33/pricing/onetouch.h"

namespace ito33 
{
 
namespace pricing
{


/// OneTouch params class
class OneTouchParams : public Params
{

public: 

  OneTouchParams(OneTouch& oneTouch,
                 const finance::SessionData& sessionData,
                 const shared_ptr<numeric::NumParams>& pNumParams,
                 const shared_ptr<numeric::MeshParams>& pMeshParams)
               : Params(oneTouch, sessionData, pNumParams, pMeshParams),
                 m_oneTouch(oneTouch)
  {
  }

  // Default dtor is ok

  void Init();

  /// Get a reference to the contract
  OneTouch& GetOneTouch() const { return m_oneTouch; }


private:
  
  OneTouch& m_oneTouch; 

  NO_COPY_CLASS(OneTouchParams);

}; // class OneTouchParams


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_ONETOUCHPARAMS_H_
