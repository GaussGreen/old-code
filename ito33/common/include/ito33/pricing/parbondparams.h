/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/parbondparams.h
// Purpose:     parbond params class
// Author:      ZHANG
// Created:     2005/05/20
// RCS-ID:      $Id: parbondparams.h,v 1.5 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/parbondparams.h
    @brief parbond params class

    Implementation of the params class for parbond.
 */

#ifndef _ITO33_PRICING_ParBondPARAMS_H_
#define _ITO33_PRICING_ParBondPARAMS_H_

#include "ito33/autoptr.h"
#include "ito33/list.h"

#include "ito33/pricing/params.h"
#include "ito33/pricing/parbond.h"

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
class ParBondParams : public Params
{

public: 

  /**
     ParBondParams constructor.

     @param parbond pricing level parbond object
   */
  ParBondParams(ParBond& parbond) : Params(parbond), m_parbond(parbond)
  { 
  }

  /**
     Constructor.

     @param parbond pricing parbond object
     @param sessionData finance session data
     @param pNumParams numerical parameters
     @param pMeshParams mesh parameters
   */
  ParBondParams(ParBond& parbond,
            const finance::SessionData& sessionData,
            const shared_ptr<numeric::NumParams>& pNumParams,
            const shared_ptr<numeric::MeshParams>& pMeshParams)
    : Params(parbond, sessionData, pNumParams, pMeshParams),
      m_parbond(parbond)
  {
  }

  // Default dtor is ok

  virtual void Init();

  /// Get a reference to the parbond contract
  ParBond& GetParBond() const { return m_parbond; }


private:

  ParBond& m_parbond;

  NO_COPY_CLASS(ParBondParams);

}; // class ParBondParams


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_ParBondPARAMS_H_
