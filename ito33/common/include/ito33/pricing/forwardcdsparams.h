/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/forwardcdsparams.h
// Purpose:     cds params class for forward pricer
// Author:      David
// Created:     2004/03/31
// RCS-ID:      $Id: forwardcdsparams.h,v 1.6 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/forwardcdsparams.h
    @brief cds params class for forward pricer

    Implementation of the params class for forward cds pricing.
 */

#ifndef _ITO33_PRICING_FORWARDCDSPARAMS_H_
#define _ITO33_PRICING_FORWARDCDSPARAMS_H_

#include "ito33/common.h"
#include "ito33/list.h"

#include "ito33/pricing/cdses.h"
#include "ito33/pricing/params.h"

namespace ito33 
{
 
namespace finance
{
  class Payoff;

  class ITO33_DLLDECL SessionData;
}

namespace pricing
{


/**
    CDS params class for forward pricer.

    Pricing parameters for cds contracts using forward equation.
 */
class ForwardCDSParams : public Params
{

public: 

  ForwardCDSParams(CDSes& cdses) : Params(cdses), m_cdses(cdses)                                 
  { 
  }

  ForwardCDSParams(CDSes& cdses,
                   const finance::SessionData& sessionData,
                   const shared_ptr<numeric::NumParams>& pNumParams,
                   const shared_ptr<numeric::MeshParams>& pMeshParams)
    : Params(cdses, sessionData, pNumParams, pMeshParams),
      m_cdses(cdses)
  {
  }

  virtual void Init();

  /// Get a reference to the cds contracts
  CDSes& GetCDSes() const { return m_cdses; }

  /// Get the payoff pointer
  const finance::Payoff *GetPayoff() const { return m_pPayoff.get(); }


protected:

  void ConstructPayoff();

  /// Helper for initial values and constraints
  shared_ptr<finance::Payoff> m_pPayoff;

  CDSes& m_cdses;


private:

  NO_COPY_CLASS(ForwardCDSParams);

}; // class ForwardCDSParams


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_FORWARDCDSPARAMS_H_
