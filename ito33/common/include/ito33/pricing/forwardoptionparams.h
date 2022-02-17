/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/forwardoptionparams.h
// Purpose:     option params class for forward pricer
// Author:      Wang
// Created:     2004/03/04
// RCS-ID:      $Id: forwardoptionparams.h,v 1.6 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/forwardoptionparams.h
    @brief option params class for forward pricer

    Implementation of the params class for options.
 */

#ifndef _ITO33_PRICING_FORWARDOPTIONPARAMS_H_
#define _ITO33_PRICING_FORWARDOPTIONPARAMS_H_

#include "ito33/pricing/options.h"
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
   Option params class for forward pricer

   Pricing parameters for European options using forward equation
 */
class ForwardOptionParams : public Params
{

public: 

  ForwardOptionParams(Options& options) : Params(options),
                                          m_options(options)                                 
  { 
  }

  ForwardOptionParams(Options& options,
                      const finance::SessionData& sessionData,
                      const shared_ptr<numeric::NumParams>& pNumParams,
                      const shared_ptr<numeric::MeshParams>& pMeshParams)
    : Params(options, sessionData, pNumParams, pMeshParams),
      m_options(options)
  {
  }

  virtual ~ForwardOptionParams() { }

  virtual void Init();

  /// Get a reference to the options contract
  Options& GetOptions() const { return m_options; }

  /// Get the payoff pointer
  const finance::Payoff* GetPayoff() const { return m_pPayoff.get(); }


protected:

  void ConstructPayoff();

  /// Helper for initial values and constraints
  shared_ptr<finance::Payoff> m_pPayoff;

  Options& m_options;


private:

  NO_COPY_CLASS(ForwardOptionParams);

}; // class ForwardOptionParams


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_FORWARDOPTIONPARAMS_H_
