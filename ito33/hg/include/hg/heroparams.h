/////////////////////////////////////////////////////////////////////////////
// Name:        hg/heroparams.h
// Purpose:     HERO params class
// Created:     2005/09/26
// RCS-ID:      $Id: heroparams.h,v 1.4 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/heroparams.h
    @brief HERO params class

    Implementation of the params class for HERO.
 */

#ifndef _HG_HEROPARAMS_H_
#define _HG_HEROPARAMS_H_

#include "ito33/pricing/params.h"

#include "hg/hedgingdata.h"

namespace ito33 
{
 
namespace finance
{
  class Payoff;
}

namespace hg
{


/**
   HERO params class

   Pricing parameters for HERO.
 */
class HeroParams : public pricing::Params
{

public: 

  /**
     Ctor for HERO params. 

     @param hedgingData reference to hedge data container (based on target)
   */
  HeroParams(HedgingData& hedgingData) : Params(hedgingData),
                                         m_hedgingData(hedgingData)
  { 
  }

  /**
     Ctor by HedgingData object and common financial and numerical datas. 

     @param hedgingData reference to target hero (based on target contract)
     @param sessionData reference to financial Session
     @param pNumParams numeric parameters
     @param pMeshParams parameters for mesh builder
   */
  HeroParams(HedgingData& hedgingData,
             const finance::SessionData& sessionData,
             const shared_ptr<numeric::NumParams>& pNumParams,
             const shared_ptr<numeric::MeshParams>& pMeshParams)
    : Params(hedgingData, sessionData, pNumParams, pMeshParams),
      m_hedgingData(hedgingData)
  {
  }

  virtual ~HeroParams() { }

  /// Initialize the object (payoff, dividends, etc)
  virtual void Init();

  /// Get the payoff pointer
  const finance::Payoff* GetPayoff() const 
  { 
    return m_pPayoff.get(); 
  }

  /// Get the hero object
  const HedgingData& GetHedgingData()
  {
    return m_hedgingData;
  }


protected:

  /// Helper for initial values and constraints
  shared_ptr<finance::Payoff> m_pPayoff;

  /// The HedgingData object containing the hedge surfaces, target, etc.
  HedgingData& m_hedgingData;


private:

  NO_COPY_CLASS(HeroParams);

}; // class HeroParams


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_HEROPARAMS_H_
