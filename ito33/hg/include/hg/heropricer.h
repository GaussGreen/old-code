/////////////////////////////////////////////////////////////////////////////
// Name:        hg/heropricer.h
// Purpose:     implementation of pricer class for HERO
// Created:     2005/09/26
// RCS-ID:      $Id: heropricer.h,v 1.3 2006/03/20 14:54:09 yann Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/heropricer.h
   @brief implementation of pricer class for HERO
 */

#ifndef _HG_HEROPRICER_H_ 
#define _HG_HEROPRICER_H_

#include "ito33/autoptr.h"
#include "ito33/common.h"
#include "ito33/dlldecl.h"

#include "hg/heromeshmanager.h"
#include "hg/heroinstdata.h"
#include "hg/herostepper.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL ComputationalFlags;
}

namespace hg
{
  class Model;
  class HeroNumOutput;
  class HeroParams;

/**
   HERO Pricer class

   This class simply constructs the relevant classes for computing the HERO,
   and calls engine to do the computation.
 */
class HeroPricer
{
public:

  /// Construct from option parameters
  HeroPricer(HeroParams& params, 
             Model& model,
             const finance::ComputationalFlags& flags);
  
  /** 
      Compute the HERO
      
      @return HeroNumOutput object containing all price data
   */
  AutoPtr<HeroNumOutput> Price();


protected:

  HeroParams& m_params;

  Model& m_model;

  const finance::ComputationalFlags& m_flags;

  HeroMeshManager m_meshes;
  
  HeroInstData m_instdata;
  
  HeroStepper m_stepper;
  
  AutoPtr<HeroNumOutput> m_pNumOutput;
  

private:

  NO_COPY_CLASS(HeroPricer);

}; // class HeroPricer


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_HEROPRICER_H_

