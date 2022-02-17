/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/parbondpricertimeonly.h
// Purpose:     time only parbond pricer class 
// Author:      ZHANG
// Created:     2005/05/20
// RCS-ID:      $Id: parbondpricertimeonly.h,v 1.2 2006/03/20 14:54:11 yann Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/parbondpricertimeonly.h
    @brief time only parbond pricer class
 */

#ifndef _IHG_ParBondPRICERTIMEONLY_H_
#define _IHG_ParBondPRICERTIMEONLY_H_

#include "ito33/autoptr.h"

#include "ito33/pricing/engine.h"

#include "ito33/pricing/parbondmeshmanager.h"

#include "ihg/parbondinstdatatimeonly.h"
#include "ihg/parbondsteppertimeonly.h"


namespace ito33
{

namespace ihg
{

  class ParBondNumOutputTimeOnly;

typedef pricing::Engine
        <
          pricing::ParBondParams,
          pricing::ParBondMeshManager,
          ParBondInstDataTimeOnly, 
          ParBondStepperTimeOnly, 
          ParBondNumOutputTimeOnly
        > ParBondEngineTimeOnly; 


/// parbond pricer class
class ParBondPricerTimeOnly
{
public:

  /// Constructor
  ParBondPricerTimeOnly(pricing::ParBondParams& params, 
                    Model& model,
                    const finance::ComputationalFlags& flags);
  
  AutoPtr<ParBondNumOutputTimeOnly> Price();


protected:

  pricing::ParBondParams& m_params;
  Model& m_model;

  const finance::ComputationalFlags& m_flags;

  pricing::ParBondMeshManager m_meshes;
  ParBondInstDataTimeOnly m_instdata;
  ParBondStepperTimeOnly m_stepper;
  
  AutoPtr<ParBondNumOutputTimeOnly> m_pNumOutput;
  
  ParBondEngineTimeOnly m_engine;


private:
  
  NO_COPY_CLASS(ParBondPricerTimeOnly);

}; // class ParBondPricerTimeOnly


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_ParBondPRICERTIMEONLY_H_

