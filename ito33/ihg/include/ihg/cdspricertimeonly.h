/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/cdspricertimeonly.h
// Purpose:     time only cds pricer class 
// Created:     2004/03/19
// RCS-ID:      $Id: cdspricertimeonly.h,v 1.11 2006/08/21 14:55:35 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/cdspricertimeonly.h
    @brief time only cds pricer class
 */

#ifndef _IHG_CDSPRICERTIMEONLY_H_
#define _IHG_CDSPRICERTIMEONLY_H_

#include "ito33/autoptr.h"

#include "ito33/pricing/engine.h"

#include "ito33/pricing/cdsmeshmanager.h"

#include "ihg/cdsinstdatatimeonly.h"
#include "ihg/cdssteppertimeonly.h"


namespace ito33
{

namespace ihg
{

  class CDSNumOutputTimeOnly;

typedef pricing::Engine
        <
          pricing::CDSParams,
          pricing::CDSMeshManager,
          CDSInstDataTimeOnly, 
          CDSStepperTimeOnly, 
          CDSNumOutputTimeOnly
        > CDSEngineTimeOnly; 


/// cds pricer class
class CDSPricerTimeOnly
{
public:

  /// Constructor
  CDSPricerTimeOnly(pricing::CDSParams& params, 
                    Model& model,
                    const finance::ComputationalFlags& flags);
  
  AutoPtr<CDSNumOutputTimeOnly> Price();


protected:

  pricing::CDSParams& m_params;
  Model& m_model;

  const finance::ComputationalFlags& m_flags;

  pricing::CDSMeshManager m_meshes;
  CDSInstDataTimeOnly m_instdata;
  CDSStepperTimeOnly m_stepper;
  
  AutoPtr<CDSNumOutputTimeOnly> m_pNumOutput;
  
  CDSEngineTimeOnly m_engine;


private:
  
  NO_COPY_CLASS(CDSPricerTimeOnly);

}; // class CDSPricerTimeOnly


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_CDSPRICERTIMEONLY_H_

