/////////////////////////////////////////////////////////////////////////////
// Name:        hg/cdspricer.h
// Purpose:     CDS pricer class 
// Created:     2005/02/16
// RCS-ID:      $Id: cdspricer.h,v 1.4 2006/03/20 14:54:09 yann Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/cdspricer.h
   @brief CDS pricer class
 */

#ifndef _HG_CDSPRICER_H_
#define _HG_CDSPRICER_H_

#include "ito33/autoptr.h"
#include "ito33/dlldecl.h"

#include "ito33/pricing/engine.h"

#include "ito33/pricing/cdsmeshmanager.h"

#include "hg/cdsinstdata.h"
#include "hg/steppertimeonly.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL ComputationalFlags;
}

namespace hg
{

  class Model;
  class NumOutputTimeOnly;

/// CDS pricer class
class CDSPricer
{
public:

  /// Constructor
  CDSPricer(pricing::CDSParams& params, 
            Model& model, 
            const finance::ComputationalFlags& flags);
  
  AutoPtr<NumOutputTimeOnly> Price();


protected:

  pricing::CDSParams& m_params;
  Model& m_model;

  const finance::ComputationalFlags& m_flags;

  pricing::CDSMeshManager m_meshes;
  
  CDSInstData m_instdata;
  StepperTimeOnly m_stepper;
  
  AutoPtr<NumOutputTimeOnly> m_pNumOutput;

private:
  
  NO_COPY_CLASS(CDSPricer);

}; // class CDSPricer


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_CDSPRICER_H_
