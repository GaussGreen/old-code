/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/cdspricer.h
// Purpose:     cds pricer class 
// Author:      Nabil
// Created:     2003/10/31
// RCS-ID:      $Id: cdspricer.h,v 1.10 2006/03/20 14:54:11 yann Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/cdspricer.h
    @brief cds pricer class

    Implementation of the Pricer and PricerRunner class for cds.
 */

#ifndef _IHG_CDSPRICER_H_
#define _IHG_CDSPRICER_H_

#include "ito33/autoptr.h"
#include "ito33/dlldecl.h"

#include "ito33/pricing/engine.h"

#include "ito33/pricing/cdsmeshmanager.h"

#include "ihg/cdsinstdata.h"
#include "ihg/cdsstepper.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL ComputationalFlags;
}

namespace ihg
{

  class Model;
  class CDSNumOutput;

/// cds pricer class
class CDSPricer
{
public:

  /// Constructor
  CDSPricer(pricing::CDSParams& params, 
            Model& model, 
            const finance::ComputationalFlags& flags);
  
  AutoPtr<CDSNumOutput> Price();


protected:

  pricing::CDSParams& m_params;
  Model& m_model;

  const finance::ComputationalFlags& m_flags;

  pricing::CDSMeshManager m_meshes;
  
  CDSInstData m_instdata;
  CDSStepper m_stepper;
  
  AutoPtr<CDSNumOutput> m_pNumOutput;

private:
  
  NO_COPY_CLASS(CDSPricer);

}; // class CDSPricer


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_CDSPRICER_H_

