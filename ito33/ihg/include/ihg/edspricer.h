/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/edspricer.h
// Purpose:     EDS pricer class 
// Created:     2005/01/26
// RCS-ID:      $Id: edspricer.h,v 1.2 2006/03/20 14:54:11 yann Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ihg/edspricer.h
   @brief EDS pricer class

   @todo This is nearly the same with CDSPricer, we should find a way to share 
         the code between them. After all, This class doesn't need special 
         information on the specific params except during the construction of 
         the objets.
 */

#ifndef _IHG_EDSPRICER_H_
#define _IHG_EDSPRICER_H_

#include "ito33/autoptr.h"
#include "ito33/dlldecl.h"

#include "ito33/pricing/engine.h"

#include "ito33/pricing/edsmeshmanager.h"

#include "ihg/edsinstdata.h"
#include "ihg/edsstepper.h"

namespace ito33
{

namespace finance
{ 
  class ITO33_DLLDECL ComputationalFlags;
}

namespace ihg
{

  class Model;
  class EDSNumOutput;

/// EDS pricer class
class EDSPricer
{
public:

  /// Constructor
  EDSPricer(pricing::EDSParams& params, 
            Model& model, 
            const finance::ComputationalFlags& flags);
  
  AutoPtr<EDSNumOutput> Price();


protected:

  pricing::EDSParams& m_params;
  Model& m_model;

  const finance::ComputationalFlags& m_flags;

  pricing::EDSMeshManager m_meshes;
  
  EDSInstData m_instdata;
  EDSStepper m_stepper;
  
  AutoPtr<EDSNumOutput> m_pNumOutput;

private:
  
  NO_COPY_CLASS(EDSPricer);

}; // class EDSPricer


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_EDSPRICER_H_
