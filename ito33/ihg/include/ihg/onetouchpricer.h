/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/onetouchpricer.h
// Purpose:     OneTouch pricer class 
// Created:     2006/08/11
// RCS-ID:      $Id: onetouchpricer.h,v 1.1 2006/08/10 23:10:42 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ihg/onetouchpricer.h
   @brief OneTouch pricer class

   @todo This is nearly the same with EDSPricer, we should find a way to share 
         the code between them. After all, This class doesn't need special 
         information on the specific params except during the construction of 
         the objets. See HG.
 */

#ifndef _IHG_ONETOUCHPRICER_H_
#define _IHG_ONETOUCHPRICER_H_

#include "ito33/autoptr.h"
#include "ito33/dlldecl.h"

#include "ito33/pricing/engine.h"

#include "ito33/pricing/onetouchmeshmanager.h"

#include "ihg/onetouchinstdata.h"
#include "ihg/onetouchstepper.h"

namespace ito33
{

namespace finance
{ 
  class ITO33_DLLDECL ComputationalFlags;
}

namespace ihg
{

  class Model;
  class OneTouchNumOutput;

/// OneTouch pricer class
class OneTouchPricer
{
public:

  /// Constructor
  OneTouchPricer(pricing::OneTouchParams& params, 
                 Model& model, 
                 const finance::ComputationalFlags& flags);
  
  AutoPtr<OneTouchNumOutput> Price();


protected:

  pricing::OneTouchParams& m_params;
  Model& m_model;

  const finance::ComputationalFlags& m_flags;

  pricing::OneTouchMeshManager m_meshes;
  
  OneTouchInstData m_instdata;
  OneTouchStepper m_stepper;
  
  AutoPtr<OneTouchNumOutput> m_pNumOutput;

private:
  
  NO_COPY_CLASS(OneTouchPricer);

}; // class OneTouchPricer


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_ONETOUCHPRICER_H_
