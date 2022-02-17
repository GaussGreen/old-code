/****************************************************************************
 * Name:        ihg/forwardoptionpricer.h
 * Purpose:     forward pricer class for European option
 * Author:      Wang
 * Created:     2004/03/10
 * RCS-ID:      $Id: forwardoptionpricer.h,v 1.8 2006/06/13 15:34:41 wang Exp $
 * Copyright:   (c) 2003-2003 Trilemma LLP
 ****************************************************************************/

#ifndef _IHG_FORWARDOPTIONPRICER_H_
#define _IHG_FORWARDOPTIONPRICER_H_

#include "ito33/autoptr.h" 

#include "ito33/pricing/engine.h"
#include "ito33/pricing/forwardoptionparams.h"
#include "ito33/pricing/forwardoptionmeshmanager.h"

#include "ihg/model.h"
#include "ihg/forwardoptioninstdata.h"
#include "ihg/forwardoptionstepper.h"
#include "ihg/forwardoptionnumoutput.h"

namespace ito33
{

namespace ihg
{

typedef pricing::Engine
       <
         pricing::ForwardOptionParams, 
         pricing::ForwardOptionMeshManager,
         ForwardOptionInstData, 
         ForwardOptionStepper, 
         ForwardOptionNumOutput
       > ForwardOptionEngine; 

/// European option Pricer using forward PDE 
class ForwardOptionPricer
{
public:

  /// Constructor
  ForwardOptionPricer(pricing::ForwardOptionParams& params, Model& model,
                      const finance::ComputationalFlags& flags)
    : m_params(params),
	    m_model(model), 
      m_meshes(m_params, m_model),
      m_instdata(m_params, model, m_meshes),
      m_stepper(m_instdata, flags),
      m_pNumOutput(new ForwardOptionNumOutput(params)),
      m_engine(m_params, m_meshes, m_instdata, m_stepper, *m_pNumOutput)
  {
    m_pNumOutput->GetComputationalFlags().SetComputeSurface
                                          (flags.GetComputeSurface()); 
  }
  
  AutoPtr<ForwardOptionNumOutput> Price()
  {
    m_params.Init();

    m_meshes.SetupMe();

    m_engine.Run();

    return m_pNumOutput;
  }


protected:

  pricing::ForwardOptionParams &m_params;

  Model &m_model;

  pricing::ForwardOptionMeshManager m_meshes;
  
  ForwardOptionInstData m_instdata;
  
  ForwardOptionStepper m_stepper;
  
  AutoPtr<ForwardOptionNumOutput> m_pNumOutput;
  
  ForwardOptionEngine m_engine;
  

private:

  NO_COPY_CLASS(ForwardOptionPricer);

}; // class ForwardOptionPricer


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_FORWARDOPTIONPRICER_H_

