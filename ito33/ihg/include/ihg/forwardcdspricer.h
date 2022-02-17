/****************************************************************************
 * Name:        ihg/forwardcdspricer.h
 * Purpose:     forward pricer class for CDS contracts
 * Author:      David
 * Created:     2004/03/29
 * RCS-ID:      $Id: forwardcdspricer.h,v 1.7 2006/06/13 15:34:41 wang Exp $
 * Copyright:   (c) 2004- Trilemma LLP
 ****************************************************************************/

#ifndef _IHG_FORWARDCDSPRICER_H_
#define _IHG_FORWARDCDSPRICER_H_

#include "ito33/autoptr.h" 

#include "ito33/pricing/engine.h"
#include "ito33/pricing/forwardcdsparams.h"
#include "ito33/pricing/forwardcdsmeshmanager.h"

#include "ihg/model.h"
#include "ihg/forwardcdsinstdata.h"
#include "ihg/forwardcdsstepper.h"
#include "ihg/forwardcdsnumoutput.h"

namespace ito33
{

namespace ihg
{

typedef pricing::Engine
       <
         pricing::ForwardCDSParams, 
         pricing::ForwardCDSMeshManager,
         ForwardCDSInstData, 
         ForwardCDSStepper, 
         ForwardCDSNumOutput
       > ForwardCDSEngine; 

/// CDS contract pricer using forward PDE 
class ForwardCDSPricer
{
public:

  /// Constructor
  ForwardCDSPricer(pricing::ForwardCDSParams& params, Model& model,
                   const finance::ComputationalFlags& flags)
    : m_params(params),
	    m_model(model), 
      m_meshes(m_params, m_model),
      m_instdata(m_params, model, m_meshes),
      m_stepper(m_instdata, flags),
      m_pNumOutput( new ForwardCDSNumOutput(params) ),
      m_engine(m_params, m_meshes, m_instdata, m_stepper, *m_pNumOutput)
  {
  }
  
  AutoPtr<ForwardCDSNumOutput> Price()
  {
    m_params.Init();

    m_meshes.SetupMe();

    m_engine.Run();

    return m_pNumOutput;
  }


protected:

  pricing::ForwardCDSParams& m_params;

  Model& m_model;

  pricing::ForwardCDSMeshManager m_meshes;
  
  ForwardCDSInstData m_instdata;
  
  ForwardCDSStepper m_stepper;
  
  AutoPtr<ForwardCDSNumOutput> m_pNumOutput;
  
  ForwardCDSEngine m_engine;
  

private:

  NO_COPY_CLASS(ForwardCDSPricer);

}; // class ForwardCDSPricer


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_FORWARDCDSPRICER_H_
