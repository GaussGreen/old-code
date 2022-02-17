/////////////////////////////////////////////////////////////////////////////
// Name:        hg/onetouchpricer.h
// Purpose:     OneTouch pricer class 
// Created:     2005/07/04
// RCS-ID:      $Id: onetouchpricer.h,v 1.3 2006/04/17 17:58:53 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/onetouchpricer.h
   @brief OneTouch pricer class
 */

#ifndef _HG_ONETOUCHPRICER_H_
#define _HG_ONETOUCHPRICER_H_

#include "ito33/autoptr.h"

#include "ito33/pricing/onetouchmeshmanager.h"

#include "hg/onetouchinstdata.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL ComputationalFlags;
}

namespace hg
{

  class Model;
  class BackwardNumOutput;

/// OneTouch pricer class
class OneTouchPricer
{
public:

  /// Constructor
  OneTouchPricer(pricing::OneTouchParams& params, 
                 Model& model, 
                 const finance::ComputationalFlags& flags);
  
  AutoPtr<BackwardNumOutput> Price();


protected:

  pricing::OneTouchParams& m_params;
  Model& m_model;

  const finance::ComputationalFlags& m_flags;

  pricing::OneTouchMeshManager m_meshes;
  
  OneTouchInstData m_instdata;
  
  AutoPtr<BackwardNumOutput> m_pNumOutput;


private:
  
  NO_COPY_CLASS(OneTouchPricer);

}; // class OneTouchPricer


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_ONETOUCHPRICER_H_
