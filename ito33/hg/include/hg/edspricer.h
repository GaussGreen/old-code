/////////////////////////////////////////////////////////////////////////////
// Name:        hg/edspricer.h
// Purpose:     EDS pricer class 
// Created:     2005/01/31
// RCS-ID:      $Id: edspricer.h,v 1.7 2006/04/17 17:58:53 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/edspricer.h
   @brief EDS pricer class

   @todo Implement a base pricer, and derive EDSPricer from it.
 */

#ifndef _HG_EDSPRICER_H_
#define _HG_EDSPRICER_H_

#include "ito33/autoptr.h"
#include "ito33/dlldecl.h"

#include "ito33/pricing/engine.h"

#include "ito33/pricing/edsmeshmanager.h"

#include "hg/edsinstdata.h"

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

/// EDS pricer class
class EDSPricer
{
public:

  /// Constructor
  EDSPricer(pricing::EDSParams& params, 
            Model& model, 
            const finance::ComputationalFlags& flags);
  
  AutoPtr<BackwardNumOutput> Price();


protected:

  pricing::EDSParams& m_params;
  Model& m_model;

  const finance::ComputationalFlags& m_flags;

  pricing::EDSMeshManager m_meshes;
  
  EDSInstData m_instdata;
  
  AutoPtr<BackwardNumOutput> m_pNumOutput;

private:
  
  NO_COPY_CLASS(EDSPricer);

}; // class EDSPricer


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_EDSPRICER_H_
