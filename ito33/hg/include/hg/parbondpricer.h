/////////////////////////////////////////////////////////////////////////////
// Name:        hg/parbondpricer.h
// Purpose:     ParBond pricer class 
// Created:     2005/02/16
// RCS-ID:      $Id: parbondpricer.h,v 1.2 2006/03/20 14:54:09 yann Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/cdspricer.h
   @brief ParBond pricer class
 */

#ifndef _HG_PARBONDPRICER_H_
#define _HG_PARBONDPRICER_H_

#include "ito33/autoptr.h"
#include "ito33/dlldecl.h"

#include "ito33/pricing/engine.h"

#include "ito33/pricing/parbondmeshmanager.h"

#include "hg/parbondinstdata.h"
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

/// ParBond pricer class
class ParBondPricer
{
public:

  /// Constructor
  ParBondPricer(pricing::ParBondParams& params, 
                Model& model, 
                const finance::ComputationalFlags& flags);
  
  AutoPtr<NumOutputTimeOnly> Price();


protected:

  pricing::ParBondParams& m_params;
  Model& m_model;

  const finance::ComputationalFlags& m_flags;

  pricing::ParBondMeshManager m_meshes;
  
  ParBondInstData m_instdata;
  StepperTimeOnly m_stepper;
  
  AutoPtr<NumOutputTimeOnly> m_pNumOutput;

private:
  
  NO_COPY_CLASS(ParBondPricer);

}; // class ParBondPricer


} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_PARBONDPRICER_H_
