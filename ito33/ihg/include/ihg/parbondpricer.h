/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/parbondpricer.h
// Purpose:     parbond pricer class 
// Author:      Nabil
// Created:     2005/05/20
// RCS-ID:      $Id: parbondpricer.h,v 1.2 2006/03/20 14:54:11 yann Exp $
// Copyright:   (c) 2003-2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/parbondpricer.h
    @brief parbond pricer class

    Implementation of the Pricer and PricerRunner class for parbond.
 */

#ifndef _IHG_ParBondPRICER_H_
#define _IHG_ParBondPRICER_H_

#include "ito33/autoptr.h"
#include "ito33/dlldecl.h"

#include "ito33/pricing/engine.h"

#include "ito33/pricing/parbondmeshmanager.h"

#include "ihg/parbondinstdata.h"
#include "ihg/parbondstepper.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL ComputationalFlags;
}

namespace ihg
{

  class Model;
  class ParBondNumOutput;

/// parbond pricer class
class ParBondPricer
{
public:

  /// Constructor
  ParBondPricer(pricing::ParBondParams& params, 
            Model& model, 
            const finance::ComputationalFlags& flags);
  
  AutoPtr<ParBondNumOutput> Price();


protected:

  pricing::ParBondParams& m_params;
  Model& m_model;

  const finance::ComputationalFlags& m_flags;

  pricing::ParBondMeshManager m_meshes;
  
  ParBondInstData m_instdata;
  ParBondStepper m_stepper;
  
  AutoPtr<ParBondNumOutput> m_pNumOutput;

private:
  
  NO_COPY_CLASS(ParBondPricer);

}; // class ParBondPricer


} // namespace ihg

} // namespace ito33

#endif // #ifndef _IHG_ParBondPRICER_H_

