/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/stepper_fix.h
// Purpose:     stepper class for fix mesh
// Author:      ZHANG Yunzhi
// Created:     2004/01/09
// RCS-ID:      $Id: stepper_fix.h,v 1.4 2004/10/04 18:04:04 pedro Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
  @file ihg/stepper_fix.h
  @brief stepper class for fix mesh

*/

#ifndef _IHG_STEPPER_FIX_H_
#define _IHG_STEPPER_FIX_H_

#include "ito33/sharedptr.h"

#include "ihg/params.h"
#include "ihg/instdata_fix.h"
#include "ihg/stepper.h"

namespace ito33
{

namespace ihg
{

class StepperFix : public ihg::Stepper
{
public:
  StepperFix(Params &params, 
             InstDataFix &instdata) 
    : Stepper(params, instdata),
      m_params(params),
      m_instdata(instdata)
  {
  }

  virtual ~StepperFix() { }

  /**
    Initialization from m_meshes and m_params

    The derived class must still define the computational grid m_pdX, and
    initialize any other contract specific data. 

    REQUIRE : m_meshes must have been done
  */
  virtual void Init()
  {
    m_nNbX = m_instdata.m_nNbS;

    Alloc(m_nNbX);

    CalculateAreaArrays(m_instdata.m_pdS, m_nNbX);
  }


  virtual void Run() = 0;

protected:

  /// In this class, we need the type information of the params
  Params &m_params;

  InstDataFix &m_instdata;

private:

  NO_COPY_CLASS(StepperFix);
};

} // namespace ihg

} // namespace ito33 


#endif // #ifndef _IHG_STEPPER_FIX_H_

