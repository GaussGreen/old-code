/////////////////////////////////////////////////////////////////////////////
// Name:        hg/steppertimeonly.h
// Purpose:     Dummy stepper class for time only problem using HG model
// Created:     2005/06/09
// RCS-ID:      $Id: steppertimeonly.h,v 1.1 2005/06/09 14:14:57 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/steppertimeonly.h
   @brief Dummy stepper class for time only problem using HG model
 */

#ifndef _HG_STEPPERTIMEONLY_H_
#define _HG_STEPPERTIMEONLY_H_

#include "hg/instdatatimeonly.h"

namespace ito33
{

namespace hg
{


/// CDS Stepper class for HG model, nearly empty because of time only
class StepperTimeOnly
{
public:
  
  StepperTimeOnly(InstDataTimeOnly& instdata) : m_instdata(instdata)
  {
  }

  ~StepperTimeOnly()
  {
  }

  /// Nothing to do here
  void Init() { };

  /**
     Make a single timestep
      
     It's more convient in this case to just call instdata to do the real work
   */
  void Run() { m_instdata.Run(); } 


private:
  
  InstDataTimeOnly& m_instdata;

  NO_COPY_CLASS(StepperTimeOnly);

}; // class StepperTimeOnly


} // namespace hg

} // namespace ito33 

#endif // #ifndef _HG_STEPPERTIMEONLY_H_
