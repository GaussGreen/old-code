/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/parbondsteppertimeonly.h
// Purpose:     parbond stepper class
// Author:      ZHANG
// Created:     2005/05/20
// RCS-ID:      $Id: parbondsteppertimeonly.h,v 1.1 2005/06/08 15:53:35 zhang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ihg/parbondsteppertimeonly.h
   @brief time only parbond stepper class
 */

#ifndef _IHG_ParBondSTEPPERTIMEONLY_H_
#define _IHG_ParBondSTEPPERTIMEONLY_H_

#include "ihg/parbondinstdatatimeonly.h"

namespace ito33
{

namespace ihg
{


/**
   Time only stepper for parbond with HazardRateTimeOnly
 */
class ParBondStepperTimeOnly 
{
public:
  
  ParBondStepperTimeOnly(ParBondInstDataTimeOnly &instdata) 
    : m_instdata(instdata)
  {
  }

  ~ParBondStepperTimeOnly()
  {
  }

  /// Nothing to do when hazard rate is time only
  void Init() { };

  /**
     Make a single timestep
      
     It's more convient in this case to just call instdata to do the real work
   */
  void Run() { m_instdata.Run(); } 


private:
  
  ParBondInstDataTimeOnly &m_instdata;

  NO_COPY_CLASS(ParBondStepperTimeOnly);

}; // class ParBondStepperTimeOnly


} // namespace ihg

} // namespace ito33 

#endif // #ifndef _IHG_ParBondSTEPPERTIMEONLY_H_

