/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/cdssteppertimeonly.h
// Purpose:     cds stepper class
// Author:      Wang
// Created:     2004/03/18
// RCS-ID:      $Id: cdssteppertimeonly.h,v 1.4 2004/10/04 18:04:04 pedro Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ihg/cdssteppertimeonly.h
   @brief time only cds stepper class
 */

#ifndef _IHG_CDSSTEPPERTIMEONLY_H_
#define _IHG_CDSSTEPPERTIMEONLY_H_

#include "ihg/cdsinstdatatimeonly.h"

namespace ito33
{

namespace ihg
{


/**
   Time only stepper for cds with HazardRateTimeOnly
 */
class CDSStepperTimeOnly 
{
public:
  
  CDSStepperTimeOnly(CDSInstDataTimeOnly &instdata) 
    : m_instdata(instdata)
  {
  }

  ~CDSStepperTimeOnly()
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
  
  CDSInstDataTimeOnly &m_instdata;

  NO_COPY_CLASS(CDSStepperTimeOnly);

}; // class CDSStepperTimeOnly


} // namespace ihg

} // namespace ito33 

#endif // #ifndef _IHG_CDSSTEPPERTIMEONLY_H_

