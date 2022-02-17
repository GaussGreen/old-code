/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/voltanhcalibrator.h
// Purpose:     calibrate tanh volatility
// Created:     2005/02/07
// RCS-ID:      $Id: voltanhcalibrator.h,v 1.1 2005/02/07 15:47:52 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file  ihg/voltanhcalibrator.h
   @brief calibrate tanh volatility

   For now, only newton iteration is implemented.
 */
#ifndef _IHG_VOLTANHCALIBRATOR_H_
#define _IHG_VOLTANHCALIBRATOR_H_

// Define one of the following macros to choose the method 
#define NEWTONITERATION


#ifdef NEWTONITERATION
#include "ihg/voltanhcalibrator_newton.h"
namespace ito33 { namespace ihg {

typedef VolTanhCalibratorNewton 
        VolTanhCalibrator;

}  } 
#endif // #ifdef NEWTONITERATION

#endif // #ifndef _IHG_VOLTANHCALIBRATOR_H_

