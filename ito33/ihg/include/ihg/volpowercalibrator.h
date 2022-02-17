/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/volpowercalibrator.h
// Purpose:     calibrate power volatility
// Author:      Ito33
// Created:     2004/11/23
// RCS-ID:      $Id: volpowercalibrator.h,v 1.2 2005/01/28 19:47:16 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file  ihg/volpowercalibrator.h
   @brief calibrate power volatility
 */
#ifndef _IHG_VOLPOWERCALIBRATOR_H_
#define _IHG_VOLPOWERCALIBRATOR_H_

// Define one of the following macros to choose the method 
#define NEWTONITERATION


#ifdef NEWTONITERATION
#include "ihg/volpowercalibrator_newton.h"
namespace ito33 { namespace ihg {

typedef VolPowerCalibratorNewton 
        VolPowerCalibrator;

}  } 
#endif // #ifdef NEWTONITERATION

#endif // #ifndef _IHG_VOLPOWERCALIBRATOR_H_

