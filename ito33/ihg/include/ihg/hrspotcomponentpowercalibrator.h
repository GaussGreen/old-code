/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/hrspotcomponentpower.h
// Purpose:     calibrate spot component power of HRSpotComponentPower
// Author:      Ito33
// Created:     2004/11/23
// RCS-ID:      $Id: hrspotcomponentpowercalibrator.h,v 1.3 2006/01/27 15:57:29 dave Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file  ihg/hrspotcomponentpower.h
   @brief calibrate spot component power of HRSpotComponentPower
 */
#ifndef _IHG_HRSPOTCOMPONENTPOWERCALIBRATOR_H_
#define _IHG_HRSPOTCOMPONENTPOWERCALIBRATOR_H_

// Define one of the following macros to choose the method 

//#define BETAITERATION
//#define NAGITERATION
//#define NEWTONITERATION
#define MIXEDITERATION

#ifdef BETAITERATION
#include "ihg/hrspotcomponentpowercalibrator_beta.h"
namespace ito33 { namespace ihg {

typedef HazardRateSpotComponentPowerCalibratorBeta 
        HazardRateSpotComponentPowerCalibrator;

}  } 
#endif // #ifdef BETAITERATION

#ifdef NAGITERATION
#include "ihg/hrspotcomponentpowercalibrator_nag.h"

#ifdef _MSC_VER
  #pragma comment(lib, "nagc.lib")
#endif

namespace ito33 { namespace ihg {

typedef HazardRateSpotComponentPowerCalibratorNAG
        HazardRateSpotComponentPowerCalibrator;

}  } 
#endif // #ifdef NAGITERATION

#ifdef NEWTONITERATION
#include "ihg/hrspotcomponentpowercalibrator_newton.h"
namespace ito33 { namespace ihg {

typedef HazardRateSpotComponentPowerCalibratorNewton 
        HazardRateSpotComponentPowerCalibrator;

}  } 
#endif // #ifdef NEWTONITERATION

#ifdef MIXEDITERATION
#include "ihg/hrspotcomponentpowercalibrator_mixed.h"
#include "ihg/hrspotcomponentpowercalibrator_beta.h"
#include "ihg/hrspotcomponentpowercalibrator_newton.h"

namespace ito33 { namespace ihg {

typedef HazardRateSpotComponentPowerCalibratorMixed 
        HazardRateSpotComponentPowerCalibrator;

}  } 
#endif // #ifdef MIXEDITERATION

#endif // #ifndef _IHG_HRSPOTCOMPONENTPOWERCALIBRATOR_H_

