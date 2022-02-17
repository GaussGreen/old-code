/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/hrspotcomponentpowercalibrator.cpp
// Purpose:     calibrate with a spot component power Hazard rate
// Author:      Ito33
// Created:     2004/22/11
// RCS-ID:      $Id: hrspotcomponentpowercalibrator.cpp,v 1.2 2005/01/28 19:47:17 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ihg/src/calibration/hrspotcomponentpowercalibrator.cpp
   @brief calibration of time component of HRSpotComponentPower
 */
#include "ihg/hrspotcomponentpowercalibrator.h"

#ifdef BETAITERATION
#include "hrspotcomponentpowercalibrator_beta.cpp"
#endif // #ifdef BETAITERATION

#ifdef NAGITERATION
#include "hrspotcomponentpowercalibrator_nag.cpp"
#endif // #ifdef NAGITERATION

#ifdef NEWTONITERATION
#include "hrspotcomponentpowercalibrator_newton.cpp"
#endif // #ifdef NEWTONITERATION

#ifdef MIXEDITERATION
#include "hrspotcomponentpowercalibrator_mixed.cpp"
#include "hrspotcomponentpowercalibrator_newton.cpp"
#include "hrspotcomponentpowercalibrator_beta.cpp"
#endif // #ifdef MIXEDITERATION
