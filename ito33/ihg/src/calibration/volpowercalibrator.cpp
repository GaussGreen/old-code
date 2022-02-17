/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/volpowercalibrator.cpp
// Purpose:     calibrate power volatility
// Author:      Ito33
// Created:     2004/22/11
// RCS-ID:      $Id: volpowercalibrator.cpp,v 1.2 2005/01/28 19:47:17 wang Exp $
// Copyright:   (c) 2004-  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ihg/src/calibration/volpowercalibrator.cpp
   @brief calibration of power volatility
 */
#include "ihg/volpowercalibrator.h"

#ifdef NEWTONITERATION
#include "volpowercalibrator_newton.cpp"
#endif // #ifdef NEWTONITERATION
