/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: argconvdefault.h,v $
 * Revision 1.1  2004/06/08 16:45:06  ebenhamou
 * Initial revision
 *
 */

    
/*----------------------------------------------------------------------------*/

/*! \file argconvdefault.h
 *
 *  \brief default table for interface reading
 *
 *	\author  Richard GUILLEMOT
 *	\version 1.0
 *	\date October 2004
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPCALCULATOR_ARGCONVDEFAUIT_H
#define _INGPCALCULATOR_ARGCONVDEFAUIT_H

/// use our macro for namespace²
#include "gpbase/port.h"
#include <gpbase/argconv.h>

CC_BEGIN_NAMESPACE( ARM )

extern const  ARM_ArgConv ARM_ArgConv_GenCalculatorCcyType;
extern const  ARM_ArgConvReverse ARM_ArgConvReverse_GenCalculatorCcyType;

extern const  ARM_ArgConv ARM_ArgConv_MRSCalibType;
extern const  ARM_ArgConvReverse ARM_ArgConvReverse_MRSCalibType;

extern const  ARM_ArgConv ARM_ArgConv_MRSStrikeCalibType;
extern const  ARM_ArgConvReverse ARM_ArgConvReverse_MRSStrikeCalibType;

extern const ARM_ArgConv ARM_ArgConv_SigmaCalibType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_SigmaCalibType;

extern const ARM_ArgConv ARM_ArgConv_TARNCalibMode;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_TARNCalibMode;

extern const ARM_ArgConv ARM_ArgConv_CaptionCalibMode;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_CaptionCalibMode;

extern const ARM_ArgConv ARM_ArgConv_CSBCalibMode;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_CSBCalibMode;

extern const ARM_ArgConv ARM_ArgConv_CSBControlVariableMode;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_CSBControlVariableMode;

extern const ARM_ArgConv ARM_ArgConv_CSOCalibMode;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_CSOCalibMode;

extern const ARM_ArgConv ARM_ArgConv_CSBTriggerMode;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_CSBTriggerMode;

extern const ARM_ArgConv ARM_ArgConv_CRASpreadCalibrationType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_CRASpreadCalibrationType;

extern const ARM_ArgConv ARM_ArgConv_CRASpreadCalibrationType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_CRASpreadCalibrationType;

extern const ARM_ArgConv ARM_ArgConv_CRASpreadCalibStrikeType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_CRASpreadCalibStrikeType;

extern const ARM_ArgConv ARM_ArgConv_PRDCCalibMode;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_PRDCCalibMode;

extern const ARM_ArgConv ARM_ArgConv_PRCSRedemptionType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_PRCSRedemptionType;

extern const ARM_ArgConv ARM_ArgConv_PRCSBasisType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_PRCSBasisType;

extern const ARM_ArgConv ARM_ArgConv_CSOCalibModeHWM2F;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_CSOCalibModeHWM2F;

extern const ARM_ArgConv ARM_ArgConv_CSOCalibModeHWM1F;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_CSOCalibModeHWM1F;

extern const ARM_ArgConv ARM_ArgConv_CSOStrike1Type;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_CSOStrike1Type;

extern const ARM_ArgConv ARM_ArgConv_CSOStrike2Type;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_CSOStrike2Type;

extern const ARM_ArgConv ARM_ArgConv_CSOStrike3Type;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_CSOStrike3Type;

extern const ARM_ArgConv ARM_ArgConv_TARNFXModelType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_TARNFXModelType;

extern const ARM_ArgConv ARM_ArgConv_FXVanillaType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_FXVanillaType;

extern const ARM_ArgConv ARM_ArgConv_FXBasketType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_FXBasketType;

extern const ARM_ArgConv ARM_ArgConv_MixPRDNoticeType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_MixPRDNoticeType;

extern const ARM_ArgConv ARM_ArgConv_TARNFXPayoffType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_TARNFXPayoffType;

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
