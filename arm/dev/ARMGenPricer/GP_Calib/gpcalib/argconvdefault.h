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
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date June 2004
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPCALIB_ARGCONVDEFAUIT_H
#define _INGPCALIB_ARGCONVDEFAUIT_H

/// use our macro for namespace
#include "gpbase/port.h"
#include <gpbase/argconv.h>

CC_BEGIN_NAMESPACE( ARM )

extern const ARM_ArgConv ARM_ArgConv_CalibMethod;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_CalibMethod;
extern const ARM_ArgConv ARM_ArgConv_TargetFuncMethod;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_TargetFuncMethod;
extern const ARM_ArgConv ARM_ArgConv_CalibDirectionMethod;
extern const ARM_ArgConv ARM_ArgConv_CalibDirection2DMethod;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_CalibDirectionMethod;
extern const ARM_ArgConv ARM_ArgConv_SolverTypeMethod;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_SolverTypeMethod;
extern const ARM_ArgConv ARM_ArgConv_OptimizerTypeMethod;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_OptimizerTypeMethod;
extern const ARM_ArgConv ARM_ArgConv_BasketCalibrationType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_BasketCalibrationType;
extern const ARM_ArgConv ARM_ArgConv_BasketCalibrationStrike;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_BasketCalibrationStrike;
extern const ARM_ArgConv ARM_ArgConv_MoneyType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_MoneyType;



CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
