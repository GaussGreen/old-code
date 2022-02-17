/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: argconvdefault.h,v $
 * Revision 1.1  2005/06/30 16:45:06  emezzine
 * Initial revision
 *
 */

    
/*----------------------------------------------------------------------------*/

/*! \file argconvdefault.h
 *
 *  \brief default table for interface reading
 *
 *	\author  El Mostafa EZZINE
 *	\version 1.0
 *	\date June 2005
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPCALIB_ARGCONVDEFAUIT_H
#define _INGPCALIB_ARGCONVDEFAUIT_H

/// use our macro for namespace
#include "gpbase/port.h"
#include <gpbase/argconv.h>

CC_BEGIN_NAMESPACE( ARM )

extern const ARM_ArgConv ARM_ArgConv_PricingModelType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_PricingModelType;

extern const ARM_ArgConv ARM_ArgConv_MultiAssetsType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_MultiAssetsType;

extern const ARM_ArgConv ARM_ArgConv_VnsPricingMethod;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_VnsPricingMethod;

extern const ARM_ArgConv ARM_ArgConv_MMCalibProxy;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_MMCalibProxy;

extern const ARM_ArgConv ARM_ArgConv_MMCorrelType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_MMCorrelType;

extern const ARM_ArgConv ARM_ArgConv_MMCalibPattern;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_MMCalibPattern;

extern const ARM_ArgConv ARM_ArgConv_HWSVFormula;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_HWSVFormula;

extern const ARM_ArgConv ARM_ArgConv_HestonMCScheme;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_HestonMCScheme;

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
