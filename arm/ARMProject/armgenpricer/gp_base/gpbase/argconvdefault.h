/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file argconvdefault.h
 *
 *  \brief default table for interface reading
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPBASE_ARGCONVDEFAULT_H
#define _INGPBASE_ARGCONVDEFAULT_H

/// use our macro for namespace
#include "port.h"
#include "argconv.h"

CC_BEGIN_NAMESPACE( ARM )

extern const ARM_ArgConv ARM_ArgConv_LgNameDayCount;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_LgNameDayCount;
extern const ARM_ArgConv ARM_ArgConv_LgNameFrequency;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_LgNameFrequency;
extern const ARM_ArgConv ARM_ArgConv_FwdRules;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_FwdRules;
extern const ARM_ArgConv ARM_ArgConv_LgTimingMod;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_LgTimingMod;
extern const ARM_ArgConv ARM_ArgConv_InterestRules;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_InterestRules;
extern const ARM_ArgConv ARM_ArgConv_StubRules;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_StubRules;
extern const ARM_ArgConv ARM_ArgConv_RcvOrPay;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_RcvOrPay;
extern const ARM_ArgConv ARM_ArgConv_InterpolType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_InterpolType;
extern const ARM_ArgConv ARM_ArgConv_CompoundingType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_CompoundingType;
extern const ARM_ArgConv ARM_ArgConv_CompoundingFrequency;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_CompoundingFrequency;


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
