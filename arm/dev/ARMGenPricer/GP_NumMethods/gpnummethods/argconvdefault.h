/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: argconvdefault.h,v $
 * Revision 1.1  2004/03/23 16:45:06  ebenhamou
 * Initial revision
 *
 */

    
/*----------------------------------------------------------------------------*/

/*! \file argconvdefault.h
 *
 *  \brief default table for interface reading
 *
 *	\author  J-M Prié
 *	\version 1.0
 *	\date November 2004
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPNUMMETHODS_ARGCONVDEFAULT_H
#define _INGPNUMMETHODS_ARGCONVDEFAULT_H

/// use our macro for namespace
#include "gpbase/port.h"
#include "gpbase/argconv.h"

CC_BEGIN_NAMESPACE( ARM )

extern const ARM_ArgConv ARM_ArgConv_SamplerType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_SamplerType;

extern const ARM_ArgConv ARM_ArgConv_TruncatorType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_TruncatorType;

extern const ARM_ArgConv ARM_ArgConv_SchedulerType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_SchedulerType;

extern const ARM_ArgConv ARM_ArgConv_ReconnectorType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_ReconnectorType;

extern const ARM_ArgConv ARM_ArgConv_SmootherType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_SmootherType;

extern const ARM_ArgConv ARM_ArgConv_PDENumSchemeType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_PDENumSchemeType;

extern const ARM_ArgConv ARM_ArgConv_PDEBoundConditionType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_PDEBoundConditionType;

extern const ARM_ArgConv ARM_ArgConv_ImpSamplerType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_ImpSamplerType;

extern const ARM_ArgConv ARM_ArgConv_PathSchemeType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_PathSchemeType;

extern const ARM_ArgConv ARM_ArgConv_PDEGridType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_PDEGridType;

extern const ARM_ArgConv ARM_ArgConv_CFmethodType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_CFmethodType;

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
