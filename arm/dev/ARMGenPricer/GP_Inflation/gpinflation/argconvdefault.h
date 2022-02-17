/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: argconvdefault.h,v $
 * Revision 1.1  2004/09/09 16:45:06  ebenhamou
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
 *	\date September 2004
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPINFLATION_ARGCONVDEFAUIT_H
#define _INGPINFLATION_ARGCONVDEFAUIT_H

/// use our macro for namespace
#include "gpbase/port.h"
#include <gpbase/argconv.h>
#include <gpinflation/typedef.h>

CC_BEGIN_NAMESPACE( ARM )

extern const ARM_ArgConv ARM_ArgConv_InfSwoptComputationMethod;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_InfSwoptComputationMethod;;

extern const ARM_ArgConv ARM_ArgConv_Copula;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_Copula;

extern const ARM_ArgConv ARM_ArgConv_InfIndex;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_InfIndex;

extern const ARM_ArgConv ARM_ArgConv_SubordIndex;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_SubordIndex;


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
