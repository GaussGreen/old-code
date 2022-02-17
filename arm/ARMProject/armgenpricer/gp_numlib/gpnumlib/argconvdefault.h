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
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date March 2004
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPNUMLIB_ARGCONVDEFAUIT_H
#define _INGPNUMLIB_ARGCONVDEFAUIT_H

/// use our macro for namespace
#include "gpbase/port.h"
#include "gpbase/argconv.h"

CC_BEGIN_NAMESPACE( ARM )

/// part for the random nb generators!
extern const ARM_ArgConv ARM_ArgConv_BaseGenAlgoType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_BaseGenAlgoType;

extern const ARM_ArgConv ARM_ArgConv_TransformAlgoType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_TransformAlgoType;

extern const ARM_ArgConv ARM_ArgConv_RandGenOrder;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_RandGenOrder;

extern const ARM_ArgConv ARM_ArgConv_ODESolverType;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_ODESolverType;

extern const ARM_ArgConv ARM_ArgConv_RegMode;
extern const ARM_ArgConvReverse ARM_ArgConvReverse_RegMode;

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
