/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: env.h,v $
 * Revision 1.1  2003/10/16 07:52:17  ebenhamou
 * Initial revision
 *
 *
 */



/*! \file env.h
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPBASE_ENV_H
#define _INGPBASE_ENV_H

#include "port.h"
#include <string>
CC_USING_NS( std, string )


CC_BEGIN_NAMESPACE( ARM )

extern string ARM_USERNAME;

extern string ARM_COMPUTERNAME;

/// function to initialize environment variables
extern void LOCALARM_InitEnvVariables();

/// in debug mode, we do a strict validation
/// strict validation can be redundant validation
/// to give safeguard when changing, inserting, creating code....
#ifdef _DEBUG
	#define __GP_STRICT_VALIDATION
	#define __GP_CHECK_SORTED_VECTOR
	#define __GP_SHOW_SHARED_NODE_COORDINATES		/// to show cordinates
	#define __GP_SFRM_STOCH_TERMS
#endif


/// in no debug mode (release ) mode, we avoid superflous checking in the kernel
/// object like ARM_Vector and ARM_Matrix
#ifdef NDEBUG
	#define __ARM_MATRIX_NO_RANGE_CHECK
	#define __ARM_VECTOR_NO_RANGE_CHECK
	#define __ARM_LINALG_NO_RANGE_CHECK
	#define __SET_PTR_TO_NULLL
	#define __GP_SFRM_STOCH_TERMS
#endif



CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

