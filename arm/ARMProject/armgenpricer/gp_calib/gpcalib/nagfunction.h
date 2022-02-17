/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 * 
 * \file nagfunction.h
 *  \brief nagfunction provides some simple function for
 *		nag optimisation
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2004
 */


#ifndef _INGPCALIB_NAGFUNCTION_H
#define _INGPCALIB_NAGFUNCTION_H

#include "gpbase/env.h"
#include "gpbase/port.h"

/// nag headers
#include "gpbase/removenagwarning.h"
#include "nag.h"
#include "nage04.h"

CC_BEGIN_NAMESPACE( ARM )

void NAG_CALL NagConfunThrowExceptionIfUsed(Integer n, Integer m, Integer needc[], double x[],
	double conf[], double cjac[], Nag_Comm* comm);

void NAG_CALL NagErrorhandler( char *strng, int code,  char *name);


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
