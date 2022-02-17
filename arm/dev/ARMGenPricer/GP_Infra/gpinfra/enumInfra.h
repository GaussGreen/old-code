/*
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 * $Log: enumInfra.h,v $
 * Revision 1.1  2006/11/14 14:51:19  emezzine
 * Initial revision
 *
 *
 */

/*! \file enumInfra.h
 *
 *  \brief 
 *
 *	\author  E.M Ezzine 
 *	\version 1.0
 *	\date NOvember 2006
 */


#ifndef _INGPINFRA_ENUMINFRA_H
#define _INGPINFRA_ENUMINFRA_H

#include "gpbase/port.h"

CC_BEGIN_NAMESPACE( ARM )

struct ARM_ModelParamBump
{
	enum BumpType
		{
            isCumulative = 0,
		    isPerturbative,
		    isNothing,
		};
};

struct ARM_FXDigitType
{
	enum DigitType
	{
		analytic=0,		/// Pure analytic digital when it is possible
		centred,		/// centred call spread ( (Call(K-eps)-Call(K+eps))/(2eps) )
		backward,		/// backward call spread ( (Call(K-eps)-Call(K))/eps )
		forward,		/// forward call spread ( (Call(K)-Call(K+eps))/eps )
	};
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/