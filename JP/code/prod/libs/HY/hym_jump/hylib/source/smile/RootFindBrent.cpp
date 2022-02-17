// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 1999 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// 1.0 10/7/99 Afshin Bayrooti
// $Header$
//

//#include "General/General.h"

#include "RootFindBrent.h"


int dummyRootFindFunction(double x, void  *data, double *result)
{
	*result =( *(AbstractFunctorWrapper*)(data) )(x) ;
	return 0;
}



