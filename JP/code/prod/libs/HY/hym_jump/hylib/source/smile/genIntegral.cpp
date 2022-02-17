// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 1999 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// 1.0 10/19/99 Neil Yang
// 
//

//#include "General/General.h"

#include "genIntegral.h"
#include "gtobf.h"

int dummyFunction(double x, void  *data, double *result)
{
	*result =( *(AbstractFunctorWrapper*)(data) )(x) ;
	return 0;
}



 



