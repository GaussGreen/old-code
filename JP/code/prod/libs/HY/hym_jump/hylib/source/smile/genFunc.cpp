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

#include "kplatform.h"

#include "genfunc.h"

double BaseFunction::deriv(double x) const
{
	double temp, dx=0.001;
	temp = (operator()(x+dx)-operator()(x))/dx;
	return temp;
}


double BaseFunction::inverse(double y) const
{
	SoverFunc  sfunc(this,y);

	double x = RootFindBrent(sfunc,
				 0,
				 ROOTLO,
				 ROOTHI);
	return x;
} 




 



