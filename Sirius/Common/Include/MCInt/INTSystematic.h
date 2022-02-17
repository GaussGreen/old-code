#ifndef __INTYSTEMATIC_H__
#define __INTSYSTEMATIC_H__


#include "INTSTL.h"
#include "INTUtilities.h"


void  systematic(const vector<CIntegrand*>& f,	const CDomain& D, STLDoubleVector& I, 
									 const STLLongVector& N,				const vector<CAlpha*>& a, 
									 const vector<CWeights*>& w,		const vector<CPeriodization*>& p, 
									 const STLIntegerVectorVector& cat,		const STLIntegerVectorVector& levels);



#endif
