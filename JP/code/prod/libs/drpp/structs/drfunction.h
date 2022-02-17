// drfunction.h: interface for the drfunction class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DRFUNCTION_H__51A01113_8EF0_11D1_8372_9A1D73000000__INCLUDED_)
#define AFX_DRFUNCTION_H__51A01113_8EF0_11D1_8372_9A1D73000000__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include "drptr.h"
#include <map>

//k  Represents a function from doubles to doubles
//k  Constructed from two valarrays: (x1, x2, ...) and (f(x1), f(x2), ...)
//k  Performs either linear or geometric interpolation for values
//k  inbetween.  However, it does not extrapolate, but rather
//k  use the value of the minimum or maximum x-value.
//k  
//k  To invoke the function, use operator()
//k		eg.  DRFunction foo ( ..);
//k			 double answer = foo (10);
//k
//k  Example:
//k     DRValarray x(2), fx(2);
//k		x[0] = 1; fx[0] = 10;
//k     x[1] = 2; fx[1] = 20;
//k		DRFunction f(x, fx);
//k     Then,
//k			f(1) = 10, f(2) = 20;  f(1.5) = 15;
//k
//k     If f used geometric interp, f(1.5) = sqrt(200)

class DRFunction : 
	protected map <double, double, less<double>, MYALLOC (double) > 

{
public:
	enum FuncInterpType { LINEAR, GEOMETRIC };
	DRFunction(DArray&, DArray&, FuncInterpType interp = LINEAR); //
	
	double operator() (double); //
	
	friend ostream& operator<<(ostream&, const DRFunction&); //
	
protected:
	double (*Func) (double, double, double, double, double);
};

typedef DRPtr<DRFunction> DRFunctionPtr;

#endif // !defined(AFX_DRFUNCTION_H__51A01113_8EF0_11D1_8372_9A1D73000000__INCLUDED_)
