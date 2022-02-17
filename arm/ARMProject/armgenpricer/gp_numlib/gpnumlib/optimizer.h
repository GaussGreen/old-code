/*!
 *
 * Copyright (c) IXIS CIB November 2004 Paris
 *	\file solver.h
 *
 *  \brief 
 *	\author  A. Schauly
 *	\version 1.0
 *	\date November 2004
 */

#ifndef _INGPNUMLIB_OPTIMIZER_H
#define _INGPNUMLIB_OPTIMIZER_H

/// gpbase
#include "gpbase/env.h"
#include "gpbase/ostringstream.h"
#include "gpbase/port.h"
#include "gpbase/utilityport.h"

#include "mrgk5.h"

/// ARM Kernel
#include "expt.h"

/// Standard Library
#include <cmath>


CC_BEGIN_NAMESPACE( ARM )

/// general constant
struct OptimizerConstant
{
	static const double DefaultXTolerance	;	
	static const double DefaultFxTolerance	;	
	static const size_t DefaultMax_Iter		;	
	static const double DefaultStepMax		;     
	static const double UpperInfiniteBound	;
	static const double LowerInfiniteBound	;
	static const double	Zero				;				
	static const double DoubleTolerance		;	
	static const bool   DefaultLocalSearch	;
	static const bool	DefaultPrint		;
};
						

 

///////////////////////////////////////////
/// \brief
///	in this file, you can find the following solvers:
///		-dichotomy optimizer
///		-brent optimizer
///////////////////////////////////////////


template <typename Func, typename T=double> class T_1DOptimizerWithInitialGuess
{
public:
	virtual T Optimize()												= 0;
	virtual void setInitialGuess(T initialGuess, 
		T lowerBoundary = OptimizerConstant::LowerInfiniteBound,
		T middlePoint = OptimizerConstant::Zero,
		T upperBoundary = OptimizerConstant::UpperInfiniteBound )		= 0;
	virtual T_1DOptimizerWithInitialGuess<Func,T>* Clone() const	= 0;
};

//////////////////////////
/// GoldenSection optimizer
/// a function should support T operator(T);
//////////////////////////

template <typename F, typename T=double> class T_GoldenSectionOptimizer : public T_1DOptimizerWithInitialGuess<F,T>
{
public:
	T_GoldenSectionOptimizer( const F& f,T lowerBoundary,T middlePoint,T upperBoundary, T precision, size_t maxIters );
	
	/// virtual functions
	~T_GoldenSectionOptimizer() { delete xmin;}
	virtual T Optimize();
	virtual T_1DOptimizerWithInitialGuess<F,T>* Clone() const
	{ return new T_GoldenSectionOptimizer<F,T>(*this);	}

	virtual void setInitialGuess( T initialGuess, 
		T lowerBoundary = OptimizerConstant::LowerInfiniteBound,
		T middlePoint = OptimizerConstant::Zero, 
		T upperBoundary = OptimizerConstant::UpperInfiniteBound )
	{ 
		itsLowerBoundary=lowerBoundary; 
		itsMiddlePoint=middlePoint; 
		itsUpperBoundary=upperBoundary;
	}

private:
	const F& itsF;
	T itsLowerBoundary;
	T itsMiddlePoint;
	T itsUpperBoundary;
	T itsPrecision;
	T *xmin;
	size_t itsMaxIters;
};


//////////////////////////
/// Brent solver
/// a function should support T operator(T);
//////////////////////////

template <typename F, typename T=double> class T_BrentOptimizer : public T_1DOptimizerWithInitialGuess<F,T>
{
public:
	T_BrentOptimizer( const F& f,T lowerBoundary,T middlePoint,T upperBoundary,T precision,	size_t maxIters );

	~T_BrentOptimizer() { delete xmin; }

	/// virtual functions
	virtual T Optimize();
	virtual T_1DOptimizerWithInitialGuess<F,T>* Clone() const
	{ return new T_BrentOptimizer<F,T>(*this); }

	virtual void setInitialGuess( T initialGuess, 
		T lowerBoundary = OptimizerConstant::LowerInfiniteBound,
		T middlepoint = OptimizerConstant::Zero,
		T upperBoundary = OptimizerConstant::UpperInfiniteBound )
	{	
		itsLowerBoundary=lowerBoundary;
		itsMiddlePoint=middlepoint;
		itsUpperBoundary=upperBoundary; 
	}

private:
	const F& itsF;
	T itsLowerBoundary;
	T itsMiddlePoint;
	T itsUpperBoundary;
	T itsPrecision;
	T *xmin;
	size_t itsMaxIters;
};

///////////////////////////////////
//// implementation
///////////////////////////////////

///////////////////////////////////
//// GoldenSection optimizer
//// Taken and adapted from numerical recipies in C
///////////////////////////////////

#define SHFT2(a,b,c) (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define R 0.61803399
#define C (1.0-R)

template <typename F, typename T> T T_GoldenSectionOptimizer<F,T>::Optimize( )
{
	T f1,f2,x0,x1,x2,x3;
	T ax,bx,cx;
	ax = itsLowerBoundary;
	bx = itsMiddlePoint;
	cx = itsUpperBoundary;
	x0=ax;
	x3=cx;

	if (fabs(cx-bx) > fabs(bx-ax)) 
	{
		x1=bx;
		x2=bx+C*(cx-bx);
	} else {
		x2=bx;
		x1=bx-C*(bx-ax);
	}
	f1=itsF(x1);
	f2=itsF(x2);

	while (fabs(x3-x0) > itsPrecision*(fabs(x1)+fabs(x2))) 
	{
		if (f2 < f1) 
		{
			SHFT3(x0,x1,x2,R*x1+C*x3)
			SHFT2(f1,f2,itsF(x2))
		} else 
		{
			SHFT3(x3,x2,x1,R*x2+C*x0)
			SHFT2(f2,f1,itsF(x1))
		}
	}
	if (f1 < f2) 
	{
		if( itsF( itsLowerBoundary ) < f1 )
			x1 = itsLowerBoundary;
		else if ( itsF( itsUpperBoundary ) < f1 )
			x1 = itsUpperBoundary;

		*xmin=x1;
		return *xmin; //f1 in the original algorithm
	} else {

		if( itsF( itsLowerBoundary ) < f2 )
			x2 = itsLowerBoundary;
		else if ( itsF( itsUpperBoundary ) < f2 )
			x2 = itsUpperBoundary;

		*xmin=x2;
		return *xmin; //f2 in the original algorithm
	}

	return *xmin;
}

#undef C
#undef R
#undef SHFT2
#undef SHFT3

template <typename F, typename T>
T_GoldenSectionOptimizer<F,T>::T_GoldenSectionOptimizer<F,T>( const F& f,T lowerBoundary,T middlePoint,T upperBoundary,
		T precision		= OptimizerConstant::DefaultFxTolerance,size_t maxIters	= OptimizerConstant::DefaultMax_Iter ):	
		itsF(f), 
		itsLowerBoundary(lowerBoundary), 
		itsMiddlePoint(middlePoint),
		itsUpperBoundary(upperBoundary), 
		itsPrecision(precision),
		itsMaxIters(maxIters)
{
	xmin = new T();
}

//////////////////////////////////////
/// Brent optimzer
//////////////////////////////////////

#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

template <typename F, typename T> T T_BrentOptimizer<F,T>::Optimize( )
{
	size_t iter;
	T a,b,d,etemp,fu,fv,fw,fx,p,q,r,itsPrecision1,itsPrecision2,u,v,w,x,xm;
	T ax,bx,cx;
	T e=0.0;

	ax = itsLowerBoundary;
	bx = itsMiddlePoint;
	cx = itsUpperBoundary;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=itsF(x);
	for (iter=1;iter<=itsMaxIters;iter++) {
		xm=0.5*(a+b);
		itsPrecision2=2.0*(itsPrecision1=itsPrecision*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (itsPrecision2-0.5*(b-a))) {

			if( itsF( itsLowerBoundary ) < itsF( x ) )
				x = itsLowerBoundary;
			else if ( itsF( itsUpperBoundary ) < itsF( x ) )
				x = itsUpperBoundary;

			*xmin=x;
			return *xmin; // fx in the original algorithm
		}
		if (fabs(e) > itsPrecision1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < itsPrecision2 || b-u < itsPrecision2)
					d=CC_SignMult(itsPrecision1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= itsPrecision1 ? x+d : x+CC_SignMult(itsPrecision1,d));
		fu=itsF(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}

	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
	"Too many iterations in <BrentOptimizer>");
}
#undef CGOLD
#undef ZEPS
#undef SHFT


template <typename F, typename T>
T_BrentOptimizer<F,T>::T_BrentOptimizer<F,T>( const F& f,	T lowerBoundary,T middlePoint,T upperBoundary,
		T precision		= OptimizerConstant::DefaultFxTolerance,
		size_t maxIters	= OptimizerConstant::DefaultMax_Iter ):	
		itsF(f), 
		itsLowerBoundary(lowerBoundary), 
		itsMiddlePoint(middlePoint),
		itsUpperBoundary(upperBoundary), 
		itsPrecision(precision),
		itsMaxIters(maxIters)
{ 
	xmin = new T; 
}

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
