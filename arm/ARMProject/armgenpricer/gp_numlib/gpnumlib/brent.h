/*!
 *
 * Copyright (c) CDC IXIS CM November 2004 Paris
 *	\file solver.h
 *
 *  \brief 
 *	\author  A. TRIKI
 *	\version 1.0
 *	\date November 2004
 */

#ifndef _INGPNUMLIB_BRENT_H
#define _INGPNUMLIB_BRENT_H

#include "solver.h"
#include "gpbase/env.h"
#include "gpbase/ostringstream.h"
#include "gpbase/port.h"
#include "mrgk5.h"
#include "expt.h"
#include <cmath>
#include <iomanip> /// for setprecision()
using namespace std;

CC_BEGIN_NAMESPACE( ARM )


//////////////////////////
/// Brent solver
/// a function should support T operator(T);
//////////////////////////

template <typename Func, typename T=double> class T_BrentSolver : public T_SolverWithInitialGuess<Func,T>
{
private:
		bool	itsKeepInitialInterval;
		bool	KeepInitialInterval(){return itsKeepInitialInterval; };
public:
		T_BrentSolver( const Func& f,
		const T target	= 0,
		const T ftolerance	= SolverConstant::DefaultFxTolerance,
		const T xtolerance	= SolverConstant::DefaultFxTolerance,
		const size_t max_iter	= SolverConstant::DefaultMax_Iter,
		const bool printLevel	= false,
		ARM_WarningKeeper* wKeeper = 0 ,
		const bool keepInitialInterval = false)
	:
		T_SolverWithInitialGuess<Func,T>(f,
			target,
			max_iter,
			ftolerance,
			xtolerance,
			printLevel,
			wKeeper), itsKeepInitialInterval(keepInitialInterval)
	{}
	~T_BrentSolver() {}
	/// virtual functions
	virtual T Solve();
	virtual T_SolverWithInitialGuess<Func,T>* Clone() const { return new T_BrentSolver<Func,T>(*this); }
};

///////////////////////////////////
//// implementation
///////////////////////////////////

//////////////////////////////////////
/// Brent Solver
//////////////////////////////////////

template <typename Func, typename T> T T_BrentSolver<Func,T>::Solve()
{
	T a=GetLowerBound(),b=GetUpperBound(),c,d,e,min1,min2;
	T fa=GetF()(a)-GetTarget(),fb=GetF()(b)-GetTarget(),fc,p,q,r,s,tol1,xm;
	
    size_t iter=0;
    T df;
	if (KeepInitialInterval())
	{
		if((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) 
			if( GetWarningKeeper() )
				GetWarningKeeper()->AddWarningMessage( "solution not found in initial interval in <BrentSolver>\n" );
	}
	else
	{
		while( ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) && iter<GetMaxIters())
		{
			/// Locate the solution
			df=(fb-fa)/(b-a);
			if( (df>SolverConstant::DoubleTolerance && fa<0) || (df<-SolverConstant::DoubleTolerance && fa>0) )
			{
				/// Shift interval towards positive values
				c=b;
				fc=fb;
				b+=(b-a);
				fb=GetF()(b)-GetTarget();
				a=c;
				fa=fc;
			}
			else if( (df>SolverConstant::DoubleTolerance && fa>0) || (df<-SolverConstant::DoubleTolerance && fa<0) )
			{
				/// Shift interval towards negative values
				c=a;
				fc=fa;
				a -= (b-a);
				fa=GetF()(a)-GetTarget();
				b=c;
				fb=fc;
			}
			else
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Slope to low to do the bracketing in <BrentSolver>");
			}
			++iter;
		}
		
		if(iter>=GetMaxIters())
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Can'y locate a valid initial solution in <BrentSolver>");
		}
	}	
    c=b;
    fc=fb;
	T bestx,bestfx;
	if(fabs(fa)<fabs(fb))
	{
		if(fabs(fc)<fabs(fa)) {bestx=c; bestfx=fc;}
		else {bestx=a; bestfx=fa;}
	}
	else
	{
		if(fabs(fc)<fabs(fb)) {bestx=c; bestfx=fc;}
		else {bestx=b; bestfx=fb;}
	}
	for (iter=1;iter<=GetMaxIters();++iter)
    {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*3.0e-8*fabs(b)+0.5*GetFxTolerance();
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)){
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
            b += (xm>0 ? tol1 : -tol1);
		fb=GetF()(b)-GetTarget();

		if(fabs(fb) < bestfx) {bestx=b; bestfx=fb;} 
	}

	/// No more throw error but add a warning message and return the best solution
   //throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
   //     "Not enough iterations (GetMaxIters()) to find root in <BrentSolver>");
	if( GetWarningKeeper() )
		GetWarningKeeper()->AddWarningMessage( "Not enough iterations to find root in <BrentSolver>\n" );
	
	return bestx;
}

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
