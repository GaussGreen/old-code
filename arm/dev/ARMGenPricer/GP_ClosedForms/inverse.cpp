/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file inverse.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include <cmath>

#include "gpclosedforms/long_double.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/inverse.h"

#include "gpbase/numericconstant.h"





CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_EPS 1.0e-13

inline double SIGN(const double &a, const double &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}


// #define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
template<class T>
inline void SWAP(T &a, T &b)
	{T dum=a; a=b; b=dum;} 

/////////////////////////////////////////////////////////////////////////
///
/// Inversion of a regular function through a zbrent algorithm
/// by default:
///			typical_x = 0.0
///			deltat_x  = 0.1
///			tol       = 1e-6
///
/////////////////////////////////////////////////////////////////////////

double Inverse::operator() (double ygoal, double typical_x0, double delta_x0, double tol)
{

	///////////////////////////////////////////////////////////////////////////////////
	/// adaptation of the functor and other variables to to the constraints
	///////////////////////////////////////////////////////////////////////////////////
	double typical_x;
	double delta_x;

	
	struct interpret_ptr_funDoubleToDouble : public DoubleToDoubleFunc
	{
		interpret_ptr_funDoubleToDouble(DoubleToDoubleFunc* f ,int imageC)
			:itsPFunc(f),itsImageConst(imageC) {}
		virtual double operator()(double x) const 
		{
			switch (itsImageConst) 
			{
			case 	Inverse::REAL :
				{
					return (*itsPFunc)(x);
						break;
				}
			case  Inverse::ALWAYSPOSITIVE :
				{
					return (*itsPFunc)(exp(x));
					break;
				}
			case  Inverse::CORRELATION :
				{
					return (*itsPFunc)(tanh(x));
					break;
				}
			case  Inverse::BOUNDEDBY0AND1 :
				{
					return (*itsPFunc)((tanh(x)+1.)/2.);
					break;
				}
			default:
				{
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Inverse::operator(): unknown image constraint" );
				}	
			}
		}
	private:
		DoubleToDoubleFunc *itsPFunc;
		int itsImageConst;
	};
	

	interpret_ptr_funDoubleToDouble pfunc(itspFunc,imageConstraints);
	
	switch (imageConstraints) 
	{
	case 	REAL :
		{
			typical_x=typical_x0;
			delta_x=delta_x0;
			break;
		}
	case	ALWAYSPOSITIVE :
		{
			typical_x=log(typical_x0);
			delta_x=log((typical_x0+delta_x0)/typical_x0);
			break;
		}
	case 	CORRELATION :
		{
			typical_x=atanh(typical_x0);
			delta_x=atanh(typical_x0+delta_x0)-atanh(typical_x0);
			break;
		}
	case 	BOUNDEDBY0AND1 :
		{
			typical_x=atanh(2.*typical_x0-1.) ;
			delta_x=atanh(2.*(typical_x0+delta_x0)-1.)-atanh(2.*typical_x0-1.) ;
			break;
		}
	default:
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"Inverse::operator(): unknown image constraint" );
		}
	}

	
	///  start of bracketing the solution
	double x1=typical_x,x2=typical_x,res1,res2;
	double i=0.0;
	if(pfunc(x1)-ygoal > 0)
	{
		x2=x1;res1=1.0;res2=1.0;
		while((res1 >=0.0) && (res2 >=0.0))
		{	
			i=2*i+1;
			if(i*delta_x>1e50) 			
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Inverse::operator():Root impossible to  bracket");

			res1=pfunc(typical_x+i*delta_x)-ygoal;
			res2=pfunc(typical_x-i*delta_x)-ygoal;

			if( res1== res2)
			{
				double u=pfunc(typical_x+i*delta_x)-ygoal;;
			}
		}
		if (res1 <0.0) x1=typical_x+i*delta_x;
		if (res2 <0.0) x1=typical_x-i*delta_x;
	}
	else
	{
		res1=-1.0;res2=-1.0;
		while((res1 <=0.0) && (res2 <=0.0))
		{	
			i=2*i+1;
			if(i>1e50) 
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Inverse::operator(): Root impossible to  bracket ");

			res1=pfunc(typical_x+i*delta_x)-ygoal;
			res2=pfunc(typical_x-i*delta_x)-ygoal;
		}
		
		if (res1 >0.0) x2=typical_x+i*delta_x;
		if (res2 >0.0) x2=typical_x-i*delta_x;
	}
	
	///   start of zbrent algorithm
	return brentSolve(pfunc,ygoal,x1,x2,tol,ARM_CF_MAXIT,imageConstraints);

}

double returnImage(double x, int ImageConst)
{
	switch (ImageConst) 
	{
	case Inverse::REAL:
		{
			return x;
		}

	case  Inverse::ALWAYSPOSITIVE :
		{
			return exp(x);
				break;
		}
	case  Inverse::CORRELATION :
		{
			return tanh(x);
			break;
		}
	case  Inverse::BOUNDEDBY0AND1 :
		{
			return (tanh(x)+1.)/2.;
			break;
		}
	default:
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Inverse::operator(): unknown image constraint" );
		}	
	}
}


double brentSolve(const DoubleToDoubleFunc& func, double ygoal, double lowerBound, double upperBound, double tol, int MAXITER, int ImageConst, double * best)
{
	int iter;
	double a = lowerBound,b = upperBound, c = upperBound,d,e,min1,min2;
	double fa = func(a)-ygoal, fb = func(b)-ygoal,fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
	{
		if(best != NULL)
		{
			*best = fabs(fa) < fabs(fb) ? a : b;
			return returnImage(*best,ImageConst);
		}
		else
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Inverse::operator(): Root is not bracketed !");
	}

	fc=fb;
	for (iter=0;iter<MAXITER;iter++) 
	{
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
		{
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb))
		{
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*ARM_CF_EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) 
		{
			return returnImage(b,ImageConst);
		}
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
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
			b += SIGN(tol1,xm);
			fb=func(b)-ygoal;
	}

	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Inverse::operator(): Maximum number of iterations exceeded in zbrenting!");
	return 0.0;

}

CC_END_NAMESPACE()


#undef ARM_CF_EPS
#undef ARM_CF_MAXIT


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

