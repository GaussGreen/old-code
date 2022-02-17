/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		BRENTSOLVER.CXX
	PROJECT:	UTIL
	
	DESCRIPTION:	this class provides methods to use a BRENT Solver
					implemented following the numerical recepies.

  -----------------------------------------------------------------

 	CREATION:	October 8, 2004

	LAST MODIF:	October 8, 2004

  -----------------------------------------------------------------
   
	ICM Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

# include "ICMKernel\util\icm_brentsolver.h"
# include "ICMKernel\util\icm_macro.h"

// # include <stdlib.h>


void BrentSolver :: Init()
{
	Reset();
}

void BrentSolver :: Reset()
{
	NbMaxIter	=	50;
	Target		=	0.0;	
}

BrentSolver* BrentSolver :: Clone() const
{
	return new BrentSolver(*this);
}

void BrentSolver::BitwiseCopy(const BrentSolver* src)
{
	NbMaxIter	=	src->NbMaxIter;
	Target		=	src->Target;
}

void BrentSolver :: Copy(const BrentSolver* src)
{
	if (src != this)
		BitwiseCopy(src);
}


double BrentSolver :: ZBrent(double a, double fa, double b, double fb, double c, double fc) const
{
	double z;

	if ((fa == 0.0) || (fb == 0.0) || (fc == 0.0))
	{
		return( 0.0 );
	}
	else if ((fa == fc) || (fb == fc) || (fb == fa))
	{
		return( b );
	}
	else
	{
		z = (fb / fa) * ((fa / fc) * ((fb / fc) - (fa / fc)) * (c - b) - (1.0 - (fb / fc)) * (b - a));
		z /= ((fa / fc) - 1.0) * ((fb / fc) - 1.0) * ((fb / fa) - 1.0);
		z += b;
		if ((fabs(z-a) < 1e-12) || (fabs(z-b) < 1e-12) )
			z = (a + b) * 0.5;
		return(z);
	}
}

//----------------------------------------------------------------------------
// Looking Outward
ReturnCode BrentSolver :: ZeroBracket(	double x1,
										double x2,
										double&	fx1,
										double&	fx2,
										unsigned int NbIterMax,
										double factor,
										ReturnCode (*fct)(void*, double, double&),
										void* param,
										bool& Result) const
{	
	ReturnCode	Err;

	Err = (*fct)(param, x1, fx1);
	if (Err != RetOk) return Err;

	Err = (*fct)(param, x2, fx2);
	if (Err != RetOk) return Err;

	unsigned int j;

	for (j=0; j<NbIterMax; j++)
	{
		fx1 -= Target;
		fx2 -= Target;

		if (fx1*fx2 < 0.)
		{
			Result = true;
			return RetOk;
		}
		else if (fabs(fx1) < fabs(fx2))
		{
			x1 += factor * (x1-x2);
			Err = (*fct)(param, x1, fx1);
			if (Err != RetOk) return Err;
		}
		else
		{
			x2 += factor * (x2-x1);
			Err = (*fct)(param, x2, fx2);
			if (Err != RetOk) return Err;
		}
	}

	Result = false;

	return RetNotOk;
}


//----------------------------------------------------------------------------
// Looking Inward
ReturnCode BrentSolver :: ZeroBracket_In(	double x1,
										double x2,
										ReturnCode (*fct)(void*, double, double&),
										void* param,
										int	n,
										double	xb1[],
										double	xb2[],
										int*		nb,
										double	fx1[],
										double	fx2[]
										) const
{	
// Given a function fx defined on the interval from x1-x2 subdivide the interval into n equally
// spaced segments, and search for zero crossings of the function. nb is input as the maximum number
// of roots sought, and is reset to the number of bracketing pairs xb1[1..nb], xb2[1..nb]
// that are found.

	ReturnCode	Err;

	int		nbb,i;
	double	x,fp,fc,dx;

	nbb=0;

	dx	=	(x2-x1)/n; // Determine the spacing appropriate to the mesh.
	
	x=x1;
	Err = (*fct)(param, x, fp);
	if (Err != RetOk) return Err;

	fp	-=	Target;

	for (i=1;i<=n;i++)
	{
		// Loop over all intervals
		x += dx;
		Err = (*fct)(param, x, fc);
		if (Err != RetOk) return Err;

		fc	-=	Target;

		if (fc*fp <= 0.0)
		{
			// If a sign change occurs then record values for the bounds
			xb1[++nbb]=x-dx;
			xb2[nbb]=x;

			fx1[nbb]=fp;
			fx2[nbb]=fc;
			if (*nb == nbb) return RetOk;
		}
		fp=fc;
	}
	*nb = nbb;

	return RetOk;
}

//----------------------------------------------------------------------------

ReturnCode BrentSolver :: Solve(
	double x1, 
	double x2, 
	double tol,
	int&	nb_iter,	
	ReturnCode (*fct)(void*, double, double&), 
	void* param, 
	double& root,
	bool	NeedZeroBracket,
	int		ZeroBracketChoice,
	int		NbIntervals) const
{
	ReturnCode	Err;
	double fx1, fx2;

	if ((fct == NULL) || (tol <= 0.0)) return RetNotOk;

	if (NeedZeroBracket)
	{
		if (ZeroBracketChoice == 0)
		{
			bool	Result;
			unsigned int NbIterMax=100;
			double factor=1.6;

			Err = ZeroBracket(x1, x2, fx1, fx2, NbIterMax, factor, fct, param, Result);
			if (Err != RetOk) return Err;
		}
		else
		{
			int	n	=	NbIntervals;
			int	nb	=	fabs(ZeroBracketChoice);	// maximum number of roots expected

			double*	xb1;
			double*	xb2;
			double*	fb1;
			double*	fb2;

			xb1	=	new double[nb+1];
			xb2	=	new double[nb+1];
			fb1	=	new double[nb+1];
			fb2	=	new double[nb+1];

			Err = ZeroBracket_In(x1, x2, fct, param, n, xb1, xb2, &nb, fb1, fb2); 
			if (Err != RetOk) return Err;

			if (nb == 0)
			{
				ICMLOG("No roots found!")
				
				delete [] xb1;
				delete [] xb2;
				delete [] fb1;
				delete [] fb2;

				return RetNotOk;
			}
			else
			{
				ICMLOG("Number of roots found: " << nb);

				// just use one root...
				x1	=	xb1[1];
				x2	=	xb2[1];

				fx1	=	fb1[1];
				fx2	=	fb2[1];
			}

			delete [] xb1
				;
			delete [] xb2;
			delete [] fb1;
			delete [] fb2;

		}
	}
	else
	{
		Err = (*fct)(param, x1, fx1);
		if (Err != RetOk) return Err;

		fx1	-=	Target;

		Err = (*fct)(param, x2, fx2);
		if (Err != RetOk) return Err;

		fx2	-=	Target;
	}

	return Solve(x1, x2, fx1, fx2, tol, nb_iter, fct, param, root);
}

ReturnCode BrentSolver :: Solve(
	double x1,
	double x2,
	double fx1,
	double fx2,
	double tol,
	int&	nb_iter,
	ReturnCode (*fct)(void*, double, double&),
	void* param,
	double& root ) const
{
	ReturnCode	Err;

	double a,b,c,fa,fb,fc;
	double Newb;

	if ((fct == NULL) || (tol <= 0.0)) return RetNotOk;

	a = x1;
	c = x2;

	fa = fx1;
	fc = fx2;
//	fa = fx1 - Target;
//	fc = fx2 - Target;

	if (fabs(fa) < tol) 
	{
		root = a; 
		return RetOk;
	}
	else
  	{
		if (fabs(fc) < tol) 
		{
			root = c; 
			return RetOk;
		}
  		else 
			b = (a * fc - c * fa) / (fc - fa);
		if (fabs(b - a) < 1e-10 || fabs( b - c)< 1e-10)
			b = (a + c) * .5;
	}  	  		
	
	for (nb_iter = 1; nb_iter <= NbMaxIter; nb_iter ++)
	{
		Err = (*fct)(param, b, fb);
		if (Err != RetOk) return Err;
		fb -= Target;
		if (fabs(fb) >= tol)
		{ 
			Newb = ZBrent(a, fa, b, fb, c, fc);                                             
			if (Newb == b) 
	    	{ 
				Newb = 0.5 * (a + c);
			}
			else if (((Newb < a) && (Newb < c)) || ((Newb > a) && (Newb > c)))
			{
					Newb = 0.5 * (a + c);
					if (Newb == b) return RetNotOk;
			}
			if (fa * fb > 0.0) 
			{ 
				a = b; 
				fa = fb;
			}
			else 
			{ 
				c = b; 
				fc = fb;
			}
			b = Newb;
		}
		else 
	  	{ 
			root = b; 
			return RetOk;
		}
	}
	return RetNotOk;
}

