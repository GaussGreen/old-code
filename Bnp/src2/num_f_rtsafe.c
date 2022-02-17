/* ==========================================================
   FILENAME :   num_f_rtsafe.c

   PURPOSE:     Find the root of a function using Newton-Rhapson method 
				or bissection
   
   DESCRIPTION: Using a combination of Newton-Rhapson and bissection, 
				find the root of a function bracketed between x1 and x2.
				The root, returned as the function value rtsafe, will be
				refined until its accuracy is known within + or - xacc.
				funcd is a user-supplied routine that returns both
				the function value and the first derivative of the function
                
   See NUMERICAL RECIPES IN C page 366								
   ========================================================== */

#include "UTALLHDR.H>
#include "math.h"

/*
 modified this version not to complain if goes outside the boundary
 */

SrtErr rtsafe(Err (*func)(double, double *, double *), 
              double x1, double x2, double xacc, int num_iter,
			  double *answer)
{
Err err = NULL;
int j;
double df,dx,dxold,f,fh,fl; 
double temp,xh,xl,rts;
                        
	err = (*func)(x1,&fl,&df);
	if(err) return err;
	err  = (*func)(x2,&fh,&df);  
	if(err) return err;

	if ((fabs(fl) < DBL_EPSILON) && (fabs(fh) < DBL_EPSILON)) 
		return serror("RTSafe - Too many solution in range");
	else
	if (fabs(fl) < DBL_EPSILON) 
	{
		*answer = x1; 
		return NULL;
	}
	else
	if (fabs(fh) < DBL_EPSILON) 
	{
		*answer = x2;
		return NULL;
	}
		
	if ( (fl>0.0 && fh>0.0) || (fl<0.0 && fh<0.0) )
	{
		if(fabs(fl) < fabs(fh)) 
		{ 
			*answer = x1; 
			return NULL;
		}
		else 
		{ 
			*answer = x2; 
			return NULL;
		}
	}

	/* Orient the search so that fl < 0 and fh > 0 */
	if (fl < 0.0)
	{ 
		xl = x1;
		xh = x2;

	} 
	else 
	{
		xh = x1;
		xl = x2;
	}                       

	rts = 0.5 * (x1 + x2);
	dxold = fabs(x2- x1);
	dx = dxold;

	err = (*func)(rts,&f,&df); 
	if(err) return err;

	for (j=1 ; j <= num_iter ; j++)
	{
		if (((	(rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0)
			||	(fabs(2.0*f) > fabs(dxold*df)))
		{
			/* bisect if Newton is out of range, or not decreasing fast enough. */
			dxold = dx;
			dx = 0.5*(xh-xl);
			rts = xl+dx;
			if (xl == rts) 
			{ 
				*answer = rts; /* the change is negligible */
				return NULL;
			}
		}
		else
		{
			/* the change in Newton is acceptable; take it */
			dxold = dx;
			dx = f/df;
			temp = rts;
			rts -= dx;
			if (temp == rts)
			{
				*answer = rts; 
				return NULL;
			}
		}


		if (fabs(dx) < xacc)
		{
			*answer = rts;
			return NULL;
		}

		err = (*func)(rts,&f,&df);
		if(err) return err;

		/* maintain the bracket on the root */
		if (f < 0.0)
			xl = rts;
		else
			xh = rts;
	}

	return serror("rtsafe hit max iterations");
}



/**************************************************************************************************************************************/
//
// rtsafe_with_par
//
// Modified version of rtsafe that passes a double * to the function being zeroed
//
/**************************************************************************************************************************************/
SrtErr rtsafe_with_par(Err (*func)(double, double *, double *, double *), 
              double x1, double x2, double xacc, int num_iter,
			  double *answer, double *par)
{
Err err = NULL;
int j;
double df,dx,dxold,f,fh,fl; 
double temp,xh,xl,rts;
                        
	err = (*func)(x1,&fl,&df,par);
	if(err) return err;
	err  = (*func)(x2,&fh,&df,par);  
	if(err) return err;

	if ((fabs(fl) < DBL_EPSILON) && (fabs(fh) < DBL_EPSILON)) 
		return serror("RTSafe - Too many solution in range");
	else
	if (fabs(fl) < DBL_EPSILON) 
	{
		*answer = x1; 
		return NULL;
	}
	else
	if (fabs(fh) < DBL_EPSILON) 
	{
		*answer = x2;
		return NULL;
	}
		
	if ( (fl>0.0 && fh>0.0) || (fl<0.0 && fh<0.0) )
	{
		if(fabs(fl) < fabs(fh)) 
		{ 
			*answer = x1; 
			return NULL;
		}
		else 
		{ 
			*answer = x2; 
			return NULL;
		}
	}

	/* Orient the search so that fl < 0 and fh > 0 */
	if (fl < 0.0)
	{ 
		xl = x1;
		xh = x2;

	} 
	else 
	{
		xh = x1;
		xl = x2;
	}                       

	rts = 0.5 * (x1 + x2);
	dxold = fabs(x2- x1);
	dx = dxold;

	err = (*func)(rts,&f,&df,par); 
	if(err) return err;

	for (j=1 ; j <= num_iter ; j++)
	{
		if (((	(rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0)
			||	(fabs(2.0*f) > fabs(dxold*df)))
		{
			/* bisect if Newton is out of range, or not decreasing fast enough. */
			dxold = dx;
			dx = 0.5*(xh-xl);
			rts = xl+dx;
			if (xl == rts) 
			{ 
				*answer = rts; /* the change is negligible */
				return NULL;
			}
		}
		else
		{
			/* the change in Newton is acceptable; take it */
			dxold = dx;
			dx = f/df;
			temp = rts;
			rts -= dx;
			if (temp == rts)
			{
				*answer = rts; 
				return NULL;
			}
		}


		if (fabs(dx) < xacc)
		{
			*answer = rts;
			return NULL;
		}

		err = (*func)(rts,&f,&df,par);
		if(err) return err;

		/* maintain the bracket on the root */
		if (f < 0.0)
			xl = rts;
		else
			xh = rts;
	}

	return serror("rtsafe_with_par hit max iterations");
}
