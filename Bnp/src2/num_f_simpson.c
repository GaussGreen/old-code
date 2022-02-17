/* ==========================================================================
   FILE_NAME:	num_f_simpson.c

   PURPOSE:     Numerical Integration using Simpson's rule    
   ========================================================================== */

#include        "utallhdr.h"
#include        <num_h_simpso.h"
#include		<math.h"


#define		EPS_SIMPSO 	(1e-7)
#define 	JMAX 		30

/* =======================================================================
   FUNC: sm_trapzd
   DESC: Trapezoid rule driver based on function in Numerical Recipes in C. 
         Slightly modified from the bits operations
   INPUTS: Migth be used with a va_list (Undefined number of parameters)
   MODIFIES:
   DECLARATION:
   ======================================================================= */

double sm_trapzd(	double (*function)(double), 
			double a, 
			double b, 
			int n, 
			double accum)
{
	double x, tnm, sum, del,ret;
	long it, j;

	if (n==1) 
	{
		accum = 0.5*(b-a)*( function(a) + function(b) );
	}
	else
	{
		it = 1;
		for (j=1 ;j < n-1; j++) 
				it *= 2;
		tnm = (double) it;
		del = (b-a)/tnm;
		x = a + 0.5 * del;
		sum = 0.0;
		for (j=1; j<=it; j++) 
		{
			sum += function(x);
			x += del;
		}
		accum = 0.5 * ( accum + (b-a)*sum/tnm );
	}
	ret = accum;
	return ret;
} /* END sm_trapzd */


/*******************************************************************************
*                                                                           
* FUNCTION     	: sm_qsimp(...)                                              
*                                                                           
* PURPOSE      	: Simpson's Rule for numerical integration        	    
*                                                                           
* DESCRIPTION  	: Based on function in "Numerical Recipes in C"		    
*                                                                           
* CALLS		: sm_trapzd(...)                                           
*                                                                           
* PARAMETERS   	: (*func)	- pointer to function returning a double    
*              	: a        	- ??                                	    
*              	: b      	- ??                                  	    
*               : precision
*                                                                           
* RETURNS      	: sum        	- ??                                        
*                                                                           
*******************************************************************************/

double 	sm_qsimp	(
			double 	(*func)(double),
			double 	a, 
			double 	b,
			double precision
			)
{       
	int 	j;
    	double 	sum;
	double	oldsum;
	double	sum2N;
	double	sumN;
			

	sumN  = sm_trapzd(func, a, b, 1, 0.0);
	sum2N = oldsum = -1.0e30;
    
	for (j = 2; j < JMAX; j++)
	{
		sum2N  = sm_trapzd(func, a, b, j, sumN);

		sum    = ( (4.0 * sum2N) - sumN ) / 3.0;
		if ( fabs(sum - oldsum) <= (precision * fabs(oldsum)) ) 
			return sum;
		oldsum = sum;
		sumN   = sum2N;
	}

	return sum;

} /* END sm_qsimp(...) */	




/* ========================================================================
   FUNC: sm_qsimp_list                                                     
   DESC: Simpson's rule for numerical integration, 
         based on function in Numerical Recipes in C.
         Undefined number of parameters. 
   MODIFIES:
   DECLARATION:
   ======================================================================== */

double sm_qsimp_list(	double (*function)(double,va_list),
			double a, 
			double b,...)

/*additional parameters may be : SwapDP p,SwapDP pf,double coupon, double clean_price,double volatility,
			 double myu, int T,
                         double redemption, double first_coupon*/
{
    va_list argptr;
    int j;
    double sum,oldsum,sum2N,sumN,ret;
    va_start(argptr,b);


    sumN = sm_trapzd_list(function,a,b,1,0.0,argptr);

    sum2N = oldsum = -1.0e30;
    for(j=2;j<JMAX;j++)
    {
	sum2N = sm_trapzd_list(function,a,b,j,sumN,argptr);
	sum = (4.0*sum2N - sumN)/3.0;
	if(fabs(sum - oldsum) < EPS_SIMPSO*fabs(oldsum)) {ret = sum; break;}
	oldsum = sum;
 	sumN = sum2N;
    }
    if(j>=JMAX)	ret = 0;
    va_end(argptr);
    return ret;
}	

/* ======================================================================= */



/* =======================================================================
   FUNC: sm_trapzd_list
   DESC: Trapezoid rule driver based on function in Numerical Recipes in C. 
         Slightly modified from the bits operations
   INPUTS: Migth be used with a va_list (Undefined number of parameters)
   MODIFIES:
   DECLARATION:
   ======================================================================= */

double sm_trapzd_list(	double (*function)(double, va_list), 
			double a, 
			double b, 
			int n, 
			double accum, 
			va_list argptr)
{
	double x, tnm, sum, del,ret;
	int it, j;

	if (n==1) {
		ret = 0.5*(b-a)*(function(a,argptr)+function(b,argptr));
	}
	else
	{
		for (it=1,j=1;j<n-1;j++) 
				it <<= 1;
		tnm = it;
		del = (b-a)/tnm;
		x = a + 0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) 
			sum += function(x,argptr);
		accum = 0.5*(accum+(b-a)*sum/tnm);
		ret = accum;
	}
	return ret;
}

#undef	EPS_SIMPSO 
#undef 	JMAX 
