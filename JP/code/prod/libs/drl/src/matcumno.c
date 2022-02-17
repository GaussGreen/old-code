/************************************************************************
 * Module:	DRL - MATH
 * Function:	Mathematical Functions
 * Function:	Cumulative Normal Density
 * Author:	C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */

#include <math.h>
#include <ctype.h>
#include <string.h>

#if defined(DRL_CLIB)
# include "gtobf.h"
# include "nvarnorm.h"	/* C Analytics */
#endif

#include "drlmath.h"	/* prototype consistency */


static	double	_Gauss1(double) ;

/*f--------------------------------------------------------------
 * Math : normal density
 *                                 
 * <br><br>
 * Returns the normal density  {1/ sqrt{2 * pi}} exp(-{x^2 / 2}).
 */

double
DrlDenNorm(double x)
{
	return exp(-0.5*x*x)*0.3989422804e0;
}


/*f--------------------------------------------------------------
 * Math : cumulative normal density
 *                         
 * <br><br>
 * Returns the cumulative normal density
 * 1/ sqrt(2 * pi) integral (-infty, x) exp(-{y^2 / 2})dy.
 */

double
DrlCumNorm(double x)
{
	double	rr;
	if (x <= -21) return 0.;
	rr = 7.071067811865475E-1;
	return(0.5 * _Gauss1(-x * rr));
}


static
double	_Gauss1(double x)
{
	double t, xup, xdown, y, inut;

	xdown = -6.25E0;
	xup = 2.66E1;

	if (x >= xup) return(0.0) ;
	if (x <= xdown) return(2.0) ;

	t = 1.0E0 - 7.5E0/(fabs(x)+3.75E0);

	y = (((((((((((((((3.328130055126039E-10 
		*t-5.718639670776992E-10)*t-4.066088879757269E-9) 
		*t+7.532536116142436E-9)*t+3.026547320064576E-8) 
		*t-7.043998994397452E-8)*t-1.822565715362025E-7) 
		*t+6.575825478226343E-7)*t+7.478317101785790E-7) 
		*t-6.182369348098529E-6)*t+3.584014089915968E-6) 
		*t+4.789838226695987E-5)*t-1.524627476123466E-4) 
		*t-2.553523453642242E-5)*t+1.802962431316418E-3) 
		*t-8.220621168415435E-3)*t+2.414322397093253E-2; 
     
		y = (((((y*t-5.480232669380236E-2)*t+1.026043120322792E-1) 
		*t-1.635718955239687E-1)*t+2.260080669166197E-1) 
		*t-2.734219314954260E-1)*t + 1.455897212750385E-1;

	inut = exp(-x*x)*y ;
	if (x <= 0.0) inut = 2.0 - inut ;
	return(inut) ; 
}

/*f--------------------------------------------------------------
 * Math : cumulative normal density
 *                         
 * <br><br>
 * Returns the cumulative normal density
 * {1/ sqrt{2 * pi}}
 *     integral_{-infty}^x exp(-{y^2 / 2})dy.
 */

double
DrlCumNormInv(double p, double *x)
{
#ifdef	DRL_CLIB
	return GtoNormalCumInv(p, x);
#else
	NOT_IMPLEMENTED;
	return (DBL_MAX);
#endif
}



/*-----------------------------------------------------------------
 * J.P. Morgan Securitie, Inc.  ---  Derivatives Research
 * Module:	Derivatives Research Exotics Toolbox
 * DDate:	
 * Authors:	P. Micottis, L. Pradier
 * Function:	Normal, Bivariate Normal density and cumulative functions.
 */

/*
 *       See Handbook of Mathematical Functions 26.2.19 page 932
 */
#define         d1      .0498673470
#define         d2      .0211410061
#define         d3      .0032776263
#define         d4      .0000380036
#define         d5      .0000488906
#define         d6      .0000053830

#undef		PI
#define         PI      3.141592653
#define         gBIG    exp(15.)




/*---------------------------------------------------------------
 * Normal: Cumulative function of the normal distribution.
 */

static	double
Normal(double x)
{
	double	y;


	/*
	*       x=+l => N(x)=1
	*/
	if (x >= gBIG)
		return (1.);
	else if (x <= -gBIG)
		return (0.);

	y = fabs(x);
	y = 1.+y*(d1+y*(d2+y*(d3+y*(d4+y*(d5+y*d6)))));

	if (x >= 0.)
		return (1.-.5/pow(y,16.));
	else
		return(.5/pow(y,16.));

}

/*---------------------------------------------------------------
 * Pdf: Density function of the normal distribution.
 */

static	double
_Pdf(double x)
{

	return (exp(-.5*pow(x,2.))/sqrt(2.*PI));

}




#define         A1      0.3253030
#define         A2      0.4211071
#define         A3      0.1334425
#define         A4      0.006374323
#define         B1      0.1337764
#define         B2      0.6243247
#define         B3      1.3425378
#define         B4      2.2626645




/*---------------------------------------------------------------
 * f:
 */

static	double
_F(
	double  x,
	double  y,
	double  a,
	double  b,
	double  p)
{
	double	f,
		a_prime,
		b_prime;


	a_prime =  a/(sqrt(2.*(1. - pow(p, 2.))));
	b_prime =  b/(sqrt(2.*(1. - pow(p, 2.))));

	f = exp(a_prime*(2.*x - a_prime) + b_prime*(2.*y - b_prime) +
		      2.*p*(x - a_prime)*(y - b_prime));

	return (f);

}



/*---------------------------------------------------------------
 * M:
 */

static	double
M(
	double a,
	double b,
	double p)
{
	int	i,
		j;
	double	sum = 0.,
		ac[5],
		bc[5] ;


	ac[1] = A1;
	ac[2] = A2;
	ac[3] = A3;
	ac[4] = A4;
	bc[1] = B1;
	bc[2] = B2;
	bc[3] = B3;
	bc[4] = B4;


	for (j = 1; j <= 4; j++)
	{
		for (i = 1; i <= 4; i++)
		{
			sum = sum + ac[i]*ac[j]*_F(bc[i],bc[j],a,b,p);
		}  /* for */
	}  /* for */

	sum = sum*sqrt(1. - pow(p,2.))/PI;

	return (sum);

}  /* M */



/*---------------------------------------------------------------
 * sgn:
 */

static	double
sgn(double a)
{

	if (a < 0.)
		return -1.;
	else
		return 1.;

}  /* sgn */



/*---------------------------------------------------------------
 * delta:
 */

static	double
delta(
	double a,
	double b)
{
	double
		d;


	d = (1. - sgn(a)*sgn(b))/4. ;

	return (d);

}  /* delta */



/*f--------------------------------------------------------------
 * Math : cumulative bivariate normal density
 * 
 * <br><br>
 * Returns the cumulative bivariate normal density
 * integral_{-\infty}^a integral_{-\infty}^b 
 *    exp (-{x^2 - 2 rho x y + y^2 / 2(1-rho^2)})
 *    {dx dy / 2 pi sqrt{1-rho^2}}.
 */

double
DrlCumBiNorm(
	double a,		/* (I) */
	double b,		/* (I) */
	double p)		/* (I) */
{
	double
		p1 = 0.,
		p2 = 0.,
		c;


	if (b < a) {
		c = a;
		a = b;
		b = c;
	}

	if (p == -1.) {
		if (a+b <= 0.)
			return (0.);
		else
			return (Normal(a) + Normal(b) - 1.);
	}

	if (p == 1.)
		return (Normal(a));

	if ((b <= 0.) && (p <= 0.))
		return(M(a,b,p));

	if (a*b*p > 0.)
	{
		p1 = (p*a -b)*sgn(a)/(sqrt(pow(a,2.) - 2*a*b*p + pow(b,2.)));
		p2 = (p*b -a)*sgn(b)/(sqrt(pow(a,2.) - 2*a*b*p + pow(b,2.)));

		return (DrlCumBiNorm(a,0.,p1)+DrlCumBiNorm(b,0.,p2)-delta(a,b));
	}  /* if */

	if (a<= 0. &&  b>=0. && p>=0.)
		return (Normal(a)-M(a,-b,-p));

	if (a>=0. && p<=0.)
		return (Normal(a)+Normal(b)-1.+M(-a,-b,p));

	return 1e30;
}



/****************************************************************
 * Multivariate Normal Cumulative
 *****************************************************************/


double
DrlCumMultiNorm(int dim, double *a, double *r)
{
static	long	inf[] = {1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L};
	long	ifail;
	int	errCode;
	double	eps = 1e-7;
	double	prob,
		bound;
#if defined(DRL_CLIB)
	errCode = GtoMultiNormalCum(dim, inf, a, a, r,
		eps, &prob, &bound, &ifail);
#else
	DrlErrMsg("%s, %d: not available.\n", __FILE__, __LINE__);
#endif
	return(errCode == SUCCESS ? prob : 1e30);
}

