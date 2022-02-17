/****************************************************************************/
/*      kutility.cpp                                                               */
/****************************************************************************/


#include <math.h>
#include <ctype.h>
#include "kutility.h"
#include "kexception.h"

double	 BiNomial( int	k, int N, double	p)
{ 
	int		i;
	double	prob;
    
	if(k<0 || k > N || p < 0 || p > 1)
		throw KException("Wrong argument value!");
	if (k == 0)
		return(pow(1. - p, (double)N));
    else
	{
        prob = 1.;
        if (k < N / 2.) 
		    for (i = 0 ; i<k; i++)
                prob *=(double)(N - i) / (k - i);
        else
            for (i = 0; i<N - k; i++)
                prob *=(double)(N - i) / (N - k - i);
		prob *= pow(p, (double)k) * pow(1. - p, (double)(N-k));
		return(prob);
	}
}



/*****  Normal  *************************************************************/
/*
*       Cumulative function of the normal distribution.
*/
double  Normal(double  x)
{
	double	y, z;
    int	i;        

	double	d1 = .0498673470;
	double	d2 = .0211410061;
	double	d3 = .0032776263;
	double	d4 = .0000380036;
	double	d5 = .0000488906;
	double	d6 = .0000053830;

	if (x >= 10)                        
		return (1.);
	else if (x <= -10)
		return (0.);

	y = fabs (x);
	y = 1. + y * (d1 + y * (d2 + y * (d3 + y * (d4 + y * (d5 + y * d6)))));
        
	for (i = 1, z = y; i < 16; i++)		
                z *= y;

	if (x >= 0.)
		return (1. - .5 / z);
	else
		return (.5 / z);

}  


double  NormalInv(double  x)
{
	double	acc = 0.00000001;
	double	up, down, mid, p;

	if(x < - acc || x > 1 + acc)
		throw KException("Probability is out of range!");

	up = 10.;
	down = -10.;
	if(x <= 0.)
		return(down);
	if(x >= 1.)
		return(up);
	mid = (up + down)/2;
	p = Normal(mid);
	while(fabs(p-x) > acc)
	{
		if(p > x)
			up = mid;
		else
			down = mid;
		mid = (up + down)/2;
		p = Normal(mid);
	}
	return(mid);
}



#define CUTOFF        20.
#define ACCURACY_H    1E-18
#define ACCURACY_L    1E-10

static double TFunction(double h,                /* (I) */
                         double a);               /* (I) */

static double TSFunction(double h,               /* (I) */
                          double a);              /* (I) 0 < a < 1 */


/* --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  --
function: BiNormal
created by: 10/14/94 Neil Yang
description: This function calculates the value for cumulative bivariate normal 
distribution by untilizing a Taylor series method.  The integration ranges for 
the bivariate distribution are from -infinity to a and -infinity to b.
References: Owen, D.B., Tables for computing bivariate normal probabilities.  
        Ann. Math. Stat. 27 (1956), 1075-1090.
        Handbook of mathematical functions, Ed M. Abramowitz and I. Stegun, 
        National Bureau of Standards, 1972, P936.

--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  -- */

double	BiNormal(double a,double b, double rho)
{
    double rho1;            /* sqrt(1-rho*rho) */
    double a_prime;         /* (b/a-rho)/sqrt(1-rho*rho) */
    double b_prime;         /* (a/b-rho)/sqrt(1-rho*rho) */
    double tFunction_a, tFunction_b;  /* temperary variables */

    
    /* check to see if |a|,|b| are too large */
    if (a > CUTOFF)
        return	Normal(b);
    if (a <-CUTOFF)
        return	0.;
    if (b > CUTOFF)
        return	Normal(a);
    if (b <-CUTOFF)
        return	0;
  
    /* for rho=0. */
    if(rho == 0)
        return Normal(a) * Normal(b);
    /* for rho=-1. */
    if(rho == -1)
    {
        if ( a + b <= 0.0)
            return 0;
        else
			return Normal(a) + Normal(b) - 1.0;
    }
    /* for rho=1. */
    if(rho == 1)
        return ((b<=a) ? Normal(b) : Normal(a));
    
    /* for rho != 1, -1, 0 */
    /* for a=0 and b=0 */
    if(a == 0 && b == 0)
        return 0.25 + asin(rho)*(1./(2.*K_PI));

    rho1 = 1./sqrt(1.-rho*rho);
    /* if a=0,b!=0) */
    if(a == 0)
    {
        if (b>0.)
            tFunction_a=0.25;
        else
            tFunction_a=-0.25;
    }
    /* a!=0, b!=0 */
    else
    {
        a_prime = (b/a-rho)*rho1;
        tFunction_a=TFunction(a,a_prime);
    }
    /* a!=0, b=0 */
    if(b == 0)
    {
        if (a>0.)
            tFunction_b=0.25;
        else
            tFunction_b=-0.25;
    }
    /* a!=0,b!=0 */
    else
    {
        b_prime = (a/b-rho)*rho1;
        tFunction_b=TFunction(b,b_prime);
    }
    if (a*b > 0.0 || (a*b == 0 && (a+b > 0.0 || a + b == 0)))
        return  0.5*(Normal(a)+Normal(b))
						-tFunction_a-tFunction_b;
    else
        return 0.5*(Normal(a)+Normal(b))
            -tFunction_a-tFunction_b-0.5;
}


/*
*      Forwad price of an option using Black&Scholes.
*/
double	BlackOption(double	fwd, 
					double	strike, 
					double	t, 
					double	vol, 
					char	coP)

{
	double	P,d,v;
	bool	cop = (toupper(coP) == 'C'? true:false);

	if(fwd == 0.)
		return (cop == true? 0.:strike);
	v  = vol * sqrt (t);
	d  = log (strike / fwd) / v;
	v /= 2.;
	if(cop == true)
		P = fwd * Normal (d + v) - strike * Normal (d - v);
	else
		P  = strike * Normal (d + v) - fwd * Normal (d - v);

	return (P);

}


/*------------------------------------------------------------------------
function: TFunction

created by: 10/14/94 - Neil Yang

description: Calculates T function used by GtoBiNormalCum.  It needs to call 
TSFunction(h,a) to get the values for 0 < a < 1.

Reference: Owen, D.B., Tables for computing bivariate normal probabilities.  Ann.
       Math. Stat. 27 (1956), 1075-1090.


-------------------------------------------------------------------------- */
static double TFunction(double h, double a)
{                 
    double      ph, pha;    /* temp variables */            
    double      a_absolute; /* |a| */                  
    double      tha;        /* Function Value T(h,a) */
    
    if(h == 0) 
        return( atan(a)/(2.*K_PI) );
    if(a == 0)
        return( 0.0 );
    if(a == 1)
    {
        ph = Normal(h);
        return(0.5*ph*(1.-ph));
    }
    /* for a!=0,1 */
    a_absolute = fabs(a);     
    if(a_absolute > 1.)       /* take care of |a| > 1 */
    {
        ph = Normal(h);
        pha = Normal(a_absolute*h);
        tha = 0.5*(ph+pha)-ph*pha
            -TSFunction(a_absolute*h,1./a_absolute);
    }
    else
        tha = TSFunction(h,a_absolute);
    return(a<0.0? -tha : tha);  
}


/*------------------------------------------------------------------------
function: TSFunction

created by: 10/14/94 - Neil Yang

description: Calculates T function for 0 < a < 1 using a Taylor series method.

Reference: Owen, D.B., Tables for computing bivariate normal probabilities. 
       Ann. Math. Stat. 27 (1956), 1075-1090.


-------------------------------------------------------------------------- */
/* T(h,a)=arctan(a)/2PI-(1/2PI)sum(cj*a^2j+1)
   cj=(-1)^j(1/(2j+1))(1-e^(-0.5*h^2)sum(to j)(h^2i)/(2^i*i!)) */

static double TSFunction(double h, double a)
{ 
    double t;                 /* value of the T function */
    double term;              /* Term in the summation of the series */
    double cterm=1.;          /* Term in the summation of coefficient series */
    double csum=1.;           /* Summation in the coefficient */
    double ap;                /* a^(2j+1) */
    double h2;                /* h^2 */
    double h4;                /* (h^2)/2 */
    double ep;                /* exp(-h^2/2) */
    double a2;                /* a^2 */
    double accuracy;          /* required accuracy */
    double cnt=1.0;           /* cnt = 2j+1 */
    double temp;

    int j=0, toggle=1;

    ap=a;
    a2=a*a;
    h2=h*h;
    h4=h2*0.5;
    ep=exp(-h4);
    t=(1.-ep)*a;
    if (fabs(h)>5.)
    {
        accuracy=ACCURACY_H;      /* used when T(h,a) become small */ 
    }
    else
    {
        accuracy=ACCURACY_L;      /* used in most cases. not better accuracy
                                  is needed as that of NormalCum is 1E-6 */
    }
    do
    {
        cterm *= h4/++j;
        csum += cterm;
        ap *= a2;
        cnt += 2.;
        if (toggle)
        {
            term=(ep*csum-1.)*ap/cnt;
        }
        else
        {
            term=(1.-ep*csum)*ap/cnt;
        }
        t += term;
        toggle=!toggle;
    }
    while (fabs(term)>accuracy);
	temp = atan(a);
    t=(temp-t)*(1./(2.*K_PI));
    return(fabs(t));
}


#undef CUTOFF
#undef ACCURACY_H
#undef ACCURACY_L

