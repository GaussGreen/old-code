/*
 * $Log: gaussian.cpp,v $
 * Revision 1.6  2003/09/16 16:20:14  mab
 * ajout de CDFInvNormal
 *
 * Revision 1.5  2003/06/24 14:20:27  jmprie
 * ajout fct de calcul des moments ou des cumul des puissance de X qd
 * X est gaussienne
 *
 * Revision 1.4  2002/11/25 16:32:10  mab
 * Formatting
 *
 */

/*----------------------------------------------------------------------------*

    gaussian.cpp
 
    Routines for densities and CDFs of Gaussian variables

*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>



#include "gaussian.h"
#include "armglob.h"





/*----------------------------------------------------------------------------*
    SYNOPSIS    double dNormal(double d) 
    
    Returns the density of Gaussian distribution (mean 0, variance 1)
    at input value d

*----------------------------------------------------------------------------*/    

double dNormal(double d)
{
    return(exp(-0.5*d*d)*0.398942280402);
}




/*----------------------------------------------------------------------------*
    SYNOPSIS    double cdfNormal(double d)
    
    Returns the CDF of Gaussian distribution (mean 0, variance 1)
    at input value d
*----------------------------------------------------------------------------*/

double cdfNormal(double d)
{
    // NEW
    if ( d < -6 )
    {
       return 0;
    }
    else if ( d > 6 )
    {
        return 1;
    }
    else
    {
       double x,y,z;
       double b1=0.2316419,b2=0.319381530,b3=-0.356563782;
       double b4=1.781477937,b5=-1.821255978,b6=1.330274429;
       double sqrtdeuxpi = 2.5066282746;
 
       if ( d >= 0 )
       {
          x=1/(1+b1*d);
          y=exp(-d*d/2)/sqrtdeuxpi;
          z=x*(b2+x*(b3+x*(b4+x*(b5+x*b6))));
          z=1-y*z;
       }
       else
       {
          x=1/(1-b1*d);
          y=exp(-d*d/2)/sqrtdeuxpi;
          z=x*(b2+x*(b3+x*(b4+x*(b5+x*b6))));
          z=y*z;
       }
 
       return(z);
    }
}


// ---------------------------------------------------------
//Calculation of  the Normale Cumulative
// ---------------------------------------------------------
double qCDFNormal(const double& z__)
{
   	/* System generated locals */
    double ret_val, d__1;

    /* Local variables */
    double zabs, p, expntl;


	/*     Normal distribution probabilities accurate to 1.e-15. */
	/*     Z = no. of standard deviations from the mean. */

	/*     Based upon algorithm 5666 for the error function, from: */
	/*     Hart, J.F. et al, 'Computer Approximations', Wiley 1968 */

	/*     Programmer: Alan Miller */
	/*     Latest revision - 30 March 1986 */

    zabs = fabs(z__);

	/*     |Z| > 37 */

    if (zabs > 37.) {
	p = 0.;
    } else {

	/*     |Z| <= 37 */

	/* Computing 2nd power */
	d__1 = zabs;
	expntl = exp(-(d__1 * d__1) / 2);

	/*     |Z| < CUTOFF = 10/SQRT(2) */
	if (zabs < 7.071067811865475) {
	    p = expntl * ((((((zabs * .03526249659989109 + .7003830644436881) 
		    * zabs + 6.37396220353165) * zabs + 33.912866078383) * 
		    zabs + 112.0792914978709) * zabs + 221.2135961699311) * 
		    zabs + 220.2068679123761) / (((((((zabs * 
		    .08838834764831844 + 1.755667163182642) * zabs + 
		    16.06417757920695) * zabs + 86.78073220294608) * zabs + 
		    296.5642487796737) * zabs + 637.3336333788311) * zabs + 
		    793.8265125199484) * zabs + 440.4137358247522);

	/*     |Z| >= CUTOFF. */
	} else {
	    p = expntl / (zabs + 1 / (zabs + 2 / (zabs + 3 / (zabs + 4 / (
		    zabs + .65))))) / 2.506628274631001;
	}
    }
    if (z__ > 0.) {
	p = 1 - p;
    }
    ret_val = p;
    return (ret_val);

}


/*----------------------------------------------------------------------------*/
/* Numerical approximation of Normal(loi Normale) inverse                     */
/*----------------------------------------------------------------------------*/

double INV_PART_FUNC_NOR(double u)
{
    double x, r;

    double a[4] = {2.50662823884, -18.61500062529, 41.39119773534, 
                   -25.44106049637};

    double b[4] = {-8.47351093090, 23.08336743743, -21.06224101826, 
                   3.13082909833};

    double c[9] = {0.3374754822726147, 0.9761690190917186, 0.1607979714918209, 
                   0.0276438810333863, 0.0038405729373609, 0.0003951896511919,
                   0.0000321767881768, 0.0000002888167364, 0.0000003960315187};
 
    x = u - 0.5;

    if ( fabs(x) < 0.42 )
    {
       r = x*x;

       r = x*(((a[3]*r+a[2])*r+a[1])*r+a[0])
              /((((b[3]*r+b[2])*r+b[1])*r+b[0])*r+1.0);

    }
    else
    {
        r = u;

        if ( x > 0 )
        {
           r = 1.0 - u;
        }

        r = log(-log(r));

        r = c[0]+r*(c[1] + r*(c[2]+r*(c[3]+r*(c[4]+r*(c[5]+r*(c[6]
             +r*(c[7]+r*(c[8]))))))));

        if ( x < 0.0 )
        {
           r = -r;
        }
    }

    return(r);
}



/*----------------------------------------------------------------------------*
    SYNOPSIS    double d2Normal(double r, double x, double y) 
    
    Returns the density of 2-dim Gaussian distribution 
    (mean (0,0), variances (1,1) and correlation r) at
    input value (x,y)
*----------------------------------------------------------------------------*/

double d2Normal(double r, double x, double y)
{
    double    z, r2;

    r2 = SQR(r);
    z = 0.5 * (SQR(x) + SQR(y) - 2.0 * r * x * y) / (1.0 - r2);
    
    return( exp(-z) * 0.159154943 / sqrt(1.0 - r2) );
}




/*----------------------------------------------------------------------------*
    SYNOPSIS    double cdf2Normal(double r, double x, double y)
    
    Returns the CDF of 2-dim Gaussian distribution 
    (mean (0,0), variances (1,1) and correlation r) at 
    input value (x,y)

*----------------------------------------------------------------------------*/

double cdf2Normal(double r, double x, double y)
{
    double    mx, my, r1, r2,    cdf;


    //    this to avoid floating overflow 

    if (fabs(x) <= K_DOUBLE_TOL) 
    {
       x = 0.0;
    }

    if (fabs(y) <= K_DOUBLE_TOL) 
    {
       y = 0.0;
    }

    if (fabs(r) <= K_DOUBLE_TOL) 
    {
       r = 0.0;
    }


    if (x*y*r <= 0.0) 
    {
       if ((x<=0.0) && (y<=0.0)) 
       {
             cdf = m2Normal(r, x, y);
       }
       else 
       {
          if ((x<=0.0) && (y>=0.0)) 
          {
             cdf = cdfNormal(x) - m2Normal(-r, x, -y);
          }
          else 
          {
             if ((x>=0.0) && (y<=0.0)) 
             {
                cdf = cdfNormal(y) - m2Normal(-r, -x, y);
             }
             else 
             {
                if ((x>=0.0) && (y>=0.0)) 
                {
                   cdf = cdfNormal(x) 
                         + cdfNormal(y) - 1.0 + m2Normal(r, -x, -y);
                }
             }
          }
       }
    }
    else 
    {
        r1 = (r*x - y) * fabs(x)/(x*sqrt(x*x - 2.0*r*x*y + y*y));
        r2 = (r*y - x) * fabs(y)/(y*sqrt(x*x - 2.0*r*x*y + y*y));

        if ((x<=0.0) && (r1<=0.0)) 
        {
            mx = m2Normal(r1, x, 0.0);
        }
        else 
        {
           if ((x>=0.0) && (r1<=0.0)) 
           {
              mx = cdfNormal(x) - 0.5 + m2Normal(r1, -x, 0.0);
           }
           else 
           {
              if ((x>=0.0) && (r1>=0.0)) 
              {
                    mx = 0.5 - m2Normal(-r1, -x, 0.0);
              }
              else 
              {
                 if ((x<=0.0) && (r1>=0.0)) 
                 {
                    mx = cdfNormal(x) - m2Normal(-r1, x, 0.0);
                 }
              }
           }
        }

        if ((y<=0.0) && (r2<=0.0)) 
        {
            my = m2Normal(r2, y, 0.0);
        }
        else 
        {
            if ((y>=0.0) && (r2<=0.0)) 
            {
               my = cdfNormal(y) - 0.5 + m2Normal(r2, -y, 0.0);
            }
            else 
            {
               if ((y>=0.0) && (r2>=0.0)) 
               {
                    my = 0.5 - m2Normal(-r2, -y, 0.0);
               }
               else 
               {
                  if ((y<=0.0) && (r2>=0.0)) 
                  {
                     my = cdfNormal(y) - m2Normal(-r2, y, 0.0);
                   }
               }
            }
        }

        cdf = mx + my - (1.0 - (fabs(x) * fabs(y)/(x * y)))/4.0;
    }


    return(cdf);
}



/*----------------------------------------------------------------------------*

    SYNOPSIS
        double f2Normal(double x, double y, double r, double a, double b)
    
    Auxiliary routine for cdf2Normal

*----------------------------------------------------------------------------*/    

double f2Normal(double x, double y, double r, double a, double b)
{
    double    aa, bb, f;


    aa = a/sqrt(2.0*(1.0-r*r));
    bb = b/sqrt(2.0*(1.0-r*r));
    f = exp(aa * (2.0*x - aa) + bb * (2.0*y - bb) + 2.0 * r * (x-aa) * (y-bb));

    return(f);
}



/*----------------------------------------------------------------------------*
    SYNOPSIS    double m2Normal(double r, double a, double b)

    Auxiliary routine for cdf2Normal
    
*----------------------------------------------------------------------------*/

double m2Normal(double r, double a, double b)
{
    double c[4], d[4];
    int    i, j;
    double    sum = 0.0;

    c[0] = 0.3253030;
    c[1] = 0.4211071;
    c[2] = 0.1334425;
    c[3] = 0.006374323;
    d[0] = 0.1337764;
    d[1] = 0.6243247;
    d[2] = 1.3425378;
    d[3] = 2.262664500;
    
    for (i=0; i<4; i++)
        for (j=0; j<4; j++) sum = sum + c[i]*c[j] 
            * f2Normal(d[i], d[j], r, a, b);

    sum = sum * sqrt(1.0-r*r)/3.141592654;
    
    return(sum);
}



/*----------------------------------------------------------------------------*
  SYNOPSIS double maxLogNormalProba(double fwd,double vol,
                                    double T,double Range)

  routine for computing the strike which gives the maximal probability for B&S

*----------------------------------------------------------------------------*/

double maxLogNormalProba(double fwd, double vol, double T, double Range)
{
    if (Range <= 0.0) 
       return (fwd);

    double strike = fwd - 0.5 * Range;
    double epsilon = 1.0e-9;

    if (strike < epsilon)
        strike = 0.01;

    double d2_K = log(fwd/strike) / (vol*sqrt(T)) - 0.5*vol*sqrt(T);
    double result = (strike+Range) * dNormal(d2_K) - strike *
                    dNormal(d2_K+log(strike/(strike+Range))/(vol*sqrt(T)));

    double derive, result_eps, strike_tmp, result_tmp;
    int timer = 0;

    while ((fabs(result) > epsilon) && (timer < 40))
    {
        strike_tmp = strike;
        result_tmp = result;

        // calcul d'une derivee
        strike += epsilon;
        d2_K    = log(fwd/strike) / (vol*sqrt(T)) - 0.5*vol*sqrt(T);
        result_eps = (strike+Range) * dNormal(d2_K) - strike *
                     dNormal(d2_K + log(strike/(strike+Range)) / (vol*sqrt(T)));

        derive = (result_eps - result) / epsilon;

        strike = strike_tmp - result/derive;

        if (strike < epsilon)
            strike = 0.01;

        d2_K    = log(fwd/strike) / (vol*sqrt(T)) - 0.5*vol*sqrt(T);
        result = (strike+Range) * dNormal(d2_K) - strike *
                 dNormal(d2_K + log(strike/(strike+Range)) / (vol*sqrt(T)));

        while ((fabs(result) > fabs(result_tmp)) && (timer < 40))
        {
            strike = 0.5 * (strike + strike_tmp);

            d2_K    = log(fwd/strike) / (vol*sqrt(T)) - 0.5*vol*sqrt(T);
            result = (strike+Range) * dNormal(d2_K) - strike *
                     dNormal(d2_K + log(strike/(strike+Range)) / (vol*sqrt(T)));

            timer++;
        }

        timer++; 
    }

    if ( strike <= 0.01 )
       return(0.0);

    if ( timer == 40 )
       return (fwd-0.5*Range);

    return(strike); 
}

 
/*--------------------------------------------------------------------------*/
/*  Calcul des esperances des puissances d'une loi normale N(mean,sigma) :  */
/*      - entre ]-infini,h] pour cumulNormalMoment()                        */
/*      - entre ]-infini,+infini[ pour normalMoment()                       */
/*--------------------------------------------------------------------------*/

double cumulNormalMoment(double h,int rank,double mean,double sigma,
                         double *cumul)
{
    if ( rank < 0 )
       return 1.0;

    double invSqrtDeuxPiSigma = 0.398942280402*sigma;

    double* result;

    if (!cumul)
    {
       result = new double[rank+1];
    }
    else
    {
       // Suppose a la bonne taille sinon...
       result = cumul;
    }

    double d = (h-mean)/sigma;
    result[0]= cdfNormal(d);

    double mhalfd2 = -0.5*d*d;

    if ( rank > 0 )
    {
       result[1] = mean*result[0]-invSqrtDeuxPiSigma*exp(mhalfd2);
    }

    double hk=1.0;
    double sigma2 = sigma*sigma;
    for (int k = 1; k < rank; k++)
    {
        hk *= h;

        result[k+1] = mean*result[k]+k*sigma2*result[k-1] 
                      -invSqrtDeuxPiSigma*hk*exp(mhalfd2);
    }

    double value = result[rank];

    if (!cumul)
       delete result;

    return value;
}



double normalMoment(int rank,double mean, double sigma,double *cumul)
{
    if ( rank < 0 )
       return(1.0);

    double invSqrtDeuxPiSigma = 0.398942280402*sigma;

    double* result;

    if (!cumul)
    {
       result = new double[rank+1];
    }
    else
    {
       // Suppose a la bonne taille sinon...
       result = cumul;
    }

    result[0]=1.0;

    if ( rank > 0 )
    {
        result[1] = mean*result[0];
    }
    double hk = 1.0;
    double sigma2 = sigma*sigma;

    for (int k = 1; k < rank; k++)
    {
        result[k+1] = mean*result[k]+k*sigma2*result[k-1];
    }

    double value=result[rank];

    if (!cumul)
       delete result;

    return value;
}



// returns the inverse of cumulative normal distribution function
// Reference: The Full Monte, by Boris Moro, Union Bank of Switzerland RISK 1995(2)
double CDFInvNormal(const double& u)
{
    static double a[4] = {2.50662823884, -18.61500062529,
                          41.39119773534, -25.44106049637};

    static double b[4] = {-8.47351093090, 23.08336743743,
                          -21.06224101826, 3.13082909833};

    static double c[9] = {0.3374754822726147, 0.9761690190917186, 
                          0.1607979714918209, 0.0276438810333863,
                          0.0038405729373609, 0.0003951896511919,
                          0.0000321767881768, 0.0000002888167364,
                          0.0000003960315187};

    double x, r;

    x = u-0.5;

    if ( fabs(x) < 0.42 )
    {
       r = x*x;
       r = x*(((a[3]*r+a[2])*r+a[1])*r+a[0])
             /((((b[3]*r+b[2])*r+b[1])*r+b[0])*r+1.0);

       return(r);
    }

    r = u;

    if ( x > 0.0 )
       r = 1.0-u;
   
    r = log(-log(r));

    r = c[0]+r*(c[1]+r*(c[2]+r*(c[3]+r*(c[4]
        +r*(c[5]+r*(c[6]+r*(c[7]+r*c[8])))))));

    if ( x < 0.0 )
       r = -r;

    return(r);
}




/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
