/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_MATHS.CPP
	PROJECT:	GLOB
	
	DESCRIPTION:	some useful math tools

  -----------------------------------------------------------------

 	CREATION:	May 10, 2005

	LAST MODIF:	May 10, 2005
  -----------------------------------------------------------------
   
	ICM Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

# include "ICMKernel\glob\icm_maths.h"

# include "ARMKernel\util\gaussian.h"
#include "ARMKernel\glob\linalg.h"

#include "ICMKernel/util/icm_macro.h"

#include <nag.h>		
#include <nags.h>	
#include <nagg01.h>	
#include <nage01.h>	
#include <nage02.h>	
#include <nage04.h>	
#include <nagg05.h>		

#include <math.h>
#include "ICMKernel\util\icm_integrator.h"

#ifndef EPS_GL
#define EPS_GL 3.0e-11
#endif

#ifndef EPS_GH
#define EPS_GH	3.0e-14				// Relative precision.
#endif

#ifndef PIM4_GH
#define PIM4_GH 0.7511255444649425	//	1/Pi^1/4.
#endif

#ifndef MAXIT_GH
#define MAXIT_GH 100				//	Maximum iterations.
#endif

#include <nag_stdlib.h>

double	 BSFormula(
					bool	Call,			// true = cap, false = floor
					double	Sigma,
					double	Spot,
					double	Strike,
					double	NbD,
					bool	Hedges,			// true if one wants all greeks to be computed
					double* Price,
					double* Delta,
					double* Gamma,
					double* Vega,
					double* Teta)
{
	double d1k, d2k, nd1k, nd2k, invexd2;
	double srt, srtpi;

	*Price=0.0; *Delta=0.0; *Gamma=0.0; *Vega=0.0; *Teta=0.0;

	// -----------------------------
	// Some basic consitency tests
	if (Sigma <= 0.0)
		return -1.0;
	if (Spot < 0.0)
		return -1.0;
	if (NbD <= 0.0)
		return -1.0;
	// -----------------------------
	
	// Degenerated case
	if (Strike <= 0.0)
	{
		if (Call)
		{
			*Price = Spot - Strike;
			if (Hedges)
				*Delta = 1.0;
		}
		return 1.0;
	}

	srt = Sigma * sqrt(NbD);
	d1k = CalBs(1.0, Spot, Strike, Sigma, NbD);
	d2k = CalBs(-1.0, Spot, Strike, Sigma, NbD);
	srtpi = srt * SQRTTWOPI;

	if (d2k < -15)
	{
		nd2k = 0.0;
		invexd2 = 0.0;
	}
	else if (d2k > 15)
	{
		nd2k = 1.0;
		invexd2 = 0.0;
	}
	else
	{
		nd2k	= cdfNormal(d2k);
		invexd2 = exp(-(d2k * d2k) / 2.0);
	}
	if (d1k < -15)
		nd1k = 0.0;
	else if (d1k > 15)
		nd1k = 1.0;
	else
		nd1k = cdfNormal(d1k);

	if (Call)
	{
		*Price = Spot * nd1k - Strike * nd2k;
		if (Hedges)
		{
			*Delta = nd1k;
			*Gamma = (Strike / (Spot * Spot * srtpi) * invexd2);
			*Vega = (Strike * sqrt(NbD)) / (SQRTTWOPI) * invexd2;
			*Teta = -(Strike * Sigma) / (2.0 * SQRTTWOPI * sqrt(NbD)) * invexd2;	// exp(-d1k^2/2) = exp(-d2k^2/2) * Strike / Sigma
		}
	}
	else
	{                                 
		*Price = Strike - Spot + (Spot * nd1k) - (Strike * nd2k);
		if (Hedges)
		{
			*Delta = -1 + nd1k;
			*Gamma = (Strike / (Spot * Spot * srtpi) * invexd2);
			*Vega = (Strike * sqrt(NbD)) / (SQRTTWOPI) * invexd2;
			*Teta = -(Strike * Sigma)/(2.0 * SQRTTWOPI * sqrt(NbD)) * invexd2;
		}    
	}

	return 1.0;
}

double	 BSNnValues(
					bool	Call,			// true = cap, false = floor
					double	Sigma,
					double	Spot,
					double	Strike,
					double	NbD,
					double* Nd1,
					double* Nd2,
					double* nd1,
					double* nd2)
{
	double d1k, d2k, nd1k, nd2k, invexd2;
	double srt, srtpi;

	*Nd1=0.0; *Nd2=0.0; *nd1=0.0; *nd2=0.0;

	// -----------------------------
	// Some basic consitency tests
	if (Sigma <= 0.0)
		return -1.0;
	if (Spot < 0.0)
		return -1.0;
	if (NbD <= 0.0)
		return -1.0;
	// -----------------------------
	
	// Degenerated case
	if (Strike <= 0.0)
	{
		*Nd1 = 1.0;
		*Nd2 = 0.0;
		*nd1 = 0.0;
		*nd1 = 0.0;

		return 1.0;
	}

	srt = Sigma * sqrt(NbD);
	d1k = CalBs(1.0, Spot, Strike, Sigma, NbD);
	d2k = CalBs(-1.0, Spot, Strike, Sigma, NbD);
	srtpi = srt * SQRTTWOPI;

	if (d2k < -15)
	{
		nd2k = 0.0;
		invexd2 = 0.0;
	}
	else if (d2k > 15)
	{
		nd2k = 1.0;
		invexd2 = 0.0;
	}
	else
	{
		nd2k	= cdfNormal(d2k);
		invexd2 = exp(-(d2k * d2k) / 2.0);
	}

	if (d1k < -15)
		nd1k = 0.0;
	else if (d1k > 15)
		nd1k = 1.0;
	else
		nd1k = cdfNormal(d1k);

	if (Call)
	{
		*Nd1 = nd1k;
		*Nd2 = nd2k;
		*nd1 = d1k;
		*nd1 = d2k;
	}
	else
	{                                 
		*Nd1 = 1.0 - nd1k;
		*Nd2 = 1.0 - nd2k;
		*nd1 = -d1k;
		*nd1 = -d2k;
	}

	return 1.0;
}


double StandardGaussianDensity(double x)
{
	if (fabs(x) > 15)
		return 0.0;
	else
		return exp(- 0.5 * x * x) * ONEOVERSQRTTWOPI;
}


double	NormalFormula(
					bool	IsCall,			// true = cap, false = floor
					double	Sigma,
					double	Spot,
					double	Strike,
					double	NbD,
					bool	Hedges,			// true if one wants all greeks to be computed
					double* Price,
					double* Delta,
					double* Gamma,
					double* Vega,
					double* Teta)
{
	double	d1k;
	double	nd1k;
	double	exd2;
	double	srt;

	// RAZ
	*Price	=	0.0;
	*Delta	=	0.0;
	*Gamma	=	0.0;
	*Vega	=	0.0;
	*Teta	=	0.0;

	// -----------------------------
	// Some basic consitency tests
	if (Sigma <= 0.0)			
		return -1.0;

	if (Spot < 0.0) 
		return -1.0;

	if (Strike<0.0)
		return -1.0;

	if (NbD <= 0.0)
		return -1.0;
	// -----------------------------

	srt		=	Sigma * sqrt(NbD);
	d1k		=	(Spot - Strike) / srt;
	exd2	=	exp(-0.5 * d1k * d1k) * ONEOVERSQRTTWOPI;
	nd1k	=	cdfNormal(d1k);
 
	if (IsCall)
	{
		*Price = (Spot - Strike) * nd1k + srt * exd2;
    
		if (Hedges)
		{
			*Delta	=	nd1k;
			*Gamma	=	exd2 / srt;
			*Vega	=	sqrt(NbD) * exd2;
			*Teta	=	0.5 * Sigma * exd2 / sqrt(NbD);
		}
	}
	else
	{                                 
		*Price = Strike - Spot + (Spot - Strike) * nd1k + srt * exd2;
    
		if (Hedges)
		{
			*Delta	=	d1k - 1.0;
			*Gamma	=	exd2 / srt;
			*Vega	=	sqrt(NbD) * exd2;
			*Teta	=	0.5 * Sigma * exd2 / sqrt(NbD);
		}
	}

	return 1.0;
}


// --------------------------------------------------------------------------------------
// BIVARIATE
// --------------------------------------------------------------------------------------

// standard normal density function
double ndf(double t)
{
return 0.398942280401433*exp(-t*t/2);
}

// standard normal cumulative distribution function
double nc(double x)
{
double result;
if (x<-7.)
        result = ndf(x)/sqrt(1.+x*x);
else if (x>7.)
        result = 1. - nc(-x);
else
{
result = 0.2316419;
static double a[5] = {0.31938153,-0.356563782,1.781477937,-1.821255978,1.330274429};
result=1./(1+result*fabs(x));
result=1-ndf(x)*(result*(a[0]+result*(a[1]+result*(a[2]+result*(a[3]+result*a[4])))));
if (x<=0.) result=1.-result;
}
return result;
}


//function needed to calculate two dimensional cumulative distribution function, see Hull
double fxy(double x, double y, double a, double b, double rho)
{
    double a_s;
    double b_s;
        double result;
    a_s = a / sqrt(2 * (1 - rho * rho));
    b_s = b / sqrt(2 * (1 - rho * rho));
    result = exp(a_s * (2 * x - a_s) + b_s * (2 * y - b_s) + 2 * rho * (x - a_s) * (y - b_s));
return result;
}


// function needed to calculate two dimensional cumulative distribution function, see Hull
// this equals ND2 if a and b and rho are all nonpositive, 
// the generalization for the other cases is ND2 below
double Ntwo(double a, double b, double rho)
{
    static double aij[4]={0.325303,
                          0.4211071,
                          0.1334425,
                          0.006374323};
        static double bij[4]={0.1337764,
                          0.6243247,
                          1.3425378,
                          2.2626645};
    int i;
    int j;
        double result;
    result = 0;
        for(i=0;i<=3;i++) 
                {
                        for(j=0;j<=3;j++)
                        {
                                result+=aij[i] * aij[j] * fxy(bij[i], bij[j], a, b, rho); 
                        }
                }
    result = result * sqrt(1 - rho * rho) / MY_PI;
return result;
}
//calculates cumulative distribution function for a bivariate normal distribution
//see John Hull: Options, Futures and Other Derivatives
double ND2(double a, double b, double rho)
{
    double rho1;
    double rho2;
    double denominator;
        double result;

    if (rho > 0.9999)
        result = nc(FMIN(a, b));
    else if (rho < -0.9999)
        result = FMAX(0, nc(a) - nc(-b));
    else
        {
        if (a * b * rho <= 0) 
                {
            if (a <= 0 && b <= 0 && rho <= 0)
                result = Ntwo(a, b, rho);
            else if (a <= 0 && b * rho >= 0)
                result = nc(a) - Ntwo(a, -b, -rho);
            else if (b <= 0 && rho >= 0)
                result = nc(b) - Ntwo(-a, b, -rho);
            else
                result = nc(a) + nc(b) - 1 + Ntwo(-a, -b, rho);
                }
        else
                {
            denominator = sqrt(a * a - 2 * rho * a * b + b * b);
            rho1 = (rho * a - b) * VB_SGN(a) / denominator;
            rho2 = (rho * b - a) * VB_SGN(b) / denominator;
            result = ND2(a, 0, rho1) + ND2(b, 0, rho2) - (1 - VB_SGN(a) * VB_SGN(b)) / 4;
        }
        if (result < 0) result = 0;
        }
        return result;
}


# define MAXIT 100
# define EPS 3.0e-7
# define FPMIN 1.0e-30
# define JMAX 20

double gammln(double xx)
// Returns the value ln|gamma(xx)| for xx > 0.
{
	//Internal arithmetic will be done in double precision, a nicety that you can omit if ve-gure
	//accuracy is good enough.
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
	24.01409824083091,-1.231739572450155,
	0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;

	return -tmp+log(2.5066282746310005*ser/x);
}

double beta(double z, double w)
// Returns the value of the beta function B(z; w).
{
	double gammln(double xx);
	return exp(gammln(z)+gammln(w)-gammln(z+w));
}

double betai(double a, double b, double x)
// Returns the incomplete beta function I x (a; b).
{
	double betacf(double a, double b, double x);
	double gammln(double xx);
	double bt;

	if (x < 0.0 || x > 1.0) throw("Bad x in routine betai");
	if (x == 0.0 || x == 1.0) bt=0.0;
	else															// Factors in front of the continued fraction.
	bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0))										// Use continued fraction directly.
	return bt*betacf(a,b,x)/a;
	else return 1.0-bt*betacf(b,a,1.0-x)/b;							// Use continued fraction after making the sym-
																	// metry transformation. 		
}

double betacf(double a, double b, double x)
// Used by betai: Evaluates continued fraction for incomplete beta function by modied Lentz's method ( x 5.2).
{
	int m,m2;
	double aa,c,d,del,h,qab,qam,qap;
	qab=a+b;								// These q's will be used in factors that occur in the coe.cients (6.4.6).
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;									// First step of Lentz's method.
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=MAXIT;m++) {
	m2=2*m;
	aa=m*(b-m)*x/((qam+m2)*(a+m2));
	d=1.0+aa*d;								// One step (the even one) of the recurrence.
	if (fabs(d) < FPMIN) d=FPMIN;
	c=1.0+aa/c;
	if (fabs(c) < FPMIN) c=FPMIN;
	d=1.0/d;
	h *= d*c;
	aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
	d=1.0+aa*d;								// Next step of the recurrence (the odd one).
	if (fabs(d) < FPMIN) d=FPMIN;
	c=1.0+aa/c;
	if (fabs(c) < FPMIN) c=FPMIN;
	d=1.0/d;
	del=d*c;
	h *= del;
	if (fabs(del-1.0) < EPS) break;			// Are we done?
	}
	if (m > MAXIT) throw("a or b too big, or MAXIT too small in betacf");
	return h;
}

double StudentDensity (int NStudent, double x)
{
	double DNStudent; DNStudent = (double)NStudent;

	if (NStudent != 4)
	{
	return 1./(sqrt(DNStudent)*
			   beta(0.5,0.5*DNStudent)*
			   pow( (1.+x*x/DNStudent), 0.5*(DNStudent+1.) )   
			   );
	}
	else
	{
	return 0.375*pow((1.+0.25*x*x),-2.5);
	}

}

double StudentCum (int NStudent, double x)
{
	double xa,NStHalf;
	xa = (double)NStudent/((double)NStudent + x*x);


	if (NStudent == 4)
	{
			if (x < 0) return 0.5 - sqrt(1.-xa)*(0.5+0.25*xa);
			else	   return 0.5 + sqrt(1.-xa)*(0.5+0.25*xa);
	}
	else
	{
			NStHalf = 0.5 * NStudent;

			if (x > 0)
			{
					return  1.-0.5*betai(NStHalf,0.5,xa);
			}
			else
			{
					return 0.5*betai(NStHalf,0.5,xa);
			}
	}


}	


double InvStudentCum (int NStudent, double p)
{
double xacc = 1.e-4;
double x1 = -6., x2 = 6.;
double df,dx,f,rtn;
double pi2 = 6.283185307;
double ac, z, x;

	if (NStudent == 4)
	{
		if (p < 0.5)
		{
		ac = acos(1.-2.*p);
		z=-2.*cos(0.333333333*(ac-pi2));
		x=2.*sqrt(1./(1.-z*z)-1.);
		return -x;
		}
		else
		{
		ac = acos(2.*p-1.);
		z=-2.*cos(0.3333333333*(ac-pi2));
		x=2.*sqrt(1./(1.-z*z)-1.);
		return +x;
		}
	}
	else 
	{	
		if (p <= StudentCum(NStudent, x1) + EPS) return x1;
		if (p >= StudentCum(NStudent, x2) - EPS) return x2;

		rtn=0.5*(x1+x2);

		for (int j=1;j<=JMAX;j++) 
		{

			f  = StudentCum(NStudent, rtn) - p;
			df = StudentDensity (NStudent, rtn);

			dx=f/df;
			rtn -= dx;
			if ((x1-rtn)*(rtn-x2) < 0.0) throw("Jumped out of brackets in InvStudent");

			if (fabs(dx) < xacc) return rtn;

		}
		throw("Maximum number of iterations exceeded in rtnewt");
		return 0.0;
	}
}
#undef JMAX
#undef MAXIT 
#undef EPS 
#undef FPMIN 



double SplineInterpol(vector<double>& x,vector<double>& y,double value,double smooth,vector<double>& weights2)
{
	double fit=0.;
	
	Nag_Spline spline;

	if (smooth<0.)
// FIXMEFRED: mig.vc8 (28/05/2007 10:25:57):cast
	{e01bac(x.size(),&(*x.begin()),&(*y.begin()),&spline,NAGERR_DEFAULT);}
	else
	{
		double fp;
		Nag_Start start=Nag_Cold;
		Nag_Comm warmstartinf;
		vector<double> weights;
		weights.resize(x.size());

		if (weights2.size()==0)
		for (int i=0;i<x.size();i++){weights[i]=1.;}
		else
		for (int i=0;i<x.size();i++){weights[i]=weights2[i];}

		e02bec(start,x.size(),&(*x.begin()),&(*y.begin()),&(*weights.begin()),smooth,(Integer)(x.size()+4),&fp,&warmstartinf,&spline,NAGERR_DEFAULT);
	}

	e02bbc(value,&fit,&spline,NAGERR_DEFAULT);

	NAG_FREE(spline.lamda);
	NAG_FREE(spline.c);

	return fit;
}


double HermiteInterpol(vector<double>& x,vector<double>& y,double value)
{
	double fit=0.;
	
	vector<double> in;in.resize(1);
	in[0]=value;
	vector<double> out;out.resize(1);

	static NagError fail;
	fail.print = false;
	vector<double> d;d.resize(x.size());

// FIXMEFRED: mig.vc8 (28/05/2007 10:30:12):cast
	e01bec(x.size(),&(*x.begin()),&(*y.begin()),&(*d.begin()),NAGERR_DEFAULT);
	e01bfc(x.size(),&(*x.begin()),&(*y.begin()),&(*d.begin()),1,&(*in.begin()),&(*out.begin()),&fail);

	fit = out[0];

	return fit;
}

double DigitalPriceBS(double& yf,double& spot,double& strike,double& vol,double& vol_epsilon,double& epsilon)
{
	double C_K,C_K_Epsilon,Delta,Gamma,Vega,Teta,output;
	output = BSFormula(true,vol,spot,strike-epsilon/2.,yf,false,&C_K,&Delta,&Gamma,&Vega,&Teta);
	output = BSFormula(true,vol_epsilon,spot,strike+epsilon/2.,yf,false,&C_K_Epsilon,&Delta,&Gamma,&Vega,&Teta);

	double result = (C_K-C_K_Epsilon)/epsilon;

	return result;
}

double CorridorPriceBS(double& yf,
					   double& spot,
					   double& K1,
					   double& K2,
					   double& volK1,
					   double& vol_epsilonk1,
					   double& volK2,
					   double& vol_epsilonk2,
					   double& epsilon)
{
	double NPV = DigitalPriceBS(yf,spot,K1,volK1,vol_epsilonk1,epsilon) - 
				DigitalPriceBS(yf,spot,K2,volK2,vol_epsilonk2,epsilon);

	return NPV;
}


/** 
void GaussLegendre_ComputeAbscissasAndWeights(FILE* fOut, 
																		 const double& x1, 
																		 const double& x2, 
																		 const int& n)
//Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
//arrays x[0..n-1] and w[0..n-1] of length n, containing the abscissas and weights of the Gauss-
//Legendre n-point quadrature formula.
{
	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1;

	if (n <= 0) return;
	double*	x	=	new	double[n];
	double*	w	=	new	double[n];

	//	High precision is a good idea for this routine.	
	m=(n+1)/2;
	//	The roots are symmetric in the interval, so we only have to .nd half of them.
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	for (i=1;i<=m;i++)
	{
		//	Loop over the desired roots.
		z=cos(3.141592654*((double)i-0.25)/((double)n+0.5));
		//	Starting with the above approximation to the ith root, we enter the main loop of
		//	refinement by Newton’s method.
		do
		{
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++)
			{
				//	Loop up the recurrence relation to get the Legendre polynomial evaluated at z.
				p3=p2;
				p2=p1;
				p1=((2.0*(double)j-1.0)*z*p2-((double)j-1.0)*p3)/(double)j;
			}

			//	p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
			//	by a standard relation involving also p2, the polynomial of one lower order.

			pp=((double)n)*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
			//	Newton’s method.
		} while (fabs(z-z1) > EPS_GL);

		x[i-1]=xm-xl*z;
		//	Scale the root to the desired interval,
		x[n-i]=xm+xl*z;
		//	and put in its symmetric counterpart.
		w[i-1]=2.0*xl/((1.0-z*z)*pp*pp);
		//	Compute the weight
		w[n-i]=w[i-1];
		//	and its symmetric counterpart.
	}

	for (i=0;i<n;i++)
		fprintf(fOut,"{%.12lf, %.12lf}\n",x[i],w[i]);

	if (x)
		delete[] x;
	if (w)
		delete[] w;

}


void GaussHermite_Display()
{
	FILE* fOut;

	if ((fOut = fopen("C:\\Credit\\Dump\\GHValues.txt", "w")) == NULL) return;

	fprintf(fOut, " ----------------- Gauss Hermite Values ----------------- \n");

	int	n;

	n=20;
	fprintf(fOut, "\nIntegration Step:\t%ul\n", n);
	GaussHermite_ComputeAbscissasAndWeights(fOut, n);

	n=40;
	fprintf(fOut, "\nIntegration Step:\t%ul\n", n);
	GaussHermite_ComputeAbscissasAndWeights(fOut, n);
	
	n=60;
	fprintf(fOut, "\nIntegration Step:\t%ul\n", n);
	GaussHermite_ComputeAbscissasAndWeights(fOut, n);
	
	n=100;
	fprintf(fOut, "\nIntegration Step:\t%ul\n", n);
	GaussHermite_ComputeAbscissasAndWeights(fOut, n);

	fclose(fOut);
}


void GaussHermite_ComputeAbscissasAndWeights(FILE* fOut, int n)
//Given n, this routine returns arrays x[0..n-1] and w[0..n-1] containing the abscissas and weights
//of the n-point Gauss-Hermite quadrature formula. The largest abscissa is returned in x[1], the
//most negative in x[n].
{
	int i,its,j,m;
	double p1,p2,p3,pp,z,z1;	// High precision is a good idea for this routine.

	if (n <= 0) return;
	double*	x	=	new	double[n];
	double*	w	=	new	double[n];

	m=(n+1)/2;
	// The roots are symmetric about the origin, so we have to .nd only half of them.
	for (i=1;i<=m;i++)
	{			
		//	Loop over the desired roots.
		if (i == 1)
		{
			// Initial guess for the largest root.
			z=sqrt((double)(2.*n+1.))-1.85575*pow((double)(2.*n+1.),-0.16667);
		}
		else if (i == 2)
		{
			// Initial guess for the second largest root.
			z -= 1.14*pow((double)n,0.426)/z;
		}
		else if (i == 3)
		{
			// Initial guess for the third largest root.
			z=1.86*z-0.86*x[0];
		}
		else if (i == 4)
		{
			// Initial guess for the fourth largest root.
			z=1.91*z-0.91*x[1];
		}
		else
		{
			// Initial guess for the other roots.
			z=2.0*z-x[i-3];
		}

		for (its=1;its<=MAXIT_GH;its++)
		{
			// Refinement by Newton’s method.
			p1=PIM4_GH;
			p2=0.0;
			for (j=1;j<=n;j++)
			{
				// Loop up the recurrence relation to get
				// the Hermite polynomial evaluated at z.
				p3=p2;
				p2=p1;
				p1=z*sqrt(2.0/(double)j)*p2-sqrt(((double)(j-1))/(double)j)*p3;
			}
			// p1 is now the desired Hermite polynomial. We next compute pp, its derivative, by
			// the relation (4.5.21) using p2, the polynomial of one lower order.
			pp=sqrt(2.*(double)n)*p2;
			z1=z;
			z=z1-p1/pp;		// Newton’s formula.
			if (fabs(z-z1) <= EPS_GH) break;
		}

		if (its > MAXIT_GH) return;	// too many iterations

		x[i-1]=z;					// Store the root
		x[n-i] = -z;			// and its symmetric counterpart.
		w[i-1]=2.0/(pp*pp);		// Compute the weight
		w[n-i]=w[i-1];			// and its symmetric counterpart.
	}

	for (i=0;i<n;i++)
		fprintf(fOut,"{%.18lf, %.80lf}\n",x[i],w[i]);

	if (x)
		delete[] x;
	if (w)
		delete[] w;

}

void GaussLegendre_Display()
{
	FILE* fOut;

	if ((fOut = fopen("C:\\Credit\\Dump\\GLValues.txt", "w")) == NULL) return;

	fprintf(fOut, " ----------------- Gauss Legendre Values ----------------- \n");

	int	n;

	n=11;
	fprintf(fOut, "\nIntegration Step:\t%ul\n", n);
	GaussLegendre_ComputeAbscissasAndWeights(fOut, 0.0, 1.0, n);

	n=21;
	fprintf(fOut, "\nIntegration Step:\t%ul\n", n);
	GaussLegendre_ComputeAbscissasAndWeights(fOut, 0.0, 1.0, n);
	
	n=51;
	fprintf(fOut, "\nIntegration Step:\t%ul\n", n);
	GaussLegendre_ComputeAbscissasAndWeights(fOut, 0.0, 1.0, n);
	
	n=101;
	fprintf(fOut, "\nIntegration Step:\t%ul\n", n);
	GaussLegendre_ComputeAbscissasAndWeights(fOut, 0.0, 1.0, n);

	fclose(fOut);
}
**/ 

// -------------------------------------------------------------------------------------------
// computation of MAX(X1+X2+...+Xn - K,0.)
// -------------------------------------------------------------------------------------------

double LinearLognProcess::Price(bool call)
{
	double put=0.;
	double spot=0.;
	int k=0;

	for (k=0;k<its_spot.size();k++)
	{spot+=its_spot[k]*its_coef[k];}
	

	if (itsIsComputed) 
	{
		if (call) return itsCallPrice;
		else
		{put = itsCallPrice - (spot-its_strike);
		return (put);
		}
	}

	

	itsIsComputed=true;

	return 0.;
}

double LinearLognProcess::MultivariateDensity(const double& x)
{
	int size = its_spot.size();

	ARM_Matrix M(size,size,its_rho);
	for (int i=0;i<size;i++)
	{M.Elt(i,i)= its_vol[i];}

	double det = M.Det(); 
	double frac = 1./(pow(TWO_PI,size/2)*sqrt(det));

	ARM_Matrix invM;
	invM.Invert(&M,det);

	return 0.;

}


double LinearLognProcess::CumulativeProba(const double& x)
{
return 0.;
}

double LinearLognProcess::BivariateNormalDensity(const double& x,const double& y,
		const double& rho,const double& volx,const double& voly,const double& mx,const double& my)
{
	double density = 1/(TWO_PI*volx*voly*sqrt(1.-rho*rho));
	double internal = (x-mx)*(x-mx)/(volx*volx) + (y-my)*(y-my)/(voly*voly) - 2.*rho*(x-mx)*(y-my)/(volx*volx);
	internal *= 1/(1-rho*rho);
	density *= exp(-0.5*internal);

	return (density);
}

double LinearLognProcess::BivariateNormalDensity_(double x, void* params)
{
	return BivariateNormalDensity(x,its_IntegrationValue,
		its_rho,its_vol[0],its_vol[1],0.,0.);
}


extern "C" void 
BivariateNormalDensity2_(void* Param, double x1, double& res)
{
	// get parameters from stack
	LinearLognProcess*	TheModel;
	TheModel	=	(LinearLognProcess*)((*(AddressVector*)Param)[0]);

	double	its_yt=TheModel->its_yt;
	double	its_strike=TheModel->its_strike;
	double	its_rho=TheModel->its_rho;
	vector<double> its_coef=TheModel->its_coef;
	vector<double> its_spot=TheModel->its_spot;			  
	vector<double> its_vol=TheModel->its_vol;
	double	its_avg1= 0.,its_avg2= 0.;
	double	x2= TheModel->its_IntegrationValue;

//	if (TheModel->TheIntegrator.GetIntegrationType() ==	qGAUSS_HERMITE)
//		x2	*=	SQRT2;

	double inv_phi = x1 /* *sqrt(2.)*(its_vol[0])*/ *sqrt(1.-its_rho*its_rho)+
					its_avg1+its_rho /* * (its_vol[0]/its_vol[1]) */ *(x2-its_avg2);

	double X1 = its_coef[0]*its_spot[0]*exp(-0.5*its_vol[0]*its_vol[0]*its_yt
											+its_vol[0]*inv_phi*sqrt(its_yt));
	double X2 = its_coef[1]*its_spot[1]*exp(-0.5*its_vol[1]*its_vol[1]*its_yt
											+its_vol[1]*x2*sqrt(its_yt));

	res = 0.;
	res = its_coef[0]*X1 + its_coef[1]*X2 - its_strike;
//	res /= sqrt(PI);
}

double LinearLognProcess::IntBND(const double& borne_inf_x,
								 const double& borne_sup_x,
								 const double& fixed_value_y)
{
	double probacond = 0.;

	if CHECK_EQUAL(borne_inf_x,borne_sup_x)
		return 0.;

	//parametres d'integration
	AddressVector	TheParameterVector;
	TheParameterVector.Append(this);

	its_IntegrationValue=fixed_value_y;
	TheIntegrator.Integrate(borne_inf_x, borne_sup_x, BivariateNormalDensity2_, &TheParameterVector, probacond); 
	return (probacond);
}

extern "C" void IntBND_(void* Param, double x, double& res)
{
	// get parameters from stack
	LinearLognProcess*	TheModel;
	TheModel	=	(LinearLognProcess*)((*(AddressVector*)Param)[0]);

	double	its_yt=TheModel->its_yt;
	double	its_strike=TheModel->its_strike;
	double	its_rho=TheModel->its_rho;
	vector<double> its_coef=TheModel->its_coef;
	vector<double> its_spot=TheModel->its_spot;			  
	vector<double> its_vol=TheModel->its_vol;
	double	its_avg1= 0.,its_avg2= 0.;

	double x2 = /* sqrt(2.) *its_vol[1]* */ x + its_avg2;

//	if (TheModel->TheIntegrator.GetIntegrationType() ==	qGAUSS_HERMITE)
//		x2	*=	SQRT2;

	double Lambda = its_strike - its_coef[1]*its_spot[1]*exp(-0.5*its_vol[1]*its_vol[1]*its_yt+its_vol[1]*x2*sqrt(its_yt));
	Lambda /= its_coef[0]*its_spot[0];
	Lambda = log(Lambda) + 0.5*its_vol[0]*its_vol[0]*its_yt;
	Lambda /= its_vol[0]*sqrt(its_yt);

	double Phi_Lambda = Lambda - its_avg1 - its_rho /* *(its_vol[0]/its_vol[1]) */ *(x2-its_avg2);
	Phi_Lambda /= sqrt(2.) /* *(its_vol[0])*/*sqrt(1-its_rho*its_rho);

	double borne_sup = 10.;
	double borne_inf = MIN(MAX(Phi_Lambda,-10.),borne_sup);

	res	=	TheModel->IntBND(borne_inf,borne_sup,x2);
//	res /= sqrt(PI);
}


double LinearLognProcess::IntBND2(const double& borne_inf_y,const double& borne_sup_y)
{

//	return IntBND2_Single(borne_inf_y,borne_sup_y);

	double proba = 0.;

	//parametres d'integration
	AddressVector	TheParameterVector;
	TheParameterVector.Append(this);
	TheIntegrator.SetIntegrationStep(its_intstep);
	//TheIntegrator.SetIntegrationType(qGAUSS_LEGENDRE);
	TheIntegrator.SetIntegrationType(qGAUSS_HERMITE);

	TheIntegrator.Integrate(borne_inf_y, borne_sup_y, IntBND_, &TheParameterVector, proba); 
	return (proba);
}

// ----------------------------------------------------------------------
// Test single
// ----------------------------------------------------------------------

extern "C" void IntBND_Single(void* Param, double x, double& res)
{
	// get parameters from stack
	LinearLognProcess*	TheModel;
	TheModel	=	(LinearLognProcess*)((*(AddressVector*)Param)[0]);

	double	its_yt=TheModel->its_yt;
	double	its_strike=TheModel->its_strike;
	double	its_rho=TheModel->its_rho;
	vector<double> its_coef=TheModel->its_coef;
	vector<double> its_spot=TheModel->its_spot;			  
	vector<double> its_vol=TheModel->its_vol;
	double	its_avg1= 0.,its_avg2= 0.;

	double x2 = x;//sqrt(2.) /* *its_vol[1]*/ *x + its_avg2;

	res = its_spot[0]*exp(-0.5*its_vol[0]*its_vol[0]*its_yt+its_vol[0]*sqrt(its_yt)*x2);
	res -= its_strike;

}

double LinearLognProcess::IntBND2_Single(const double& borne_inf_y,const double& borne_sup_y)
{
	double proba = 0.;

	//parametres d'integration
	AddressVector	TheParameterVector;
	TheParameterVector.Append(this);
	TheIntegrator.SetIntegrationStep(its_intstep);
	//TheIntegrator.SetIntegrationType(qGAUSS_LEGENDRE);
	TheIntegrator.SetIntegrationType(qGAUSS_HERMITE);

	double Lambda = log(its_strike/its_spot[0])+0.5*its_vol[0]*its_vol[0]*its_yt;
	Lambda /= its_vol[0]*sqrt(its_yt);
	//Lambda /= sqrt(2.);

	double borne_sup = 100.;
	double borne_inf = Lambda; //MIN(MAX(Lambda,-10.),borne_sup);

	TheIntegrator.Integrate(borne_inf, borne_sup, IntBND_Single, &TheParameterVector, proba); 
	return (proba);
}
