

/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  Functions for   BlackSholes 
 *
 *	\file vanilla_bs.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

#include <glob/firsttoinc.h>

#include "gpnumlib/gaussiananalytics.h"
#include "gpbase/numericconstant.h"
#include "gpbase/utilityport.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/vanille_bs_formula.h"
#include "gpclosedforms/lambert_function.h"
#include "gpinfra/argconvdefault.h"
#include "gpclosedforms/normal.h"
#include "gpclosedforms/inverse.h"



#include <cmath>
#include <glob/expt.h>

#include <limits>


CC_BEGIN_NAMESPACE( ARM )

#define ARM_CF_IMPLICITVOL_ACCURACY 1e-15

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///    Utilitary functions
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double BS(double f, 
		  double k, 
		  double t, 
		  double v, 
		  int callput,
		  double* delta,
		  double* gamma,
		  double* vega,
		  double* theta)
{
	// some validation
	// maybe we should do the following tests:
	// f >  0
	// v >= 0
	// t >= 0

	// maturity <= 0 or vol <= 0 --> intrinsic value
	if (t < ARM_NumericConstants::ARM_TOLERENCE || v < ARM_NumericConstants::ARM_TOLERENCE)
	{
		if (delta)  *delta	= (double)callput;
		if (gamma)	*gamma	= 0.0;
		if (vega)	*vega	= 0.0;
		if (theta)  *theta	= 0.0;
		
		return CC_Max (callput * (f-k), 0.0);
	}

	double sqrtt	= sqrt(t);
	double vsqrt	= v * sqrtt;
	double d1		= log(f/k)/vsqrt + 0.5 * vsqrt;
	double d2		= d1 - vsqrt;
	double Nd1		= ARM_GaussianAnalytics::cdfNormal(d1);
	double Nd2		= ARM_GaussianAnalytics::cdfNormal(d2);
	double call		= f * Nd1 - k * Nd2;
	
	if (delta)	*delta = callput * Nd1 ;
	if (gamma)  *gamma = ARM_GaussianAnalytics::dNormal(d1) / (f*vsqrt);
	if (vega)	*vega  = f * sqrtt * ARM_GaussianAnalytics::dNormal(d1);
	if (theta)  *theta = -f * ARM_GaussianAnalytics::dNormal(d1) * v / (2.*sqrtt);

	if (callput == K_PUT)
		return call - f + k;
	else
		return call;
}

double VanillaImpliedVol_BS (double f, double k, double t, double target, int callput, double* guess, bool* success)
{
	// initial guess for implied vol 
	double v0;
	if (guess)
		v0 = *guess;
	else
	{		
		if (fabs(f-k)/f < 1e-2) 
			v0 = target * sqrt(2. * ARM_NumericConstants::ARM_PI / t) / f;
		else
			v0 = sqrt(2.*fabs(log(f/k))/t);
	}

	/// option too OTM
	if (fabs(target)<ARM_NumericConstants::ARM_TOLERENCE)
		return v0;

	/// option too ITM
	if (fabs(callput*(f-k) - target)<ARM_NumericConstants::ARM_TOLERENCE)
		return v0;
		
	// newton raphson
	vector<double> init (12);			
	
	init[0] = v0;	
	init[1] = 1.5*v0;
	init[2] = 3.*v0;
	init[3] = 8.*v0;
	init[4] = 0.1*v0;
	init[5] = 0.01*v0;
	init[6] = 0.10;
	init[7] = 0.01;
	init[8] = 0.20;
	init[9] = 0.50;
	init[10]= 1.0;
	init[11]= 5.0;

	double			_MIN (0.0);		
	double			_MAX (1.e12);	
	unsigned int	nbitermax (50);		
	double			tol_fct (1e-9);		
	double			tol_x   (1e-9);		
	
	//------------------------------------------------------------------------
	//--------- NEWTON RAPHSON -----------------------------------------------
	//------------------------------------------------------------------------
	double func, dfunc;
	int				p	= 0;
	double			v	= init[0];
	unsigned int	cpt	= 1;
		
	while (true) 
	{	
		func  = BS(f, k, t, v, callput, NULL, NULL, &dfunc, NULL) - target;
		
		if ( fabs(func/target)<tol_fct ) 
		{
			break ;
		}

		v = v - func/dfunc;
		
		if ( fabs(func/dfunc/v) < tol_x  &&  fabs(dfunc)>DBL_EPSILON) 
		{			
			break;  // sortie de la boucle while
		}
		
		if (v<=_MIN || v>=_MAX || fabs(v)<DBL_EPSILON || fabs(dfunc)<DBL_EPSILON || p==nbitermax )
		{	
			if(cpt<init.size())
			{
				v = init[cpt];
				cpt++;
				p = -1;
			}
			else 
			{
				if (success)
				{
					*success = false;
					return 0.0;
				}
				else
					ARM_THROW(  ERR_INVALID_ARGUMENT, "VanillaImpliedVol_BS : all initial values failed." );
			}
				
		}
		p++;
	}

	if (success)
		*success = true;
	return v;
	
}

double BSImpliedVol(double f, double k, double t, double target, int callPut, bool * success)
{
	ARM_ImpliedVolBS func(f,k,t,callPut);
	return func.vol(target,success);
}

double VanillaOptionFromCall(double forward,double k,double call,int CallOrPutFlag)
{
	switch (CallOrPutFlag) 
	{
	case K_CALL :
		{
			return call;
			break;
		}
	case K_PUT :
		{
			return call - (forward -k);
			break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"VanillaOptionFromCall : CallOrPutFlag  bad input :");
			break;
		}
	}
	
}



/*----------------------------------------------------------------------------*
   NAME    BlackSholes_Formula

   Compute value of european Vanilla option using Black & Scholes formula
   in its more general formulation

   Input:

   forward : the expected value of the underlying under the risk neutral 
     measure, assuming that the underlying has a lognormal distribution of
	 at reset time under the given probability..
	 
   totalvolatility: the volatility of underlying at reset time ..
     In cas of a process with a deterministic volatility, this is 
     the square root of the integral of the square of the 
     deterministic volatility of the process  between the start of 
	 the option an the reset time of the option.

   bondprice : the price of a bond with notionnal=1, maturiring at the payment
     time of the option.(payment time >= reset time)


   Output:

   the Black Scholes price for an equity type option,

*----------------------------------------------------------------------------*/
double BlackSholes_Formula(	double forward,
							double totalvolatility,
							double bondprice,
							double strike,
							double CallPut)
{
	if(totalvolatility<0.0) return 0.0;
	/// handle the case of negative strike
	if(strike<=K_DOUBLE_TOL)
		return CallPut == K_CALL? bondprice*(forward-strike): 0.0;

	/// standard Black Scholes case
    double d1		= (log(forward/strike))/totalvolatility+0.5*totalvolatility ;
    double  value	= CallPut*bondprice*(forward*ARM_GaussianAnalytics::cdfNormal(CallPut*d1)
		-strike*ARM_GaussianAnalytics::cdfNormal(CallPut*(d1-totalvolatility)));
    return value;
}

double BlackSholes_Formula(	double forward,
							double volatility,
							double bondprice,
							double strike,
							double maturity,
							double CallPut)
{
	if(maturity<0.0) return 0.0;
	double totalvolatility=volatility*sqrt(maturity);
	/// handle the case of negative strike
	if(strike<=K_DOUBLE_TOL)
		return CallPut == K_CALL? bondprice*(forward-strike): 0.0;

	/// standard Black Scholes case
    double d1		= (log(forward/strike))/totalvolatility+0.5*totalvolatility ;
    double  value	= CallPut*bondprice*(forward*NormalCDF(CallPut*d1)-strike*NormalCDF(CallPut*(d1-totalvolatility)));
    return value;
}


double BlackSholes_Derivative_1(	double forward,			//delta
									double totalvolatility,
									double bondprice,
									double strike,
									double CallPut)
{
	if(totalvolatility<0.0) return 0.0;
	if(strike<=K_DOUBLE_TOL)
		strike = K_DOUBLE_TOL;
    double d1	= (log(forward/strike))/totalvolatility+0.5*totalvolatility ;
    double value= bondprice*CallPut*ARM_GaussianAnalytics::cdfNormal(CallPut*d1);
    return value;
}


double BlackSholes_Derivative_2(	double forward,			//vega total
									double totalvolatility,
									double bondprice,
									double strike,
									double CallPut)
{
	if(totalvolatility<0.0) return 0.0;
	if(strike<=K_DOUBLE_TOL)
		strike = K_DOUBLE_TOL;
    double d1		= (log(forward/strike))/totalvolatility+0.5*totalvolatility ;
    double value	= bondprice*forward*exp(-d1*d1/2)*ARM_NumericConstants::ARM_INVSQRT2PI;
    return value;
}

double BlackSholes_Derivative_3(	double forward,			// bond price total
									double totalvolatility,
									double bondprice,
									double strike,
									double CallPut)
{
	if(totalvolatility<0.0) return 0.0;
	/// handle the case of negative strike
	if(strike<=K_DOUBLE_TOL)
		return CallPut == K_CALL? (forward-strike): 0.0;

	/// standard Black Scholes case
    double d1		= (log(forward/strike))/totalvolatility+0.5*totalvolatility ;
    double  value	= CallPut*(forward*ARM_GaussianAnalytics::cdfNormal(CallPut*d1)
		-strike*ARM_GaussianAnalytics::cdfNormal(CallPut*(d1-totalvolatility)));
    return value;
}

double BlackSholes_Derivative_4(	double forward,			// strike
									double totalvolatility,
									double bondprice,
									double strike,
									double CallPut)
{
	if(totalvolatility<0.0) return 0.0;
	if(strike<=K_DOUBLE_TOL)
		strike = K_DOUBLE_TOL;
    double d2	= (log(forward/strike))/totalvolatility-0.5*totalvolatility ;
    double value= -bondprice*CallPut*ARM_GaussianAnalytics::cdfNormal(CallPut*d2);
    return value;
}

double BlackSholes_Derivative_1_1(	double forward,		//gamma
									double totalvolatility,
									double bondprice,
									double strike,
									double CallPut)
{
	if(totalvolatility<0.0) return 0.0;
	if(strike<=K_DOUBLE_TOL)
		strike = K_DOUBLE_TOL;
    double d1 = (log(forward/strike))/totalvolatility+0.5*totalvolatility ;
    double value = exp(-(d1*d1)/2.)*ARM_NumericConstants::ARM_INVSQRT2PI/(forward*totalvolatility);
    return value;
}



///////////////////////////////////////////
/// Digital part
///////////////////////////////////////////

double DigitalBlackSholes_Formula(	double forward,
							double totalvolatility,
							double bondprice,
							double strike,
							double CallPut)
{
	/// handle the case of negative strike
	if(strike<=K_DOUBLE_TOL)
		return CallPut == K_CALL? bondprice: 0.0;

	/// standard Black Scholes case
    double d2		= (log(forward/strike))/totalvolatility-0.5*totalvolatility ;
    double  value	= bondprice*ARM_GaussianAnalytics::cdfNormal(CallPut*d2);
    return value;
}

double DigitalBlackSholesSmooth_Formula(	double forward,
											double totalvolatility_plus,
											double totalvolatility_moins,
											double bondprice,
											double strike,
											double CallPut,
											double spread)
{
	/// handle the case of zero spread
	if(fabs(spread)<=K_DOUBLE_TOL)
	{
		double totalvolatility = 0.5*(totalvolatility_plus+totalvolatility_moins);
		return DigitalBlackSholes_Formula(forward,totalvolatility, bondprice, strike, CallPut);
	}

	/// Call spread formula 
	double sp = fabs(spread);
	double  BS_strike_Plus_sp	= BlackSholes_Formula(forward,totalvolatility_plus, 1.0, strike+sp, CallPut);
	double  BS_strike_moins_sp	= BlackSholes_Formula(forward,totalvolatility_moins, 1.0, strike-sp, CallPut);
	double value = bondprice*(BS_strike_moins_sp-BS_strike_Plus_sp)/(2*CallPut*sp);
    return value;
}



 /////////////////////////////////////////////////////////////////////////////////////////////////////
 ///
 ///  Calcul de vol implicite
 ///
 ///
 /////////////////////////////////////////////////////////////////////////////////////////////////////
double CumulativeNormalDerivative(double x)
{ 
	return ARM_NumericConstants::ARM_INVSQRT2PI * exp( 0.5 * x * x); 
}



double ARM_CF_BS_Formula::callimplicit_totalvolatility(double forward, double discount,
													   double strike, double CallPut,double opt,double accuracy)
{
	if(CallPut==1)
	{
		return  BSphiInverse(log(forward/strike),opt/(strike*discount),accuracy) ;
	}
	else
	{
		double call=opt+(forward-strike)*discount;
		return  BSphiInverse(log(forward/strike),call/(strike*discount),accuracy) ;
	}
}

double ARM_CF_BS_Formula::callimplicit_totalvolatility1(double forward, double discount,
													   double strike, double CallPut,double opt,double accuracy)
{
	if(CallPut==1)
	{
		return  BSphiInverse1(log(forward/strike),opt/(strike*discount),accuracy) ;
	}
	else
	{
		double call=opt+(forward-strike)*discount;
		return  BSphiInverse1(log(forward/strike),call/(strike*discount),accuracy) ;
	}
}


/// Formule plus juste et plus universelle:

double ARM_CF_BS_Formula::callimplicit_totalvolatility2(double forward, double discount,
													   double strike, double CallPut,double opt,double accuracy)
{
	if(CallPut==1)
	{
		return  BSphiInverse2(log(forward/strike),opt/(strike*discount),accuracy) ;
	}
	else
	{
		double call=opt+(forward-strike)*discount;
		return  BSphiInverse2(log(forward/strike),call/(strike*discount),accuracy) ;
	}
}



double ARM_CF_BS_Formula::callimplicit_totalvolatility_DerOpt(double forward, double discount,
															  double strike, double CallPut,double opt,double accuracy)
{
	double totalvol=ARM_CF_BS_Formula::callimplicit_totalvolatility(forward,discount,strike,CallPut,opt,accuracy);
	return 1./BlackSholes_Derivative_2(forward,totalvol,discount,strike,CallPut)	;
}


double ARM_CF_BS_Formula::callimplicit_totalvolatility_DerForward(double forward, double discount,
															  double strike, double CallPut,double opt,double accuracy)
{
	double totalvol=ARM_CF_BS_Formula::callimplicit_totalvolatility(forward,discount,strike,CallPut,opt,accuracy);
	return -BlackSholes_Derivative_1(forward,totalvol,discount,strike,CallPut)/
			BlackSholes_Derivative_2(forward,totalvol,discount,strike,CallPut);
}

double ARM_CF_BS_Formula::callimplicit_totalvolatility_DerStrike(double forward, double discount,
															  double strike, double CallPut,double opt,double accuracy)
{
	double totalvol=ARM_CF_BS_Formula::callimplicit_totalvolatility(forward,discount,strike,CallPut,opt,accuracy);
	return -BlackSholes_Derivative_3(forward,totalvol,discount,strike,CallPut)/
			BlackSholes_Derivative_2(forward,totalvol,discount,strike,CallPut);
}


	/// auxiliary functions:
	
	/// Asymptotic value of phi when u tend toward 0+
	double BSphi_positive(double m,double u)
	{
		double m16=16.+4.*m*m/(u*u)+u*u;
		double sq1=sqrt(m16-4.*m);
		double sq2=sqrt(m16+4.*m);
		double argexp=(-2.*m+u*u)/u;
		argexp=-argexp*argexp/2.;

		return -1.+exp(m)+2.*exp(argexp)*ARM_NumericConstants::ARM_SQRT_2_DIVIDED_BY_SQRT_PI*u*(1./(2.*m+u*(-u+sq1))-1./(2.*m+u*(u+sq2)));
	}

	/// Asymptotic value of phi when u tend toward 0-
	double BSphi_negative(double m,double u)
	{
		double m16=16.+4.*m*m/(u*u)+u*u;
		double sq1=sqrt(m16-4.*m);
		double sq2=sqrt(m16+4.*m);
		double argexp=(-2.*m+u*u)/u;
		argexp=-argexp*argexp/2.;

		return 2.*exp(argexp)*ARM_NumericConstants::ARM_SQRT_2_DIVIDED_BY_SQRT_PI*u*(1./(2.*m-u*(u+sq1))-1./(2.*m+u*(u-sq2)));
	}

	double BSphi_Asymptotic(double m,double u)
	{
		if(m>=0) 
		{
			return BSphi_positive(m,u);
		}
		else 
		{
			return BSphi_negative(m,u);
		}
	}


	double ARM_CF_BS_Formula::BSphi(double m, double u)
	{
		if(fabs(m/u)>5.)
		{
			return BSphi_Asymptotic(m,u);
		}
		else
		{
			double up=u/2.0;
			double um=m/u;
			return ARM_GaussianAnalytics::cdfNormal(um+up)*exp(m) -ARM_GaussianAnalytics::cdfNormal(um-up);
		}
	}



	double ARM_CF_BS_Formula::BSphiDerivative(double m, double u)
	{
		double up=u/2.0;
		double um=m/u;
		double u2m=um/u;
		return CumulativeNormalDerivative(um+up)*exp(m)*(0.5-u2m)+CumulativeNormalDerivative(um-up)*(u2m+0.5);
	}

	
	double ARM_CF_BS_Formula::BSphiDerivative2(double m, double u)
	{
		double arg=-m/u+u/2.0;
		if (fabs(arg)>706.0 ) 
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_BS_Formula::BSphiDerivative: arg=-m/u+u/2 too big! ");	
		}
		return CumulativeNormalDerivative(arg);
	}



	void ARM_CF_BS_Formula::funcd(const double x,const double param,const double objective,double &fn, double &df)
	{
		double valphi=BSphi(param,x);
        fn=valphi-objective;
        df=BSphiDerivative(param,x);
	} 
	void ARM_CF_BS_Formula::funcd2(const double x,const double param,const double objective,double &fn, double &df)
	{
		double valphi=BSphi(param,x);
        fn=valphi-objective;
		if((fabs(-param/x+x/2.0))>706.0)
		{
			df=1.0e+300;
		}
		else
		{
			df=BSphiDerivative2(param,x);
		}
	} 
	

	double ARM_CF_BS_Formula::findroot(void funcd(const double, const double, const double,double &, double &), const double x1, const double x2,
		const double xacc,const double param,const double objective)
	{
		const int MAXIT=200;
		int j;
		double df,dx,dxold,f,fh,fl,temp,xh,xl,rts;
		fl=0.;fh=0.;
		
		funcd(x1,param,objective,fl,df);
		funcd(x2,param,objective,fh,df);
		if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_BS_Formula::findroot: impossible to braket the solution");
		if (fabs(fl) <= ARM_CF_IMPLICITVOL_ACCURACY) 
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_BS_Formula::findroot: volatility to low or moneyness too far from 1 ");;
		}
		if (fabs(fh) <= ARM_CF_IMPLICITVOL_ACCURACY) 
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_BS_Formula::findroot: volatility too high (> 1000 %)");;
		}
		if (fl < 0.0) {
			xl=x1;
			xh=x2;
		} else {
			xh=x1;
			xl=x2;
		}
		rts=0.5*(x1+x2);
		dxold=fabs(x2-x1);
		dx=dxold;
		funcd(rts,param,objective,f,df);
		for (j=0;j<MAXIT;j++) {
			if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)
				|| (fabs(2.0*f) > fabs(dxold*df))) {
				dxold=dx;
				dx=0.5*(xh-xl);
				rts=xl+dx;
				if (xl == rts) return rts;
			} else {
				dxold=dx;
				dx=f/df;
				temp=rts;
				rts -= dx;
				if (temp == rts) return rts;
			}
			if (fabs(dx) < xacc) return rts;
			funcd(rts,param,objective,f,df);
			if (f < 0.0)
				xl=rts;
			else
				xh=rts;
		}
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"findroot:Maximum number of iterations exceeded in rtsafe");
		return 0.0;
	}
	
	
	double ARM_CF_BS_Formula::BSphiInverse(double param,double objective, double accuracy)
	{
        return findroot(funcd,0.0001,10.0,accuracy,param,objective);
	}
	
	double ARM_CF_BS_Formula::BSphiInverse1(double param,double objective, double accuracy)
	{
        return findroot(funcd2,0.0001,10.0,accuracy,param,objective);
	}
	
	
	double ARM_CF_BS_Formula::BSphiInverse2(double param,double objective, double accuracy)
	{
        class FunctionToInverse : public DoubleToDoubleFunc
		{
		public:
			double p;
			FunctionToInverse(double p0):p(p0) {}
			double operator() (double x) const
			{
				return ARM_CF_BS_Formula::BSphi(p,x);
			}
		};
		FunctionToInverse xx(param);	
		return Inverse(xx,Inverse::ALWAYSPOSITIVE)(objective,0.2,0.05,1e-12);	
	}
	

	/////////////////////////////////////////////////////////////////////////////////////////
	///
	///				Asymptotic formula
	///
	////////////////////////////////////////////////////////////////////////////////////////

	double Asymptotic_BlackSholesTimeValue(	double forward,
										double strike,
										double totalvolatility)
	{
		if(totalvolatility<0.0) return 0.0;
		double m=log(forward/strike);
		return strike*(sqrt(forward/(strike*ARM_NumericConstants::ARM_PI*2.))*exp(-m*m/(2.*totalvolatility*totalvolatility))*
			totalvolatility*totalvolatility*totalvolatility/(m*m));
	}


	double Asymptotic_BlackSholesTimeValue_ImplicitVol(
										double forward,
										double strike,
										double optprice)
	{
			double m=fabs(log(forward/strike));
			return m/sqrt(3.* lambertfunc(pow(m/optprice,2./3.)/(3.*pow(ARM_NumericConstants::ARM_PI*2./(forward*strike),1./3.))));
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	///
	///				Ratio and Product Options
	///
	////////////////////////////////////////////////////////////////////////////////////////
	
	double 	BS_RatioOption(	double S1,double Mu1,double Sigma1,
		double S2,double Mu2,double Sigma2,double Rho,double K,double T,int CallPut)
	{
		if(T<0.0) return 0.0;
		double arg=(log(K/S1)-(Mu1-3.*Sigma1*Sigma1/2.+Sigma1*Sigma1*Rho)*T)/(Sigma1*sqrt(T));
		switch (CallPut) 
		{
		case K_CALL :
			{
				return S2/S1*exp((Mu2-Mu1-Sigma1*Sigma2*Rho)*T)*(1.-NormalCDF(arg));
				break;
			}
		case K_PUT :
			{
				return S2/S1*exp((Mu2-Mu1-Sigma1*Sigma2*Rho)*T)*(NormalCDF(arg));
				break;
			}
		default :
			{	
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BS_RatioOption : CallPut  bad input :");
				break;
			}
		}
		
	}
	
	double 	BS_ProductOption(	double S1,double Mu1,double Sigma1,
		double S2,double Mu2,double Sigma2,double Rho,double K,double T,int CallPut)
	{
		if(T<0.0) return 0.0;
		double arg=(log(K/S1)-(Mu1+Sigma1*Sigma1/2.+Sigma1*Sigma1*Rho)*T)/(Sigma1*sqrt(T));
		switch (CallPut) 
		{
		case K_CALL :
			{
				return S2*S1*exp((Mu2+Mu1+Sigma1*Sigma2*Rho)*T)*(1.-NormalCDF(arg));
				break;
			}
		case K_PUT :
			{
				return S2*S1*exp((Mu2+Mu1+Sigma1*Sigma2*Rho)*T)*(NormalCDF(arg));
				break;
			}
		default :
			{	
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BS_ProductOption : CallPut  bad input :");
				break;
			}
		}
		
	}
	/////////////////////////////////////////////////////////////////////////////////////////
	///
	///				Payoff L2 si L1>K
	///
	////////////////////////////////////////////////////////////////////////////////////////

	double 	BS_PayFloat_DigitalOption(	double S1,double Mu1,double Sigma1,
		double S2,double Mu2,double Sigma2,double Rho,double K,double T,int CallPut)
	{
		if(T<0.0) return 0.0;
		double arg=(log(K/S1)-(Mu1+Sigma1*Sigma1/2.+Sigma1*Sigma1*Rho)*T)/(Sigma1*sqrt(T));
		switch (CallPut) 
		{
		case K_CALL :
			{
				return S2*exp(Mu2*T)*(1.-NormalCDF(arg));
				break;
			}
		case K_PUT :
			{
				return S2*exp(Mu2*T)*(NormalCDF(arg));
				break;
			}
		default :
			{	
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BS_BiOption : CallPut  bad input :");
				break;
			}
		}		
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	///
	///				Payoff forward2 si forward1>K if Call 
	///
	///////////////////////////////////////////////////////////////////////////////////////	
	double BS_PayFloat_DigitalOptionSmooth_Formula(	double forward1,
													double totalvolatility1,
													double forward2,
													double totalvolatility2,
													double correl,
													double bondprice,
													double strike,
													double CallPut,
													double spread)
	{
		if((totalvolatility1<0.0)||(totalvolatility2<0.0)) return 0.0;
		/// change numeraire of forward1
		double forward1_tilde = forward1*exp(correl*totalvolatility1*totalvolatility2);

		//digitalPrice under the new numeraire 
		double digitalPrice = DigitalBlackSholesSmooth_Formula(forward1_tilde,totalvolatility1,totalvolatility1,1.0,strike,CallPut,spread);
		double value = forward2*bondprice*digitalPrice;
		return value;
	}




CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
