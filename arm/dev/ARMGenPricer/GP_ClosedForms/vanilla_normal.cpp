/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file vanilla_normal.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpbase/utilityport.h"
#include "gpclosedforms/vanilla_normal.h"

#include <cmath>

#include "gpnumlib/gaussiananalytics.h"
#include "gpbase/numericconstant.h"
#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/normal.h"
#include "gpclosedforms/lambert_function.h"
#include "gpclosedforms/vanilla_bs.h"

/// gpnumlib
#include "gpnumlib/numfunction.h"
#include "gpnumlib/dichotomy.h"
#include "gpnumlib/newtonraphson.h"


#include <glob/expt.h>

#include <limits>



CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_IMPLICITVOL_ACCURACY 1e-308


//////////////////////////////////////////////////////////////////////////////////////////
///
///   Pricing Functions :
///
///   because we are dealing with normal models, the volatility information needs
///   to be introduces as a standard deviation .
///   therfore to convert from a lognormal volatility information :
///   standardDev=vol*F
///
//////////////////////////////////////////////////////////////////////////////////////////


double VanillaDigitalOption_N(double F,
        double standardDev,
        double k,
        double t,
        int callput)
{
	double v1 = standardDev*sqrt(t);
	double d1 = (F-k)/v1;

	return (1.-callput)/2.+callput*NormalCDF(d1);
}

double VegaVanillaDigitalOption_N(double F,
        double standardDev,
        double k,
        double t ,
        int callOrPut)
{	
	double v1 = standardDev*sqrt(t);
	double d1 = (F-k)/v1;

	return callOrPut*d1/standardDev*NormalPDF(d1);
}

double VanillaIndexPaying_DigitalCall_N(double F,
        double standardDev,
        double k,
        double t )
{	
	double v1 = standardDev*sqrt(t);
	double d1 = (F-k)/v1;
	return F*(NormalCDF(d1)
		+v1*NormalPDF(d1));

}


double VanillaOption_N( double F, 
        double stdDev,
        double k, 
        double t,  
        int callOrPut )
{
	if (callOrPut==0)
	{
		if (t <= ARM_NumericConstants::ARM_TOLERENCE)
			return 0.0;

		double v1 = stdDev*sqrt(t);
		double d1 = (F-k)/v1;
		if(fabs(F-k) <ARM_NumericConstants::ARM_TOLERENCE)
		{
			return 2.*(v1*ARM_NumericConstants::ARM_INVSQRT2PI);
		}
		else if(callOrPut*(F-k) < ARM_NumericConstants::ARM_TOLERENCE) 
		{ 
			if(v1 < ARM_NumericConstants::ARM_TOLERENCE)
				return 0.0;
		}
		else if(callOrPut*(F-k) > ARM_NumericConstants::ARM_TOLERENCE)
		{ 
			if(v1 < ARM_NumericConstants::ARM_TOLERENCE)
				return 0.0;
		}

		return (F-k)*(2*NormalCDF(d1)-1)+2*v1*NormalPDF(d1);
	}
	else
	{
		if (t <= ARM_NumericConstants::ARM_TOLERENCE)
			return CC_Max(callOrPut * (F-k), 0.0);

		double v1 = stdDev*sqrt(t);
		double d1 = (F-k)/v1;
		if(fabs(F-k) <ARM_NumericConstants::ARM_TOLERENCE)
		{
			return (v1*ARM_NumericConstants::ARM_INVSQRT2PI);
		}
		else if(callOrPut*(F-k) < ARM_NumericConstants::ARM_TOLERENCE) 
		{ 
			if(v1 < ARM_NumericConstants::ARM_TOLERENCE)
				return 0.0;
		}
		else if(callOrPut*(F-k) > ARM_NumericConstants::ARM_TOLERENCE)
		{ 
			if(v1 < ARM_NumericConstants::ARM_TOLERENCE)
				return callOrPut*(F-k);
		}

		return callOrPut*(F-k)*NormalCDF(callOrPut*d1)+
				   v1*NormalPDF(callOrPut*d1);
	}
}

double DeltaVanillaOption_N( double F,
        double stdDev,
        double k,
        double t,  
        int callOrPut )
{
	if (callOrPut==0)
	{
		double v = stdDev*sqrt(t);	
		double d=(F-k)/v;

		return 2.*NormalCDF(d)-1.;
	}
	else
	{
		double v = stdDev*sqrt(t);	
		double d=(F-k)/v;

		return callOrPut*NormalCDF(callOrPut*d);
	}
}

double GammaVanillaOption_N( double F,
        double stdDev,
        double k,
        double t,  
        int callOrPut )
{
	if (callOrPut==0)
	{
		double v = stdDev*sqrt(t);	
		double d=(F-k)/v;
		
		double gamma = 2.0/v*NormalPDF(d);
		return gamma;

	}
	else
	{
		double v = stdDev*sqrt(t);	
		double d=(F-k)/v;
		
		double gamma = 1.0/v*NormalPDF(callOrPut*d);
		return gamma;
	}
}

double VegaVanillaOption_N( double F,
        double stdDev,
        double k,
        double t,  
        int callOrPut )
{
	if (callOrPut==0)
	{
		double v = stdDev*sqrt(t);	
		double h=(F-k)/v;

		return 2.*exp(-h*h/2.0)*sqrt(t)*ARM_NumericConstants::ARM_INVSQRT2PI;
	}
	else
	{
		double v = stdDev*sqrt(t);	
		double h=(F-k)/v;

		return exp(-h*h/2.0)*sqrt(t)*ARM_NumericConstants::ARM_INVSQRT2PI;
	}
}


//
//  Implied normal volatility for a vanilla call/put 
//	the problem being very specificn, we do not used a generic solver...
//
double VanillaImpliedVol_N( double f, 
							double target, 
							double k, 
							double t, 
							int callput,
							double* guess,
							bool*	success)
{
	
	// initial guess for implied vol 
	double v0;
	if (guess)
		v0 = *guess;
	else
	{
		if(fabs(f-k)/f <1.e-2) 
			v0 = target * sqrt(2. * ARM_NumericConstants::ARM_PI / t);
		else 
		{
			v0 = sqrt(2. * fabs(log(f/k)) / t);
			v0 = v0 * (f-k) / log(f/k) * (1. - 1./24. * v0 * v0 * t);
		}
	}

	
	/// option too OTM
	if (fabs(target)<ARM_NumericConstants::ARM_TOLERENCE)
		return v0;

	/// option too ITM
	if (fabs(callput*(f-k) - target)<ARM_NumericConstants::ARM_TOLERENCE)
		return v0;
	
	
	/// let's be a little tolerant
	if ( callput*(f-k) > target )
	{
		if ( (callput*(f-k) - target) / target < 1.0e-4 )
			return v0;
		else
			ARM_THROW(  ERR_INVALID_ARGUMENT, "VanillaImpliedVol_N : invalid target price." );
	}

		
	// newton raphson
	vector<double> init (12);			
	
	init[0] = v0;	
	init[1] = 2.*v0;
	init[2] = 5.*v0;
	init[3] = 10.*v0;
	init[4] = 100.*v0;
	init[5] = 0.1*v0;
	init[6] = f*0.1;
	init[7] = f;
	init[8] = f*10;
	init[9] = 1000.*v0;
	init[10]= target * 0.1;
	init[11]= target;

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
		func  = VanillaOption_N(f, v, k, t, callput) - target;
		dfunc = VegaVanillaOption_N(f, v, k, t, callput);

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
					ARM_THROW(  ERR_INVALID_ARGUMENT, "VanillaImpliedVol_N : all initial values failed." );
			}
		}
		p++;
	}

	if (success)
		*success = true;
	return v;

}


////////////////////////////////////////////////////////////////////////
///
///
///			Nouvelle Implementation de Normal Vanilla Implied vol
///
///
////////////////////////////////////////////////////////////////////////

double VanillaCall_N_Phi( double m,
        double TotalstdDev)
{
	double d1 = m/TotalstdDev;
	return m*NormalCDF(d1)+TotalstdDev*NormalPDF(d1);
}

double VanillaCall_N_Phi_Derivative( double m,
        double TotalstdDev)
{
	double d1 = m/TotalstdDev;
	return NormalPDF(d1);
}

	void VanillaCall_N_funcd(const double x,const double param,const double objective,double &fn, double &df)
	{
		double valphi=VanillaCall_N_Phi(param,x);
        fn=valphi-objective;
        df=VanillaCall_N_Phi_Derivative(param,x);
	} 
	

	double VanillaCall_N_findroot(void funcd(const double, const double, const double,double &, double &), const double x1, const double x2,
		const double xacc,const double param,const double objective)
	{
		const int MAXIT=200;
		int j;
		double df,dx,dxold,f,fh,fl,temp,xh,xl,rts;
		fl=0.0;fh=0.0;
		
		funcd(x1,param,objective,fl,df);
		funcd(x2,param,objective,fh,df);
		if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"VanillaCall_N_findroot: impossible to braket the solution");
		if (fabs(fl) <= ARM_CF_IMPLICITVOL_ACCURACY) 
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"VanillaCall_N_findroot: volatility to low or moneyness too far from 1 ");;
		}
		if (fabs(fh) <= ARM_CF_IMPLICITVOL_ACCURACY) 
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"VanillaCall_N_findroot: volatility too high (> 1000 %)");;
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
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"VanillaCall_N_findroot:Maximum number of iterations exceeded in rtsafe");
		return 0.0;
	}
	
	double Asymptotic_VanillaCall_N_BSphi_Inverse(double m,double objective, double accuracy)
	{
		return m/sqrt(3.* lambertfunc(pow(m/objective,2./3.)/(3.*pow(ARM_NumericConstants::ARM_PI*2.,1./3.))));
		
	}



	double VanillaCall_N_BSphi_Inverse(double param,double objective, double accuracy)
	{
        return VanillaCall_N_findroot(VanillaCall_N_funcd,0.00001,10.0,accuracy,param,objective);
	}


////////////////////////////////////////////////////////////////////////
///
///
///			Nouvelle Implementation de Normal Digital Implied vol
///
///
////////////////////////////////////////////////////////////////////////




double DigitalCall_N_ImpliedVol(double f,double k, double opt,int callput)
{
	return (f-k)/NormalCDFInverse((opt-(1.-callput)/2.)/callput);
}




///********************************************************************
///////////////////////////////////////////////////////////////////////
///
///   Analytics for Double Corridors pricing
///   ======================================
///	
///		Payoff = K x 1(Libor<U) x 1(Spread>D)
///
///////////////////////////////////////////////////////////////////////
///********************************************************************

/// -----------------------------------------------------------------
///  BivariateExpectations_N
/// ------------------------
///		This function computes E{ X^p Y^q 1(X>x)1(Y>y) }
///		where:
///			p, q in {0,1}
///			X = N(EX, sigmaX)
///			Y = N(EY, sigmaY)
///			cor(X,Y) = rho
/// -----------------------------------------------------------------
void BivariateExpectations_N
					 	  (	double x, double y, 
							double EX, double sigmaX, 
							double EY, double sigmaY,
							double rho,
							/// results
							double& E00, double& E10, double& E01, double& E11)
{
	double std_x = - (x - EX) / sigmaX;
	double std_y = - (y - EY) / sigmaY;

	double std_E00 = NormalCDF(std_x, std_y, rho);
	double std_E10 = Normal_X_Expectation(std_x, std_y, rho);
	double std_E01 = Normal_Y_Expectation(std_x, std_y, rho);
	double std_E11 = Normal_XY_Expectation(std_x, std_y, rho);

	E00 = std_E00;
	E10 = EX * std_E00 - sigmaX * std_E10;
	E01 = EY * std_E00 - sigmaY * std_E01;
	E11 = EX * EY * std_E00 
		- EY * sigmaX * std_E10
		- EX * sigmaY * std_E01
		+ sigmaX * sigmaY * std_E11;
}

/// -----------------------------------------------------------------
///  DoubleCall_N
/// ---------------
///		This function computes E{ (X(T) - Kx)^+ (Y(T)- Ky)^+}
///		where:
///			X(T) = N(X0, volX*sqrt(T))
///			Y(T) = N(Y0, volY*sqrt(T))
///			cor(X,Y) = correl
/// -----------------------------------------------------------------

/// -----------------------------------
/// direct integration of double call
double DoubleCall_N (double maturity, 
					 double Kx, double Ky, 
					 double X0, double volX, double Y0, double volY, double correl)
{
	double sqrtT  = sqrt(maturity);
	double sigmaX = volX * sqrtT;
	double sigmaY = volY * sqrtT;

	int N		= 250; // nb of intervals
	double xmin = (Kx - X0) / sigmaX ;
	double xmax = 7.0;
	double dx   = (xmax - xmin) / double(N);
	double x    = xmin + 0.5 * dx;
	
	double dblcall (0.0);
	double sqrtunmrho2 = sqrt(1.0 - correl*correl);
	double func;
	
	for (int i(0); i<N; i++)
	{
		func = (X0 + sigmaX * x - Kx);
		func *= exp(-0.5 * x * x);
		func *= VanillaOption_N(Y0 + correl * x * sigmaY, sigmaY * sqrtunmrho2, Ky, 1.0, 1);
		func *= dx;
		
		dblcall += func;

		x += dx;
	}

	dblcall *= ARM_NumericConstants::ARM_INVSQRT2PI;

	return dblcall;
}

/// old method using O.Croissant analytics
double DoubleCall_N_old (double maturity, 
					 double Kx, double Ky, 
					 double X0, double volX, double Y0, double volY, double correl)
{
	double sqrtT  = sqrt(maturity);
	double sigmaX = volX * sqrtT;
	double sigmaY = volY * sqrtT;

	double E00, E10, E01, E11;
	
	BivariateExpectations_N (Kx, Ky, X0, sigmaX, Y0, sigmaY, correl, E00, E10, E01, E11);

	return E11 - Ky * E10 - Kx * E01 + Kx * Ky * E00;
}

/// -----------------------------------------------------------------
///  DoubleDigital_N
/// ------------------
///		This function computes E{ 1(X< Kx or X > Kx) 1(Y<Ky or Y > Ky) } using call spread
///		decomposition
///
/// -----------------------------------------------------------------
double DoubleDigital_N_1 (double maturity, 
						double Kx, double spread_x,
						double Ky, double spread_y,
						double X0, double volXplus, double volXminus,
						double Y0, double volYplus, double volYminus, 
						double correl, int XCap, int YCap)
{
	double KxMinus = Kx - spread_x;
	double KxPlus  = Kx + spread_x;
	double KyMinus = Ky - spread_y;
	double KyPlus  = Ky + spread_y;

	double CallYminus = VanillaOption_N(Y0, volYminus, KyMinus, maturity, K_CALL);
	double CallYplus  = VanillaOption_N(Y0, volYplus,  KyPlus,  maturity, K_CALL);

	double CallXminus = VanillaOption_N(X0, volXminus, KxMinus, maturity, K_CALL);
	double CallXplus  = VanillaOption_N(X0, volXplus,  KxPlus,  maturity, K_CALL);
	
	double den = 4.0 * spread_x * spread_y ;

	/// digital 1(X>Kx) 1(Y>Ky)
	double result =  DoubleCall_N(maturity, KxMinus, KyMinus, X0, volXminus, Y0, volYminus,  correl) / den
				   - DoubleCall_N(maturity, KxMinus, KyPlus,  X0, volXminus, Y0, volYplus,   correl) / den
				   - DoubleCall_N(maturity, KxPlus,  KyMinus, X0, volXplus,  Y0, volYminus,  correl) / den
				   + DoubleCall_N(maturity, KxPlus,  KyPlus,  X0, volXplus,  Y0, volYplus,   correl) / den;

	if( (XCap != K_CAP) && (YCap != K_CAP) )
		result += 1 - (CallYminus - CallYplus) / (2.*spread_y) - (CallXminus - CallXplus) / (2.*spread_x);
	else if( XCap != K_CAP )
		result = (CallYminus - CallYplus) / (2.*spread_y) - result;
	else if( YCap != K_CAP )
		result = (CallXminus - CallXplus) / (2.*spread_x) - result;
	
	return result;
}

/// --------------------------------------
/// New model with vol up / down states
///
class DigitalAsFuncOfStdevUp : public ARM_GP::UnaryFunc<double,double> 
{
private:
	double	itsFwd, itsStrike, itsStdevDown, itsProbaDown;
public: 
		DigitalAsFuncOfStdevUp(double fwd, double strike, double stdevDown, double probaDown) 
			:	itsFwd(fwd),
				itsStrike(strike),
				itsStdevDown(stdevDown),
				itsProbaDown(probaDown) {}
		
		inline virtual double operator() (double stdevUp) const
		{
			return	(1.0-itsProbaDown) * NormalCDF( (itsStrike - itsFwd) / stdevUp ) 
				+   itsProbaDown	   * NormalCDF( (itsStrike - itsFwd) / itsStdevDown ) ;
		}
};

double DoubleDigital_N_2 (double maturity, 
						double Kx, double spread_x,
						double Ky, double spread_y,
						double X0, double volXplus, double volXminus,
						double Y0, double volYplus, double volYminus, 
						double correl, int XCap, int YCap)
{
	double sqrtT   = sqrt(maturity);
	double KxMinus = Kx - spread_x;
	double KxPlus  = Kx + spread_x;
	double KyMinus = Ky - spread_y;
	double KyPlus  = Ky + spread_y;

	double CallXminus = VanillaOption_N(X0, volXminus, KxMinus, maturity, K_CALL);
	double CallXplus  = VanillaOption_N(X0, volXplus,  KxPlus,  maturity, K_CALL);

	double CallYminus = VanillaOption_N(Y0, volYminus, KyMinus, maturity, K_CALL);
	double CallYplus  = VanillaOption_N(Y0, volYplus,  KyPlus,  maturity, K_CALL);

	double DigitalX = 1.0 - (CallXminus - CallXplus) /  (2.0 * spread_x);
	double DigitalY = 1.0 - (CallYminus - CallYplus) /  (2.0 * spread_y);

	double stdevXdown = 0.5 * (volXplus + volXminus) * sqrtT;
	double stdevYdown = 0.5 * (volYplus + volYminus) * sqrtT;

	/// --------------------------
	/// Idea : find stdevXup and stdevYup so that:
	/// DigitalX = 0.5 * ( NormalCDF( (Kx - X0) / stdevXup ) +  NormalCDF( (Kx - X0) / stdevXdown ) ) ;
	/// DigitalY = 0.5 * ( NormalCDF( (Ky - Y0) / stdevYup ) +  NormalCDF( (Ky - Y0) / stdevYdown ) ) ;
	double stdevXup, stdevYup;
	
	double DEFAULT_PRECISION = 1.0e-10;
	int	   MAX_ITER		     = 100;
	double NSTDEV_INTERP	 = 0.20;
	double tol = 1.0e-8;

	bool KxTooCloseToATM = ( fabs(X0 - Kx) < NSTDEV_INTERP * stdevXdown ) ;
	bool KyTooCloseToATM = ( fabs(Y0 - Ky) < NSTDEV_INTERP * stdevYdown ) ;

	/// a bit of recursivity
	if ( KxTooCloseToATM || KyTooCloseToATM )
	{
		double pourEtrePeinard = 1e-10;

		/// Kx is too close from ATM
		if (KxTooCloseToATM && !KyTooCloseToATM)
		{
			double KxDown    = X0 - NSTDEV_INTERP * stdevXdown - pourEtrePeinard;
			double KxUp      = X0 + NSTDEV_INTERP * stdevXdown + pourEtrePeinard;
			double priceDown = DoubleDigital_N(maturity, KxDown, spread_x, Ky, spread_y, X0, volXplus, volXminus, Y0, volYplus, volYminus, correl, XCap, YCap);
			double priceUp   = DoubleDigital_N(maturity, KxUp,   spread_x, Ky, spread_y, X0, volXplus, volXminus, Y0, volYplus, volYminus, correl, XCap, YCap);
			double price     = ( (Kx - KxDown) * priceUp + (KxUp - Kx) * priceDown ) / (KxUp - KxDown);
			return price;
		}
		/// Ky is too close from ATM
		else if (!KxTooCloseToATM && KyTooCloseToATM)
		{
			double KyDown    = Y0 - NSTDEV_INTERP * stdevYdown - pourEtrePeinard;
			double KyUp      = Y0 + NSTDEV_INTERP * stdevYdown + pourEtrePeinard;
			double priceDown = DoubleDigital_N(maturity, Kx, spread_x, KyDown, spread_y, X0, volXplus, volXminus, Y0, volYplus, volYminus, correl, XCap, YCap);
			double priceUp   = DoubleDigital_N(maturity, Kx, spread_x, KyUp,   spread_y, X0, volXplus, volXminus, Y0, volYplus, volYminus, correl, XCap, YCap);
			double price     = ( (Ky - KyDown) * priceUp + (KyUp - Ky) * priceDown ) / (KyUp - KyDown);
			return price;
		}
		/// both Kx & Ky are too close from ATM 
		else
		{
			double KxDown = X0 - NSTDEV_INTERP * stdevXdown - pourEtrePeinard;
			double KxUp   = X0 + NSTDEV_INTERP * stdevXdown + pourEtrePeinard;
			double KyDown = Y0 - NSTDEV_INTERP * stdevYdown - pourEtrePeinard;
			double KyUp   = Y0 + NSTDEV_INTERP * stdevYdown + pourEtrePeinard;
			double priceDownDown = DoubleDigital_N(maturity, KxDown, spread_x, KyDown, spread_y, X0, volXplus, volXminus, Y0, volYplus, volYminus, correl, XCap, YCap);
			double priceUpUp     = DoubleDigital_N(maturity, KxUp,   spread_x, KyUp,   spread_y, X0, volXplus, volXminus, Y0, volYplus, volYminus, correl, XCap, YCap);
			double priceDownUp   = DoubleDigital_N(maturity, KxDown, spread_x, KyUp,   spread_y, X0, volXplus, volXminus, Y0, volYplus, volYminus, correl, XCap, YCap);
			double priceUpDown   = DoubleDigital_N(maturity, KxUp,   spread_x, KyDown, spread_y, X0, volXplus, volXminus, Y0, volYplus, volYminus, correl, XCap, YCap);
			double priceDown = ( (Kx - KxDown) * priceUpDown + (KxUp - Kx) * priceDownDown ) / (KxUp - KxDown); ;
			double priceUp   = ( (Kx - KxDown) * priceUpUp   + (KxUp - Kx) * priceDownUp )   / (KxUp - KxDown); ;
			double price     = ( (Ky - KyDown) * priceUp + (KyUp - Ky) * priceDown) / (KyUp - KyDown);
		}
	}


	/// -- find stdevXup
	if (fabs(DigitalX)>tol && fabs(1.0-DigitalX)>tol)
	{
		DigitalAsFuncOfStdevUp funcX(X0, Kx, stdevXdown, 0.5);
		UnaryFuncWithNumDerivative<double> funcWithDerivativeX(funcX);
		T_SmoothNewtonRaphsonSolver<UnaryFuncWithNumDerivative<double> > solverX(funcWithDerivativeX, DigitalX, DEFAULT_PRECISION, DEFAULT_PRECISION, MAX_ITER);
		solverX.setInitialGuess(stdevXdown);
		stdevXup = solverX.Solve();
	}
	else
		stdevXup = stdevXdown;

	
	/// -- find stdevYup
	if (fabs(DigitalY)>tol && fabs(1.0-DigitalY)>tol)
	{
		DigitalAsFuncOfStdevUp funcY(Y0, Ky, stdevYdown, 0.5);
		UnaryFuncWithNumDerivative<double> funcWithDerivativeY(funcY);
		T_SmoothNewtonRaphsonSolver<UnaryFuncWithNumDerivative<double> > solverY(funcWithDerivativeY, DigitalY, DEFAULT_PRECISION, DEFAULT_PRECISION, MAX_ITER);
		solverY.setInitialGuess(stdevYdown);
		stdevYup = solverY.Solve();
	}
	else
		stdevYup = stdevYdown;
	
	/// -- define probas for the 4 states of the world
	double p_upup     = 0.25; /// to be interfaced if necessary
	double p_downdown = p_upup;
	double p_updown   = 0.5 * (1.0 - p_upup - p_downdown);
	double p_downup   = p_updown;

	double price = 0.0;
	
	/// up up case
	price += p_upup     * NormalCDF((Kx - X0)/stdevXup, (Ky - Y0)/stdevYup, correl);
	/// down down case
	price += p_downdown * NormalCDF((Kx - X0)/stdevXdown, (Ky - Y0)/stdevYdown, correl);
	/// up down case
	price += p_updown   * NormalCDF((Kx - X0)/stdevXup, (Ky - Y0)/stdevYdown, correl);
	/// down up case
	price += p_downup   * NormalCDF((Kx - X0)/stdevXdown, (Ky - Y0)/stdevYup, correl);

	if( XCap == K_CAP && YCap == K_CAP)
		price += 1 - DigitalX - DigitalY;
	else if ( XCap == K_CAP && YCap == K_FLOOR)
		price = DigitalY - price;
	else if ( XCap == K_FLOOR && YCap == K_CAP)
		price = DigitalX - price;

	return price;
}





/// --------------------------------------
/// New model with vol up / down states
///
class DigitalAsFuncOfLognormStdevUp : public ARM_GP::UnaryFunc<double,double> 
{
private:
	double	itsFwd, itsShift, itsStrike, itsStdevDown, itsProbaDown, itsDdown;
public: 
		DigitalAsFuncOfLognormStdevUp(double fwd, double shift, double strike, double stdevDown, double probaDown) 
			:	itsFwd(fwd),
				itsShift(shift),
				itsStrike(strike),
				itsStdevDown(stdevDown),
				itsProbaDown(probaDown) 
		{
			itsDdown  = log( (itsStrike+itsShift)/(itsFwd+itsShift) ) + 0.5 * itsStdevDown * itsStdevDown;
			itsDdown /= itsStdevDown;
		}
		
		inline virtual double operator() (double stdevUp) const
		{
			double d_up;
			d_up  = log( (itsStrike+itsShift)/(itsFwd+itsShift) ) + 0.5 * stdevUp * stdevUp;
			d_up /= stdevUp;
			
			return	(1.0-itsProbaDown) * NormalCDF( d_up ) 
				+   itsProbaDown	   * NormalCDF( itsDdown ) ;
		}
};

double DoubleDigital_N_3 (double maturity, 
						double Kx, double spread_x,
						double Ky, double spread_y,
						double X0, double volXplus, double volXminus,
						double Y0, double volYplus, double volYminus, 
						double correl, int XCap, int YCap)
{
	double sqrtT   = sqrt(maturity);
	double KxMinus = Kx - spread_x;
	double KxPlus  = Kx + spread_x;
	double KyMinus = Ky - spread_y;
	double KyPlus  = Ky + spread_y;

	double CallXminus = VanillaOption_N(X0, volXminus, KxMinus, maturity, K_CALL);
	double CallXplus  = VanillaOption_N(X0, volXplus,  KxPlus,  maturity, K_CALL);

	double CallYminus = VanillaOption_N(Y0, volYminus, KyMinus, maturity, K_CALL);
	double CallYplus  = VanillaOption_N(Y0, volYplus,  KyPlus,  maturity, K_CALL);

	double DigitalX = 1.0 - (CallXminus - CallXplus) /  (2.0 * spread_x);
	double DigitalY = 1.0 - (CallYminus - CallYplus) /  (2.0 * spread_y);

	double normVolXdown = 0.5 * (volXplus + volXminus) ;
	double normVolYdown = 0.5 * (volYplus + volYminus) ;

	double CallX = VanillaOption_N(X0, normVolXdown, Kx, maturity, K_CALL);
	double CallY = VanillaOption_N(Y0, normVolYdown, Ky, maturity, K_CALL);

	double nstdevShift = 2.0;
	double shiftX = - Kx + nstdevShift * normVolXdown * sqrtT;
	double shiftY = - Ky + nstdevShift * normVolYdown * sqrtT;

	double guessX = normVolXdown / (X0 + shiftX);
	double guessY = normVolYdown / (Y0 + shiftY);
	double stdevXdown = VanillaImpliedVol_BS(X0+shiftX, Kx+shiftX, maturity, CallX, K_CALL, &guessX);
	double stdevYdown = VanillaImpliedVol_BS(Y0+shiftY, Ky+shiftY, maturity, CallY, K_CALL, &guessY);


	/// --------------------------
	/// Idea : find stdevXup and stdevYup so that:
	/// DigitalX = 0.5 * ( NormalCDF( (Kx - X0) / stdevXup ) +  NormalCDF( (Kx - X0) / stdevXdown ) ) ;
	/// DigitalY = 0.5 * ( NormalCDF( (Ky - Y0) / stdevYup ) +  NormalCDF( (Ky - Y0) / stdevYdown ) ) ;
	double stdevXup, stdevYup;
	
	double DEFAULT_PRECISION = 1.0e-10;
	int	   MAX_ITER		     = 100;
	double tol = 1.0e-8;
	

	/// -- find stdevXup
	if (fabs(DigitalX)>tol && fabs(1.0-DigitalX)>tol)
	{
		DigitalAsFuncOfLognormStdevUp funcX(X0, shiftX, Kx, stdevXdown, 0.5);
		UnaryFuncWithNumDerivative<double> funcWithDerivativeX(funcX);
		T_SmoothNewtonRaphsonSolver<UnaryFuncWithNumDerivative<double> > solverX(funcWithDerivativeX, DigitalX, DEFAULT_PRECISION, DEFAULT_PRECISION, MAX_ITER);
		solverX.setInitialGuess(stdevXdown);
		stdevXup = solverX.Solve();
	}
	else
		stdevXup = stdevXdown;

	
	/// -- find stdevYup
	if (fabs(DigitalY)>tol && fabs(1.0-DigitalY)>tol)
	{
		DigitalAsFuncOfLognormStdevUp funcY(Y0, shiftY, Ky, stdevYdown, 0.5);
		UnaryFuncWithNumDerivative<double> funcWithDerivativeY(funcY);
		T_SmoothNewtonRaphsonSolver<UnaryFuncWithNumDerivative<double> > solverY(funcWithDerivativeY, DigitalY, DEFAULT_PRECISION, DEFAULT_PRECISION, MAX_ITER);
		solverY.setInitialGuess(stdevYdown);
		stdevYup = solverY.Solve();
	}
	else
		stdevYup = stdevYdown;
	
	/// -- define probas for the 4 states of the world
	double p_upup     = 0.25; /// to be interfaced if necessary
	double p_downdown = p_upup;
	double p_updown   = 0.5 * (1.0 - p_upup - p_downdown);
	double p_downup   = p_updown;

	double dxup   = ( log( (Kx+shiftX)/(X0+shiftX) ) + 0.5 * stdevXup * stdevXup ) / stdevXup;
	double dxdown = ( log( (Kx+shiftX)/(X0+shiftX) ) + 0.5 * stdevXdown * stdevXdown ) / stdevXdown;
	double dyup   = ( log( (Ky+shiftY)/(Y0+shiftY) ) + 0.5 * stdevYup * stdevYup ) / stdevYup;
	double dydown = ( log( (Ky+shiftY)/(Y0+shiftY) ) + 0.5 * stdevYdown * stdevYdown ) / stdevYdown;

	double price = 0.0;
	
	/// up up case
	price += p_upup     * NormalCDF(dxup,   dyup,   correl);
	/// down down case
	price += p_downdown * NormalCDF(dxdown, dydown, correl);
	/// up down case
	price += p_updown   * NormalCDF(dxup,   dydown, correl);
	/// down up case
	price += p_downup   * NormalCDF(dxdown, dyup,   correl);

	if( XCap == K_CAP && YCap == K_CAP)
		price += 1 - DigitalX - DigitalY;
	else if ( XCap == K_CAP && YCap == K_FLOOR)
		price = DigitalY - price;
	else if ( XCap == K_FLOOR && YCap == K_CAP)
		price = DigitalX - price;

	return price;
}



/// --------------------------------------
///
/// Barrier Shift Method ...
///
///
class DigitalAsFuncOfStrike : public ARM_GP::UnaryFunc<double,double> 
{
private:
	double	itsFwd, itsStdev;
	bool itsIsLognormal;
public: 
		DigitalAsFuncOfStrike(double fwd, double stdev, bool isLognormal = false) 
			:	itsFwd(fwd),
				itsStdev(stdev),
				itsIsLognormal(isLognormal) {}
					
		inline virtual double operator() (double strike) const
		{
			if (!itsIsLognormal)
				return NormalCDF( (strike - itsFwd) / itsStdev ) ;
			else
			{
				double d; 
				d  = log( strike/itsFwd ) + 0.5 * itsStdev * itsStdev;
				d /= itsStdev;
				return NormalCDF( d ) ;
			}
		}
};



double DoubleDigital_N_4 (double maturity, 
						double Kx, double spread_x,
						double Ky, double spread_y,
						double X0, double volXplus, double volXminus,
						double Y0, double volYplus, double volYminus, 
						double correl, int XCap, int YCap,
						double* volX,double* volY,
						double* KxAdj,double* KyAdj,
						double* probaXY,double* probaX,double* probaY)
{
	double sqrtT   = sqrt(maturity);
	double KxMinus = Kx - spread_x;
	double KxPlus  = Kx + spread_x;
	double KyMinus = Ky - spread_y;
	double KyPlus  = Ky + spread_y;

	double CallXminus = VanillaOption_N(X0, volXminus, KxMinus, maturity, K_CALL);
	double CallXplus  = VanillaOption_N(X0, volXplus,  KxPlus,  maturity, K_CALL);

	double CallYminus = VanillaOption_N(Y0, volYminus, KyMinus, maturity, K_CALL);
	double CallYplus  = VanillaOption_N(Y0, volYplus,  KyPlus,  maturity, K_CALL);

	double DigitalX = 1.0 - (CallXminus - CallXplus) /  (2.0 * spread_x);
	double DigitalY = 1.0 - (CallYminus - CallYplus) /  (2.0 * spread_y);

	double meanVolX = 0.5 * (volXplus + volXminus);
	double meanVolY = 0.5 * (volYplus + volYminus);
	double stdevX = meanVolX * sqrtT;
	double stdevY = meanVolY * sqrtT;
	
	double KxShifted,KyShifted;
	
	double DEFAULT_PRECISION = 1.0e-10;
	int	   MAX_ITER		     = 100;
	double TOL				 = 1.0e-8;

	/// -- find stdevXup
	if (DigitalX>TOL && 1.0-DigitalX>TOL)
	{
		KxShifted = X0 + stdevX * NormalCDFInverse(DigitalX);
/*		DigitalAsFuncOfStrike funcX(X0, stdevX);
		UnaryFuncWithNumDerivative<double> funcWithDerivativeX(funcX);
		T_NewtonRaphsonSolver<UnaryFuncWithNumDerivative<double> > solverX(funcWithDerivativeX, DigitalX, DEFAULT_PRECISION, DEFAULT_PRECISION, MAX_ITER);
		solverX.setInitialGuess(Kx);
		KxShifted = solverX.Solve();
*/	}
	else
		KxShifted = Kx;

	
	/// -- find stdevYup
	if (DigitalY>TOL && 1.0-DigitalY>TOL)
	{
		KyShifted = Y0 + stdevY * NormalCDFInverse(DigitalY);
/*		DigitalAsFuncOfStrike funcY(Y0, stdevY);
		UnaryFuncWithNumDerivative<double> funcWithDerivativeY(funcY);
		T_NewtonRaphsonSolver<UnaryFuncWithNumDerivative<double> > solverY(funcWithDerivativeY, DigitalY, DEFAULT_PRECISION, DEFAULT_PRECISION, MAX_ITER);
		solverY.setInitialGuess(Ky);
		KyShifted = solverY.Solve();
*/	}
	else
		KyShifted = Ky;

	double prXY = NormalCDF((KxShifted - X0)/stdevX, (KyShifted - Y0)/stdevY, correl);

	double price = prXY; // FLOOR & FLOOR case

	if( XCap == K_CAP && YCap == K_CAP)
		price += 1 - DigitalX - DigitalY;
	else if ( XCap == K_CAP && YCap == K_FLOOR)
		price = DigitalY - price;
	else if ( XCap == K_FLOOR && YCap == K_CAP)
		price = DigitalX - price;

	if(volX && volY)
	{
		*volX = meanVolX;
		*volY = meanVolY;
	}
	if(KxAdj && KyAdj)
	{
		*KxAdj = KxShifted;
		*KyAdj = KyShifted;
	}
	if(probaX)
		*probaX = DigitalX;
	if(probaY)
		*probaY = DigitalY;
	if(probaXY)
		*probaXY = prXY;

	return price;
}


///
/// switch from one method to another (default = Barrier Shift)
///
int DoubleDigitalMethod = 4;

double DoubleDigital_N(double maturity, 
						double Kx, double spread_x,
						double Ky, double spread_y,
						double X0, double volXplus, double volXminus,
						double Y0, double volYplus, double volYminus, 
						double correl, int XCap, int YCap,
						double* volX, double* volY,
						double* KxShifted,double* KyShifted,
						double* probaXY,double* probaX,double* probaY)
{
	if (DoubleDigitalMethod == 1)
		return DoubleDigital_N_1(maturity, Kx, spread_x, Ky, spread_y, X0, volXplus, volXminus, Y0, volYplus, volYminus, correl, XCap, YCap);
	else if (DoubleDigitalMethod == 2)
		return DoubleDigital_N_2(maturity, Kx, spread_x, Ky, spread_y, X0, volXplus, volXminus, Y0, volYplus, volYminus, correl, XCap, YCap);
	else if (DoubleDigitalMethod == 3)
		return DoubleDigital_N_3(maturity, Kx, spread_x, Ky, spread_y, X0, volXplus, volXminus, Y0, volYplus, volYminus, correl, XCap, YCap);
	else if (DoubleDigitalMethod == 4)
		return DoubleDigital_N_4(maturity, Kx, spread_x, Ky, spread_y, X0, volXplus, volXminus, Y0, volYplus, volYminus, correl, XCap, YCap,
								 volX, volY, KxShifted, KyShifted, probaXY, probaX, probaY);

	return 0;
}

CC_END_NAMESPACE()
 

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


