
/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 11/15/2005
 *
 *  Formulaes of SABR Model in case Beta<>1
 *
 *	\file SABR_Analytics.cpp
 *
 *  \brief
 *
 *	\author  M.Abdelmoumni
 *	\version 1.0
 *	\date November 2005
 */

#ifndef unix
#include "gpclosedforms/sabrbdiff1.h"
#else
#include "sabrbdiff1.h"
#endif
#include <glob/expt.h>

#include <math.h>

#ifdef unix
#include <ieeefp.h>
#else
#include <float.h>
#define finite _finite
#endif

#define ARM_GP_CF_MAX_STRIKE 10000
#define ARM_GP_CF_MAX_VOL 10000



#define ARM_SABR_SERIE_USAGE_LIMIT 0.1 


#define SABR_SQR(x) ((x)*(x))


static int IS_ONE_TEST(double x) // at 1e-6
{
    long truncated_x = long(floor(x*100000.0));

    return( truncated_x == 100000L );
}


static int IS_ZERO_TEST(double x) // at 1e-8
{
		return fabs(x) < 1e-8  ?  1 : 0 ;
}


double CptSABR_implicit_vol_direct_Series(double f,double K, double T,
                                          double alpha, double beta,
                                          double rho, double nu,
                                          int flag)
{
    double res;


#ifdef SABR_ORDER1

    double a0;
    double a1;


    a0 = 2.0*alpha*pow(3.0, 0.5)*pow(K, beta)
                *pow(-(T*alpha*alpha*pow(-1.0+beta, 2.0)
                *pow(K, 2.0*beta))-6.0*alpha*beta*nu*rho*T*pow(K, 1.0+beta)
                +K*K*(12.0+T*nu*nu*(-2.0+3.0*rho*rho)),
                -0.5);

    a1 = pow(3.0, 0.5)*(-6.0*(-1.0+beta)*beta*nu*rho*T
                *alpha*alpha*pow(K, 2.0*beta)
                +nu*rho*K*K*(-12.0+T*nu*nu*(5.0-6.0*rho*rho))
                +alpha*(-1.0+beta)*pow(K, 1.0+beta)
                *(12.0+T*nu*nu*(-2.0+3.0*rho*rho)))
                *pow(-(T*alpha*alpha*pow(-1.0+beta, 2.0)
                *pow(K, 2.0*beta))-6.0*alpha*beta*nu*rho*T*pow(K, 1.0+beta)
                +K*K*(12.0+T*nu*nu*(-2.0+3.0*rho*rho)), -1.5);


    res = a0+(f-K)*a1;

#else

    double a0,a1,a2 = 0.0;
	
    switch(flag) 
	{
	    case K_SABR_DIRECTEXACT :
	    case K_SABR_IMPLNVOL :  /// Directe/Exacte / Arithmetic
		{
			a0 = alpha*pow(f,-1 + beta)*pow(1 - (T*pow(f,-2)*(pow(alpha,2)*pow(-1 + beta,2)*pow(f,2*beta) + 6*alpha*beta*nu*rho*pow(f,1 + beta) + 
				pow(f,2)*pow(nu,2)*(2 - 3*pow(rho,2))))/12.,-0.5);
			
            a1 = pow(3,0.5)*(-3*(-1 + beta)*beta*nu*rho*T*pow(alpha,2)*pow(f,2*beta) + alpha*(-1 + beta)*pow(f,1 + beta)*(12 + T*pow(nu,2)*(-2 + 3*pow(rho,2))) + 
				nu*rho*pow(f,2)*(12 + T*pow(nu,2)*(-5 + 6*pow(rho,2))))*pow(-(T*pow(alpha,2)*pow(-1 + beta,2)*pow(f,2*beta)) - 6*alpha*beta*nu*rho*T*pow(f,1 + beta) + 
				pow(f,2)*(12 + T*pow(nu,2)*(-2 + 3*pow(rho,2))),-1.5);
			
            a2 = (pow(3,-0.5)*pow(alpha,-1)*pow(f,-2 - beta)*(3*(-1 + beta)*T*pow(alpha,4)*pow(f,2 + 4*beta)*
				(-36 + pow(beta,2)*(-588 + T*pow(nu,2)*(98 - 657*pow(rho,2))) + 3*T*pow(nu,2)*(2 - 3*pow(rho,2)) + 29*beta*(12 + T*pow(nu,2)*(-2 + 3*pow(rho,2))) + 
				pow(beta,3)*(276 + T*pow(nu,2)*(-46 + 39*pow(rho,2)))) + 
				30*(-1 + beta)*nu*rho*T*pow(alpha,3)*pow(f,3 + 3*beta)*(10*beta*(12 + T*pow(nu,2)) + 4*(12 + T*pow(nu,2)*(-5 + 6*pow(rho,2))) + 
				pow(beta,2)*(84 + T*pow(nu,2)*(-32 + 39*pow(rho,2)))) + pow(alpha,6)*pow(-1 + beta,6)*pow(f,6*beta)*pow(T,2) - 
				18*beta*(7 + 8*beta)*nu*rho*pow(alpha,5)*pow(-1 + beta,3)*pow(f,1 + 5*beta)*pow(T,2) + 
				6*alpha*nu*rho*pow(f,5 + beta)*(-1440 + 120*T*pow(nu,2)*(7 - 9*pow(rho,2) + beta*(-4 + 15*pow(rho,2))) + 
				pow(nu,4)*(-10*(10 - 27*pow(rho,2) + 18*pow(rho,4)) + beta*(184 - 1380*pow(rho,2) + 1305*pow(rho,4)))*pow(T,2)) + 
				pow(f,6)*pow(nu,2)*(2880*(2 - 3*pow(rho,2)) - 36*T*pow(nu,2)*(88 - 300*pow(rho,2) + 225*pow(rho,4)) + 
				pow(nu,4)*(368 - 1062*pow(rho,2) + 1350*pow(rho,4) - 675*pow(rho,6))*pow(T,2)) + 
				3*pow(alpha,2)*pow(f,4 + 2*beta)*(3*(1280 + 120*T*pow(nu,2)*(-4 + 7*pow(rho,2)) + pow(nu,4)*(56 - 260*pow(rho,2) + 225*pow(rho,4))*pow(T,2)) - 
				2*beta*(2400 + 240*T*pow(nu,2)*(-4 + 9*pow(rho,2)) + pow(nu,4)*(128 - 810*pow(rho,2) + 765*pow(rho,4))*pow(T,2)) + 
				pow(beta,2)*(960 + 120*T*pow(nu,2)*(-4 + 27*pow(rho,2)) + pow(nu,4)*(88 - 1440*pow(rho,2) + 1575*pow(rho,4))*pow(T,2))))*
				pow(-(T*pow(alpha,2)*pow(-1 + beta,2)*pow(f,2*beta)) - 6*alpha*beta*nu*rho*T*pow(f,1 + beta) + pow(f,2)*(12 + T*pow(nu,2)*(-2 + 3*pow(rho,2))),-2.5))/40.;
			
			
			break;
		}

	    case K_SABR_DIRECTGEOMETRIC :
		{
			a0 = alpha*pow(K,-1 + beta)*pow(1 - 
				(alpha*beta*nu*rho*T*pow(K,-1 + beta))/2. + 
				(T*nu*nu*(-2 + 3*rho*rho))/12.,-0.5);
		
            a1 = 1./K/K*(-6*(-1 + beta)*beta*nu*rho*T*alpha*alpha*
				pow(K,2*beta) + nu*rho*pow(K,2)*
				(-12 + T*nu*nu*(5 - 6*rho*rho)) + 
				alpha*(-1 + beta)*pow(K,1 + beta)*
				(12 + T*nu*nu*(-2 + 3*rho*rho)))*
				pow(4 - 2*alpha*beta*nu*rho*T*pow(K,-1 + beta) + 
				T*nu*nu*(-0.6666666666666666 + rho*rho),-0.5)*
				pow(-6*alpha*beta*nu*rho*T*pow(K,beta) + 
				K*(12 + T*nu*nu*(-2 + 3*rho*rho)),-1);
			
			break;
		}

	    case K_SABR_DIRECTARITHMETIC :
		{
			a0 = 2*alpha*pow(3,0.5)*pow(K,-1 + beta)*
				pow(1./K/K*((-1 + 3*beta)*T*alpha*alpha*pow(K,2*beta) - 
				6*alpha*beta*nu*rho*T*pow(K,1 + beta) + 
				K*K*(12 + T*nu*nu*(-2 + 3*rho*rho))),-0.5);
			
            a1 = pow(3,0.5)/K/K/K*
				(-6*(-1 + beta)*beta*nu*rho*T*alpha*alpha*pow(K,2*beta) + 
				nu*rho*pow(K,2)*(-12 + T*nu*nu*(5 - 6*rho*rho)) + 
				alpha*(-1 + beta)*pow(K,1 + beta)*
				(12 + T*nu*nu*(-2 + 3*rho*rho)))*
				pow(1./K/K*((-1 + 3*beta)*T*alpha*alpha*pow(K,2*beta) - 
				6*alpha*beta*nu*rho*T*pow(K,1 + beta) + 
				pow(K,2)*(12 + T*nu*nu*(-2 + 3*rho*rho))),-1.5);
			
			
			break;
		}

	    case K_SABR_IMPLNVOL2 :  
	    case K_SABR_DIRECTEXACTSTRIKE :  
		{
			a0 = alpha*pow(f,-1 + beta)*pow(1 - (T*pow(f,-2)*(pow(alpha,2)*pow(-1 + beta,2)*pow(f,2*beta) + 6*alpha*beta*nu*rho*pow(f,1 + beta) + 
				pow(f,2)*pow(nu,2)*(2 - 3*pow(rho,2))))/12.,-0.5);
			
            a1 = pow(3,0.5)*pow(f,-2)*(alpha*(-1 + beta)*pow(f,beta)*(12 + T*pow(nu,2)*(-2 + 3*pow(rho,2))) + f*nu*rho*(12 + T*pow(nu,2)*(-5 + 6*pow(rho,2))))*
				pow(pow(f,-2)*(-(T*pow(alpha,2)*pow(-1 + beta,2)*pow(f,2*beta)) - 6*alpha*beta*nu*rho*T*pow(f,1 + beta) + pow(f,2)*(12 + T*pow(nu,2)*(-2 + 3*pow(rho,2)))),-1.5);
			
            a2 = (pow(3,-0.5)*pow(alpha,-1)*pow(f,-3 - beta)*(60*(-1 + beta)*nu*rho*T*pow(alpha,3)*pow(f,3 + 3*beta)*
				(beta*(-84 + T*pow(nu,2)*(29 - 36*pow(rho,2))) + 2*(12 + T*pow(nu,2)*(-5 + 6*pow(rho,2))) + pow(beta,2)*(132 + T*pow(nu,2)*(-31 + 42*pow(rho,2)))) - 
				3*T*pow(alpha,4)*pow(-1 + beta,2)*pow(f,2 + 4*beta)*(-36 + 3*T*pow(nu,2)*(2 - 3*pow(rho,2)) + 26*beta*(12 + T*pow(nu,2)*(-2 + 3*pow(rho,2))) + 
				pow(beta,2)*(-276 + T*pow(nu,2)*(46 + 51*pow(rho,2)))) + pow(alpha,6)*pow(-1 + beta,6)*pow(f,6*beta)*pow(T,2) - 
				54*beta*nu*rho*pow(alpha,5)*pow(-1 + beta,4)*pow(f,1 + 5*beta)*pow(T,2) + 
				6*alpha*nu*rho*pow(f,5 + beta)*(-1440 + 120*T*pow(nu,2)*(7 - 9*pow(rho,2) + beta*(-4 + 15*pow(rho,2))) + 
				pow(nu,4)*(-10*(10 - 27*pow(rho,2) + 18*pow(rho,4)) + beta*(184 - 1380*pow(rho,2) + 1305*pow(rho,4)))*pow(T,2)) + 
				pow(f,6)*pow(nu,2)*(2880*(2 - 3*pow(rho,2)) - 36*T*pow(nu,2)*(88 - 300*pow(rho,2) + 225*pow(rho,4)) + 
				pow(nu,4)*(368 - 1062*pow(rho,2) + 1350*pow(rho,4) - 675*pow(rho,6))*pow(T,2)) + 
				3*pow(alpha,2)*pow(f,4 + 2*beta)*(3*(1280 + 120*T*pow(nu,2)*(-4 + 7*pow(rho,2)) + pow(nu,4)*(56 - 260*pow(rho,2) + 225*pow(rho,4))*pow(T,2)) - 
				2*beta*(2400 + 120*T*pow(nu,2)*(-8 + 27*pow(rho,2)) + pow(nu,4)*(128 - 1260*pow(rho,2) + 1305*pow(rho,4))*pow(T,2)) + 
				pow(beta,2)*(960 + 120*T*pow(nu,2)*(-4 + 45*pow(rho,2)) + pow(nu,4)*(88 - 2340*pow(rho,2) + 2655*pow(rho,4))*pow(T,2))))*
				pow(T*pow(alpha,2)*pow(-1 + beta,2)*pow(f,2*beta) + 6*alpha*beta*nu*rho*T*pow(f,1 + beta) + pow(f,2)*(-12 + T*pow(nu,2)*(2 - 3*pow(rho,2))),-2)*
				pow(pow(f,-2)*(-(T*pow(alpha,2)*pow(-1 + beta,2)*pow(f,2*beta)) - 6*alpha*beta*nu*rho*T*pow(f,1 + beta) + pow(f,2)*(12 + T*pow(nu,2)*(-2 + 3*pow(rho,2)))),-0.5)
				)/40.;
			break;
		}
		
	    default :
		{	
			ARM_THROW( ERR_INVALID_ARGUMENT,  "CptSABR_implicit_vol_direct_Series : flag:  wrong input!");
			break;
		}
	}
	

	res = a0+(K-f)*a1+(K-f)*(K-f)*a2;

#endif

    return(res);
}



double CptSABR_implicit_vol_direct(double f, double K, double tex,
                                   double alpha, double beta,
                                   double rho, double nu,
                                   int flag)
{
    double vol;


    if ( fabs(K-f) < alpha*ARM_SABR_SERIE_USAGE_LIMIT*f )
    {
       vol = CptSABR_implicit_vol_direct_Series(f, K, tex, alpha,
                                                beta, rho, nu,
                                                flag);
    }
	else if (K > ARM_GP_CF_MAX_STRIKE*f)
		return ARM_GP_CF_MAX_VOL;
    else
    {
		/*
		if ((beta<1.00001)&&(beta>0.99999))
		{
			return CptSABR_BetEqOne_ImplVol( f, K, tex,  alpha, rho,  nu);
		}
		*/
		double z;
		double b1;

		switch(flag) 
		{
		    case K_SABR_DIRECTEXACT :
		    case K_SABR_IMPLNVOL :
			{
				z  = (pow(f,(1.-beta))-pow(K,(1.-beta)))/(alpha*(1.-beta));
				b1 = beta*pow((f+K)/2.,beta-1.);
				
                break;
			}

		    case K_SABR_DIRECTGEOMETRIC :
			{
				z  = log(f/K)/alpha*pow(f*K,(1. - beta)/2.);
				b1 = beta*pow((f+K)/2.,beta-1.);
				
                break;
			}

		    case K_SABR_DIRECTARITHMETIC :
			{
				z  = (f - K)*pow(2./(f + K),beta)/alpha;
				b1 = beta*pow((f+K)/2.,beta-1.);
				
                break;
			}

		    case K_SABR_IMPLNVOL2 :  
		    case K_SABR_DIRECTEXACTSTRIKE :  
			{
				z  = (pow(f,(1.-beta))-pow(K,(1.-beta)))/(alpha*(1.-beta));
				b1 = beta*pow(K,beta-1.);
				
                break;
			}

		    default :
			{	
				ARM_THROW( ERR_INVALID_ARGUMENT,  "CptSABR_implicit_vol_direct : flag: Wrong input :");

                break;
			}
		}

		double z1     = sqrt(1.-2.*nu*rho*z+nu*nu*z*z);
		
        double x     = log((-rho+nu*z+z1)/(1.-rho))/nu;
		
        double theta = 0.25*alpha*b1*nu*rho*z*z
			+log(alpha*pow((f*K),beta/2.)*z/(f-K))
			+log(x*sqrt(sqrt(1.-2.*nu*rho*z+nu*nu*z*z))/z);
	
		vol = log(f/K)/(x*sqrt(1.-2.*tex*theta/(x*x)+2.*tex/(x*x)*log(sqrt(f*K)*log(f/K)/(f-K))));
    }

    return(vol);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///						zs beta=1 (SABR_A in ARM) (CompatibleKernel)
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////

double CptSABR_BetEqOne_ImplVol(double f,double K,double tex, double alpha,double rho, double nu)
{
    double zeta=2.*(nu/alpha)*(f-K)/(f+K);
	double res,x,factor;
	
    factor=(1.+(rho*nu*alpha/4.+(2.-3.*rho*rho)*nu*nu/24.)*tex); // stupar
	if (fabs((f-K)/K)<1e-5)
	{
		res = (1-zeta/12.*(6*rho+zeta*(-2.+3.*rho*rho)))*alpha * factor; // stupar
	}
	else
	{
		if((rho<1)&&(rho>-1))
		{
			x=log((sqrt(1.-2.*rho*zeta+zeta*zeta)+zeta-rho)/(1-rho))/nu;
		}
		else
		{
			x=(1.+zeta/fabs(1-zeta))/nu;
		}
		res = zeta*alpha/(nu*x) * factor;
	}

	return res;
}



double CptSABR_implicit_vol_normal_Series(double f, double K, double tex,
                                          double alpha, double beta,
                                          double rho, double nu,
                                          int flag)
{
    double a0;
    double a1;


    a0 = (alpha*pow(K,-3 + beta)*(tex*alpha*alpha*pow(-1 + beta,2)*
		pow(K,2*beta) + 6*alpha*beta*nu*rho*tex*pow(K,1 + beta) + 
		pow(K,2)*(24 + tex*nu*nu*(2 - 3*rho*rho))))/24.;

	a1 = (pow(K,-4)*(3*tex*alpha*alpha*alpha*pow(-1 + beta,3)*pow(K,3*beta) + 
		nu*rho*tex*alpha*alpha*(-1 - 10*beta + 11*beta*beta)*
		pow(K,1 + 2*beta) + nu*rho*pow(K,3)*
		(-24 + tex*nu*nu*(-2 + 3*rho*rho)) - 
		alpha*pow(K,2 + beta)*
		(24 + tex*nu*nu*(2 - 3*rho*rho) + 
		beta*(-24 + tex*nu*nu*(-2 + 9*rho*rho)))))/48. ;

    double res = a0+(f-K)*a1;

    return(res);
}




double CptSABR_implicit_vol_normal(double f, double K, double tex,
                                   double alpha, double beta, double rho,
                                   double nu,
                                   int flag)
{
    if ( fabs(f-K) < (f*ARM_SABR_SERIE_USAGE_LIMIT*alpha ) )
    {
       double res;

       res = CptSABR_implicit_vol_normal_Series(f, K, tex, alpha, beta, rho, nu, flag);

       return(res);
    }
    else
    {
       double res;

       double z;
		
       switch (flag) 
       {
            case K_SABR_GEO :
		    case K_SABR_NORMALEXACT :
			{
				z = (pow(f,(1.-beta))-pow(K,(1.-beta)))/(alpha*(1.-beta));
			
                break;
			}

		    case K_SABR_NORMALGEOMETRIC :
			{
				z = log(f/K)/alpha*pow(f*K,(1. - beta)/2.);
				
                break;
			}

		    case K_SABR_NORMALARITHMETIC :
			{
				z = (f - K)*pow(2./(f + K),beta)/alpha;
				
                break;
			}
			
		    default :
			{	
				ARM_THROW( ERR_INVALID_ARGUMENT,  "CptSABR_implicit_vol_normal : flag: Wrong input :");
				
                break;
			}
		}
		
		double x = log(pow(1. - rho,-1)*(-rho + nu*z + 
			       pow(1. - 2.*nu*rho*z + nu*nu*z*z,0.5)))/nu;
		
        res = alpha*z*pow(f*K,(-1. + beta)/2.)*
			(1. + tex*((alpha*beta*nu*rho*pow(f*K,(-1. + beta)/2.))/4. + 
			(alpha*alpha*pow(1. - beta,2)*pow(f*K,-1. + beta))/24. + 
			(nu*nu*(2 - 3*rho*rho))/24.))/x*
			pow(1. + (pow(1. - beta,2)*pow(log(f/K),2))/24. + 
			(pow(1. - beta,4)*pow(log(f/K),4))/1920.,-1);
  

       return(res);
    }
}



double PhiAlpha(double f, double K,
                double maturity,
                double sigmaATM,
                double alpha, double rho, double nu, double beta,
                double weight, int type)
{
    double phi = 0.0;

#ifdef SABR_APPROX_NORMAL

    switch(type)
    {
        case K_SABR_ARITH :
        {
            double ratio = 1.0;

            double oneMoinsBeta = 1.0-beta;

            double SQROneMoinsBeta = oneMoinsBeta*oneMoinsBeta;

            double sig = alpha/(pow(f*K, oneMoinsBeta/2.0)*(1+SQROneMoinsBeta
                         *SQR(log(f/K))/24.0+SQR(SQROneMoinsBeta)
                         *SQR(log(f/K))*SQR(log(f/K))/1920.0))
                         *(ratio)
                         *(1.0+(SQR(oneMoinsBeta)*SQR(alpha)/(24*pow(f*K, oneMoinsBeta))
                         +rho*beta*nu*alpha/(4*pow(f*K, oneMoinsBeta/2.0))
                         +(2.0-3.0*SQR(rho))*SQR(nu)/24)*maturity);

            double res = sig-sigmaATM;

            return(res);

        };
        break;

        case K_SABR_IMPLNVOL :
        {
            double sig;

            sig = alpha/(pow(f*K, 0.5*(1-beta))
                  *sqrt(1.0+2.0*maturity*(-1.0/24.0*pow(alpha*(beta-1.0), 2.0)
                  /pow(f, 2.0*(1.0-beta))-1.0/4.0*alpha*beta*rho*nu/pow(f,1.0-beta)
                  +nu*nu*(-1.0/12.0+1.0/8.0*rho*rho))));

            double res = sig-sigmaATM;

            return(res);
        };
        break;

        case K_SABR_WEIGHT :
        {
            double sig;

            double ratio = 1.0;

            double oneMoinsBeta = 1.0-beta;

            double SQROneMoinsBeta = oneMoinsBeta*oneMoinsBeta;

            double sig1 = alpha/(pow(f*K, oneMoinsBeta/2.0)*(1+SQROneMoinsBeta
                         *SQR(log(f/K))/24.0+SQR(SQROneMoinsBeta)
                         *SQR(log(f/K))*SQR(log(f/K))/1920.0))
                         *(ratio)
                         *(1.0+(SQR(oneMoinsBeta)*SQR(alpha)/(24*pow(f*K, oneMoinsBeta))
                         +rho*beta*nu*alpha/(4*pow(f*K, oneMoinsBeta/2.0))
                         +(2.0-3.0*SQR(rho))*SQR(nu)/24)*maturity);

            double sig2 = alpha/(pow(f*K, 0.5*(1-beta))
                          *sqrt(1.0+2.0*maturity*(-1.0/24.0*pow(alpha*(beta-1.0), 2.0)
                          /pow(f, 2.0*(1.0-beta))-1.0/4.0*alpha*beta*rho
                          *nu/pow(f,1.0-beta)
                          +nu*nu*(-1.0/12.0+1.0/8.0*rho*rho))));

            sig = weight*sig1+(1.0-weight)*sig2;

            double res = sig-sigmaATM;

            return(res);
        };
    }

#else

    switch(type)
    {
        case K_SABR_GEO :
        {
            double sig = CptSABR_implicit_vol_normal(f, K, maturity, alpha, beta,
                                                     rho, nu, K_SABR_GEO);

            double res = sig-sigmaATM;

            return(res);

        };
        break;

        case K_SABR_IMPLNVOL :
        {
            double sig;

            sig = CptSABR_implicit_vol_direct(f, K, maturity, alpha, beta,
                                              rho, nu, K_SABR_IMPLNVOL);

            double res = sig-sigmaATM;

            return(res);
        };
        break;

        case K_SABR_WEIGHT :
        {
            double sig;

            double sig1 = CptSABR_implicit_vol_normal(f, K, maturity, alpha, beta,
                                                      rho, nu, K_SABR_GEO);

            double sig2 = CptSABR_implicit_vol_direct(f, K, maturity, alpha, beta,
                                                      rho, nu, K_SABR_IMPLNVOL);

            sig = weight*sig1+(1.0-weight)*sig2;

            double res = sig-sigmaATM;

            return(res);
        };
    }

#endif

    return(phi);
}



double CptSABR_implicit_vol_direct_DerAlpha(double f, double K, double tex,
                                            double alpha, double beta,
                                            double rho, double nu,
                                            int flag)
{
    double deriv;

		if ( fabs(K-f) < alpha*ARM_SABR_SERIE_USAGE_LIMIT*f )
    {
		switch (flag) 
		{
		    case K_SABR_DIRECTEXACT :
		    case K_SABR_IMPLNVOL :
			{
				deriv = pow(3,0.5)*pow(K,beta)*(-2*alpha*tex*pow(K,beta)*(-3*beta*K*nu*rho - alpha*pow(-1 + beta,2)*pow(K,beta))*
					(-(tex*alpha*alpha*pow(-1 + beta,2)*pow(K,2*beta)) - 6*alpha*beta*nu*rho*tex*pow(K,1 + beta) + pow(K,2)*(12 + tex*nu*nu*(-2 + 3*rho*rho))) + 
					(-1 + beta)*(f - K)*(-12*alpha*beta*nu*rho*tex*pow(K,beta) + K*(12 + tex*nu*nu*(-2 + 3*rho*rho)))*
					(-(tex*alpha*alpha*pow(-1 + beta,2)*pow(K,2*beta)) - 6*alpha*beta*nu*rho*tex*pow(K,1 + beta) + pow(K,2)*(12 + tex*nu*nu*(-2 + 3*rho*rho))) - 
					3*(f - K)*tex*(-3*beta*K*nu*rho - alpha*pow(-1 + beta,2)*pow(K,beta))*
					(-6*(-1 + beta)*beta*nu*rho*tex*alpha*alpha*pow(K,2*beta) + nu*rho*pow(K,2)*(-12 + tex*nu*nu*(5 - 6*rho*rho)) + 
					alpha*(-1 + beta)*pow(K,1 + beta)*(12 + tex*nu*nu*(-2 + 3*rho*rho))) + 
					2*pow(tex*alpha*alpha*pow(-1 + beta,2)*pow(K,2*beta) + 6*alpha*beta*nu*rho*tex*pow(K,1 + beta) + pow(K,2)*(-12 + tex*nu*nu*(2 - 3*rho*rho)),2))*
					pow(-(tex*alpha*alpha*pow(-1 + beta,2)*pow(K,2*beta)) - 6*alpha*beta*nu*rho*tex*pow(K,1 + beta) + pow(K,2)*(12 + tex*nu*nu*(-2 + 3*rho*rho)),-2.5);
				
                return(deriv);

                break;
			}

		    case K_SABR_DIRECTEXACTSTRIKE :
		    case K_SABR_IMPLNVOL2 :
			{
				deriv = pow(3,0.5)*pow(K,beta)*(3*(-f + K)*tex*(3*beta*K*nu*rho + alpha*pow(-1 + beta,2)*pow(K,beta))*
					(6*(-1 + beta)*beta*nu*rho*tex*pow(alpha,2)*pow(K,2*beta) - alpha*(-1 + beta)*pow(K,1 + beta)*(12 + tex*pow(nu,2)*(-2 + 3*pow(rho,2))) + 
					nu*rho*pow(K,2)*(12 + tex*pow(nu,2)*(-5 + 6*pow(rho,2))))*pow(-(tex*pow(alpha,2)*pow(-1 + beta,2)*pow(K,2*beta)) - 6*alpha*beta*nu*rho*tex*pow(K,1 + beta) + 
					pow(K,2)*(12 + tex*pow(nu,2)*(-2 + 3*pow(rho,2))),-2.5) + (-1 + beta)*(-f + K)*(12*alpha*beta*nu*rho*tex*pow(K,beta) + K*(-12 + tex*pow(nu,2)*(2 - 3*pow(rho,2))))*
					pow(-(tex*pow(alpha,2)*pow(-1 + beta,2)*pow(K,2*beta)) - 6*alpha*beta*nu*rho*tex*pow(K,1 + beta) + pow(K,2)*(12 + tex*pow(nu,2)*(-2 + 3*pow(rho,2))),-1.5) + 
					2*alpha*tex*pow(K,-3 + beta)*(3*beta*K*nu*rho + alpha*pow(-1 + beta,2)*pow(K,beta))*
					pow(pow(K,-2)*(-(tex*pow(alpha,2)*pow(-1 + beta,2)*pow(K,2*beta)) - 6*alpha*beta*nu*rho*tex*pow(K,1 + beta) + pow(K,2)*(12 + tex*pow(nu,2)*(-2 + 3*pow(rho,2)))),-1.5) + 
					2*pow(K,-1)*pow(pow(K,-2)*(-(tex*pow(alpha,2)*pow(-1 + beta,2)*pow(K,2*beta)) - 6*alpha*beta*nu*rho*tex*pow(K,1 + beta) + pow(K,2)*(12 + tex*pow(nu,2)*(-2 + 3*pow(rho,2)))),-0.5));
			
                return(deriv);

                break;
			}

		    case K_SABR_DIRECTGEOMETRIC :
			{
			    deriv = pow(K,-2 + beta)*pow(4 - 2*alpha*beta*nu*rho*tex*pow(K,-1 + beta) + tex*nu*nu*(-0.6666666666666666 + rho*rho),-0.5)*
					(-(f*(9*alpha*(-1 + beta)*beta*nu*rho*tex*pow(K,1 + beta)*(12 + tex*nu*nu*(-2 + 3*rho*rho)) - 
					18*(-1 + beta)*alpha*alpha*pow(beta,2)*pow(K,2*beta)*nu*nu*rho*rho*tex*tex + 
					pow(K,2)*(beta*(-144 + 12*tex*nu*nu*(4 + 3*rho*rho) + pow(nu,4)*(-4 - 33*rho*rho + 45*pow(rho,4))*tex*tex) + pow(12 + tex*nu*nu*(-2 + 3*rho*rho),2))
					)) + K*(9*alpha*(-3 + beta)*beta*nu*rho*tex*pow(K,1 + beta)*(12 + tex*nu*nu*(-2 + 3*rho*rho)) - 
					18*(-3 + beta)*alpha*alpha*pow(beta,2)*pow(K,2*beta)*nu*nu*rho*rho*tex*tex + 
					pow(K,2)*(beta*(-144 + 12*tex*nu*nu*(4 + 3*rho*rho) + pow(nu,4)*(-4 - 33*rho*rho + 45*pow(rho,4))*tex*tex) + 3*pow(12 + tex*nu*nu*(-2 + 3*rho*rho),2))
					))*pow(-6*alpha*beta*nu*rho*tex*pow(K,beta) + K*(12 + tex*nu*nu*(-2 + 3*rho*rho)),-2);
				
                return(deriv);

                break;
			}

		    case K_SABR_DIRECTARITHMETIC :
			{	
				deriv = pow(3,0.5)*pow(K,-5 + beta)*(-2*alpha*tex*pow(K,beta)*(-3*beta*K*nu*rho + alpha*(-1 + 3*beta)*pow(K,beta))*
					((-1 + 3*beta)*tex*alpha*alpha*pow(K,2*beta) - 6*alpha*beta*nu*rho*tex*pow(K,1 + beta) + pow(K,2)*(12 + tex*nu*nu*(-2 + 3*rho*rho))) + 
					(-1 + beta)*(f - K)*(-12*alpha*beta*nu*rho*tex*pow(K,beta) + K*(12 + tex*nu*nu*(-2 + 3*rho*rho)))*
					((-1 + 3*beta)*tex*alpha*alpha*pow(K,2*beta) - 6*alpha*beta*nu*rho*tex*pow(K,1 + beta) + pow(K,2)*(12 + tex*nu*nu*(-2 + 3*rho*rho))) - 
					3*(f - K)*tex*(-3*beta*K*nu*rho + alpha*(-1 + 3*beta)*pow(K,beta))*
					(-6*(-1 + beta)*beta*nu*rho*tex*alpha*alpha*pow(K,2*beta) + nu*rho*pow(K,2)*(-12 + tex*nu*nu*(5 - 6*rho*rho)) + 
					alpha*(-1 + beta)*pow(K,1 + beta)*(12 + tex*nu*nu*(-2 + 3*rho*rho))) + 
					2*pow((-1 + 3*beta)*tex*alpha*alpha*pow(K,2*beta) - 6*alpha*beta*nu*rho*tex*pow(K,1 + beta) + pow(K,2)*(12 + tex*nu*nu*(-2 + 3*rho*rho)),2))*
					pow(pow(K,-2)*((-1 + 3*beta)*tex*alpha*alpha*pow(K,2*beta) - 6*alpha*beta*nu*rho*tex*pow(K,1 + beta) + pow(K,2)*(12 + tex*nu*nu*(-2 + 3*rho*rho))),-2.5);
				
                return(deriv);

                break;
			}
			
		    default :
			{	
				ARM_THROW( ERR_INVALID_ARGUMENT,  "CptSABR_implicit_vol_direct_DerAlpha : flag: Wrong input :");
				
                break;
			}
		}
    }
    else
    {
		double z,x,b1,theta,Derzalpha;
		
		switch (flag) 
		{
		    case K_SABR_DIRECTEXACT :
		    case K_SABR_IMPLNVOL :
			{
				z = 1./(alpha*(1.-beta))*
					(pow(f,1 - beta) - pow(K,1 - beta));
				
				Derzalpha = -z/alpha;
				
                b1 = beta/pow((f+K)/2.,1 - beta);
			
                break;
			}
		
		    case K_SABR_DIRECTEXACTSTRIKE :
		    case K_SABR_IMPLNVOL2 :
			{
				z = 1./(alpha*(1.-beta))*
					(pow(f,1 - beta) - pow(K,1 - beta));
				
				Derzalpha = -z/alpha;
				b1 = beta/pow(K,1 - beta);
				break;
			}

		    case K_SABR_DIRECTGEOMETRIC :
			{
				
				z=log(f/K)/alpha*pow(f*K,(1. - beta)/2.);
				
				Derzalpha=-z/alpha;
				b1=beta/pow(K,1 - beta);
				break;
			}

		    case K_SABR_DIRECTARITHMETIC :
			{	
				z = (f - K)*pow(2./(f + K),beta)/alpha;
				
				Derzalpha = -z/alpha;
			
                b1 = beta/pow(K,1 - beta);
				
				break;
			}
			
		    default :
			{	
				ARM_THROW( ERR_INVALID_ARGUMENT,  "CptSABR_implicit_vol_direct_DerAlpha : flag: Wrong input :");
				break;
			}
		}
		
		x = log(pow(1 - rho,-1)*(-rho + nu*z + 
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			0.5)))*pow(nu,-1);
		
		theta = log(alpha*z*pow(f - K,-1)*pow(f*K,beta/2.)) + 
			log(x*pow(z,-1)*pow(1 - 2*nu*rho*z + 
			nu*nu*pow(z,2),0.25)) + 
			(alpha*b1*nu*rho*pow(z,2))/4.;
		
		
		
		double DerImpVoltheta = tex*log(f*pow(K,-1))*pow(x,-3)*pow(pow(x,-2)*
			(-2*tex*theta + 2*tex*log(log(f*pow(K,-1))*pow(f - K,-1)*pow(f*K,0.5)) + 
			pow(x,2)),-1.5);
		
		
		double DerImpVolx = -(log(f*pow(K,-1))*(4*tex*theta*pow(x,-3) - 
			4*tex*log(log(f*pow(K,-1))*pow(f - K,-1)*pow(f*K,0.5))*pow(x,-3))*
			pow(x,-1)*pow(1 - 2*tex*theta*pow(x,-2) + 
			2*tex*log(log(f*pow(K,-1))*pow(f - K,-1)*pow(f*K,0.5))*pow(x,-2),-1.5)
			)/2. - log(f*pow(K,-1))*pow(x,-2)*
			pow(1 - 2*tex*theta*pow(x,-2) + 
			2*tex*log(log(f*pow(K,-1))*pow(f - K,-1)*pow(f*K,0.5))*pow(x,-2),-0.5);
		
		double Derxz = pow(nu,-1)*(nu + ((-2*nu*rho + 2*z*nu*nu)*
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			-0.5))/2.)*pow(-rho + nu*z + 
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),0.5),
			-1);
		
		double Derthetab1 = (alpha*nu*rho*pow(z,2))/4.;
		
		double Derthetax = 1.0/x;
		
		double Derthetaalpha = pow(alpha,-1)+(b1*nu*rho*pow(z,2))/4.;
				
		
		double Derthetaz = (alpha*b1*nu*rho*z)/2. + pow(z,-1) + 
			z*pow(x,-1)*pow(1 - 2*nu*rho*z + 
			nu*nu*pow(z,2),-0.25)*
			((x*(-2*nu*rho + 2*z*nu*nu)*pow(z,-1)*
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			-0.75))/4. - 
			x*pow(z,-2)*pow(1 - 2*nu*rho*z + 
			nu*nu*pow(z,2),0.25));
		
		deriv = DerImpVoltheta*(Derthetaalpha + Derthetaz*Derzalpha + 
			    Derthetax*Derxz*Derzalpha) + DerImpVolx*Derxz*Derzalpha;
       
        return(deriv);
    }
}



double CptSABR_implicit_vol_normal_DerAlpha(double f, double K, double tex,
                                            double alpha, double beta, double rho,
                                            double nu,
                                            int flag)

{
    double deriv;

  if ( fabs(K-f) < alpha*ARM_SABR_SERIE_USAGE_LIMIT*f )
    {
       deriv = (pow(K,-4 + beta)*(f*(9*tex*alpha*alpha*pow(-1 + beta,3)*
			pow(K,2*beta) + 2*alpha*nu*rho*tex*(-1 - 10*beta + 11*pow(beta,2))*pow(K,1 + beta) + 
			pow(K,2)*(-24 + beta*(24 + tex*nu*nu*(2 - 9*rho*rho)) +
			tex*nu*nu*(-2 + 3*rho*rho))) + 
			K*(-3*(-5 + 3*beta)*tex*alpha*alpha*pow(-1 + beta,2)*pow(K,2*beta) -
			2*alpha*nu*rho*tex*(-1 - 22*beta + 11*pow(beta,2))*pow(K,1 + beta) + 
			pow(K,2)*(72 + tex*nu*nu*(6 - 9*rho*rho) + beta*(-24 +
			tex*nu*nu*(-2 + 9*rho*rho))))))/48.;

       return(deriv);
    }
    else
    {
		double z, x, b1, Derzalpha, dersigmaalpha, dersigmax, dersigmaz, Derxz;

	    switch(flag) 
        {
            case K_SABR_GEO :
	        case K_SABR_NORMALEXACT :
            {
			    z = 1./(alpha*(1.-beta))*
			       (pow(f,1 - beta) - pow(K,1 - beta));
			
			    Derzalpha = -z/alpha;
			    b1 = beta/pow(K,1 - beta);
			
                break;
            }

	        case K_SABR_NORMALGEOMETRIC :
            {
			    z = log(f/K)/alpha*pow(f*K,(1. - beta)/2.);
			
			    Derzalpha = -z/alpha;

			    b1 = beta/pow(K,1 - beta);
			
                break;
            }

	        case K_SABR_NORMALARITHMETIC :
            {
			    z = (f - K)*pow(2./(f + K),beta)/alpha;
			
			    Derzalpha = -z/alpha;

			    b1 = beta/pow(K,1 - beta);
			
			    break;
            }
		
	        default :
            {	
			    throw "CptSABR_implicit_vol_normal_DerAlpha : flag: Wrong input :";
			
                break;
            }
        }
		
		x = log(pow(1 - rho,-1)*(-rho + nu*z + 
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			0.5)))*pow(nu,-1);

		Derxz = pow(nu,-1)*(nu + ((-2*nu*rho + 2*z*nu*nu)*
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			-0.5))/2.)*pow(-rho + nu*z + 
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),0.5),
			-1);
		
		dersigmaalpha = 80*z*pow(f,-2)*pow(K,-2)*pow(f*K,beta/2.)*(3*tex*alpha*alpha*
			pow(-1 + beta,2)*pow(f*K,0.5 + beta) + 
			f*K*(24*pow(f*K,0.5) + 12*alpha*beta*nu*rho*tex*pow(f*K,beta/2.) + 
			tex*pow(f*K,0.5)*nu*nu*(2 - 3*rho*rho)))*pow(x,-1)*
			pow(1920 + 80*pow(-1 + beta,2)*pow(log(f*pow(K,-1)),2) + pow(-1 + beta,4)*
			pow(log(f*pow(K,-1)),4),-1);
		
		dersigmax = -(alpha*z*pow(f*K,(-1 + beta)/2.)*(1 + tex*((alpha*beta*nu*rho*
			pow(f*K,(-1 + beta)/2.))/4. + (alpha*alpha*pow(-1 + beta,2)*
			pow(f*K,-1 + beta))/24. + 
			(nu*nu*(2 - 3*rho*rho))/24.))*pow(x,-2)*pow(1 + (pow(-1 + beta,2)*
			pow(log(f*pow(K,-1)),2))/24. + (pow(-1 + beta,4)*pow(log(f*pow(K,-1)),4))/1920.,-1));
		
		dersigmaz = alpha*pow(f*K,(-1 + beta)/2.)*(1 + tex*((alpha*beta*nu*rho*pow(f*K,(-1 + beta)/2.))/4. +
			(alpha*alpha*pow(-1 + beta,2)*pow(f*K,-1 + beta))/24. + 
			(nu*nu*(2 - 3*rho*rho))/24.))*pow(x,-1)*pow(1 + (pow(-1 + beta,2)*
			pow(log(f*pow(K,-1)),2))/24. + (pow(-1 + beta,4)*pow(log(f*pow(K,-1)),4))/1920.,-1);

	    deriv = dersigmaalpha+dersigmax*Derxz*Derzalpha+dersigmaz*Derzalpha;

        return(deriv);
    }
}



double DerivPhiAlpha(double f, double K,
                     double mat,
                     double alpha, double rho, double nu, double beta,
                     double weight, int type,
                     double sigATM)
{
    double deriv;


#ifdef SABR_APPROX_NORMAL

    switch(type)
    {
        case K_SABR_ARITH :
        {
            deriv = dsigBSappATM(f, K, mat, alpha, rho, nu, beta);
        };
        break;

        case K_SABR_IMPLNVOL :
        {
            deriv = dsigBS2(f, K, mat, alpha, rho, nu, beta);
        };
        break;

        case K_SABR_WEIGHT :
        {
            double deriv1 = dsigBSappATM(f, K, mat,
                                         alpha, rho, nu, beta);

            double deriv2 = dsigBS2(f, K, mat,
                                    alpha, rho, nu, beta);

            deriv = weight*deriv1+(1.0-weight)*deriv2;
        };
    }

#else

    switch(type)
    {
        case K_SABR_GEO :
        {
            deriv = CptSABR_implicit_vol_normal_DerAlpha(f, K, mat,
                                                         alpha,
                                                         beta,
                                                         rho, nu, K_SABR_GEO);
        };
        break;

        case K_SABR_IMPLNVOL :
        {
            deriv = CptSABR_implicit_vol_direct_DerAlpha(f, K, mat, alpha,
                                                         beta, rho, nu, K_SABR_IMPLNVOL);
        };
        break;

        case K_SABR_WEIGHT :
        {
            double deriv1 = CptSABR_implicit_vol_normal_DerAlpha(f, K, mat,
                                                                 alpha,
                                                                 beta,
                                                                 rho, nu, K_SABR_GEO);

            double deriv2 = CptSABR_implicit_vol_direct_DerAlpha(f, K, mat,
                                                                 alpha,
                                                                 beta,
                                                                 rho, nu,
                                                                 K_SABR_IMPLNVOL);

            deriv = weight*deriv1+(1.0-weight)*deriv2;
        };
    }

#endif

    return(deriv);
}



double ComputeAlphaBetaEqOne(double f, double K,
                             double maturity, 
                             double sigmaATM,
                             double Rho_t, 
							 double Nu_t,
                             int type)
{
    double spot   = f;
    double strike = K;
	double Alpha_t = 0.0;

    // Beta == 1

    if ( type != K_LD ) // SABR type
    {
       double ATMVol= sigmaATM;    
       
       if (IS_ONE_TEST(Rho_t))
       {
          Rho_t =  0.99 ;
       }
       
       if (IS_ZERO_TEST(Rho_t))
       {
          Alpha_t = ATMVol/(1.0+(SABR_SQR(Nu_t)*maturity/12.0)); // SQR replaced explicitely
         
       }
       else
       {
          double sqrtTerm;
            
          double rho2_3_24;
            
          rho2_3_24 = ((2.0-3.0*SABR_SQR(Rho_t))/24.0)*SABR_SQR(Nu_t) * maturity;
            
          sqrtTerm = sqrt(SABR_SQR(1.0+rho2_3_24) + (ATMVol*Rho_t*Nu_t*maturity));
            
          Alpha_t = (2.0/(Rho_t*Nu_t*maturity))*(sqrtTerm-(1.0+rho2_3_24));

       }  
    }
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "ComputeAlphaBetaEqOne() unavailable for LD type" );
	
	return Alpha_t;
}



/*----------------------------------------------------------------------------*
    To get Alpha from the current Rho and nu, beta ...  by Newton raphson
*----------------------------------------------------------------------------*/
double ComputeAlpha(double f, 
					double K,
                    double mat,
                    double sigmaATM,
                    double rho, 
					double nu, 
					double beta,
                    double weight, 
					int type)
{
    double alpha0 = sigmaATM*pow(f, 1.0-beta);

    double alpha_i = alpha0;
    double alpha_iplus1 = alpha0;

    double prec = 1e-5;
    int nbIterMax = 100;

    if (IS_ONE_TEST(beta))
    {
       double alphaRes = ComputeAlphaBetaEqOne(f, K, mat, sigmaATM, rho, nu,type);

       return alphaRes;
    }

    double ATMVal = f;

	int i = 0;
    do
    {
        i++;

        alpha_i = alpha_iplus1;

        double denomDeriv = DerivPhiAlpha(ATMVal, ATMVal, mat,
                                          alpha_i, rho, nu, beta,
                                          weight, type, sigmaATM);

        if ((!finite(denomDeriv)) || IS_ZERO_TEST(denomDeriv))
        {
           ARM_THROW( ERR_INVALID_ARGUMENT, "ComputeAlpha: impossible to compute Alpha : Deriv==0");
        }

        alpha_iplus1 = alpha_i-PhiAlpha(ATMVal, ATMVal, mat,
                                        sigmaATM,
                                        alpha_i, rho, nu, beta,
                                        weight, type)/denomDeriv;
    }
    while (( fabs(alpha_iplus1-alpha_i) > prec ) && ( i < nbIterMax ));

    if ( i == nbIterMax )
    {
       ARM_THROW( ERR_INVALID_ARGUMENT, "problem in ComputeAlpha: reach max iteration");
    }
	
	return alpha_i ;
}


double ARM_CptSABRVolBetaDiffOneNormalExactSABR_GEO(double maturity,
                                                    double volOrAlpha,
                                                    double f,
                                                    double K,
                                                    double rho,
                                                    double nu,
                                                    double beta,
                                                    int sigmaOrAlphaFLAG)
{
    double alpha = volOrAlpha;

    double theWeight = 0.0;
    if (sigmaOrAlphaFLAG) // The input is ATM Vol
    {
       alpha = ComputeAlpha(f, K, maturity, volOrAlpha, rho, nu, beta, theWeight, K_SABR_GEO);
    }

    double SABRVol = CptSABR_implicit_vol_normal(f, K, maturity, alpha, beta, rho, nu, K_SABR_GEO);

    return SABRVol;
}


double ARM_CptSABRVolBetaDiffOneDirectExactSABR_IMPLNVOL(double maturity,
                                                         double volOrAlpha,
                                                         double f,
                                                         double K,
                                                         double rho,
                                                         double nu,
                                                         double beta,
                                                         int sigmaOrAlphaFLAG)
{
    double alpha = volOrAlpha;

    double theWeight = 0.0;
    if (sigmaOrAlphaFLAG) // The input is ATM Vol
    {
       alpha = ComputeAlpha(f, K, maturity, volOrAlpha, rho, nu, beta, theWeight, K_SABR_IMPLNVOL);
    }

    double SABRVol = CptSABR_implicit_vol_direct(f, K, maturity, alpha, beta, rho, nu, K_SABR_IMPLNVOL);

    return SABRVol;
}




/************************************************************/
/*                                                          */
/*                 Principal SABR Call                      */
/*                                                          */
/************************************************************/

double ARM_CalcSABRVolBetaDiffOne(double maturity,
                                  double volOrAlpha,
                                  double f,
                                  double K,
                                  double rho,
                                  double nu,
                                  double beta,
                                  int sigmaOrAlphaFLAG,
                                  int SABR_TYPE,
                                  double SABRWeight)
{
    double SABRVol = 0.0;


    switch(SABR_TYPE)
    {
        case K_SABR_GEO :
        {
            SABRVol = ARM_CptSABRVolBetaDiffOneNormalExactSABR_GEO(maturity,
                                                                   volOrAlpha,
                                                                   f,
                                                                   K,
                                                                   rho,
                                                                   nu,
                                                                   beta,
                                                                   sigmaOrAlphaFLAG);
            return(SABRVol);
        };
        break;

        case K_SABR_IMPLNVOL :
        {
            // case "Direct exact"

            SABRVol = ARM_CptSABRVolBetaDiffOneDirectExactSABR_IMPLNVOL(maturity,
                                                                        volOrAlpha,
                                                                        f,
                                                                        K,
                                                                        rho,
                                                                        nu,
                                                                        beta,
																		sigmaOrAlphaFLAG);


            return(SABRVol);

        };
        break;

        case K_SABR_WEIGHT :
        {
            double sig1 = ARM_CptSABRVolBetaDiffOneNormalExactSABR_GEO(maturity,
                                                                   volOrAlpha,
                                                                   f,
                                                                   K,
                                                                   rho,
                                                                   nu,
                                                                   beta,
                                                                   sigmaOrAlphaFLAG);

            double sig2 = ARM_CptSABRVolBetaDiffOneDirectExactSABR_IMPLNVOL(maturity,
                                                                        volOrAlpha,
                                                                        f,
                                                                        K,
                                                                        rho,
                                                                        nu,
                                                                        beta,
																		sigmaOrAlphaFLAG); 

            double w = SABRWeight;

            SABRVol = w*sig1+(1.0-w)*sig2;

            return(SABRVol);

        };
        break;

        default :
        {
            ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_CalcSABRVolBetaDiffOne: (Beta<>1): invalid type!?");
        }
        break;
    };

    return(SABRVol);
}




/*----------------------------------------------------------------------------*/
/*---- End Of File ----*/
