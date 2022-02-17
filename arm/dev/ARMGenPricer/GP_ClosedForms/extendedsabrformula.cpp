/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file extendedsabrformula.cpp
 *
 *  \brief
 *
 *	\author  O.Croissant & E.Ezzine
 *	\version 1.0
 *	\date April 2004
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpbase/gpvector.h"


#include <cmath>
#include <complex>
#include "gpbase/gpmatrix.h"
#include "gpbase/gpvector.h"
#include "gpbase/removenagwarning.h"
#include "nag.h"
#include "nage04.h"

#include "gpnumlib/gaussiananalytics.h"
#include "gpbase/numericconstant.h"

#include "gpclosedforms/sabrvanilla.h"
#include "gpclosedforms/sabrimpliedvol.h"
#include "gpclosedforms/smile_sabr.h"
#include "gpclosedforms/extendedsabrformula.h"

#include "gpclosedforms/optimization1.h"


#include <glob/expt.h>   // for the exceptions

using namespace std;

CC_BEGIN_NAMESPACE(ARM)

#define ARM_GP_CF_SABR_SERIE_USAGE_LIMIT 0.1

double SABR_ComputeImpliedVol( double f,
       double K,
       double T,
       double alpha,
       double beta,
       double rho, 
       double nu,
       int SABRflag )
{
     double value = 0.0;
    ///Test if we ar aroud the money
    if (fabs(K-f)<ARM_GP_CF_SABR_SERIE_USAGE_LIMIT*f*alpha) 
	{
        value = SABR_ComputeImpliedVolAroundATM(f,
            K,
            T,
            alpha,
            beta,
            rho,
            nu,
            SABRflag);

        return value;
    }
   value =SABR_implicit_vol_direct(f,
            K,
            T,
            alpha,
            beta,
            rho,
            nu,
            SABRflag);
   return value;
}

double SABR_ComputeImpliedVolAroundATM( double f,
       double K,
       double T,
       double alpha,
       double beta,
       double rho, 
       double nu,
       int SABRflag )
{
    double value=0.0;
	double alpha2=alpha*alpha;
    double alpha3=alpha2*alpha;
    double alpha4=alpha3*alpha;
    double alpha5=alpha4*alpha;
    double alpha6=alpha5*alpha;
	double beta2=beta*beta;
	double rho2=rho*rho;
    double rho3=rho2*rho;
    double rho4=rho3*rho;
    double rho5=rho4*rho;
    double rho6=rho5*rho;
	double nu2=nu*nu;
    double nu3=nu2*nu;
    double nu4=nu3*nu;
	double f2=f*f;double f3=f2*f;double f4=f3*f;double f5=f4*f;double f6=f5*f;
	double T2=T*T;
	switch (SABRflag) 
	{
	case ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACT :
	case ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL :
		{
			value = (pow(f,-1 + beta)*alpha)/sqrt(1 - (T*(pow(f,2*beta)*alpha2*pow(-1 + \
				beta,2) + 6*pow(f,1 + beta)*alpha*beta*nu*rho + f2*nu2*(2 - \
				3*rho2)))/(12.*f2)) + (sqrt(3.)*f*(-f + K)*(pow(f,beta)*alpha*(-1 + beta)*(12 \
				+ T*nu2*(-2 + 3*rho2)) + f*nu*rho*(12 + T*nu2*(-5 + \
				6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),1.5) + (pow(f,-2 - \
				beta)*pow(-f + K,2)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 54*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 - 3*rho2) - \
				36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + 1350*rho4 - \
				675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + beta)*nu*rho*(beta*(-84 + \
				T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + \
				T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + 4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 \
				+ 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + \
				T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 \
				- 9*rho2 + beta*(-4 + 15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + \
				beta*(184 - 1380*rho2 + 1305*rho4))) + 3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + \
				120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + \
				120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 - 1260*rho2 + 1305*rho4)) + beta2*(960 \
				+ 120*T*nu2*(-4 + 45*rho2) + T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(40.*sqrt(3.)*alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5));
            break;
			
		}
	case ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL2 :  
	case ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACTSTRIKE :  
		{
			value = (pow(f,-1 + beta)*alpha)/sqrt(1 - (T*(pow(f,2*beta)*alpha2*pow(-1 + \
				beta,2) + 6*pow(f,1 + beta)*alpha*beta*nu*rho + f2*nu2*(2 - \
				3*rho2)))/(12.*f2)) + (sqrt(3.)*f*(-f + K)*(pow(f,beta)*alpha*(-1 + beta)*(12 \
				+ T*nu2*(-2 + 3*rho2)) + f*nu*rho*(12 + T*nu2*(-5 + \
				6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),1.5) + (pow(f,-2 - \
				beta)*pow(-f + K,2)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 54*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 - 3*rho2) - \
				36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + 1350*rho4 - \
				675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + beta)*nu*rho*(beta*(-84 + \
				T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + \
				T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + 4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 \
				+ 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + \
				T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 \
				- 9*rho2 + beta*(-4 + 15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + \
				beta*(184 - 1380*rho2 + 1305*rho4))) + 3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + \
				120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + \
				120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 - 1260*rho2 + 1305*rho4)) + beta2*(960 \
				+ 120*T*nu2*(-4 + 45*rho2) + T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(40.*sqrt(3.)*alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5));
            break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"SABR_DetermineDerivatives : flag  bad input :");
			break;
		}
	}
    return value;
}

double SABR_ComputeAlphaFromSigmaATM(double f,
		                             double K,
		                             double tex,
		                             double sigma, 
		                             double beta,
		                             double rho, 
		                             double nu,
		                             int SABRflag)

{
	double alpha;
	if(fabs(beta -1.0)< ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE)
		alpha = SABR_BetaEqOne_CompatibleKernel_FromATMsigma_to_alpha(f,K,tex,sigma,rho,nu);
	else
		alpha = SABR_Direct_FromATMsigma_to_alpha( f,K,tex,sigma, beta, rho, nu, SABRflag);

	return alpha;
}


double SABR_ComputePartialDerivative( double f,
       double K,
       double T,
       double alpha,
       double beta,
       double rho, 
       double nu,
       const string& ParamStr,
       int SABRflag )
{
    double value = 0.0;

    ///Test if we ar aroud the money
    if (fabs(K-f)<ARM_GP_CF_SABR_SERIE_USAGE_LIMIT*f*alpha) 
	{
        value = SABR_ComputePartialDerivativeAroundATM(f,
            K,
            T,
            alpha,
            beta,rho,
            nu,ParamStr,
            SABRflag);

        return value;
    }
    
	double z        = (pow(f,1 - beta) - pow(K,1 - beta))/(alpha*(1 - beta));
	double x        = log((z*nu - rho + sqrt(1 + pow(z,2)*pow(nu,2) - 2*z*nu*rho))/(1 - rho))/nu;
	double theta    = pow(2,-1 - beta)*pow(f + K,-1 + beta)*pow(z,2)*alpha*beta*nu*rho + log((pow(f*K,beta/2.)*z*alpha)/(f - K)) + log((x*pow(1 + pow(z,2)*pow(nu,2) - 2*z*nu*rho,0.25))/z);
	double sqrtTheta=1 - (2*T*(theta - log((sqrt(f*K)*log(f/K))/(f - K))))/pow(x,2);
	if(sqrtTheta<=0)
		return value;
    else
	{
		//sig = log(f/K)/(x*sqrt(sqrtTheta));
		switch (SABRflag) 
		{
		case ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACT :
		case ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL :
			{
                double thetaDerx = 1/x;
                if(ParamStr == "Alpha")
                {
                    
                    double sigDerx          = -(log(f/K)/(sqrt(1 - (2*T*(theta - log((sqrt(f*K)*log(f/K))/(f - K))))/ \
                                                pow(x,2))*(pow(x,2) - 2*T*theta + 2*T*log((sqrt(f*K)*log(f/K))/(f - K)))));
                    double xDerz            = 1/sqrt(1 + pow(z,2)*pow(nu,2) - 2*z*nu*rho);
                    double zDeralpha        = (pow(f,1 - beta) - pow(K,1 - beta))/(pow(alpha,2)*(-1 + beta));
                    double sigDertheta      = (T*log(f/K))/(pow(x,3)*pow(1 - (2*T*(theta - \
				                                log((sqrt(f*K)*log(f/K))/(f - K))))/pow(x,2),1.5)); 
                    double thetaDerz        = (pow(2,-1 - beta)*nu*(pow(2,beta)*f*(z*nu - rho) + \
			                                    pow(2,beta)*K*(z*nu - rho) + 2*pow(f + K,beta)*z*alpha*beta*rho*(1 + \
			                                    pow(z,2)*pow(nu,2) - 2*z*nu*rho)))/((f + K)*(1 + pow(z,2)*pow(nu,2) - \
			                                    2*z*nu*rho));
                     double thetaDeralpha   = 1/alpha + pow(2,-1 - beta)*pow(f + K,-1 + \
				                                beta)*pow(z,2)*beta*nu*rho;

                    value = xDerz*zDeralpha*(sigDerx + sigDertheta*thetaDerx) + \
					        sigDertheta*(zDeralpha*thetaDerz + thetaDeralpha);
                    return value;
                }
                if(ParamStr == "Beta")
                {
                    double thetaDerx    = 1/x;
                    double xDerz        = 1/sqrt(1 + pow(z,2)*pow(nu,2) - 2*z*nu*rho);
                    double zDerbeta     = (-(pow(f,beta)*K) + f*pow(K,beta) + f*pow(K,beta)*(-1 + \
		                                     beta)*log(f) - pow(f,beta)*K*(-1 + \
		                                     beta)*log(K))/(pow(f,beta)*pow(K,beta)*alpha*pow(-1 + beta,2));
                    double sigDerx      = -(log(f/K)/(sqrt(1 - (2*T*(theta - log((sqrt(f*K)*log(f/K))/(f - K))))/ \
                                             pow(x,2))*(pow(x,2) - 2*T*theta + 2*T*log((sqrt(f*K)*log(f/K))/(f - K)))));
                    double sigDertheta  = (T*log(f/K))/(pow(x,3)*pow(1 - (2*T*(theta - \
				                            log((sqrt(f*K)*log(f/K))/(f - K))))/pow(x,2),1.5)); 
                    double thetaDerz    = (pow(2,-1 - beta)*nu*(pow(2,beta)*f*(z*nu - rho) + \
			                                pow(2,beta)*K*(z*nu - rho) + 2*pow(f + K,beta)*z*alpha*beta*rho*(1 + \
			                                pow(z,2)*pow(nu,2) - 2*z*nu*rho)))/((f + K)*(1 + pow(z,2)*pow(nu,2) - \
			                                2*z*nu*rho));
                    double thetaDerbeta = (pow(2.,-1 - beta)*(pow(2.,beta)*(f + K)*log(f*K) - pow(f + \
			                                K,beta)*pow(z,2)*alpha*nu*rho*(-1 + beta*log(2.) - beta*log(f + K))))/(f + K);

                    value = xDerz*zDerbeta*(sigDerx + sigDertheta*thetaDerx) + \
                                    sigDertheta*(zDerbeta*thetaDerz + thetaDerbeta);
                    return value;
                }
                if(ParamStr == "Correlation")
                {
                    double sigDertheta = (T*log(f/K))/(pow(x,3)*pow(1 - (2*T*(theta - \
                                            log((sqrt(f*K)*log(f/K))/(f - K))))/pow(x,2),1.5)); 
                    double sigDerx      = -(log(f/K)/(sqrt(1 - (2*T*(theta - log((sqrt(f*K)*log(f/K))/(f - K))))/ \
                                            pow(x,2))*(pow(x,2) - 2*T*theta + 2*T*log((sqrt(f*K)*log(f/K))/(f - K)))));
                    double thetaDerrho  = (z*nu*((pow(f + K,-1 + beta)*z*alpha*beta)/pow(2,beta) + 1/(-1 \
                                             - pow(z,2)*pow(nu,2) + 2*z*nu*rho)))/2.;
                    double xDerrho      = ((1 - rho)*(z*nu - rho + sqrt(1 + pow(z,2)*pow(nu,2) - 2*z*nu*rho) \
                                            + (1 - rho)*(-1 - (z*nu)/sqrt(1 + pow(z,2)*pow(nu,2) - \
                                            2*z*nu*rho))))/(nu*pow(-1 + rho,2)*(z*nu - rho + sqrt(1 + pow(z,2)*pow(nu,2) \
                                            - 2*z*nu*rho)));
                    value = xDerrho*(sigDerx + sigDertheta*thetaDerx) + \
                            sigDertheta*thetaDerrho;
                    return value;
                }
                if(ParamStr == "VolOfVol")
                {
                    double xDernu       = ((z*nu)/sqrt(1 + pow(z,2)*pow(nu,2) - 2*z*nu*rho) - log((z*nu - rho \
                                            + sqrt(1 + pow(z,2)*pow(nu,2) - 2*z*nu*rho))/(1 - rho)))/pow(nu,2);
                    double sigDerx      = -(log(f/K)/(sqrt(1 - (2*T*(theta - log((sqrt(f*K)*log(f/K))/(f - K))))/ \
                                            pow(x,2))*(pow(x,2) - 2*T*theta + 2*T*log((sqrt(f*K)*log(f/K))/(f - K)))));
                    double sigDertheta  = (T*log(f/K))/(pow(x,3)*pow(1 - (2*T*(theta - \
                                            log((sqrt(f*K)*log(f/K))/(f - K))))/pow(x,2),1.5)); 
                    double thetaDernu   = pow(2,-1 - beta)*pow(f + K,-1 + beta)*pow(z,2)*alpha*beta*rho + \
                                            (z*(z*nu - rho))/(2 + 2*pow(z,2)*pow(nu,2) - 4*z*nu*rho);
                    value = xDernu*(sigDerx + sigDertheta*thetaDerx) + \
                                sigDertheta*thetaDernu;
                    return value;
                }
                if(ParamStr == "Forward")
                {
                    double thetaDerz    = (pow(2,-1 - beta)*nu*(pow(2,beta)*f*(z*nu - rho) + \
                                            pow(2,beta)*K*(z*nu - rho) + 2*pow(f + K,beta)*z*alpha*beta*rho*(1 + \
                                            pow(z,2)*pow(nu,2) - 2*z*nu*rho)))/((f + K)*(1 + pow(z,2)*pow(nu,2) - \
                                            2*z*nu*rho));
                    double sigDertheta  = (T*log(f/K))/(pow(x,3)*pow(1 - (2*T*(theta - \
                                            log((sqrt(f*K)*log(f/K))/(f - K))))/pow(x,2),1.5)); 
                    double zDerf        = 1/(pow(f,beta)*alpha);
                    double thetaDerf    = (2*f - f*beta + K*beta)/(-2*pow(f,2) + 2*f*K) + pow(2,-1 - \
                                            beta)*pow(f + K,-2 + beta)*pow(z,2)*alpha*(-1 + beta)*beta*nu*rho;
                    double xDerz        = 1/sqrt(1 + pow(z,2)*pow(nu,2) - 2*z*nu*rho);
                    double sigDerf      = ((f + K)*T*log(f/K) - 2*(f - K)*(T - pow(x,2) + 2*T*theta - \
                                            2*T*log((sqrt(f*K)*log(f/K))/(f - K))))/(2.*f*(f - K)*x*sqrt(1 - (2*T*(theta \
                                            - log((sqrt(f*K)*log(f/K))/(f - K))))/pow(x,2))*(pow(x,2) - 2*T*theta + \
                                            2*T*log((sqrt(f*K)*log(f/K))/(f - K))));

                    double sigDerx      = -(log(f/K)/(sqrt(1 - (2*T*(theta - log((sqrt(f*K)*log(f/K))/(f - K))))/ \
                                            pow(x,2))*(pow(x,2) - 2*T*theta + 2*T*log((sqrt(f*K)*log(f/K))/(f - K)))));

                    value = sigDerf + xDerz*zDerf*(sigDerx + sigDertheta*thetaDerx) + \
                                            sigDertheta*(thetaDerf + zDerf*thetaDerz);
                    return value;
                }
                else
                        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"SABR_DetermineDerivatives : Param Str  bad input :");
                }
				
				break;
				
		case ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL2 :  
		case ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACTSTRIKE :  
			{
                double thetaDerx = 1/x;
				double sigDerx = -(log(f/K)/(sqrt(1 - (2*T*(theta - log((sqrt(f*K)*log(f/K))/(f - \
					                    K))))/pow(x,2))*(pow(x,2) - 2*T*theta + 2*T*log((sqrt(f*K)*log(f/K))/(f - \
					                    K)))));
				double sigDertheta = (T*log(f/K))/(pow(x,3)*pow(1 - (2*T*(theta - \
					                    log((sqrt(f*K)*log(f/K))/(f - K))))/pow(x,2),1.5));
				double sigDerf = ((f + K)*T*log(f/K) - 2*(f - K)*(T - pow(x,2) + 2*T*theta - \
					                    2*T*log((sqrt(f*K)*log(f/K))/(f - K))))/(2.*f*(f - K)*x*sqrt(1 - (2*T*(theta \
					                    - log((sqrt(f*K)*log(f/K))/(f - K))))/pow(x,2))*(pow(x,2) - 2*T*theta + \
					                    2*T*log((sqrt(f*K)*log(f/K))/(f - K))));
				
				double thetaDerz = (nu*(K*(z*nu - rho) + pow(K,beta)*z*alpha*beta*rho*(1 + \
					                    pow(z,2)*pow(nu,2) - 2*z*nu*rho)))/(2.*K*(1 + pow(z,2)*pow(nu,2) - \
					                    2*z*nu*rho));
				double thetaDeralpha = 1/alpha + (pow(K,-1 + beta)*pow(z,2)*beta*nu*rho)/4.;
				double thetaDerbeta = (pow(K,beta)*pow(z,2)*alpha*nu*rho + \
				                    	pow(K,beta)*pow(z,2)*alpha*beta*nu*rho*log(K) + 2*K*log(f*K))/(4.*K);
				double thetaDerrho = (z*nu*(pow(K,-1 + beta)*z*alpha*beta - 2/(1 + \
					                    pow(z,2)*pow(nu,2) - 2*z*nu*rho)))/4.;
				double thetaDernu = (z*(pow(K,-1 + beta)*z*alpha*beta*rho + (2*z*nu - 2*rho)/(1 + \
					                    pow(z,2)*pow(nu,2) - 2*z*nu*rho)))/4.;
				double thetaDerf = (2*f - f*beta + K*beta)/(-2*pow(f,2) + 2*f*K);
				double xDerz = 1/sqrt(1 + pow(z,2)*pow(nu,2) - 2*z*nu*rho);
			    double xDerrho = ((1 - rho)*(z*nu - rho + sqrt(1 + pow(z,2)*pow(nu,2) - 2*z*nu*rho) \
					                    + (1 - rho)*(-1 - (z*nu)/sqrt(1 + pow(z,2)*pow(nu,2) - \
					                    2*z*nu*rho))))/(nu*pow(-1 + rho,2)*(z*nu - rho + sqrt(1 + pow(z,2)*pow(nu,2) \
				                    	- 2*z*nu*rho)));
				double xDernu = ((z*nu)/sqrt(1 + pow(z,2)*pow(nu,2) - 2*z*nu*rho) - log((z*nu - rho \
					                    + sqrt(1 + pow(z,2)*pow(nu,2) - 2*z*nu*rho))/(1 - rho)))/pow(nu,2);
				double zDerf = 1/(pow(f,beta)*alpha);
				double zDeralpha = (pow(f,1 - beta) - pow(K,1 - beta))/(pow(alpha,2)*(-1 + beta));
				double zDerbeta = (-(pow(f,beta)*K) + f*pow(K,beta) + f*pow(K,beta)*(-1 + \
					                beta)*log(f) - pow(f,beta)*K*(-1 + \
					                beta)*log(K))/(pow(f,beta)*pow(K,beta)*alpha*pow(-1 + beta,2));

                if(ParamStr == "Alpha")
                {
                    value = xDerz*zDeralpha*(sigDerx + sigDertheta*thetaDerx) + \
					                sigDertheta*(zDeralpha*thetaDerz + thetaDeralpha);
                }
                if(ParamStr == "Beta")
                {
                    value = xDerz*zDerbeta*(sigDerx + sigDertheta*thetaDerx) + \
					                sigDertheta*(zDerbeta*thetaDerz + thetaDerbeta);
                }
                if(ParamStr == "Correlation")
                {
                    value = xDerrho*(sigDerx + sigDertheta*thetaDerx) + \
					                    sigDertheta*thetaDerrho;
                }
                if(ParamStr == "VolOfVol")
                {
                    value = xDernu*(sigDerx + sigDertheta*thetaDerx) + \
				                sigDertheta*thetaDernu;
                }
                if(ParamStr == "Forward")
                {
                    value = sigDerf + xDerz*zDerf*(sigDerx + sigDertheta*thetaDerx) + \
					                sigDertheta*(thetaDerf + zDerf*thetaDerz);
                }
                else
                   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"SABR_DetermineDerivatives : Param Str  bad input :");

                }
				break;
        default :
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"SABR_DetermineDerivatives : flag  bad input :");
            break;
        }
    }
    return value;
}

double SABR_ComputePartialDerivativeAroundATM( double f,
       double K,
       double T,
       double alpha,
       double beta,
       double rho, 
       double nu,
       const string& ParamStr,
       int SABRflag )
{
    double value = 0.0;
	double alpha2=alpha*alpha;
    double alpha3=alpha2*alpha;
    double alpha4=alpha3*alpha;
    double alpha5=alpha4*alpha;
    double alpha6=alpha5*alpha;
	double beta2=beta*beta;
    double beta3=beta2*beta;
    double beta4=beta3*beta;
    double beta5=beta4*beta;
    double beta6=beta5*beta;
	double rho2=rho*rho;
    double rho3=rho2*rho;
    double rho4=rho3*rho;
    double rho5=rho4*rho;
    double rho6=rho5*rho;
	double nu2=nu*nu;
    double nu3=nu2*nu;
    double nu4=nu3*nu;
    double nu5=nu4*nu;
    double nu6=nu5*nu;
	double f2=f*f;
    double f3=f2*f;
    double f4=f3*f;
    double f5=f4*f;
    double f6=f5*f;
	double T2=T*T;
    double T3=T2*T;
    double T4=T3*T;
    double T5=T4*T;
    double T6=T5*T;
    
    switch (SABRflag) 
	{
	case ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACT :
	case ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL :
		{
			if (ParamStr == "Alpha")
            {
			    value = ((240*pow(f,beta)*(f - K)*(-1 + \
				beta)*(6*pow(f,beta)*T*alpha*beta*nu*rho + f*(-12 + T*nu2*(2 - \
				3*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),1.5) + \
				(480*pow(f,-3 + 2*beta)*T*alpha*(pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho))/pow((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 \
				+ beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2,1.5) + \
				(480*pow(f,-1 + beta))/sqrt((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2) + \
				(720*pow(f,beta)*(f - K)*T*(pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho)*(3*pow(f,2*beta)*T*alpha2*(-1 + beta)*beta*nu*rho + \
				f2*nu*rho*(-12 + T*nu2*(5 - 6*rho2)) - pow(f,1 + beta)*alpha*(-1 + beta)*(12 \
				+ T*nu2*(-2 + 3*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5) + \
				(12*pow(f - K,2)*(pow(f,5*beta)*T2*alpha5*pow(-1 + beta,6) - 15*pow(f,1 + \
				4*beta)*T2*alpha4*pow(-1 + beta,3)*beta*(7 + 8*beta)*nu*rho + 2*pow(f,2 + \
				3*beta)*T*alpha3*(-1 + beta)*(-36 + 3*T*nu2*(2 - 3*rho2) + beta2*(-588 + \
				T*nu2*(98 - 657*rho2)) + 29*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta3*(276 + \
				T*nu2*(-46 + 39*rho2))) + 15*pow(f,3 + 2*beta)*T*alpha2*(-1 + \
				beta)*nu*rho*(10*beta*(12 + T*nu2) + 4*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(84 \
				+ T*nu2*(-32 + 39*rho2))) + f5*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + \
				beta*(-4 + 15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - \
				1380*rho2 + 1305*rho4))) + pow(f,4 + beta)*alpha*(3*(1280 + 120*T*nu2*(-4 + \
				7*rho2) + T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + 240*T*nu2*(-4 + \
				9*rho2) + T2*nu4*(128 - 810*rho2 + 765*rho4)) + beta2*(960 + 120*T*nu2*(-4 + \
				27*rho2) + T2*nu4*(88 - 1440*rho2 + \
				1575*rho4)))))/(f2*alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5)) - \
				(10*pow(f - K,2)*T*(-(pow(f,beta)*alpha*pow(-1 + beta,2)) - \
				3*f*beta*nu*rho)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 18*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,3)*beta*(7 + 8*beta)*nu*rho + f6*nu2*(2880*(2 \
				- 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + \
				1350*rho4 - 675*rho6)) + 3*pow(f,2 + 4*beta)*T*alpha4*(-1 + beta)*(-36 + \
				3*T*nu2*(2 - 3*rho2) + beta2*(-588 + T*nu2*(98 - 657*rho2)) + 29*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta3*(276 + T*nu2*(-46 + 39*rho2))) + 30*pow(f,3 + \
				3*beta)*T*alpha3*(-1 + beta)*nu*rho*(10*beta*(12 + T*nu2) + 4*(12 + T*nu2*(-5 \
				+ 6*rho2)) + beta2*(84 + T*nu2*(-32 + 39*rho2))) + 6*pow(f,5 + \
				beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 240*T*nu2*(-4 + 9*rho2) + T2*nu4*(128 \
				- 810*rho2 + 765*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 27*rho2) + T2*nu4*(88 \
				- 1440*rho2 + 1575*rho4)))))/(f2*alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),3.5)) - (2*pow(f,-2 - beta)*pow(f - \
				K,2)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 18*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,3)*beta*(7 + 8*beta)*nu*rho + f6*nu2*(2880*(2 \
				- 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + \
				1350*rho4 - 675*rho6)) + 3*pow(f,2 + 4*beta)*T*alpha4*(-1 + beta)*(-36 + \
				3*T*nu2*(2 - 3*rho2) + beta2*(-588 + T*nu2*(98 - 657*rho2)) + 29*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta3*(276 + T*nu2*(-46 + 39*rho2))) + 30*pow(f,3 + \
				3*beta)*T*alpha3*(-1 + beta)*nu*rho*(10*beta*(12 + T*nu2) + 4*(12 + T*nu2*(-5 \
				+ 6*rho2)) + beta2*(84 + T*nu2*(-32 + 39*rho2))) + 6*pow(f,5 + \
				beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 240*T*nu2*(-4 + 9*rho2) + T2*nu4*(128 \
				- 810*rho2 + 765*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 27*rho2) + T2*nu4*(88 \
				- 1440*rho2 + 1575*rho4)))))/(alpha2*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5)))/(80.*sqrt(3.));
                return value;
            }

			if (ParamStr == "Beta")
            {
			    value = ((480*pow(f,-1 + \
				beta)*alpha*log(f))/sqrt((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2) - \
				(2*pow(f,-2 - beta)*pow(f - K,2)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - \
				18*pow(f,1 + 5*beta)*T2*alpha5*pow(-1 + beta,3)*beta*(7 + 8*beta)*nu*rho + \
				f6*nu2*(2880*(2 - 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 \
				- 1062*rho2 + 1350*rho4 - 675*rho6)) + 3*pow(f,2 + 4*beta)*T*alpha4*(-1 + \
				beta)*(-36 + 3*T*nu2*(2 - 3*rho2) + beta2*(-588 + T*nu2*(98 - 657*rho2)) + \
				29*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta3*(276 + T*nu2*(-46 + 39*rho2))) + \
				30*pow(f,3 + 3*beta)*T*alpha3*(-1 + beta)*nu*rho*(10*beta*(12 + T*nu2) + \
				4*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(84 + T*nu2*(-32 + 39*rho2))) + \
				6*pow(f,5 + beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + \
				15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + \
				1305*rho4))) + 3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) \
				+ T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + 240*T*nu2*(-4 + 9*rho2) \
				+ T2*nu4*(128 - 810*rho2 + 765*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 27*rho2) \
				+ T2*nu4*(88 - 1440*rho2 + \
				1575*rho4))))*log(f))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5)) + \
				(480*pow(f,-3 + 2*beta)*T*alpha2*(pow(f,beta)*alpha*(-1 + beta) + 3*f*nu*rho \
				+ (pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho)*log(f)))/pow((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)))/f2,1.5) - (720*pow(f,beta)*(f - \
				K)*T*alpha*(-3*pow(f,2*beta)*T*alpha2*(-1 + beta)*beta*nu*rho + pow(f,1 + \
				beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 + 3*rho2)) + f2*nu*rho*(12 + \
				T*nu2*(-5 + 6*rho2)))*(pow(f,beta)*alpha*(-1 + beta) + 3*f*nu*rho + \
				(pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho)*log(f)))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5) + \
				(10*pow(f - K,2)*T*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 18*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,3)*beta*(7 + 8*beta)*nu*rho + f6*nu2*(2880*(2 \
				- 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + \
				1350*rho4 - 675*rho6)) + 3*pow(f,2 + 4*beta)*T*alpha4*(-1 + beta)*(-36 + \
				3*T*nu2*(2 - 3*rho2) + beta2*(-588 + T*nu2*(98 - 657*rho2)) + 29*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta3*(276 + T*nu2*(-46 + 39*rho2))) + 30*pow(f,3 + \
				3*beta)*T*alpha3*(-1 + beta)*nu*rho*(10*beta*(12 + T*nu2) + 4*(12 + T*nu2*(-5 \
				+ 6*rho2)) + beta2*(84 + T*nu2*(-32 + 39*rho2))) + 6*pow(f,5 + \
				beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 240*T*nu2*(-4 + 9*rho2) + T2*nu4*(128 \
				- 810*rho2 + 765*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 27*rho2) + T2*nu4*(88 \
				- 1440*rho2 + 1575*rho4))))*(pow(f,beta)*alpha*(-1 + beta) + 3*f*nu*rho + \
				(pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho)*log(f)))/(f2*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) \
				- 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),3.5)) \
				+ (240*pow(f,beta)*(f - K)*alpha*(3*pow(f,beta)*T*alpha*(-1 + 2*beta)*nu*rho \
				+ f*(-12 + T*nu2*(2 - 3*rho2)) + (-1 + \
				beta)*(6*pow(f,beta)*T*alpha*beta*nu*rho + f*(-12 + T*nu2*(2 - \
				3*rho2)))*log(f)))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 \
				+ beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),1.5) + (6*pow(f - \
				K,2)*(2*pow(f,5*beta)*T2*alpha5*pow(-1 + beta,5) - 48*pow(f,1 + \
				4*beta)*T2*alpha4*pow(-1 + beta,3)*beta*nu*rho - 6*pow(f,1 + \
				4*beta)*T2*alpha4*pow(-1 + beta,3)*(7 + 8*beta)*nu*rho - 18*pow(f,1 + \
				4*beta)*T2*alpha4*pow(-1 + beta,2)*beta*(7 + 8*beta)*nu*rho + \
				2*f5*T*nu3*rho*(120*(-4 + 15*rho2) + T*nu2*(184 - 1380*rho2 + 1305*rho4)) + \
				pow(f,2 + 3*beta)*T*alpha3*(-36 + 3*T*nu2*(2 - 3*rho2) + beta2*(-588 + \
				T*nu2*(98 - 657*rho2)) + 29*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta3*(276 + \
				T*nu2*(-46 + 39*rho2))) + 20*pow(f,3 + 2*beta)*T*alpha2*(-1 + \
				beta)*nu*rho*(60 + 5*T*nu2 + beta*(84 + T*nu2*(-32 + 39*rho2))) + 10*pow(f,3 \
				+ 2*beta)*T*alpha2*nu*rho*(10*beta*(12 + T*nu2) + 4*(12 + T*nu2*(-5 + \
				6*rho2)) + beta2*(84 + T*nu2*(-32 + 39*rho2))) + pow(f,2 + \
				3*beta)*T*alpha3*(-1 + beta)*(29*(12 + T*nu2*(-2 + 3*rho2)) + 3*beta2*(276 + \
				T*nu2*(-46 + 39*rho2)) - 2*beta*(588 + T*nu2*(-98 + 657*rho2))) + pow(f,4 + \
				beta)*alpha*(-2*(2400 + 240*T*nu2*(-4 + 9*rho2) + T2*nu4*(128 - 810*rho2 + \
				765*rho4)) + 2*beta*(960 + 120*T*nu2*(-4 + 27*rho2) + T2*nu4*(88 - 1440*rho2 \
				+ 1575*rho4))) + 2*pow(f,5*beta)*T2*alpha5*pow(-1 + beta,6)*log(f) - \
				30*pow(f,1 + 4*beta)*T2*alpha4*pow(-1 + beta,3)*beta*(7 + \
				8*beta)*nu*rho*log(f) + 4*pow(f,2 + 3*beta)*T*alpha3*(-1 + beta)*(-36 + \
				3*T*nu2*(2 - 3*rho2) + beta2*(-588 + T*nu2*(98 - 657*rho2)) + 29*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta3*(276 + T*nu2*(-46 + 39*rho2)))*log(f) + \
				30*pow(f,3 + 2*beta)*T*alpha2*(-1 + beta)*nu*rho*(10*beta*(12 + T*nu2) + \
				4*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(84 + T*nu2*(-32 + 39*rho2)))*log(f) + \
				2*f5*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + \
				1305*rho4)))*log(f) + 2*pow(f,4 + beta)*alpha*(3*(1280 + 120*T*nu2*(-4 + \
				7*rho2) + T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + 240*T*nu2*(-4 + \
				9*rho2) + T2*nu4*(128 - 810*rho2 + 765*rho4)) + beta2*(960 + 120*T*nu2*(-4 + \
				27*rho2) + T2*nu4*(88 - 1440*rho2 + \
				1575*rho4)))*log(f)))/(f2*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5)))/(80.*sqrt(3.));
                return value;
            }

			if (ParamStr == "Correlation")
            {
			    value = (sqrt(3.)*nu*((240*pow(f,-2 + \
				beta)*T*alpha*(pow(f,beta)*alpha*beta - \
				f*nu*rho))/pow((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2,1.5) + \
				(360*f*(f - K)*T*(-(pow(f,beta)*alpha*beta) + \
				f*nu*rho)*(-3*pow(f,2*beta)*T*alpha2*(-1 + beta)*beta*nu*rho + pow(f,1 + \
				beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 + 3*rho2)) + f2*nu*rho*(12 + \
				T*nu2*(-5 + 6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5) - \
				(40*(f - K)*(-3*pow(f,2*beta)*T*alpha2*(-1 + beta)*beta + 6*pow(f,1 + \
				beta)*T*alpha*(-1 + beta)*nu*rho + f2*(12 + T*nu2*(-5 + \
				18*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),1.5) - (5*pow(f,-1 \
				- beta)*pow(f - K,2)*T*(-(pow(f,beta)*alpha*beta) + \
				f*nu*rho)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 18*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,3)*beta*(7 + 8*beta)*nu*rho + f6*nu2*(2880*(2 \
				- 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + \
				1350*rho4 - 675*rho6)) + 3*pow(f,2 + 4*beta)*T*alpha4*(-1 + beta)*(-36 + \
				3*T*nu2*(2 - 3*rho2) + beta2*(-588 + T*nu2*(98 - 657*rho2)) + 29*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta3*(276 + T*nu2*(-46 + 39*rho2))) + 30*pow(f,3 + \
				3*beta)*T*alpha3*(-1 + beta)*nu*rho*(10*beta*(12 + T*nu2) + 4*(12 + T*nu2*(-5 \
				+ 6*rho2)) + beta2*(84 + T*nu2*(-32 + 39*rho2))) + 6*pow(f,5 + \
				beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 240*T*nu2*(-4 + 9*rho2) + T2*nu4*(128 \
				- 810*rho2 + 765*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 27*rho2) + T2*nu4*(88 \
				- 1440*rho2 + 1575*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),3.5)) - (2*pow(f,-1 - beta)*pow(f - \
				K,2)*(3*pow(f,5*beta)*T2*alpha5*pow(-1 + beta,3)*beta*(7 + 8*beta) - \
				3*pow(f,1 + 4*beta)*T2*alpha4*(3 - 32*beta + 248*beta2 - 232*beta3 + \
				13*beta4)*nu*rho + 3*f5*nu*rho*(960 + 600*T*nu2*(-2 + 3*rho2) + T2*nu4*(118 - \
				300*rho2 + 225*rho4)) - 30*pow(f,3 + 2*beta)*T*alpha2*nu*rho*(84 + T*nu2*(-26 \
				+ 45*rho2) - 6*beta*(24 + T*nu2*(-9 + 17*rho2)) + 3*beta2*(36 + T*nu2*(-16 + \
				35*rho2))) - 5*pow(f,2 + 3*beta)*T*alpha3*(-1 + beta)*(48 + 10*beta*(12 + \
				T*nu2) + 4*T*nu2*(-5 + 18*rho2) + beta2*(84 + T*nu2*(-32 + 117*rho2))) - \
				pow(f,4 + beta)*alpha*(-1440 + 120*T*nu2*(7 - 27*rho2 + beta*(-4 + 45*rho2)) \
				+ T2*nu4*(-10*(10 - 81*rho2 + 90*rho4) + beta*(184 - 4140*rho2 + \
				6525*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5))))/40.;
                return value;
            }

			if (ParamStr == "VolOfVol")
            {
			    value = ((480*pow(f,-2 + beta)*T*alpha*(3*pow(f,beta)*alpha*beta*rho + \
				f*nu*(2 - 3*rho2)))/pow((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)))/f2,1.5) - (240*(f - K)*(-3*pow(f,2*beta)*T*alpha2*(-1 + \
				beta)*beta*rho + 2*pow(f,1 + beta)*T*alpha*(-1 + beta)*nu*(-2 + 3*rho2) + \
				3*f2*rho*(4 + T*nu2*(-5 + 6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),1.5) + (720*f*(f - K)*T*(-3*pow(f,beta)*alpha*beta*rho + f*nu*(-2 + \
				3*rho2))*(-3*pow(f,2*beta)*T*alpha2*(-1 + beta)*beta*nu*rho + pow(f,1 + \
				beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 + 3*rho2)) + f2*nu*rho*(12 + \
				T*nu2*(-5 + 6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5) - \
				(12*pow(f,-1 - beta)*pow(f - K,2)*(3*pow(f,5*beta)*T2*alpha5*pow(-1 + \
				beta,3)*beta*(7 + 8*beta)*rho - pow(f,1 + 4*beta)*T2*alpha4*(-1 + beta)*nu*(6 \
				- 9*rho2 + beta2*(98 - 657*rho2) + 29*beta*(-2 + 3*rho2) + beta3*(-46 + \
				39*rho2)) + f5*nu*(960*(-2 + 3*rho2) + 24*T*nu2*(88 - 300*rho2 + 225*rho4) + \
				T2*nu4*(-368 + 1062*rho2 - 1350*rho4 + 675*rho6)) - 15*pow(f,2 + \
				3*beta)*T*alpha3*(-1 + beta)*rho*(10*beta*(4 + T*nu2) + 4*(4 + T*nu2*(-5 + \
				6*rho2)) + beta2*(28 + T*nu2*(-32 + 39*rho2))) - 5*pow(f,4 + \
				beta)*alpha*rho*(-288 + 72*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) - \
				2*pow(f,3 + 2*beta)*T*alpha2*nu*(3*(60*(-4 + 7*rho2) + T*nu2*(56 - 260*rho2 + \
				225*rho4)) - 2*beta*(120*(-4 + 9*rho2) + T*nu2*(128 - 810*rho2 + 765*rho4)) + \
				beta2*(60*(-4 + 27*rho2) + T*nu2*(88 - 1440*rho2 + \
				1575*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5)) - \
				(10*pow(f,-1 - beta)*pow(f - K,2)*T*(-3*pow(f,beta)*alpha*beta*rho + f*nu*(-2 \
				+ 3*rho2))*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 18*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,3)*beta*(7 + 8*beta)*nu*rho + f6*nu2*(2880*(2 \
				- 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + \
				1350*rho4 - 675*rho6)) + 3*pow(f,2 + 4*beta)*T*alpha4*(-1 + beta)*(-36 + \
				3*T*nu2*(2 - 3*rho2) + beta2*(-588 + T*nu2*(98 - 657*rho2)) + 29*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta3*(276 + T*nu2*(-46 + 39*rho2))) + 30*pow(f,3 + \
				3*beta)*T*alpha3*(-1 + beta)*nu*rho*(10*beta*(12 + T*nu2) + 4*(12 + T*nu2*(-5 \
				+ 6*rho2)) + beta2*(84 + T*nu2*(-32 + 39*rho2))) + 6*pow(f,5 + \
				beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 240*T*nu2*(-4 + 9*rho2) + T2*nu4*(128 \
				- 810*rho2 + 765*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 27*rho2) + T2*nu4*(88 \
				- 1440*rho2 + 1575*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),3.5)))/(80.*sqrt(3.));
                return value;
            }

			if (ParamStr == "Forward")
            {
			    value = ((480*pow(f,-4 + 2*beta)*T*alpha2*(-1 + \
				beta)*(pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho))/pow((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 \
				+ beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2,1.5) + \
				(480*pow(f,-2 + beta)*alpha*(-1 + \
				beta))/sqrt((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2) - \
				(240*f*(pow(f,beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 + 3*rho2)) + \
				f*nu*rho*(12 + T*nu2*(-5 + 6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),1.5) + (240*(-f + K)*(pow(f,beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 \
				+ 3*rho2)) + f*nu*rho*(12 + T*nu2*(-5 + \
				6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),1.5) + (720*(f - \
				K)*(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)*beta) - 3*pow(f,1 + \
				beta)*T*alpha*beta*(1 + beta)*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)))*(pow(f,beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 + 3*rho2)) + \
				f*nu*rho*(12 + T*nu2*(-5 + 6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5) - (240*(f - K)*(pow(f,beta)*alpha*(-1 + beta)*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + f*nu*rho*(12 + T*nu2*(-5 + \
				6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),1.5) - (4*pow(f,-2 \
				- beta)*(-f + K)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 54*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 - 3*rho2) - \
				36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + 1350*rho4 - \
				675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + beta)*nu*rho*(beta*(-84 + \
				T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + \
				T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + 4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 \
				+ 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + \
				T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 \
				- 9*rho2 + beta*(-4 + 15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + \
				beta*(184 - 1380*rho2 + 1305*rho4))) + 3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + \
				120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + \
				120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 - 1260*rho2 + 1305*rho4)) + beta2*(960 \
				+ 120*T*nu2*(-4 + 45*rho2) + T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5)) + \
				(2*pow(f,-3 - beta)*pow(f - K,2)*(-2 - beta)*(pow(f,6*beta)*T2*alpha6*pow(-1 \
				+ beta,6) - 54*pow(f,1 + 5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + \
				f6*nu2*(2880*(2 - 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 \
				- 1062*rho2 + 1350*rho4 - 675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + \
				beta)*nu*rho*(beta*(-84 + T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + \
				6*rho2)) + beta2*(132 + T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + \
				4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 + 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta2*(-276 + T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + \
				beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 \
				- 1260*rho2 + 1305*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 45*rho2) + \
				T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5)) - \
				(5*pow(f,-2 - beta)*pow(f - K,2)*(-2*pow(f,-1 + 2*beta)*T*alpha2*pow(-1 + \
				beta,2)*beta - 6*pow(f,beta)*T*alpha*beta*(1 + beta)*nu*rho + f*(24 + \
				2*T*nu2*(-2 + 3*rho2)))*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - \
				54*pow(f,1 + 5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 \
				- 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + \
				1350*rho4 - 675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + \
				beta)*nu*rho*(beta*(-84 + T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + \
				6*rho2)) + beta2*(132 + T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + \
				4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 + 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta2*(-276 + T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + \
				beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 \
				- 1260*rho2 + 1305*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 45*rho2) + \
				T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),3.5)) + \
				(2*pow(f,-2 - beta)*pow(f - K,2)*(6*pow(f,-1 + 6*beta)*T2*alpha6*pow(-1 + \
				beta,6)*beta - 54*pow(f,5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*(1 + \
				5*beta)*nu*rho + 6*f5*nu2*(2880*(2 - 3*rho2) - 36*T*nu2*(88 - 300*rho2 + \
				225*rho4) + T2*nu4*(368 - 1062*rho2 + 1350*rho4 - 675*rho6)) + 180*pow(f,2 + \
				3*beta)*T*alpha3*(-1 + beta)*(1 + beta)*nu*rho*(beta*(-84 + T*nu2*(29 - \
				36*rho2)) + 2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + T*nu2*(-31 + \
				42*rho2))) - 3*pow(f,1 + 4*beta)*T*alpha4*pow(-1 + beta,2)*(2 + 4*beta)*(-36 \
				+ 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + \
				T*nu2*(46 + 51*rho2))) + 6*pow(f,4 + beta)*alpha*(5 + beta)*nu*rho*(-1440 + \
				120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + \
				18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + 6*pow(f,3 + \
				2*beta)*alpha2*(2 + beta)*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 \
				- 1260*rho2 + 1305*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 45*rho2) + \
				T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5)))/(80.*sqrt(3.));
                return value;
            }
            else
                throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"SABR_DetermineDerivatives : Param Str  bad input :");
            break;
			
		}
	case ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL2 :  
	case ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACTSTRIKE :  
		{
            if (ParamStr == "Forward")
            {
			    value = ((480*pow(f,-4 + 2*beta)*T*alpha2*(-1 + \
				beta)*(pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho))/pow((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 \
				+ beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2,1.5) + \
				(480*pow(f,-2 + beta)*alpha*(-1 + \
				beta))/sqrt((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2) - \
				(240*f*(pow(f,beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 + 3*rho2)) + \
				f*nu*rho*(12 + T*nu2*(-5 + 6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),1.5) + (240*(-f + K)*(pow(f,beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 \
				+ 3*rho2)) + f*nu*rho*(12 + T*nu2*(-5 + \
				6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),1.5) + (720*(f - \
				K)*(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)*beta) - 3*pow(f,1 + \
				beta)*T*alpha*beta*(1 + beta)*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)))*(pow(f,beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 + 3*rho2)) + \
				f*nu*rho*(12 + T*nu2*(-5 + 6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5) - (240*(f - K)*(pow(f,beta)*alpha*(-1 + beta)*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + f*nu*rho*(12 + T*nu2*(-5 + \
				6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),1.5) - (4*pow(f,-2 \
				- beta)*(-f + K)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 54*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 - 3*rho2) - \
				36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + 1350*rho4 - \
				675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + beta)*nu*rho*(beta*(-84 + \
				T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + \
				T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + 4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 \
				+ 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + \
				T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 \
				- 9*rho2 + beta*(-4 + 15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + \
				beta*(184 - 1380*rho2 + 1305*rho4))) + 3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + \
				120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + \
				120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 - 1260*rho2 + 1305*rho4)) + beta2*(960 \
				+ 120*T*nu2*(-4 + 45*rho2) + T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5)) + \
				(2*pow(f,-3 - beta)*pow(f - K,2)*(-2 - beta)*(pow(f,6*beta)*T2*alpha6*pow(-1 \
				+ beta,6) - 54*pow(f,1 + 5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + \
				f6*nu2*(2880*(2 - 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 \
				- 1062*rho2 + 1350*rho4 - 675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + \
				beta)*nu*rho*(beta*(-84 + T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + \
				6*rho2)) + beta2*(132 + T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + \
				4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 + 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta2*(-276 + T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + \
				beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 \
				- 1260*rho2 + 1305*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 45*rho2) + \
				T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5)) - \
				(5*pow(f,-2 - beta)*pow(f - K,2)*(-2*pow(f,-1 + 2*beta)*T*alpha2*pow(-1 + \
				beta,2)*beta - 6*pow(f,beta)*T*alpha*beta*(1 + beta)*nu*rho + f*(24 + \
				2*T*nu2*(-2 + 3*rho2)))*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - \
				54*pow(f,1 + 5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 \
				- 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + \
				1350*rho4 - 675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + \
				beta)*nu*rho*(beta*(-84 + T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + \
				6*rho2)) + beta2*(132 + T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + \
				4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 + 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta2*(-276 + T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + \
				beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 \
				- 1260*rho2 + 1305*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 45*rho2) + \
				T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),3.5)) + \
				(2*pow(f,-2 - beta)*pow(f - K,2)*(6*pow(f,-1 + 6*beta)*T2*alpha6*pow(-1 + \
				beta,6)*beta - 54*pow(f,5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*(1 + \
				5*beta)*nu*rho + 6*f5*nu2*(2880*(2 - 3*rho2) - 36*T*nu2*(88 - 300*rho2 + \
				225*rho4) + T2*nu4*(368 - 1062*rho2 + 1350*rho4 - 675*rho6)) + 180*pow(f,2 + \
				3*beta)*T*alpha3*(-1 + beta)*(1 + beta)*nu*rho*(beta*(-84 + T*nu2*(29 - \
				36*rho2)) + 2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + T*nu2*(-31 + \
				42*rho2))) - 3*pow(f,1 + 4*beta)*T*alpha4*pow(-1 + beta,2)*(2 + 4*beta)*(-36 \
				+ 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + \
				T*nu2*(46 + 51*rho2))) + 6*pow(f,4 + beta)*alpha*(5 + beta)*nu*rho*(-1440 + \
				120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + \
				18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + 6*pow(f,3 + \
				2*beta)*alpha2*(2 + beta)*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 \
				- 1260*rho2 + 1305*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 45*rho2) + \
				T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5)))/(80.*sqrt(3.));
                return value;
            }

			if (ParamStr == "Alpha")
            {
			    value = ((240*pow(f,4 + beta)*(-f + K)*(-1 + beta)*(12 + T*nu2*(-2 + \
				3*rho2)))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),1.5) + \
				(480*pow(f,2*beta)*T*alpha*(pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho))/pow((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 \
				+ beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2,1.5) + \
				(480*pow(f,2 + beta))/sqrt((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2) - \
				(720*pow(f,4 + beta)*(f - K)*T*(pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho)*(pow(f,beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 + 3*rho2)) + \
				f*nu*rho*(12 + T*nu2*(-5 + 6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5) + (12*f*pow(f - K,2)*(pow(f,5*beta)*T2*alpha5*pow(-1 + beta,6) \
				- 45*pow(f,1 + 4*beta)*T2*alpha4*pow(-1 + beta,4)*beta*nu*rho + 30*pow(f,3 + \
				2*beta)*T*alpha2*(-1 + beta)*nu*rho*(beta*(-84 + T*nu2*(29 - 36*rho2)) + \
				2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + T*nu2*(-31 + 42*rho2))) - \
				2*pow(f,2 + 3*beta)*T*alpha3*pow(-1 + beta,2)*(-36 + 3*T*nu2*(2 - 3*rho2) + \
				26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + T*nu2*(46 + 51*rho2))) + \
				f5*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				pow(f,4 + beta)*alpha*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 \
				- 1260*rho2 + 1305*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 45*rho2) + \
				T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5)) - \
				(10*f*pow(f - K,2)*T*(-(pow(f,beta)*alpha*pow(-1 + beta,2)) - \
				3*f*beta*nu*rho)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 54*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 - 3*rho2) - \
				36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + 1350*rho4 - \
				675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + beta)*nu*rho*(beta*(-84 + \
				T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + \
				T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + 4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 \
				+ 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + \
				T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 \
				- 9*rho2 + beta*(-4 + 15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + \
				beta*(184 - 1380*rho2 + 1305*rho4))) + 3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + \
				120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + \
				120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 - 1260*rho2 + 1305*rho4)) + beta2*(960 \
				+ 120*T*nu2*(-4 + 45*rho2) + T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),3.5)) - \
				(2*pow(f,1 - beta)*pow(f - K,2)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - \
				54*pow(f,1 + 5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 \
				- 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + \
				1350*rho4 - 675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + \
				beta)*nu*rho*(beta*(-84 + T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + \
				6*rho2)) + beta2*(132 + T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + \
				4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 + 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta2*(-276 + T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + \
				beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 \
				- 1260*rho2 + 1305*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 45*rho2) + \
				T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha2*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5)))/(80.*sqrt(3.)*f3);
                return value;
            }

			if (ParamStr == "Beta")
            {
			    value = ((480*pow(f,2 + \
				beta)*alpha*log(f))/sqrt((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2) - \
				(2*pow(f,1 - beta)*pow(f - K,2)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - \
				54*pow(f,1 + 5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 \
				- 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + \
				1350*rho4 - 675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + \
				beta)*nu*rho*(beta*(-84 + T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + \
				6*rho2)) + beta2*(132 + T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + \
				4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 + 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta2*(-276 + T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + \
				beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 \
				- 1260*rho2 + 1305*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 45*rho2) + \
				T2*nu4*(88 - 2340*rho2 + \
				2655*rho4))))*log(f))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5)) - \
				(240*pow(f,4 + beta)*(f - K)*alpha*(12 + T*nu2*(-2 + 3*rho2))*(1 + (-1 + \
				beta)*log(f)))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),1.5) + \
				(480*pow(f,2*beta)*T*alpha2*(pow(f,beta)*alpha*(-1 + beta) + 3*f*nu*rho + \
				(pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho)*log(f)))/pow((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)))/f2,1.5) - (720*pow(f,4 + beta)*(f - \
				K)*T*alpha*(pow(f,beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 + 3*rho2)) + \
				f*nu*rho*(12 + T*nu2*(-5 + 6*rho2)))*(pow(f,beta)*alpha*(-1 + beta) + \
				3*f*nu*rho + (pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho)*log(f)))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5) + \
				(10*f*pow(f - K,2)*T*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 54*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 - 3*rho2) - \
				36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + 1350*rho4 - \
				675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + beta)*nu*rho*(beta*(-84 + \
				T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + \
				T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + 4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 \
				+ 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + \
				T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 \
				- 9*rho2 + beta*(-4 + 15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + \
				beta*(184 - 1380*rho2 + 1305*rho4))) + 3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + \
				120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + \
				120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 - 1260*rho2 + 1305*rho4)) + beta2*(960 \
				+ 120*T*nu2*(-4 + 45*rho2) + T2*nu4*(88 - 2340*rho2 + \
				2655*rho4))))*(pow(f,beta)*alpha*(-1 + beta) + 3*f*nu*rho + \
				(pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho)*log(f)))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),3.5) + \
				(6*f*pow(f - K,2)*(2*pow(f,5*beta)*T2*alpha5*pow(-1 + beta,5) - 18*pow(f,1 + \
				4*beta)*T2*alpha4*pow(-1 + beta,4)*nu*rho - 72*pow(f,1 + \
				4*beta)*T2*alpha4*pow(-1 + beta,3)*beta*nu*rho + 2*f5*T*nu3*rho*(120*(-4 + \
				15*rho2) + T*nu2*(184 - 1380*rho2 + 1305*rho4)) + 20*pow(f,3 + \
				2*beta)*T*alpha2*nu*rho*(beta*(-84 + T*nu2*(29 - 36*rho2)) + 2*(12 + \
				T*nu2*(-5 + 6*rho2)) + beta2*(132 + T*nu2*(-31 + 42*rho2))) + 20*pow(f,3 + \
				2*beta)*T*alpha2*(-1 + beta)*nu*rho*(-84 + T*nu2*(29 - 36*rho2) + beta*(264 + \
				2*T*nu2*(-31 + 42*rho2))) - pow(f,2 + 3*beta)*T*alpha3*pow(-1 + \
				beta,2)*(26*(12 + T*nu2*(-2 + 3*rho2)) + 2*beta*(-276 + T*nu2*(46 + \
				51*rho2))) - 2*pow(f,2 + 3*beta)*T*alpha3*(-1 + beta)*(-36 + 3*T*nu2*(2 - \
				3*rho2) + 26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + T*nu2*(46 + \
				51*rho2))) + pow(f,4 + beta)*alpha*(-2*(2400 + 120*T*nu2*(-8 + 27*rho2) + \
				T2*nu4*(128 - 1260*rho2 + 1305*rho4)) + 2*beta*(960 + 120*T*nu2*(-4 + \
				45*rho2) + T2*nu4*(88 - 2340*rho2 + 2655*rho4))) + \
				2*pow(f,5*beta)*T2*alpha5*pow(-1 + beta,6)*log(f) - 90*pow(f,1 + \
				4*beta)*T2*alpha4*pow(-1 + beta,4)*beta*nu*rho*log(f) + 60*pow(f,3 + \
				2*beta)*T*alpha2*(-1 + beta)*nu*rho*(beta*(-84 + T*nu2*(29 - 36*rho2)) + \
				2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + T*nu2*(-31 + 42*rho2)))*log(f) - \
				4*pow(f,2 + 3*beta)*T*alpha3*pow(-1 + beta,2)*(-36 + 3*T*nu2*(2 - 3*rho2) + \
				26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + T*nu2*(46 + \
				51*rho2)))*log(f) + 2*f5*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + \
				15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + \
				1305*rho4)))*log(f) + 2*pow(f,4 + beta)*alpha*(3*(1280 + 120*T*nu2*(-4 + \
				7*rho2) + T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + 120*T*nu2*(-8 + \
				27*rho2) + T2*nu4*(128 - 1260*rho2 + 1305*rho4)) + beta2*(960 + 120*T*nu2*(-4 \
				+ 45*rho2) + T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))*log(f)))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5))/(80.*sqrt(3.)*f3);
                return value;
            }

			if (ParamStr == "Correlation")
            {
			    value = (sqrt(3.)*nu*((240*pow(f,beta)*T*alpha*(pow(f,beta)*alpha*beta - \
				f*nu*rho))/pow((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2,1.5) + \
				(360*f4*(f - K)*T*(-(pow(f,beta)*alpha*beta) + \
				f*nu*rho)*(pow(f,beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 + 3*rho2)) + \
				f*nu*rho*(12 + T*nu2*(-5 + 6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5) - (40*f3*(f - K)*(6*pow(f,beta)*T*alpha*(-1 + beta)*nu*rho + \
				f*(12 + T*nu2*(-5 + 18*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),1.5) - (5*pow(f,1 - beta)*pow(f - K,2)*T*(-(pow(f,beta)*alpha*beta) \
				+ f*nu*rho)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 54*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 - 3*rho2) - \
				36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + 1350*rho4 - \
				675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + beta)*nu*rho*(beta*(-84 + \
				T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + \
				T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + 4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 \
				+ 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + \
				T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 \
				- 9*rho2 + beta*(-4 + 15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + \
				beta*(184 - 1380*rho2 + 1305*rho4))) + 3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + \
				120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + \
				120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 - 1260*rho2 + 1305*rho4)) + beta2*(960 \
				+ 120*T*nu2*(-4 + 45*rho2) + T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),3.5)) - \
				(2*pow(f,1 - beta)*pow(f - K,2)*(9*pow(f,5*beta)*T2*alpha5*pow(-1 + \
				beta,4)*beta + 3*pow(f,1 + 4*beta)*T2*alpha4*pow(-1 + beta,2)*(-3 + 26*beta + \
				17*beta2)*nu*rho + 3*f5*nu*rho*(960 + 600*T*nu2*(-2 + 3*rho2) + T2*nu4*(118 - \
				300*rho2 + 225*rho4)) - 30*pow(f,3 + 2*beta)*T*alpha2*nu*rho*(84 + T*nu2*(-26 \
				+ 45*rho2) - 6*beta*(36 + T*nu2*(-14 + 29*rho2)) + 3*beta2*(60 + T*nu2*(-26 + \
				59*rho2))) - 10*pow(f,2 + 3*beta)*T*alpha3*(-1 + beta)*(24 + 2*T*nu2*(-5 + \
				18*rho2) + beta*(-84 + T*nu2*(29 - 108*rho2)) + beta2*(132 + T*nu2*(-31 + \
				126*rho2))) - pow(f,4 + beta)*alpha*(-1440 + 120*T*nu2*(7 - 27*rho2 + \
				beta*(-4 + 45*rho2)) + T2*nu4*(-10*(10 - 81*rho2 + 90*rho4) + beta*(184 - \
				4140*rho2 + 6525*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5))))/(40.*f2);
                return value;
            }

		if (ParamStr == "VolOfVol")
            {
			   value = ((240*pow(f,beta)*T*alpha*(3*pow(f,beta)*alpha*beta*rho + f*nu*(2 \
				- 3*rho2)))/pow((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2,1.5) - \
				(120*f3*(f - K)*(2*pow(f,beta)*T*alpha*(-1 + beta)*nu*(-2 + 3*rho2) + \
				3*f*rho*(4 + T*nu2*(-5 + 6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),1.5) + (360*f4*(f - K)*T*(-3*pow(f,beta)*alpha*beta*rho + f*nu*(-2 + \
				3*rho2))*(pow(f,beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 + 3*rho2)) + \
				f*nu*rho*(12 + T*nu2*(-5 + 6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5) - (6*pow(f,1 - beta)*pow(f - \
				K,2)*(9*pow(f,5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*rho + pow(f,1 + \
				4*beta)*T2*alpha4*pow(-1 + beta,2)*nu*(6 - 9*rho2 + 26*beta*(-2 + 3*rho2) + \
				beta2*(46 + 51*rho2)) + f5*nu*(960*(-2 + 3*rho2) + 24*T*nu2*(88 - 300*rho2 + \
				225*rho4) + T2*nu4*(-368 + 1062*rho2 - 1350*rho4 + 675*rho6)) - 30*pow(f,2 + \
				3*beta)*T*alpha3*(-1 + beta)*rho*(beta*(-28 + T*nu2*(29 - 36*rho2)) + 2*(4 + \
				T*nu2*(-5 + 6*rho2)) + beta2*(44 + T*nu2*(-31 + 42*rho2))) - 5*pow(f,4 + \
				beta)*alpha*rho*(-288 + 72*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) - \
				2*pow(f,3 + 2*beta)*T*alpha2*nu*(3*(60*(-4 + 7*rho2) + T*nu2*(56 - 260*rho2 + \
				225*rho4)) + beta*(960 - 3240*rho2 - 2*T*nu2*(128 - 1260*rho2 + 1305*rho4)) + \
				beta2*(60*(-4 + 45*rho2) + T*nu2*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5)) - \
				(5*pow(f,1 - beta)*pow(f - K,2)*T*(-3*pow(f,beta)*alpha*beta*rho + f*nu*(-2 + \
				3*rho2))*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 54*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 - 3*rho2) - \
				36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + 1350*rho4 - \
				675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + beta)*nu*rho*(beta*(-84 + \
				T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + \
				T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + 4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 \
				+ 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + \
				T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 \
				- 9*rho2 + beta*(-4 + 15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + \
				beta*(184 - 1380*rho2 + 1305*rho4))) + 3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + \
				120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + \
				120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 - 1260*rho2 + 1305*rho4)) + beta2*(960 \
				+ 120*T*nu2*(-4 + 45*rho2) + T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),3.5)))/(40.*sqrt(3.)*f2);
               return value;
        }
        else
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"SABR_DetermineDerivatives : ParamStr  bad input :");
        break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"SABR_DetermineDerivatives : flag  bad input :");
			break;
		}
	}
    return value;
}

void SABR_DetermineGenDerivatives( double f,
	double K,
	double T, 
	double alpha,
	double beta, 
	double rho, 
	double nu,
	int flag,
	double* implvol,
	double* der_alpha, 
	double* der_beta, 
	double* der_rho, 
	double* der_nu, 
	double* der_f)
{
	double z,x,theta,sig;
	double sigDerx,sigDertheta,sigDerf,thetaDerx,thetaDerz,thetaDeralpha;
	double thetaDerbeta,thetaDerrho,thetaDernu,thetaDerf,xDerz,xDerrho,
            xDernu,zDerf,zDeralpha,zDerbeta;
	double ZTotalSigDeralpha,ZTotalSigDerbeta,ZTotalSigDerrho,ZTotalSigDernu,ZTotalSigDerf;
	z = (pow(f,1 - beta) - pow(K,1 - beta))/(alpha*(1 - beta));
	x = log((z*nu - rho + sqrt(1 + pow(z,2)*pow(nu,2) - 2*z*nu*rho))/(1 - rho))/nu;
	theta = pow(2,-1 - beta)*pow(f + K,-1 + beta)*pow(z,2)*alpha*beta*nu*rho + log((pow(f*K,beta/2.)*z*alpha)/(f - K)) + log((x*pow(1 + pow(z,2)*pow(nu,2) - 2*z*nu*rho,0.25))/z);
	double sqrtTheta=1 - (2*T*(theta - log((sqrt(f*K)*log(f/K))/(f - K))))/pow(x,2);
	if(sqrtTheta<=0)
	{
		(*der_alpha)=0;
		(*der_beta)=0;
		(*der_rho)=0;
		(*der_nu)=0;
		(*der_f)=0;
		(*implvol)=0;
	}
	else
	{
		sig = log(f/K)/(x*sqrt(sqrtTheta));
		switch (flag) 
		{
		case ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACT :
		case ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL :
			{
				sigDerx = -(log(f/K)/(sqrt(1 - (2*T*(theta - log((sqrt(f*K)*log(f/K))/(f - \
					K))))/pow(x,2))*(pow(x,2) - 2*T*theta + 2*T*log((sqrt(f*K)*log(f/K))/(f - \
					K)))));
				sigDertheta = (T*log(f/K))/(pow(x,3)*pow(1 - (2*T*(theta - \
					log((sqrt(f*K)*log(f/K))/(f - K))))/pow(x,2),1.5));
				sigDerf = ((f + K)*T*log(f/K) - 2*(f - K)*(T - pow(x,2) + 2*T*theta - \
					2*T*log((sqrt(f*K)*log(f/K))/(f - K))))/(2.*f*(f - K)*x*sqrt(1 - (2*T*(theta \
					- log((sqrt(f*K)*log(f/K))/(f - K))))/pow(x,2))*(pow(x,2) - 2*T*theta + \
					2*T*log((sqrt(f*K)*log(f/K))/(f - K))));
				thetaDerx = 1/x;
				thetaDerz = (pow(2,-1 - beta)*nu*(pow(2,beta)*f*(z*nu - rho) + \
					pow(2,beta)*K*(z*nu - rho) + 2*pow(f + K,beta)*z*alpha*beta*rho*(1 + \
					pow(z,2)*pow(nu,2) - 2*z*nu*rho)))/((f + K)*(1 + pow(z,2)*pow(nu,2) - \
					2*z*nu*rho));
				thetaDeralpha = 1/alpha + pow(2,-1 - beta)*pow(f + K,-1 + \
					beta)*pow(z,2)*beta*nu*rho;
				thetaDerbeta = (pow(2,-1 - beta)*(pow(2,beta)*(f + K)*log(f*K) - pow(f + \
					K,beta)*pow(z,2)*alpha*nu*rho*(-1 + beta*log(2.) - beta*log(f + K))))/(f + K);
				thetaDerrho = (z*nu*((pow(f + K,-1 + beta)*z*alpha*beta)/pow(2,beta) + 1/(-1 \
					- pow(z,2)*pow(nu,2) + 2*z*nu*rho)))/2.;
				thetaDernu = pow(2,-1 - beta)*pow(f + K,-1 + beta)*pow(z,2)*alpha*beta*rho + \
					(z*(z*nu - rho))/(2 + 2*pow(z,2)*pow(nu,2) - 4*z*nu*rho);
				thetaDerf = (2*f - f*beta + K*beta)/(-2*pow(f,2) + 2*f*K) + pow(2,-1 - \
					beta)*pow(f + K,-2 + beta)*pow(z,2)*alpha*(-1 + beta)*beta*nu*rho;
				xDerz = 1/sqrt(1 + pow(z,2)*pow(nu,2) - 2*z*nu*rho);
				xDerrho = ((1 - rho)*(z*nu - rho + sqrt(1 + pow(z,2)*pow(nu,2) - 2*z*nu*rho) \
					+ (1 - rho)*(-1 - (z*nu)/sqrt(1 + pow(z,2)*pow(nu,2) - \
					2*z*nu*rho))))/(nu*pow(-1 + rho,2)*(z*nu - rho + sqrt(1 + pow(z,2)*pow(nu,2) \
					- 2*z*nu*rho)));
				xDernu = ((z*nu)/sqrt(1 + pow(z,2)*pow(nu,2) - 2*z*nu*rho) - log((z*nu - rho \
					+ sqrt(1 + pow(z,2)*pow(nu,2) - 2*z*nu*rho))/(1 - rho)))/pow(nu,2);
				zDerf = 1/(pow(f,beta)*alpha);
				zDeralpha = (pow(f,1 - beta) - pow(K,1 - beta))/(pow(alpha,2)*(-1 + beta));
				zDerbeta = (-(pow(f,beta)*K) + f*pow(K,beta) + f*pow(K,beta)*(-1 + \
					beta)*log(f) - pow(f,beta)*K*(-1 + \
					beta)*log(K))/(pow(f,beta)*pow(K,beta)*alpha*pow(-1 + beta,2));
				ZTotalSigDeralpha = xDerz*zDeralpha*(sigDerx + sigDertheta*thetaDerx) + \
					sigDertheta*(zDeralpha*thetaDerz + thetaDeralpha);
				ZTotalSigDerbeta = xDerz*zDerbeta*(sigDerx + sigDertheta*thetaDerx) + \
					sigDertheta*(zDerbeta*thetaDerz + thetaDerbeta);
				ZTotalSigDerrho = xDerrho*(sigDerx + sigDertheta*thetaDerx) + \
					sigDertheta*thetaDerrho;
				ZTotalSigDernu = xDernu*(sigDerx + sigDertheta*thetaDerx) + \
					sigDertheta*thetaDernu;
				ZTotalSigDerf = sigDerf + xDerz*zDerf*(sigDerx + sigDertheta*thetaDerx) + \
					sigDertheta*(thetaDerf + zDerf*thetaDerz);
				break;
				
			}
		case ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL2 :  
		case ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACTSTRIKE :  
			{
				sigDerx = -(log(f/K)/(sqrt(1 - (2*T*(theta - log((sqrt(f*K)*log(f/K))/(f - \
					K))))/pow(x,2))*(pow(x,2) - 2*T*theta + 2*T*log((sqrt(f*K)*log(f/K))/(f - \
					K)))));
				sigDertheta = (T*log(f/K))/(pow(x,3)*pow(1 - (2*T*(theta - \
					log((sqrt(f*K)*log(f/K))/(f - K))))/pow(x,2),1.5));
				sigDerf = ((f + K)*T*log(f/K) - 2*(f - K)*(T - pow(x,2) + 2*T*theta - \
					2*T*log((sqrt(f*K)*log(f/K))/(f - K))))/(2.*f*(f - K)*x*sqrt(1 - (2*T*(theta \
					- log((sqrt(f*K)*log(f/K))/(f - K))))/pow(x,2))*(pow(x,2) - 2*T*theta + \
					2*T*log((sqrt(f*K)*log(f/K))/(f - K))));
				thetaDerx = 1/x;
				thetaDerz = (nu*(K*(z*nu - rho) + pow(K,beta)*z*alpha*beta*rho*(1 + \
					pow(z,2)*pow(nu,2) - 2*z*nu*rho)))/(2.*K*(1 + pow(z,2)*pow(nu,2) - \
					2*z*nu*rho));
				thetaDeralpha = 1/alpha + (pow(K,-1 + beta)*pow(z,2)*beta*nu*rho)/4.;
				thetaDerbeta = (pow(K,beta)*pow(z,2)*alpha*nu*rho + \
					pow(K,beta)*pow(z,2)*alpha*beta*nu*rho*log(K) + 2*K*log(f*K))/(4.*K);
				thetaDerrho = (z*nu*(pow(K,-1 + beta)*z*alpha*beta - 2/(1 + \
					pow(z,2)*pow(nu,2) - 2*z*nu*rho)))/4.;
				thetaDernu = (z*(pow(K,-1 + beta)*z*alpha*beta*rho + (2*z*nu - 2*rho)/(1 + \
					pow(z,2)*pow(nu,2) - 2*z*nu*rho)))/4.;
				thetaDerf = (2*f - f*beta + K*beta)/(-2*pow(f,2) + 2*f*K);
				xDerz = 1/sqrt(1 + pow(z,2)*pow(nu,2) - 2*z*nu*rho);
				xDerrho = ((1 - rho)*(z*nu - rho + sqrt(1 + pow(z,2)*pow(nu,2) - 2*z*nu*rho) \
					+ (1 - rho)*(-1 - (z*nu)/sqrt(1 + pow(z,2)*pow(nu,2) - \
					2*z*nu*rho))))/(nu*pow(-1 + rho,2)*(z*nu - rho + sqrt(1 + pow(z,2)*pow(nu,2) \
					- 2*z*nu*rho)));
				xDernu = ((z*nu)/sqrt(1 + pow(z,2)*pow(nu,2) - 2*z*nu*rho) - log((z*nu - rho \
					+ sqrt(1 + pow(z,2)*pow(nu,2) - 2*z*nu*rho))/(1 - rho)))/pow(nu,2);
				zDerf = 1/(pow(f,beta)*alpha);
				zDeralpha = (pow(f,1 - beta) - pow(K,1 - beta))/(pow(alpha,2)*(-1 + beta));
				zDerbeta = (-(pow(f,beta)*K) + f*pow(K,beta) + f*pow(K,beta)*(-1 + \
					beta)*log(f) - pow(f,beta)*K*(-1 + \
					beta)*log(K))/(pow(f,beta)*pow(K,beta)*alpha*pow(-1 + beta,2));
				ZTotalSigDeralpha = xDerz*zDeralpha*(sigDerx + sigDertheta*thetaDerx) + \
					sigDertheta*(zDeralpha*thetaDerz + thetaDeralpha);
				ZTotalSigDerbeta = xDerz*zDerbeta*(sigDerx + sigDertheta*thetaDerx) + \
					sigDertheta*(zDerbeta*thetaDerz + thetaDerbeta);
				ZTotalSigDerrho = xDerrho*(sigDerx + sigDertheta*thetaDerx) + \
					sigDertheta*thetaDerrho;
				ZTotalSigDernu = xDernu*(sigDerx + sigDertheta*thetaDerx) + \
					sigDertheta*thetaDernu;
				ZTotalSigDerf = sigDerf + xDerz*zDerf*(sigDerx + sigDertheta*thetaDerx) + \
					sigDertheta*(thetaDerf + zDerf*thetaDerz);
				break;
			}
		default :
			{	
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"SABR_DetermineDerivatives : flag  bad input :");
				break;
			}
	}
	(*der_alpha)=ZTotalSigDeralpha;
	(*der_beta)=ZTotalSigDerbeta;
	(*der_rho)=ZTotalSigDerrho;
	(*der_nu)=ZTotalSigDernu;
	(*der_f)=ZTotalSigDerf;
	(*implvol)=sig;
	}
}

/// development in series of the derivatives of SABR option formula around the money

void SABR_AroundAtTheMoney_DetermineDerivatives( double f,double K,double T, double alpha, double beta, double rho, double nu,int flag,
			double* implvol,double* der_alpha, double* der_beta, double* der_rho, double* der_nu, double* der_f,
			bool der_alphaflag, bool der_betaflag, bool der_rhoflag, bool der_nuflag, bool der_fflag)
{
	double deralphaATM=0,derbetaATM=0,derrhoATM=0,dernuATM=0,derfATM=0,ATMopt;
	double alpha2=alpha*alpha;double alpha3=alpha2*alpha;double \
		alpha4=alpha3*alpha;double alpha5=alpha4*alpha;double alpha6=alpha5*alpha;
	double beta2=beta*beta;double beta3=beta2*beta;double beta4=beta3*beta;double \
		beta5=beta4*beta;double beta6=beta5*beta;
	double rho2=rho*rho;double rho3=rho2*rho;double rho4=rho3*rho;double \
		rho5=rho4*rho;double rho6=rho5*rho;
	double nu2=nu*nu;double nu3=nu2*nu;double nu4=nu3*nu;double nu5=nu4*nu;double \
		nu6=nu5*nu;
	double f2=f*f;double f3=f2*f;double f4=f3*f;double f5=f4*f;double f6=f5*f;
	double T2=T*T;double T3=T2*T;double T4=T3*T;double T5=T4*T;double T6=T5*T;
		switch (flag) 
	{
	case ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACT :
	case ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL :
		{
			if (der_alphaflag)
			deralphaATM = ((240*pow(f,beta)*(f - K)*(-1 + \
				beta)*(6*pow(f,beta)*T*alpha*beta*nu*rho + f*(-12 + T*nu2*(2 - \
				3*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),1.5) + \
				(480*pow(f,-3 + 2*beta)*T*alpha*(pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho))/pow((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 \
				+ beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2,1.5) + \
				(480*pow(f,-1 + beta))/sqrt((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2) + \
				(720*pow(f,beta)*(f - K)*T*(pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho)*(3*pow(f,2*beta)*T*alpha2*(-1 + beta)*beta*nu*rho + \
				f2*nu*rho*(-12 + T*nu2*(5 - 6*rho2)) - pow(f,1 + beta)*alpha*(-1 + beta)*(12 \
				+ T*nu2*(-2 + 3*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5) + \
				(12*pow(f - K,2)*(pow(f,5*beta)*T2*alpha5*pow(-1 + beta,6) - 15*pow(f,1 + \
				4*beta)*T2*alpha4*pow(-1 + beta,3)*beta*(7 + 8*beta)*nu*rho + 2*pow(f,2 + \
				3*beta)*T*alpha3*(-1 + beta)*(-36 + 3*T*nu2*(2 - 3*rho2) + beta2*(-588 + \
				T*nu2*(98 - 657*rho2)) + 29*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta3*(276 + \
				T*nu2*(-46 + 39*rho2))) + 15*pow(f,3 + 2*beta)*T*alpha2*(-1 + \
				beta)*nu*rho*(10*beta*(12 + T*nu2) + 4*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(84 \
				+ T*nu2*(-32 + 39*rho2))) + f5*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + \
				beta*(-4 + 15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - \
				1380*rho2 + 1305*rho4))) + pow(f,4 + beta)*alpha*(3*(1280 + 120*T*nu2*(-4 + \
				7*rho2) + T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + 240*T*nu2*(-4 + \
				9*rho2) + T2*nu4*(128 - 810*rho2 + 765*rho4)) + beta2*(960 + 120*T*nu2*(-4 + \
				27*rho2) + T2*nu4*(88 - 1440*rho2 + \
				1575*rho4)))))/(f2*alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5)) - \
				(10*pow(f - K,2)*T*(-(pow(f,beta)*alpha*pow(-1 + beta,2)) - \
				3*f*beta*nu*rho)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 18*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,3)*beta*(7 + 8*beta)*nu*rho + f6*nu2*(2880*(2 \
				- 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + \
				1350*rho4 - 675*rho6)) + 3*pow(f,2 + 4*beta)*T*alpha4*(-1 + beta)*(-36 + \
				3*T*nu2*(2 - 3*rho2) + beta2*(-588 + T*nu2*(98 - 657*rho2)) + 29*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta3*(276 + T*nu2*(-46 + 39*rho2))) + 30*pow(f,3 + \
				3*beta)*T*alpha3*(-1 + beta)*nu*rho*(10*beta*(12 + T*nu2) + 4*(12 + T*nu2*(-5 \
				+ 6*rho2)) + beta2*(84 + T*nu2*(-32 + 39*rho2))) + 6*pow(f,5 + \
				beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 240*T*nu2*(-4 + 9*rho2) + T2*nu4*(128 \
				- 810*rho2 + 765*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 27*rho2) + T2*nu4*(88 \
				- 1440*rho2 + 1575*rho4)))))/(f2*alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),3.5)) - (2*pow(f,-2 - beta)*pow(f - \
				K,2)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 18*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,3)*beta*(7 + 8*beta)*nu*rho + f6*nu2*(2880*(2 \
				- 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + \
				1350*rho4 - 675*rho6)) + 3*pow(f,2 + 4*beta)*T*alpha4*(-1 + beta)*(-36 + \
				3*T*nu2*(2 - 3*rho2) + beta2*(-588 + T*nu2*(98 - 657*rho2)) + 29*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta3*(276 + T*nu2*(-46 + 39*rho2))) + 30*pow(f,3 + \
				3*beta)*T*alpha3*(-1 + beta)*nu*rho*(10*beta*(12 + T*nu2) + 4*(12 + T*nu2*(-5 \
				+ 6*rho2)) + beta2*(84 + T*nu2*(-32 + 39*rho2))) + 6*pow(f,5 + \
				beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 240*T*nu2*(-4 + 9*rho2) + T2*nu4*(128 \
				- 810*rho2 + 765*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 27*rho2) + T2*nu4*(88 \
				- 1440*rho2 + 1575*rho4)))))/(alpha2*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5)))/(80.*sqrt(3.));
			if (der_betaflag)
			derbetaATM = ((480*pow(f,-1 + \
				beta)*alpha*log(f))/sqrt((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2) - \
				(2*pow(f,-2 - beta)*pow(f - K,2)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - \
				18*pow(f,1 + 5*beta)*T2*alpha5*pow(-1 + beta,3)*beta*(7 + 8*beta)*nu*rho + \
				f6*nu2*(2880*(2 - 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 \
				- 1062*rho2 + 1350*rho4 - 675*rho6)) + 3*pow(f,2 + 4*beta)*T*alpha4*(-1 + \
				beta)*(-36 + 3*T*nu2*(2 - 3*rho2) + beta2*(-588 + T*nu2*(98 - 657*rho2)) + \
				29*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta3*(276 + T*nu2*(-46 + 39*rho2))) + \
				30*pow(f,3 + 3*beta)*T*alpha3*(-1 + beta)*nu*rho*(10*beta*(12 + T*nu2) + \
				4*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(84 + T*nu2*(-32 + 39*rho2))) + \
				6*pow(f,5 + beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + \
				15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + \
				1305*rho4))) + 3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) \
				+ T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + 240*T*nu2*(-4 + 9*rho2) \
				+ T2*nu4*(128 - 810*rho2 + 765*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 27*rho2) \
				+ T2*nu4*(88 - 1440*rho2 + \
				1575*rho4))))*log(f))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5)) + \
				(480*pow(f,-3 + 2*beta)*T*alpha2*(pow(f,beta)*alpha*(-1 + beta) + 3*f*nu*rho \
				+ (pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho)*log(f)))/pow((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)))/f2,1.5) - (720*pow(f,beta)*(f - \
				K)*T*alpha*(-3*pow(f,2*beta)*T*alpha2*(-1 + beta)*beta*nu*rho + pow(f,1 + \
				beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 + 3*rho2)) + f2*nu*rho*(12 + \
				T*nu2*(-5 + 6*rho2)))*(pow(f,beta)*alpha*(-1 + beta) + 3*f*nu*rho + \
				(pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho)*log(f)))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5) + \
				(10*pow(f - K,2)*T*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 18*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,3)*beta*(7 + 8*beta)*nu*rho + f6*nu2*(2880*(2 \
				- 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + \
				1350*rho4 - 675*rho6)) + 3*pow(f,2 + 4*beta)*T*alpha4*(-1 + beta)*(-36 + \
				3*T*nu2*(2 - 3*rho2) + beta2*(-588 + T*nu2*(98 - 657*rho2)) + 29*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta3*(276 + T*nu2*(-46 + 39*rho2))) + 30*pow(f,3 + \
				3*beta)*T*alpha3*(-1 + beta)*nu*rho*(10*beta*(12 + T*nu2) + 4*(12 + T*nu2*(-5 \
				+ 6*rho2)) + beta2*(84 + T*nu2*(-32 + 39*rho2))) + 6*pow(f,5 + \
				beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 240*T*nu2*(-4 + 9*rho2) + T2*nu4*(128 \
				- 810*rho2 + 765*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 27*rho2) + T2*nu4*(88 \
				- 1440*rho2 + 1575*rho4))))*(pow(f,beta)*alpha*(-1 + beta) + 3*f*nu*rho + \
				(pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho)*log(f)))/(f2*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) \
				- 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),3.5)) \
				+ (240*pow(f,beta)*(f - K)*alpha*(3*pow(f,beta)*T*alpha*(-1 + 2*beta)*nu*rho \
				+ f*(-12 + T*nu2*(2 - 3*rho2)) + (-1 + \
				beta)*(6*pow(f,beta)*T*alpha*beta*nu*rho + f*(-12 + T*nu2*(2 - \
				3*rho2)))*log(f)))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 \
				+ beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),1.5) + (6*pow(f - \
				K,2)*(2*pow(f,5*beta)*T2*alpha5*pow(-1 + beta,5) - 48*pow(f,1 + \
				4*beta)*T2*alpha4*pow(-1 + beta,3)*beta*nu*rho - 6*pow(f,1 + \
				4*beta)*T2*alpha4*pow(-1 + beta,3)*(7 + 8*beta)*nu*rho - 18*pow(f,1 + \
				4*beta)*T2*alpha4*pow(-1 + beta,2)*beta*(7 + 8*beta)*nu*rho + \
				2*f5*T*nu3*rho*(120*(-4 + 15*rho2) + T*nu2*(184 - 1380*rho2 + 1305*rho4)) + \
				pow(f,2 + 3*beta)*T*alpha3*(-36 + 3*T*nu2*(2 - 3*rho2) + beta2*(-588 + \
				T*nu2*(98 - 657*rho2)) + 29*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta3*(276 + \
				T*nu2*(-46 + 39*rho2))) + 20*pow(f,3 + 2*beta)*T*alpha2*(-1 + \
				beta)*nu*rho*(60 + 5*T*nu2 + beta*(84 + T*nu2*(-32 + 39*rho2))) + 10*pow(f,3 \
				+ 2*beta)*T*alpha2*nu*rho*(10*beta*(12 + T*nu2) + 4*(12 + T*nu2*(-5 + \
				6*rho2)) + beta2*(84 + T*nu2*(-32 + 39*rho2))) + pow(f,2 + \
				3*beta)*T*alpha3*(-1 + beta)*(29*(12 + T*nu2*(-2 + 3*rho2)) + 3*beta2*(276 + \
				T*nu2*(-46 + 39*rho2)) - 2*beta*(588 + T*nu2*(-98 + 657*rho2))) + pow(f,4 + \
				beta)*alpha*(-2*(2400 + 240*T*nu2*(-4 + 9*rho2) + T2*nu4*(128 - 810*rho2 + \
				765*rho4)) + 2*beta*(960 + 120*T*nu2*(-4 + 27*rho2) + T2*nu4*(88 - 1440*rho2 \
				+ 1575*rho4))) + 2*pow(f,5*beta)*T2*alpha5*pow(-1 + beta,6)*log(f) - \
				30*pow(f,1 + 4*beta)*T2*alpha4*pow(-1 + beta,3)*beta*(7 + \
				8*beta)*nu*rho*log(f) + 4*pow(f,2 + 3*beta)*T*alpha3*(-1 + beta)*(-36 + \
				3*T*nu2*(2 - 3*rho2) + beta2*(-588 + T*nu2*(98 - 657*rho2)) + 29*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta3*(276 + T*nu2*(-46 + 39*rho2)))*log(f) + \
				30*pow(f,3 + 2*beta)*T*alpha2*(-1 + beta)*nu*rho*(10*beta*(12 + T*nu2) + \
				4*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(84 + T*nu2*(-32 + 39*rho2)))*log(f) + \
				2*f5*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + \
				1305*rho4)))*log(f) + 2*pow(f,4 + beta)*alpha*(3*(1280 + 120*T*nu2*(-4 + \
				7*rho2) + T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + 240*T*nu2*(-4 + \
				9*rho2) + T2*nu4*(128 - 810*rho2 + 765*rho4)) + beta2*(960 + 120*T*nu2*(-4 + \
				27*rho2) + T2*nu4*(88 - 1440*rho2 + \
				1575*rho4)))*log(f)))/(f2*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5)))/(80.*sqrt(3.));
			if(der_rhoflag)
			derrhoATM = (sqrt(3.)*nu*((240*pow(f,-2 + \
				beta)*T*alpha*(pow(f,beta)*alpha*beta - \
				f*nu*rho))/pow((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2,1.5) + \
				(360*f*(f - K)*T*(-(pow(f,beta)*alpha*beta) + \
				f*nu*rho)*(-3*pow(f,2*beta)*T*alpha2*(-1 + beta)*beta*nu*rho + pow(f,1 + \
				beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 + 3*rho2)) + f2*nu*rho*(12 + \
				T*nu2*(-5 + 6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5) - \
				(40*(f - K)*(-3*pow(f,2*beta)*T*alpha2*(-1 + beta)*beta + 6*pow(f,1 + \
				beta)*T*alpha*(-1 + beta)*nu*rho + f2*(12 + T*nu2*(-5 + \
				18*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),1.5) - (5*pow(f,-1 \
				- beta)*pow(f - K,2)*T*(-(pow(f,beta)*alpha*beta) + \
				f*nu*rho)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 18*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,3)*beta*(7 + 8*beta)*nu*rho + f6*nu2*(2880*(2 \
				- 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + \
				1350*rho4 - 675*rho6)) + 3*pow(f,2 + 4*beta)*T*alpha4*(-1 + beta)*(-36 + \
				3*T*nu2*(2 - 3*rho2) + beta2*(-588 + T*nu2*(98 - 657*rho2)) + 29*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta3*(276 + T*nu2*(-46 + 39*rho2))) + 30*pow(f,3 + \
				3*beta)*T*alpha3*(-1 + beta)*nu*rho*(10*beta*(12 + T*nu2) + 4*(12 + T*nu2*(-5 \
				+ 6*rho2)) + beta2*(84 + T*nu2*(-32 + 39*rho2))) + 6*pow(f,5 + \
				beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 240*T*nu2*(-4 + 9*rho2) + T2*nu4*(128 \
				- 810*rho2 + 765*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 27*rho2) + T2*nu4*(88 \
				- 1440*rho2 + 1575*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),3.5)) - (2*pow(f,-1 - beta)*pow(f - \
				K,2)*(3*pow(f,5*beta)*T2*alpha5*pow(-1 + beta,3)*beta*(7 + 8*beta) - \
				3*pow(f,1 + 4*beta)*T2*alpha4*(3 - 32*beta + 248*beta2 - 232*beta3 + \
				13*beta4)*nu*rho + 3*f5*nu*rho*(960 + 600*T*nu2*(-2 + 3*rho2) + T2*nu4*(118 - \
				300*rho2 + 225*rho4)) - 30*pow(f,3 + 2*beta)*T*alpha2*nu*rho*(84 + T*nu2*(-26 \
				+ 45*rho2) - 6*beta*(24 + T*nu2*(-9 + 17*rho2)) + 3*beta2*(36 + T*nu2*(-16 + \
				35*rho2))) - 5*pow(f,2 + 3*beta)*T*alpha3*(-1 + beta)*(48 + 10*beta*(12 + \
				T*nu2) + 4*T*nu2*(-5 + 18*rho2) + beta2*(84 + T*nu2*(-32 + 117*rho2))) - \
				pow(f,4 + beta)*alpha*(-1440 + 120*T*nu2*(7 - 27*rho2 + beta*(-4 + 45*rho2)) \
				+ T2*nu4*(-10*(10 - 81*rho2 + 90*rho4) + beta*(184 - 4140*rho2 + \
				6525*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5))))/40.;
			if(der_nuflag)
			dernuATM = ((480*pow(f,-2 + beta)*T*alpha*(3*pow(f,beta)*alpha*beta*rho + \
				f*nu*(2 - 3*rho2)))/pow((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)))/f2,1.5) - (240*(f - K)*(-3*pow(f,2*beta)*T*alpha2*(-1 + \
				beta)*beta*rho + 2*pow(f,1 + beta)*T*alpha*(-1 + beta)*nu*(-2 + 3*rho2) + \
				3*f2*rho*(4 + T*nu2*(-5 + 6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),1.5) + (720*f*(f - K)*T*(-3*pow(f,beta)*alpha*beta*rho + f*nu*(-2 + \
				3*rho2))*(-3*pow(f,2*beta)*T*alpha2*(-1 + beta)*beta*nu*rho + pow(f,1 + \
				beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 + 3*rho2)) + f2*nu*rho*(12 + \
				T*nu2*(-5 + 6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5) - \
				(12*pow(f,-1 - beta)*pow(f - K,2)*(3*pow(f,5*beta)*T2*alpha5*pow(-1 + \
				beta,3)*beta*(7 + 8*beta)*rho - pow(f,1 + 4*beta)*T2*alpha4*(-1 + beta)*nu*(6 \
				- 9*rho2 + beta2*(98 - 657*rho2) + 29*beta*(-2 + 3*rho2) + beta3*(-46 + \
				39*rho2)) + f5*nu*(960*(-2 + 3*rho2) + 24*T*nu2*(88 - 300*rho2 + 225*rho4) + \
				T2*nu4*(-368 + 1062*rho2 - 1350*rho4 + 675*rho6)) - 15*pow(f,2 + \
				3*beta)*T*alpha3*(-1 + beta)*rho*(10*beta*(4 + T*nu2) + 4*(4 + T*nu2*(-5 + \
				6*rho2)) + beta2*(28 + T*nu2*(-32 + 39*rho2))) - 5*pow(f,4 + \
				beta)*alpha*rho*(-288 + 72*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) - \
				2*pow(f,3 + 2*beta)*T*alpha2*nu*(3*(60*(-4 + 7*rho2) + T*nu2*(56 - 260*rho2 + \
				225*rho4)) - 2*beta*(120*(-4 + 9*rho2) + T*nu2*(128 - 810*rho2 + 765*rho4)) + \
				beta2*(60*(-4 + 27*rho2) + T*nu2*(88 - 1440*rho2 + \
				1575*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5)) - \
				(10*pow(f,-1 - beta)*pow(f - K,2)*T*(-3*pow(f,beta)*alpha*beta*rho + f*nu*(-2 \
				+ 3*rho2))*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 18*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,3)*beta*(7 + 8*beta)*nu*rho + f6*nu2*(2880*(2 \
				- 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + \
				1350*rho4 - 675*rho6)) + 3*pow(f,2 + 4*beta)*T*alpha4*(-1 + beta)*(-36 + \
				3*T*nu2*(2 - 3*rho2) + beta2*(-588 + T*nu2*(98 - 657*rho2)) + 29*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta3*(276 + T*nu2*(-46 + 39*rho2))) + 30*pow(f,3 + \
				3*beta)*T*alpha3*(-1 + beta)*nu*rho*(10*beta*(12 + T*nu2) + 4*(12 + T*nu2*(-5 \
				+ 6*rho2)) + beta2*(84 + T*nu2*(-32 + 39*rho2))) + 6*pow(f,5 + \
				beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 240*T*nu2*(-4 + 9*rho2) + T2*nu4*(128 \
				- 810*rho2 + 765*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 27*rho2) + T2*nu4*(88 \
				- 1440*rho2 + 1575*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),3.5)))/(80.*sqrt(3.));
			ATMopt = (pow(f,-1 + beta)*alpha)/sqrt(1 - (T*(pow(f,2*beta)*alpha2*pow(-1 + \
				beta,2) + 6*pow(f,1 + beta)*alpha*beta*nu*rho + f2*nu2*(2 - \
				3*rho2)))/(12.*f2)) + (sqrt(3.)*f*(-f + K)*(pow(f,beta)*alpha*(-1 + beta)*(12 \
				+ T*nu2*(-2 + 3*rho2)) + f*nu*rho*(12 + T*nu2*(-5 + \
				6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),1.5) + (pow(f,-2 - \
				beta)*pow(-f + K,2)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 54*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 - 3*rho2) - \
				36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + 1350*rho4 - \
				675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + beta)*nu*rho*(beta*(-84 + \
				T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + \
				T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + 4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 \
				+ 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + \
				T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 \
				- 9*rho2 + beta*(-4 + 15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + \
				beta*(184 - 1380*rho2 + 1305*rho4))) + 3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + \
				120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + \
				120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 - 1260*rho2 + 1305*rho4)) + beta2*(960 \
				+ 120*T*nu2*(-4 + 45*rho2) + T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(40.*sqrt(3.)*alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5));
			if (der_fflag)
			derfATM = ((480*pow(f,-4 + 2*beta)*T*alpha2*(-1 + \
				beta)*(pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho))/pow((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 \
				+ beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2,1.5) + \
				(480*pow(f,-2 + beta)*alpha*(-1 + \
				beta))/sqrt((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2) - \
				(240*f*(pow(f,beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 + 3*rho2)) + \
				f*nu*rho*(12 + T*nu2*(-5 + 6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),1.5) + (240*(-f + K)*(pow(f,beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 \
				+ 3*rho2)) + f*nu*rho*(12 + T*nu2*(-5 + \
				6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),1.5) + (720*(f - \
				K)*(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)*beta) - 3*pow(f,1 + \
				beta)*T*alpha*beta*(1 + beta)*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)))*(pow(f,beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 + 3*rho2)) + \
				f*nu*rho*(12 + T*nu2*(-5 + 6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5) - (240*(f - K)*(pow(f,beta)*alpha*(-1 + beta)*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + f*nu*rho*(12 + T*nu2*(-5 + \
				6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),1.5) - (4*pow(f,-2 \
				- beta)*(-f + K)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 54*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 - 3*rho2) - \
				36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + 1350*rho4 - \
				675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + beta)*nu*rho*(beta*(-84 + \
				T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + \
				T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + 4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 \
				+ 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + \
				T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 \
				- 9*rho2 + beta*(-4 + 15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + \
				beta*(184 - 1380*rho2 + 1305*rho4))) + 3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + \
				120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + \
				120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 - 1260*rho2 + 1305*rho4)) + beta2*(960 \
				+ 120*T*nu2*(-4 + 45*rho2) + T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5)) + \
				(2*pow(f,-3 - beta)*pow(f - K,2)*(-2 - beta)*(pow(f,6*beta)*T2*alpha6*pow(-1 \
				+ beta,6) - 54*pow(f,1 + 5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + \
				f6*nu2*(2880*(2 - 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 \
				- 1062*rho2 + 1350*rho4 - 675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + \
				beta)*nu*rho*(beta*(-84 + T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + \
				6*rho2)) + beta2*(132 + T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + \
				4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 + 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta2*(-276 + T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + \
				beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 \
				- 1260*rho2 + 1305*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 45*rho2) + \
				T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5)) - \
				(5*pow(f,-2 - beta)*pow(f - K,2)*(-2*pow(f,-1 + 2*beta)*T*alpha2*pow(-1 + \
				beta,2)*beta - 6*pow(f,beta)*T*alpha*beta*(1 + beta)*nu*rho + f*(24 + \
				2*T*nu2*(-2 + 3*rho2)))*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - \
				54*pow(f,1 + 5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 \
				- 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + \
				1350*rho4 - 675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + \
				beta)*nu*rho*(beta*(-84 + T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + \
				6*rho2)) + beta2*(132 + T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + \
				4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 + 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta2*(-276 + T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + \
				beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 \
				- 1260*rho2 + 1305*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 45*rho2) + \
				T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),3.5)) + \
				(2*pow(f,-2 - beta)*pow(f - K,2)*(6*pow(f,-1 + 6*beta)*T2*alpha6*pow(-1 + \
				beta,6)*beta - 54*pow(f,5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*(1 + \
				5*beta)*nu*rho + 6*f5*nu2*(2880*(2 - 3*rho2) - 36*T*nu2*(88 - 300*rho2 + \
				225*rho4) + T2*nu4*(368 - 1062*rho2 + 1350*rho4 - 675*rho6)) + 180*pow(f,2 + \
				3*beta)*T*alpha3*(-1 + beta)*(1 + beta)*nu*rho*(beta*(-84 + T*nu2*(29 - \
				36*rho2)) + 2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + T*nu2*(-31 + \
				42*rho2))) - 3*pow(f,1 + 4*beta)*T*alpha4*pow(-1 + beta,2)*(2 + 4*beta)*(-36 \
				+ 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + \
				T*nu2*(46 + 51*rho2))) + 6*pow(f,4 + beta)*alpha*(5 + beta)*nu*rho*(-1440 + \
				120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + \
				18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + 6*pow(f,3 + \
				2*beta)*alpha2*(2 + beta)*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 \
				- 1260*rho2 + 1305*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 45*rho2) + \
				T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5)))/(80.*sqrt(3.));
			
				break;
			
		}
	case ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL2 :  
	case ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACTSTRIKE :  
		{
			
			ATMopt = (pow(f,-1 + beta)*alpha)/sqrt(1 - (T*(pow(f,2*beta)*alpha2*pow(-1 + \
				beta,2) + 6*pow(f,1 + beta)*alpha*beta*nu*rho + f2*nu2*(2 - \
				3*rho2)))/(12.*f2)) + (sqrt(3.)*f*(-f + K)*(pow(f,beta)*alpha*(-1 + beta)*(12 \
				+ T*nu2*(-2 + 3*rho2)) + f*nu*rho*(12 + T*nu2*(-5 + \
				6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),1.5) + (pow(f,-2 - \
				beta)*pow(-f + K,2)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 54*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 - 3*rho2) - \
				36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + 1350*rho4 - \
				675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + beta)*nu*rho*(beta*(-84 + \
				T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + \
				T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + 4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 \
				+ 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + \
				T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 \
				- 9*rho2 + beta*(-4 + 15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + \
				beta*(184 - 1380*rho2 + 1305*rho4))) + 3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + \
				120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + \
				120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 - 1260*rho2 + 1305*rho4)) + beta2*(960 \
				+ 120*T*nu2*(-4 + 45*rho2) + T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(40.*sqrt(3.)*alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5));
			if (der_fflag)
			derfATM = ((480*pow(f,-4 + 2*beta)*T*alpha2*(-1 + \
				beta)*(pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho))/pow((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 \
				+ beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2,1.5) + \
				(480*pow(f,-2 + beta)*alpha*(-1 + \
				beta))/sqrt((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2) - \
				(240*f*(pow(f,beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 + 3*rho2)) + \
				f*nu*rho*(12 + T*nu2*(-5 + 6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),1.5) + (240*(-f + K)*(pow(f,beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 \
				+ 3*rho2)) + f*nu*rho*(12 + T*nu2*(-5 + \
				6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),1.5) + (720*(f - \
				K)*(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)*beta) - 3*pow(f,1 + \
				beta)*T*alpha*beta*(1 + beta)*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)))*(pow(f,beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 + 3*rho2)) + \
				f*nu*rho*(12 + T*nu2*(-5 + 6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5) - (240*(f - K)*(pow(f,beta)*alpha*(-1 + beta)*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + f*nu*rho*(12 + T*nu2*(-5 + \
				6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),1.5) - (4*pow(f,-2 \
				- beta)*(-f + K)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 54*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 - 3*rho2) - \
				36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + 1350*rho4 - \
				675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + beta)*nu*rho*(beta*(-84 + \
				T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + \
				T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + 4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 \
				+ 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + \
				T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 \
				- 9*rho2 + beta*(-4 + 15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + \
				beta*(184 - 1380*rho2 + 1305*rho4))) + 3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + \
				120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + \
				120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 - 1260*rho2 + 1305*rho4)) + beta2*(960 \
				+ 120*T*nu2*(-4 + 45*rho2) + T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5)) + \
				(2*pow(f,-3 - beta)*pow(f - K,2)*(-2 - beta)*(pow(f,6*beta)*T2*alpha6*pow(-1 \
				+ beta,6) - 54*pow(f,1 + 5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + \
				f6*nu2*(2880*(2 - 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 \
				- 1062*rho2 + 1350*rho4 - 675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + \
				beta)*nu*rho*(beta*(-84 + T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + \
				6*rho2)) + beta2*(132 + T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + \
				4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 + 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta2*(-276 + T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + \
				beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 \
				- 1260*rho2 + 1305*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 45*rho2) + \
				T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5)) - \
				(5*pow(f,-2 - beta)*pow(f - K,2)*(-2*pow(f,-1 + 2*beta)*T*alpha2*pow(-1 + \
				beta,2)*beta - 6*pow(f,beta)*T*alpha*beta*(1 + beta)*nu*rho + f*(24 + \
				2*T*nu2*(-2 + 3*rho2)))*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - \
				54*pow(f,1 + 5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 \
				- 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + \
				1350*rho4 - 675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + \
				beta)*nu*rho*(beta*(-84 + T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + \
				6*rho2)) + beta2*(132 + T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + \
				4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 + 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta2*(-276 + T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + \
				beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 \
				- 1260*rho2 + 1305*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 45*rho2) + \
				T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),3.5)) + \
				(2*pow(f,-2 - beta)*pow(f - K,2)*(6*pow(f,-1 + 6*beta)*T2*alpha6*pow(-1 + \
				beta,6)*beta - 54*pow(f,5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*(1 + \
				5*beta)*nu*rho + 6*f5*nu2*(2880*(2 - 3*rho2) - 36*T*nu2*(88 - 300*rho2 + \
				225*rho4) + T2*nu4*(368 - 1062*rho2 + 1350*rho4 - 675*rho6)) + 180*pow(f,2 + \
				3*beta)*T*alpha3*(-1 + beta)*(1 + beta)*nu*rho*(beta*(-84 + T*nu2*(29 - \
				36*rho2)) + 2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + T*nu2*(-31 + \
				42*rho2))) - 3*pow(f,1 + 4*beta)*T*alpha4*pow(-1 + beta,2)*(2 + 4*beta)*(-36 \
				+ 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + \
				T*nu2*(46 + 51*rho2))) + 6*pow(f,4 + beta)*alpha*(5 + beta)*nu*rho*(-1440 + \
				120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + \
				18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + 6*pow(f,3 + \
				2*beta)*alpha2*(2 + beta)*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 \
				- 1260*rho2 + 1305*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 45*rho2) + \
				T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5)))/(80.*sqrt(3.));
			if (der_alphaflag)
			deralphaATM = ((240*pow(f,4 + beta)*(-f + K)*(-1 + beta)*(12 + T*nu2*(-2 + \
				3*rho2)))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),1.5) + \
				(480*pow(f,2*beta)*T*alpha*(pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho))/pow((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 \
				+ beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2,1.5) + \
				(480*pow(f,2 + beta))/sqrt((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2) - \
				(720*pow(f,4 + beta)*(f - K)*T*(pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho)*(pow(f,beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 + 3*rho2)) + \
				f*nu*rho*(12 + T*nu2*(-5 + 6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5) + (12*f*pow(f - K,2)*(pow(f,5*beta)*T2*alpha5*pow(-1 + beta,6) \
				- 45*pow(f,1 + 4*beta)*T2*alpha4*pow(-1 + beta,4)*beta*nu*rho + 30*pow(f,3 + \
				2*beta)*T*alpha2*(-1 + beta)*nu*rho*(beta*(-84 + T*nu2*(29 - 36*rho2)) + \
				2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + T*nu2*(-31 + 42*rho2))) - \
				2*pow(f,2 + 3*beta)*T*alpha3*pow(-1 + beta,2)*(-36 + 3*T*nu2*(2 - 3*rho2) + \
				26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + T*nu2*(46 + 51*rho2))) + \
				f5*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				pow(f,4 + beta)*alpha*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 \
				- 1260*rho2 + 1305*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 45*rho2) + \
				T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5)) - \
				(10*f*pow(f - K,2)*T*(-(pow(f,beta)*alpha*pow(-1 + beta,2)) - \
				3*f*beta*nu*rho)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 54*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 - 3*rho2) - \
				36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + 1350*rho4 - \
				675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + beta)*nu*rho*(beta*(-84 + \
				T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + \
				T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + 4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 \
				+ 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + \
				T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 \
				- 9*rho2 + beta*(-4 + 15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + \
				beta*(184 - 1380*rho2 + 1305*rho4))) + 3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + \
				120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + \
				120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 - 1260*rho2 + 1305*rho4)) + beta2*(960 \
				+ 120*T*nu2*(-4 + 45*rho2) + T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),3.5)) - \
				(2*pow(f,1 - beta)*pow(f - K,2)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - \
				54*pow(f,1 + 5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 \
				- 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + \
				1350*rho4 - 675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + \
				beta)*nu*rho*(beta*(-84 + T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + \
				6*rho2)) + beta2*(132 + T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + \
				4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 + 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta2*(-276 + T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + \
				beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 \
				- 1260*rho2 + 1305*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 45*rho2) + \
				T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha2*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5)))/(80.*sqrt(3.)*f3);
			if (der_betaflag)
			derbetaATM = ((480*pow(f,2 + \
				beta)*alpha*log(f))/sqrt((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2) - \
				(2*pow(f,1 - beta)*pow(f - K,2)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - \
				54*pow(f,1 + 5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 \
				- 3*rho2) - 36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + \
				1350*rho4 - 675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + \
				beta)*nu*rho*(beta*(-84 + T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + \
				6*rho2)) + beta2*(132 + T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + \
				4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 + 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + \
				T*nu2*(-2 + 3*rho2)) + beta2*(-276 + T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + \
				beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) + \
				3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + 120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - \
				260*rho2 + 225*rho4)) - 2*beta*(2400 + 120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 \
				- 1260*rho2 + 1305*rho4)) + beta2*(960 + 120*T*nu2*(-4 + 45*rho2) + \
				T2*nu4*(88 - 2340*rho2 + \
				2655*rho4))))*log(f))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5)) - \
				(240*pow(f,4 + beta)*(f - K)*alpha*(12 + T*nu2*(-2 + 3*rho2))*(1 + (-1 + \
				beta)*log(f)))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),1.5) + \
				(480*pow(f,2*beta)*T*alpha2*(pow(f,beta)*alpha*(-1 + beta) + 3*f*nu*rho + \
				(pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho)*log(f)))/pow((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)))/f2,1.5) - (720*pow(f,4 + beta)*(f - \
				K)*T*alpha*(pow(f,beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 + 3*rho2)) + \
				f*nu*rho*(12 + T*nu2*(-5 + 6*rho2)))*(pow(f,beta)*alpha*(-1 + beta) + \
				3*f*nu*rho + (pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho)*log(f)))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5) + \
				(10*f*pow(f - K,2)*T*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 54*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 - 3*rho2) - \
				36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + 1350*rho4 - \
				675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + beta)*nu*rho*(beta*(-84 + \
				T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + \
				T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + 4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 \
				+ 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + \
				T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 \
				- 9*rho2 + beta*(-4 + 15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + \
				beta*(184 - 1380*rho2 + 1305*rho4))) + 3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + \
				120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + \
				120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 - 1260*rho2 + 1305*rho4)) + beta2*(960 \
				+ 120*T*nu2*(-4 + 45*rho2) + T2*nu4*(88 - 2340*rho2 + \
				2655*rho4))))*(pow(f,beta)*alpha*(-1 + beta) + 3*f*nu*rho + \
				(pow(f,beta)*alpha*pow(-1 + beta,2) + \
				3*f*beta*nu*rho)*log(f)))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),3.5) + \
				(6*f*pow(f - K,2)*(2*pow(f,5*beta)*T2*alpha5*pow(-1 + beta,5) - 18*pow(f,1 + \
				4*beta)*T2*alpha4*pow(-1 + beta,4)*nu*rho - 72*pow(f,1 + \
				4*beta)*T2*alpha4*pow(-1 + beta,3)*beta*nu*rho + 2*f5*T*nu3*rho*(120*(-4 + \
				15*rho2) + T*nu2*(184 - 1380*rho2 + 1305*rho4)) + 20*pow(f,3 + \
				2*beta)*T*alpha2*nu*rho*(beta*(-84 + T*nu2*(29 - 36*rho2)) + 2*(12 + \
				T*nu2*(-5 + 6*rho2)) + beta2*(132 + T*nu2*(-31 + 42*rho2))) + 20*pow(f,3 + \
				2*beta)*T*alpha2*(-1 + beta)*nu*rho*(-84 + T*nu2*(29 - 36*rho2) + beta*(264 + \
				2*T*nu2*(-31 + 42*rho2))) - pow(f,2 + 3*beta)*T*alpha3*pow(-1 + \
				beta,2)*(26*(12 + T*nu2*(-2 + 3*rho2)) + 2*beta*(-276 + T*nu2*(46 + \
				51*rho2))) - 2*pow(f,2 + 3*beta)*T*alpha3*(-1 + beta)*(-36 + 3*T*nu2*(2 - \
				3*rho2) + 26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + T*nu2*(46 + \
				51*rho2))) + pow(f,4 + beta)*alpha*(-2*(2400 + 120*T*nu2*(-8 + 27*rho2) + \
				T2*nu4*(128 - 1260*rho2 + 1305*rho4)) + 2*beta*(960 + 120*T*nu2*(-4 + \
				45*rho2) + T2*nu4*(88 - 2340*rho2 + 2655*rho4))) + \
				2*pow(f,5*beta)*T2*alpha5*pow(-1 + beta,6)*log(f) - 90*pow(f,1 + \
				4*beta)*T2*alpha4*pow(-1 + beta,4)*beta*nu*rho*log(f) + 60*pow(f,3 + \
				2*beta)*T*alpha2*(-1 + beta)*nu*rho*(beta*(-84 + T*nu2*(29 - 36*rho2)) + \
				2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + T*nu2*(-31 + 42*rho2)))*log(f) - \
				4*pow(f,2 + 3*beta)*T*alpha3*pow(-1 + beta,2)*(-36 + 3*T*nu2*(2 - 3*rho2) + \
				26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + T*nu2*(46 + \
				51*rho2)))*log(f) + 2*f5*nu*rho*(-1440 + 120*T*nu2*(7 - 9*rho2 + beta*(-4 + \
				15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + \
				1305*rho4)))*log(f) + 2*pow(f,4 + beta)*alpha*(3*(1280 + 120*T*nu2*(-4 + \
				7*rho2) + T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + 120*T*nu2*(-8 + \
				27*rho2) + T2*nu4*(128 - 1260*rho2 + 1305*rho4)) + beta2*(960 + 120*T*nu2*(-4 \
				+ 45*rho2) + T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))*log(f)))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5))/(80.*sqrt(3.)*f3);
			if(der_rhoflag)
			derrhoATM = (sqrt(3.)*nu*((240*pow(f,beta)*T*alpha*(pow(f,beta)*alpha*beta - \
				f*nu*rho))/pow((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2,1.5) + \
				(360*f4*(f - K)*T*(-(pow(f,beta)*alpha*beta) + \
				f*nu*rho)*(pow(f,beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 + 3*rho2)) + \
				f*nu*rho*(12 + T*nu2*(-5 + 6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5) - (40*f3*(f - K)*(6*pow(f,beta)*T*alpha*(-1 + beta)*nu*rho + \
				f*(12 + T*nu2*(-5 + 18*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),1.5) - (5*pow(f,1 - beta)*pow(f - K,2)*T*(-(pow(f,beta)*alpha*beta) \
				+ f*nu*rho)*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 54*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 - 3*rho2) - \
				36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + 1350*rho4 - \
				675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + beta)*nu*rho*(beta*(-84 + \
				T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + \
				T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + 4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 \
				+ 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + \
				T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 \
				- 9*rho2 + beta*(-4 + 15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + \
				beta*(184 - 1380*rho2 + 1305*rho4))) + 3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + \
				120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + \
				120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 - 1260*rho2 + 1305*rho4)) + beta2*(960 \
				+ 120*T*nu2*(-4 + 45*rho2) + T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),3.5)) - \
				(2*pow(f,1 - beta)*pow(f - K,2)*(9*pow(f,5*beta)*T2*alpha5*pow(-1 + \
				beta,4)*beta + 3*pow(f,1 + 4*beta)*T2*alpha4*pow(-1 + beta,2)*(-3 + 26*beta + \
				17*beta2)*nu*rho + 3*f5*nu*rho*(960 + 600*T*nu2*(-2 + 3*rho2) + T2*nu4*(118 - \
				300*rho2 + 225*rho4)) - 30*pow(f,3 + 2*beta)*T*alpha2*nu*rho*(84 + T*nu2*(-26 \
				+ 45*rho2) - 6*beta*(36 + T*nu2*(-14 + 29*rho2)) + 3*beta2*(60 + T*nu2*(-26 + \
				59*rho2))) - 10*pow(f,2 + 3*beta)*T*alpha3*(-1 + beta)*(24 + 2*T*nu2*(-5 + \
				18*rho2) + beta*(-84 + T*nu2*(29 - 108*rho2)) + beta2*(132 + T*nu2*(-31 + \
				126*rho2))) - pow(f,4 + beta)*alpha*(-1440 + 120*T*nu2*(7 - 27*rho2 + \
				beta*(-4 + 45*rho2)) + T2*nu4*(-10*(10 - 81*rho2 + 90*rho4) + beta*(184 - \
				4140*rho2 + 6525*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5))))/(40.*f2);
			if (der_nuflag)
			dernuATM = ((240*pow(f,beta)*T*alpha*(3*pow(f,beta)*alpha*beta*rho + f*nu*(2 \
				- 3*rho2)))/pow((-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - 6*pow(f,1 + \
				beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)))/f2,1.5) - \
				(120*f3*(f - K)*(2*pow(f,beta)*T*alpha*(-1 + beta)*nu*(-2 + 3*rho2) + \
				3*f*rho*(4 + T*nu2*(-5 + 6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),1.5) + (360*f4*(f - K)*T*(-3*pow(f,beta)*alpha*beta*rho + f*nu*(-2 + \
				3*rho2))*(pow(f,beta)*alpha*(-1 + beta)*(12 + T*nu2*(-2 + 3*rho2)) + \
				f*nu*rho*(12 + T*nu2*(-5 + 6*rho2))))/pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + \
				beta,2)) - 6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),2.5) - (6*pow(f,1 - beta)*pow(f - \
				K,2)*(9*pow(f,5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*rho + pow(f,1 + \
				4*beta)*T2*alpha4*pow(-1 + beta,2)*nu*(6 - 9*rho2 + 26*beta*(-2 + 3*rho2) + \
				beta2*(46 + 51*rho2)) + f5*nu*(960*(-2 + 3*rho2) + 24*T*nu2*(88 - 300*rho2 + \
				225*rho4) + T2*nu4*(-368 + 1062*rho2 - 1350*rho4 + 675*rho6)) - 30*pow(f,2 + \
				3*beta)*T*alpha3*(-1 + beta)*rho*(beta*(-28 + T*nu2*(29 - 36*rho2)) + 2*(4 + \
				T*nu2*(-5 + 6*rho2)) + beta2*(44 + T*nu2*(-31 + 42*rho2))) - 5*pow(f,4 + \
				beta)*alpha*rho*(-288 + 72*T*nu2*(7 - 9*rho2 + beta*(-4 + 15*rho2)) + \
				T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + beta*(184 - 1380*rho2 + 1305*rho4))) - \
				2*pow(f,3 + 2*beta)*T*alpha2*nu*(3*(60*(-4 + 7*rho2) + T*nu2*(56 - 260*rho2 + \
				225*rho4)) + beta*(960 - 3240*rho2 - 2*T*nu2*(128 - 1260*rho2 + 1305*rho4)) + \
				beta2*(60*(-4 + 45*rho2) + T*nu2*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + 3*rho2)),2.5)) - \
				(5*pow(f,1 - beta)*pow(f - K,2)*T*(-3*pow(f,beta)*alpha*beta*rho + f*nu*(-2 + \
				3*rho2))*(pow(f,6*beta)*T2*alpha6*pow(-1 + beta,6) - 54*pow(f,1 + \
				5*beta)*T2*alpha5*pow(-1 + beta,4)*beta*nu*rho + f6*nu2*(2880*(2 - 3*rho2) - \
				36*T*nu2*(88 - 300*rho2 + 225*rho4) + T2*nu4*(368 - 1062*rho2 + 1350*rho4 - \
				675*rho6)) + 60*pow(f,3 + 3*beta)*T*alpha3*(-1 + beta)*nu*rho*(beta*(-84 + \
				T*nu2*(29 - 36*rho2)) + 2*(12 + T*nu2*(-5 + 6*rho2)) + beta2*(132 + \
				T*nu2*(-31 + 42*rho2))) - 3*pow(f,2 + 4*beta)*T*alpha4*pow(-1 + beta,2)*(-36 \
				+ 3*T*nu2*(2 - 3*rho2) + 26*beta*(12 + T*nu2*(-2 + 3*rho2)) + beta2*(-276 + \
				T*nu2*(46 + 51*rho2))) + 6*pow(f,5 + beta)*alpha*nu*rho*(-1440 + 120*T*nu2*(7 \
				- 9*rho2 + beta*(-4 + 15*rho2)) + T2*nu4*(-10*(10 - 27*rho2 + 18*rho4) + \
				beta*(184 - 1380*rho2 + 1305*rho4))) + 3*pow(f,4 + 2*beta)*alpha2*(3*(1280 + \
				120*T*nu2*(-4 + 7*rho2) + T2*nu4*(56 - 260*rho2 + 225*rho4)) - 2*beta*(2400 + \
				120*T*nu2*(-8 + 27*rho2) + T2*nu4*(128 - 1260*rho2 + 1305*rho4)) + beta2*(960 \
				+ 120*T*nu2*(-4 + 45*rho2) + T2*nu4*(88 - 2340*rho2 + \
				2655*rho4)))))/(alpha*pow(-(pow(f,2*beta)*T*alpha2*pow(-1 + beta,2)) - \
				6*pow(f,1 + beta)*T*alpha*beta*nu*rho + f2*(12 + T*nu2*(-2 + \
				3*rho2)),3.5)))/(40.*sqrt(3.)*f2);

						break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"SABR_DetermineDerivatives : flag  bad input :");
			break;
		}
	}
	
	(*der_alpha)=deralphaATM;
	(*der_beta)=derbetaATM;
	(*der_rho)=derrhoATM;
	(*der_nu)=dernuATM;
	(*der_f)=derfATM;
	(*implvol)=ATMopt;

}

/// Set of Derivatives for the SABR Model
void SABR_DetermineAllDerivatives( double f,
     double K,
     double T, 
     double alpha, 
     double beta, 
     double rho, 
     double nu,
     int flag,
     double* implvol,
     double* der_alpha, 
     double* der_beta,
     double* der_rho, 
     double* der_nu, 
     double* der_f,
     bool der_alphaflag, 
     bool der_betaflag, 
     bool der_rhoflag, 
     bool der_nuflag, 
     bool der_fflag)
{
	if (fabs(K-f)<ARM_GP_CF_SABR_SERIE_USAGE_LIMIT*f*alpha) 
	{
		SABR_AroundAtTheMoney_DetermineDerivatives(f,K,T,alpha,beta,rho,nu,flag,implvol,der_alpha,der_beta,der_rho,der_nu,der_f,
			der_alphaflag,der_betaflag,der_rhoflag,der_nuflag,der_fflag);
	}
	else
	{
		SABR_DetermineGenDerivatives(f,K,T,alpha,beta,rho,nu,flag,implvol,der_alpha,der_beta,der_rho,der_nu,der_f);
	}
}


CC_END_NAMESPACE()


#undef ARM_CF_EPS
#undef ARM_CF_MAXIT

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
