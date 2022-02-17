/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file SABR_Analytics.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpclosedforms/smile_sabr.h"
#include "gpclosedforms/sabrvanilla.h"
#include "gpclosedforms/sabr_calibration.h"
#include "gpclosedforms/smile_calibration.h"
#include "gpclosedforms/sabrbdiff1.h"
#include "gpclosedforms/sabr_newversion.h"



# define ARM_GP_CF_SABR_SERIE_USAGE_LIMIT 0.01


CC_BEGIN_NAMESPACE(ARM)



///////////////////////////////////////////////////////////////////////
///  
///			Begining Exportable  Pricing Functions 
///
///////////////////////////////////////////////////////////////////////


/// callput =  1  for call
/// callput = -1  for put

double Export_SABR_ImplicitVol(double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int flag,
							   int nbsteps,
							   double C_Alpha_Exp,
							   double C_Alpha_Tanh,
							   double C_Kb_Tanh)
{
	ArgumentList a(f,K,tex,alpha,beta,rho,nu,flag,nbsteps,C_Alpha_Exp,C_Alpha_Tanh,C_Kb_Tanh);

	Power_Expression<ARM_CF_SABR_ImplicitVol_Formula> y;
	return y(a);
}

double Export_SABR_ImplicitVol(double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int flag,
							   int nbsteps)
{
	ArgumentList a(f,K,tex,alpha,beta,rho,nu,flag,nbsteps,0.01,0.03,1.5);

	Power_Expression<ARM_CF_SABR_ImplicitVol_Formula> y;
	return y(a);
}
///////////////////////////////////////////////////////////////////////
///  
///			1st Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_SABR_ImplicitVol(int i,
							   double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int flag,
							   int nbsteps,
							   double C_Alpha_Exp,
							   double C_Alpha_Tanh,
							   double C_Kb_Tanh)
{
	ArgumentList a(f,K,tex,alpha,beta,rho,nu,flag,nbsteps,C_Alpha_Exp,C_Alpha_Tanh,C_Kb_Tanh);
	
	Power_Expression<ARM_CF_SABR_ImplicitVol_Formula> y;
	return y(i,a);
}

double Export_SABR_ImplicitVol(int i,
							   double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int flag,
							   int nbsteps)
{
	ArgumentList a(f,K,tex,alpha,beta,rho,nu,flag,nbsteps,0.01,0.03,1.5);
	
	Power_Expression<ARM_CF_SABR_ImplicitVol_Formula> y;
	return y(i,a);
}

///////////////////////////////////////////////////////////////////////
///  
///			2nd Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_SABR_ImplicitVol(int i,int j,
							   double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int flag,
							   int nbsteps,
							   double C_Alpha_Exp,
							   double C_Alpha_Tanh,
							   double C_Kb_Tanh)
{
	ArgumentList a(f,K,tex,alpha,beta,rho,nu,flag,nbsteps,C_Alpha_Exp,C_Alpha_Tanh,C_Kb_Tanh);
	
	Power_Expression<ARM_CF_SABR_ImplicitVol_Formula> y;
	return y(i,j,a);
}

double Export_SABR_ImplicitVol(int i,int j,
							   double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int flag,
							   int nbsteps)
{
	ArgumentList a(f,K,tex,alpha,beta,rho,nu,flag,nbsteps,0.01,0.03,1.5);
	
	Power_Expression<ARM_CF_SABR_ImplicitVol_Formula> y;
	return y(i,j,a);
}
////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
///  
///			Begining Exportable  Pricing Functions 
///
///////////////////////////////////////////////////////////////////////


/// callput =  1  for call
/// callput = -1  for put

double Export_SABR_VanillaOption(double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int callputflag,
							   int flag,
							   int nbsteps,
							   double C_Alpha_Exp,
							   double C_Alpha_Tanh,
							   double C_Kb_Tanh)
{
	ArgumentList a(f,K,tex,alpha,beta,rho,nu,callputflag,flag,nbsteps,C_Alpha_Exp,C_Alpha_Tanh,C_Kb_Tanh);

	Power_Expression<ARM_CF_SABR_VanillaOption_Formula> y;
	return y(a);
}

double Export_SABR_VanillaOption(double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int callputflag,
							   int flag,
							   int nbsteps)
{
	ArgumentList a(f,K,tex,alpha,beta,rho,nu,callputflag,flag,nbsteps,0.01,0.03,1.5);

	Power_Expression<ARM_CF_SABR_VanillaOption_Formula> y;
	return y(a);
}


double Export_SABR_FromGaussianToDistribution(double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int flag,
							   int nbsteps,double st1,double st2,double st3)
{
	ArgumentList Underlying(f,alpha,beta,rho,nu,flag,nbsteps,st1,st2,st3);
	return SABR_smile::gaussian_to_distribution(Underlying,  K,  tex);
}

double Export_SABR_FromDistributionToGaussian(double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int flag,
							   int nbsteps,double st1,double st2,double st3)
{
	ArgumentList Underlying(f,alpha,beta,rho,nu,flag,nbsteps,st1,st2,st3);
	return SABR_smile::distribution_to_gaussian(Underlying,  K,  tex);
}


///////////////////////////////////////////////////////////////////////
///  
///			1st Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_SABR_VanillaOption(int i,
							   double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int callputflag,
							   int flag,
							   int nbsteps,
							   double C_Alpha_Exp,
							   double C_Alpha_Tanh,
							   double C_Kb_Tanh)
{
	ArgumentList a(f,K,tex,alpha,beta,rho,nu,callputflag,flag,nbsteps,C_Alpha_Exp,C_Alpha_Tanh,C_Kb_Tanh);
	
	Power_Expression<ARM_CF_SABR_VanillaOption_Formula> y;
	return y(i,a);
}

double Export_SABR_VanillaOption(int i,
							   double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int callputflag,
							   int flag,
							   int nbsteps)
{
	ArgumentList a(f,K,tex,alpha,beta,rho,nu,callputflag,flag,nbsteps,0.01,0.03,1.5);
	
	Power_Expression<ARM_CF_SABR_VanillaOption_Formula> y;
	return y(i,a);
}

///////////////////////////////////////////////////////////////////////
///  
///			2nd Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_SABR_VanillaOption(int i,int j,
							   double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int callputflag,
							   int flag,
							   int nbsteps,
							   double C_Alpha_Exp,
							   double C_Alpha_Tanh,
							   double C_Kb_Tanh)
{
	ArgumentList a(f,K,tex,alpha,beta,rho,nu,callputflag,flag,nbsteps,C_Alpha_Exp,C_Alpha_Tanh,C_Kb_Tanh);
	
	Power_Expression<ARM_CF_SABR_VanillaOption_Formula> y;
	return y(i,j,a);
}

double Export_SABR_VanillaOption(int i,int j,
							   double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int callputflag,
							   int flag,
							   int nbsteps)
{
	ArgumentList a(f,K,tex,alpha,beta,rho,nu,callputflag,flag,nbsteps,0.01,0.03,1.5);
	
	Power_Expression<ARM_CF_SABR_VanillaOption_Formula> y;
	return y(i,j,a);
}

///////////////////////////////////////////////////////////////////////
///  
///			SABR Calibration  
///
///////////////////////////////////////////////////////////////////////

///  beta, alpha,rho, and nu are searched

SABR_ParameterSet*  Export_SABR_Model_Calibrate(ARM_GP_Vector* K_Vec,ARM_GP_Vector* ImpVol_Vec,ARM_GP_Vector* Weigth_Vec,double f,double tex,int flag,
												int nbsteps,int algorithm,
												double alpha0,double beta0,double rho0,double nu0)
{
	int arg1=K_Vec->size();
	int arg2=ImpVol_Vec->size();
	if (arg1!=arg2)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_SABR_Model_Calibrate::value  : nb element of K_Vec and ImpVol_Vec should be the same !");
	}
	int arg3=Weigth_Vec->size();
	if (arg1!=arg3)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_SABR_Model_Calibrate::value  : nb element of K_Vec and Weigth_Vec should be the same !");
	}
	SABR_ParameterSet* setptr= SABR_CalibrateToSmile(K_Vec,ImpVol_Vec,Weigth_Vec,f,tex,flag,nbsteps,algorithm,
		alpha0, beta0,rho0, nu0);
	return setptr;
	
}

/// the most general search : f, beta, alpha,rho, and nu are searched

SABR_ParameterSet*  Export_SABR_Model_Calibrate(ARM_GP_Vector* K_Vec,ARM_GP_Vector* ImpVol_Vec,ARM_GP_Vector* Weigth_Vec,double tex,int flag,
												int nbsteps,int algorithm,
												double alpha0,double beta0,double rho0,double nu0,double f0)
{
	int arg1=K_Vec->size();
	int arg2=ImpVol_Vec->size();
	if (arg1!=arg2)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_SABR_Model_Calibrate::value  : nb element of K_Vec and ImpVol_Vec should be the same !");
	}
	int arg3=Weigth_Vec->size();
	if (arg1!=arg3)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_SABR_Model_Calibrate::value  : nb element of K_Vec and Weigth_Vec should be the same !");
	}
	SABR_ParameterSet* setptr= SABR_CalibrateToSmile(K_Vec,ImpVol_Vec,Weigth_Vec,tex,flag,nbsteps,algorithm,
		alpha0, beta0,rho0, nu0,f0);
	return setptr;
	
}

///  f and beta fixed
SABR_ParameterSet*  Export_SABR_Model_Calibrate(ARM_GP_Vector* K_Vec,ARM_GP_Vector* ImpVol_Vec,ARM_GP_Vector* Weigth_Vec,double f,double beta,double tex,
												int flag,double alpha_exp,double alpha_tanh,double kb_tanh,int nbsteps,int algorithm,
												double alpha0,double rho0,double nu0)
{
	int arg1=K_Vec->size();
	int arg2=ImpVol_Vec->size();
	if (arg1!=arg2)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_SABR_Model_Calibrate::value  : nb element of K_Vec and ImpVol_Vec should be the same !");
	}
	int arg3=Weigth_Vec->size();
	if (arg1!=arg3)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_SABR_Model_Calibrate::value  : nb element of K_Vec and Weigth_Vec should be the same !");
	}
	SABR_ParameterSet* setptr= SABR_CalibrateToSmile_withstrikecuter(K_Vec,ImpVol_Vec,Weigth_Vec,f,beta,
		tex,flag,alpha_exp,alpha_tanh,kb_tanh,nbsteps,algorithm,
		alpha0, rho0, nu0);

	return setptr;
	
}
///  f and beta fixed
///  with weigth on alpha, rho and nu 
SABR_ParameterSet*  Export_SABR_Model_Calibrate(ARM_GP_Vector* K_Vec,ARM_GP_Vector* ImpVol_Vec,ARM_GP_Vector* Weigth_Vec,double f,double beta,double tex,int flag,
												int nbsteps,int algorithm,
												double alpha0,double rho0,double nu0,double alphap,double rhop, double nup,
												double rweight_alpha,double rweight_rho,double rweight_nu)
{
	int arg1=K_Vec->size();
	int arg2=ImpVol_Vec->size();
	if (arg1!=arg2)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_SABR_Model_Calibrate::value  : nb element of K_Vec and ImpVol_Vec should be the same !");
	}
	int arg3=Weigth_Vec->size();
	if (arg1!=arg3)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_SABR_Model_Calibrate::value  : nb element of K_Vec and Weigth_Vec should be the same !");
	}
	SABR_ParameterSet* setptr= SABR_CalibrateToSmile(K_Vec,ImpVol_Vec,Weigth_Vec,f,beta,tex,flag,nbsteps,algorithm,
		alpha0, rho0, nu0,alphap,rhop,nup,rweight_alpha,rweight_rho,rweight_nu);
	return setptr;
	
}


///////////////////////////////////////////////////////////////////////
///  
///			conversion from sigma(K) to alpha  
///
///////////////////////////////////////////////////////////////////////
double Export_SABR_From_Sigma_To_Alpha (
							   double f, 
							   double K, 
							   double tex,
							   double sigma, 
							   double beta, 
							   double rho, 
							   double nu,
							   int flag,
							   int nbsteps)
{
	switch(flag)
	{
		///		Here the limitation of the dependence in the smile and the copula is that they 
		///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
		///
		
	case ARM_CF_SABR_VanillaOption_Formula::DIRECTEXACT :
		return SABR_Direct_FromATMsigma_to_alpha( f, K, tex,  sigma,  beta,  rho,  nu, flag);
	case ARM_CF_SABR_VanillaOption_Formula::DIRECTEXACTSTRIKE :
		return SABR_Direct_FromATMsigma_to_alpha( f, K, tex,  sigma,  beta,  rho,  nu, flag);
	case ARM_CF_SABR_VanillaOption_Formula::DIRECTGEOMETRIC :
		return SABR_Direct_FromATMsigma_to_alpha( f, K, tex,  sigma,  beta,  rho,  nu, flag);
	case ARM_CF_SABR_VanillaOption_Formula::DIRECTARITHMETIC :
		return SABR_Direct_FromATMsigma_to_alpha( f, K, tex,  sigma,  beta,  rho,  nu, flag);
	case ARM_CF_SABR_VanillaOption_Formula::NORMALEXACT :
		return SABR_Normal_FromATMsigma_to_alpha( f, K, tex,  sigma,  beta,  rho,  nu, flag);
	case ARM_CF_SABR_VanillaOption_Formula::NORMALGEOMETRIC :
		return SABR_Normal_FromATMsigma_to_alpha( f, K, tex,  sigma,  beta,  rho,  nu, flag);
	case ARM_CF_SABR_VanillaOption_Formula::NORMALARITHMETIC :
		return SABR_Normal_FromATMsigma_to_alpha( f, K, tex,  sigma,  beta,  rho,  nu, flag);
	case ARM_CF_SABR_VanillaOption_Formula::ANALYTICZP0 :  
		return SABR_NumericalIntegration_FromATMsigma_to_alpha( f, K, tex,  sigma,  beta,  rho,  nu, nbsteps);
	case ARM_CF_SABR_VanillaOption_Formula::ANALYTICZP2 :  
		return SABR_NumericalIntegration_FromATMsigma_to_alpha( f, K, tex,  sigma,  beta,  rho,  nu, nbsteps);
	case ARM_CF_SABR_VanillaOption_Formula::SABR_IMPLNVOL :  
		return SABR_Direct_FromATMsigma_to_alpha( f, K, tex,  sigma,  beta,  rho,  nu, flag);
	case ARM_CF_SABR_VanillaOption_Formula::SABR_IMPLNVOL2 :  
		return SABR_Direct_FromATMsigma_to_alpha( f, K, tex,  sigma,  beta,  rho,  nu, flag);
	case ARM_CF_SABR_VanillaOption_Formula::SABR_A :  
		return SABR_BetaEqOne_CompatibleKernel_FromATMsigma_to_alpha( f, K, tex,  sigma,  rho,  nu);
	case ARM_CF_SABR_VanillaOption_Formula::SABR_G :  
		return SABR_Normal_FromATMsigma_to_alpha( f, K, tex,  sigma,  beta,  rho,  nu, flag);
		
	default :
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_SABR_From_Sigma_To_Alpha::flag : incorrect input");
	}
}

double Export_SABR2B_ImplicitVol(double fwd, 
							   double strike, 
							   double mat,
							   double alpha, 
							   double beta1, 
							   double beta2, 
							   double rho, 
							   double nu,
							   double zero,
							   double lambda)
{
	//return sabr2b_implicit_vol(fwd,strike,mat,alpha,beta1,beta2,rho,nu,zero,lambda);
	return ARM_SmileCalibration_SABR2beta::vol(fwd,strike,mat,alpha,beta1,beta2,rho,nu,zero,lambda);
}

CC_END_NAMESPACE()
 

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/