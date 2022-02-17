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

#include <cmath>
#include <algorithm>

#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/inverse.h"
#include "gpclosedforms/erf.h"
#include "gpclosedforms/gamma.h"

#include "gpbase\utilityport.h"

#include "gpclosedforms/sabr_calibration.h"
#include "gpclosedforms/extendedsabrformula.h"

#include "gpbase/numericconstant.h"

#include <glob/expt.h>

using std::sqrt;
using std::exp;


#define ARM_GP_CF_SABR_SERIE_USAGE_LIMIT 0.05
#define ARM_GP_CF_SABR_MAX_STRIKE 10000
#define ARM_GP_CF_SABR_MAX_VOL 10000

CC_BEGIN_NAMESPACE(ARM)

double SABR_BetEqOne_ImplVol(double f,double K,double tex, double alpha,double rho, double nu);

double SABR_DirectExact_ImplicitVol(
									   double f,
									   double K1,
									   double T,
									   double alpha,
									   double beta1,
									   double rho,
									   double nu
									   )
{
	double beta,K;
	if(fabs(beta1-1.0)<0.0001)
	{
		beta=0.9999;
	}
	else
	{
		beta=beta1;
	}
	if(K1<0.0001*f)
	{
		K=0.0001*f;
	}
	else
	{
		K=K1;
	}
	double vol;
    if ( fabs(K-f) < alpha*ARM_GP_CF_SABR_SERIE_USAGE_LIMIT*f )
    {
		double	a0 = alpha*pow(f,-1 + beta)*pow(1 - (T*pow(f,-2)*(pow(alpha,2)*pow(-1 + beta,2)*pow(f,2*beta) + 6*alpha*beta*nu*rho*pow(f,1 + beta) + 
			pow(f,2)*pow(nu,2)*(2 - 3*pow(rho,2))))/12.,-0.5);
		
		double    a1 = pow(3,0.5)*(-3*(-1 + beta)*beta*nu*rho*T*pow(alpha,2)*pow(f,2*beta) + alpha*(-1 + beta)*pow(f,1 + beta)*(12 + T*pow(nu,2)*(-2 + 3*pow(rho,2))) + 
			nu*rho*pow(f,2)*(12 + T*pow(nu,2)*(-5 + 6*pow(rho,2))))*pow(-(T*pow(alpha,2)*pow(-1 + beta,2)*pow(f,2*beta)) - 6*alpha*beta*nu*rho*T*pow(f,1 + beta) + 
			pow(f,2)*(12 + T*pow(nu,2)*(-2 + 3*pow(rho,2))),-1.5);
		
		double   a2 = (pow(3,-0.5)*pow(alpha,-1)*pow(f,-2 - beta)*(3*(-1 + beta)*T*pow(alpha,4)*pow(f,2 + 4*beta)*
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
		
		return a0+(K-f)*a1+(K-f)*(K-f)*a2;
    }
	else if (K > ARM_GP_CF_SABR_MAX_STRIKE*f)
		return ARM_GP_CF_SABR_MAX_VOL;
    else
    {
	/*
	if ((beta<1.00001)&&(beta>0.99999))
	{
	return CptSABR_BetEqOne_ImplVol( f, K, tex,  alpha, rho,  nu);
	}
		*/
		double	z  = (pow(f,(1.-beta))-pow(K,(1.-beta)))/(alpha*(1.-beta));
		double	b1 = beta*pow((f+K)/2.,beta-1.);
		
		double z1     = sqrt(1.-2.*nu*rho*z+nu*nu*z*z);
		
        double x     = log((-rho+nu*z+z1)/(1.-rho))/nu;
		
        double theta = 0.25*alpha*b1*nu*rho*z*z
			+log(alpha*pow((f*K),beta/2.)*z/(f-K))
			+log(x*sqrt(sqrt(1.-2.*nu*rho*z+nu*nu*z*z))/z);
		
		vol = log(f/K)/(x*sqrt(1.-2.*T*theta/(x*x)+2.*T/(x*x)*log(sqrt(f*K)*log(f/K)/(f-K))));
    }
	
    return(vol);
}

double SABR_StrikeCuterExp_ImplicitVol(
									   double f,
									   double K1,
									   double T,
									   double alpha,
									   double beta,
									   double rho,
									   double nu,
									   double alpha_exp)
{
	double K;
	if (fabs(alpha_exp)<1e-10)
	{
		K=K1;
	}
	else
	{
		K=K1+exp(-K1/alpha_exp)*alpha_exp;
	}
	return SABR_DirectExact_ImplicitVol(f,K,T,alpha,beta,rho,nu);
}


double SABR_StrikeCuterTanh_ImplicitVol(
									   double f,
									   double K1,
									   double T,
									   double alpha,
									   double beta,
									   double rho,
									   double nu,
									   double alpha_tanh,
									   double kb_tanh)
{
	double K;
	if (fabs(alpha_tanh*kb_tanh)<1e-10)
	{
		K=K1;
	}
	else
	{
		if(K1>kb_tanh)
		{
			K=K1;
		}
		else
		{
			K=kb_tanh*(1.0+alpha_tanh*tanh((K1-kb_tanh)/(alpha_tanh*kb_tanh)));
		}
		
	}
	return SABR_DirectExact_ImplicitVol(f,K,T,alpha,beta,rho,nu);
}


double SABR_StrikeCuterExp_VanillaOption(
									   double f,
									   double K,
									   double T,
									   double alpha,
									   double beta,
									   double rho,
									   double nu,
									   double alpha_exp,
									   int callput)
{
	return BlackSholes_Formula(f,
							 SABR_StrikeCuterExp_ImplicitVol(f,K,T,alpha,beta,rho,nu,alpha_exp),
							 1.0,
							 K,
							 T,
							 callput);
}

double SABR_StrikeCuterTanh_VanillaOption(
									   double f,
									   double K,
									   double T,
									   double alpha,
									   double beta,
									   double rho,
									   double nu,
									   double alpha_tanh,
									   double kb_tanh,
									   int callput)
{
	return BlackSholes_Formula(f,
							 SABR_StrikeCuterTanh_ImplicitVol(f,K,T,alpha,beta,rho,nu,alpha_tanh,kb_tanh),
							 1.0,
							 K,
							 T,
							 callput);
}


/// Set of Derivatives for the SABR Model
void SABR_DetermineAllDerivatives_withstrikecuter( double f,
     double K2,
     double T, 
     double alpha, 
     double beta1, 
     double rho, 
     double nu,
     int    flag,
	 double alpha_exp,
	 double alpha_tanh,
	 double kb_tanh,
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
	double beta,K1,K;
	if(fabs(beta1-1.0)<0.0001)
	{
		beta=0.9999;
	}
	else
	{
		beta=beta1;
	}
	if(K2<0.0001*f)
	{
		K1=0.0001*f;
	}
	else
	{
		K1=K2;
	}
	switch(flag)
	{	
	case ARM_CF_SABR_ImplicitVol_Formula::SABR_DIRECT_SC1 :
		{
			if (fabs(alpha_exp)<1e-10)
			{
				K=K1;
			}
			else
			{
				K=K1+exp(-K1/alpha_exp)*alpha_exp;
			}
			break;
		}
	case ARM_CF_SABR_ImplicitVol_Formula::SABR_DIRECT_SC2 :
		{
			if (fabs(alpha_tanh*kb_tanh)<1e-10)
			{
				K=K1;
			}
			else
			{
				if(K1>kb_tanh)
				{
					K=K1;
				}
				else
				{
					K=kb_tanh*(1.0+alpha_tanh*tanh((K1-kb_tanh)/(alpha_tanh*kb_tanh)));
				}
				
			}
			break;
		}
	default :
		{	
			K=K1;
			break;
		}
	}
//	int flag1=ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACT;
	int flag1=flag;
	SABR_DetermineAllDerivatives(f,K,T,alpha,beta,rho,nu,flag1,implvol,
		der_alpha,der_beta,der_rho,der_nu,der_f,
		der_alphaflag,der_betaflag,der_rhoflag,der_nuflag,der_fflag);
}

/// Calibration of 3 parameters : Beta Fixed
SABR_ParameterSet*  SABR_CalibrateToSmile_withstrikecuter(ARM_GP_Vector* K_Vec,ARM_GP_Vector* ImpVol_Vec,
										  ARM_GP_Vector* Weigth_Vec,
										  double f,double beta,double tex,int flag,
										  double alpha_exp,double alpha_tanh, double kb_tanh,
									      int nbsteps,int algorithm,double alpha0,double rho0,double nu0)
{
	int i;
	int n=3;
	int m=K_Vec->size();
	if(ImpVol_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"SABR_CalibrateToSmile_withstrikecuter: K_Vec and ImpVol_Vec do not have the same size!" );
	}
	ARM_GP_Vector* strikelist= new ARM_GP_Vector(m);
	ARM_GP_Vector* pricelist= new ARM_GP_Vector(m);
	ARM_GP_Vector* weightlist= new ARM_GP_Vector(m);
	ARM_GP_Vector* initialparamlist= new ARM_GP_Vector(n);
	ARM_GP_Vector* lowerboundaryparamlist= new ARM_GP_Vector(n);
	ARM_GP_Vector* upperboundaryparamlist= new ARM_GP_Vector(n);
	for(i=0;i<m;i++)
	{
		(*strikelist)[i]=(*K_Vec)[i];
		(*pricelist)[i]=(*ImpVol_Vec)[i];
		(*weightlist)[i]=1.;
	}
	(*initialparamlist)[0]=alpha0;
	(*initialparamlist)[1]=rho0;
	(*initialparamlist)[2]=nu0;
	(*lowerboundaryparamlist)[0]=0.001;
	(*lowerboundaryparamlist)[1]=-0.9999;
	(*lowerboundaryparamlist)[2]=0.00001;
	(*upperboundaryparamlist)[0]=100000000;;
	(*upperboundaryparamlist)[1]=0.9999;
	(*upperboundaryparamlist)[2]=10000000;

	class objectiveFuntion : public Optimization_ObjectiveFuntion
	{
	private:
		ARM_GP_Vector* strike_list;
		double maturity;
		int flag;
		int nbsteps;
		double forward;
		double beta;
		double alpha_exp;
		double alpha_tanh;
		double kb_tanh;
	public:
		void NAG_CALL operator() (Integer m, Integer n, 
			double x[], /// input
			double f[],	/// output (f(x))
			double fjac[],  /// output  (Df(x,i))
			Integer tdfjac, Nag_Comm *comm)
			
		{
			int i;
			int m_max=(m<strike_list->size())?m:strike_list->size();
			double K,alpha,rho,nu,implvol,der_alpha,der_beta,der_rho,der_nu,der_f;
			for(i=0;i<m_max;i++)
			{
				K=(*strike_list)[i];
				alpha=	x[0];
				rho=	x[1];
				nu=		x[2];
				double beta_eff;
				if (beta>0.9999)
				{
					beta_eff=0.9999;
				}
				else
				{
					beta_eff=beta;
				}
				
				SABR_DetermineAllDerivatives_withstrikecuter(forward,
                    K,
                    maturity,
                    alpha,
                    beta_eff,
                    rho,
                    nu,
					flag,
					alpha_exp,
					alpha_tanh,
					kb_tanh,
                    &implvol,
                    &der_alpha,
                    &der_beta,
                    &der_rho,
                    &der_nu,
                    &der_f,
					1,
                    0,
                    1,
                    1,
                    0);
					
				f[i]		=implvol;
				fjac[i*3]	=der_alpha;
				fjac[i*3+1]	=der_rho;
				fjac[i*3+2]	=der_nu;
			}
			
		}
		objectiveFuntion(ARM_GP_Vector* strike_list0,double f0,double beta0,double tex0,int flag0,
			int alpha_exp0,int alpha_tanh0,int kb_tanh0,int nbsteps0):
		strike_list(strike_list0),maturity(tex0),nbsteps(nbsteps0),flag(flag0),forward(f0),beta(beta0),
			alpha_exp(alpha_exp0),alpha_tanh(alpha_tanh0),kb_tanh(kb_tanh0)
		{}
		
	};
	objectiveFuntion func(strikelist,f,beta,tex, flag,alpha_exp,alpha_tanh,kb_tanh, nbsteps);
	string tracefile("C:\\Nag_SABR_CAlibrate.txt");

	Optimization_Result_Set* result=OptimizeWithDerivatives(
		strikelist,
		pricelist,
		Weigth_Vec,
		&func,
		initialparamlist,
		lowerboundaryparamlist,
		upperboundaryparamlist,
		algorithm,
		TRUE,
//		FALSE,
		tracefile
		);
	SABR_ParameterSet* setptr= new SABR_ParameterSet(*result,f,beta);
	return setptr;

	delete strikelist;
	delete pricelist;
	delete weightlist;
	delete initialparamlist;
	delete lowerboundaryparamlist;
	delete upperboundaryparamlist;
	

}


CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/