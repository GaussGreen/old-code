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
#include "gpclosedforms/vanille_bs_formula.h"
#include "gpclosedforms/sabrimpliedvol.h"
//#include "gpclosedforms/extended_sabr_formula.h"
#include "gpclosedforms/inverse.h"
#include "gpclosedforms/smile_sabr.h"
#include "gpclosedforms/erf.h"
#include "gpclosedforms/extendedsabrformula.h"
#include "gpclosedforms/sabrbdiff1.h"
#include "gpclosedforms/gamma.h"

#include "gpbase\utilityport.h"

#include "gpclosedforms/sabr_calibration.h"

#include "gpbase/numericconstant.h"

#include <glob/expt.h>

using std::sqrt;
using std::exp;


#define ARM_GP_CF_SABR_SERIE_USAGE_LIMIT 0.05


CC_BEGIN_NAMESPACE(ARM)

//////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///			SABR Implicit, Vol Direct Method,  
///			value of z= exact (flag=DIRECTEXACT), Geom (flag=DIRECTGEOMETRIC) , arith (flag=DIRECTARITHMETIC)
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////

double SABR_implicit_vol_direct_Series(double f,
		                               double K,
		                               double T, 
		                               double alpha,
		                               double beta,
		                               double rho, 
		                               double nu,
		                               int flag)
{
    double res = CptSABR_implicit_vol_direct_Series(f, K, T, alpha, beta, rho, nu, flag);

    return res;
}



double SABR_implicit_vol_direct(double f,
		                        double K,
		                        double tex, 
		                        double alpha,
		                        double beta, 
		                        double rho, 
		                        double nu,
		                        int flag)
{
    double res = CptSABR_implicit_vol_direct(f, K, tex, alpha, beta, rho, nu, flag);

    return res;
}



double SABR_Direct_FromATMsigma_to_alpha(double f,
		                                 double K,
		                                 double tex,
		                                 double sigma, 
		                                 double beta,
		                                 double rho, 
		                                 double nu,
		                                 int flag)
{
	class SABR_vol: public DoubleToDoubleFunc 
	{
	public: 
		double f0;
		double K0;
		double tex0;
		double beta0;
		double rho0;
		double nu0;
		int flag0;
		SABR_vol(double f1,
            double K1,
            double tex1, 
            double beta1,
            double rho1, 
            double nu1,
            int flag1)
       :
		f0(f1),
        K0(K1),
        tex0(tex1),
        beta0(beta1),
        rho0(rho1),
        nu0(nu1),
        flag0(flag1)
        {}

		virtual double operator() (double alpha0)  const
		{
			return SABR_implicit_vol_direct(f0,K0,tex0,alpha0,beta0,rho0,nu0,flag0);
		}
	};

	SABR_vol x(f,K,tex,beta,rho,nu,flag);
	double gess=sigma*pow(f,1.-beta);
	return Inverse(x,Inverse::ALWAYSPOSITIVE)(sigma,gess,gess/10,1e-12); 
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///			SABR Implicit, Vol Normal Method,  
///		 value of z= exact (flag=NORMALEXACT), Geom (flag=NORMALGEOMETRIC) , arith (flag=NORMALARITHMETIC)
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////
double SABR_implicit_vol_normal_Series(double f,
                                       double K,
                                       double tex,
                                       double alpha,
                                       double beta, 
                                       double rho, 
                                       double nu,
                                       int flag)
{
	double res = CptSABR_implicit_vol_normal_Series(f, K, tex,alpha, beta,rho, nu,flag);
	
	return res;
}



double SABR_implicit_vol_normal(double f,
                                double K,
                                double tex,
                                double alpha,
                                double beta,
                                double rho, 
                                double nu,
                                int flag)
{
    double res = CptSABR_implicit_vol_normal(f, K, tex,alpha, beta, rho,nu,flag);

    return res;
}



double SABR_Normal_FromATMsigma_to_alpha(double f,
                                         double K,
                                         double tex, 
                                         double sigma,
                                         double beta,
                                         double rho, 
                                         double nu,
                                         int flag)
{
	class SABR_vol: public DoubleToDoubleFunc 
	{
	public: 
		double f0;
		double K0;
		double tex0;
		double beta0;
		double rho0;
		double nu0;
		int flag0;
		SABR_vol(double f1,
            double K1,
            double tex1, 
            double beta1,
            double rho1, 
            double nu1,
            int flag1):
		f0(f1),
        K0(K1),
        tex0(tex1),
        beta0(beta1),
        rho0(rho1),
        nu0(nu1),
        flag0(flag1) 
        {}

		virtual double operator() (double alpha0)  const
		{
			return SABR_implicit_vol_normal(f0,K0,tex0,alpha0,beta0,rho0,nu0,flag0);
		}
	};

	SABR_vol x(f,K,tex,beta,rho,nu,flag);
	double gess=sigma*pow(f,1.-beta);
	return Inverse(x,Inverse::ALWAYSPOSITIVE)(sigma,gess,gess/10,1e-12); 


}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///				Numerical integration method; zp =0
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/// 
double SABR_NumericalIntegration_Call(double f,
                                      double K, 
                                      double tex,
                                      double alpha, 
                                      double beta, 
                                      double rho, 
                                      double nu,
                                      int nbsteps)
{
	double z,x,b1,theta,xhi,kappa,theta0,theta2,theta3,xhi0,xhi1,xhi2,xhi3;
	
	if(beta==1.)
	{
		z=log(f/K)/alpha;
	}
	else
	{
		z=(pow(f,1. - beta) - pow(K,1. - beta))/(alpha*(1.0-beta));
	}
    x=log(pow(1. - rho,-1)*(-rho + nu*z + pow(1. - 2.*nu*rho*z + nu*nu*pow(z,2),0.5)))*pow(nu,-1);
    b1=beta*pow(2,1. - beta)*pow(f + K,-1. + beta);
	if ((fabs(f-K)<0.002)&&(beta!=1.0))
	{
		theta0=log(pow(K,-beta)*pow(pow(pow(K,1. - beta),-2*beta*pow(-1. + beta,-1)),0.5));
		theta2=(pow(alpha,-2)*pow(K,-2*(1. + beta))*((-2 + beta)*beta*alpha*alpha*
			pow(K,2*beta) + 6.*alpha*beta*nu*rho*pow(K,1. + beta) + 
			pow(K,2)*nu*nu*(2. - 3.*rho*rho)))/24.;
		theta3=(pow(alpha,-3)*pow(K,-3*(1 + beta))*(-((-2 + beta)*
			beta*pow(alpha,3)*pow(K,3*beta)) - 6*nu*rho*alpha*alpha*pow(beta,2)*pow(K,1. + 2*beta) + 
			rho*pow(K,3)*pow(nu,3)*(5. - 6.*rho*rho) + alpha*beta*pow(K,2 + beta)*nu*nu*(-2. + 3.*rho*rho)))/24.;
		theta=theta0+theta2*(f-K)*(f-K)+theta3*(f-K)*(f-K)*(f-K);
		xhi0=alpha*pow(K,beta);
		xhi1=(-(nu*rho) + alpha*beta*pow(K,-1. + beta))/2.;
		xhi2=(pow(alpha,-1)*pow(K,-2 - beta)*((-2. + beta)*beta*alpha*alpha*
			pow(K,2*beta) + pow(K,2)*nu*nu*(2. - 3.*rho*rho)))/12.;
		xhi3=(pow(alpha,-2)*pow(K,-3 - 2*beta)*(-((-2. + beta)*beta*pow(alpha,3)*
			pow(K,3*beta)) + rho*pow(K,3)*pow(nu,3)*(5. - 6.*rho*rho) + 
			alpha*beta*pow(K,2 + beta)*nu*nu*(-2. + 3.*rho*rho)))/24.;
		xhi=xhi0+(f-K)*xhi1+(f-K)*(f-K)*xhi2+(f-K)*(f-K)*(f-K)*xhi3;
	}
	else
	{
		if((fabs(f-K)<1e-10)&&(beta==1.0))
		{
			theta=log(alpha*pow(f*K,beta/2.)/(K*alpha)) + 
				log(pow(1. - 2.*nu*rho*z + 
				nu*nu*pow(z,2),0.25)) + 
				(alpha*b1*nu*rho*pow(z,2))/4.;
		}
		else
		{
			theta=log(alpha*z*pow(f - K,-1)*pow(f*K,beta/2.)) + 
				log(x*pow(z,-1)*pow(1. - 2.*nu*rho*z + 
				nu*nu*pow(z,2),0.25)) + 
				(alpha*b1*nu*rho*pow(z,2))/4.;
		}
		
		
		if(fabs(x)<1e-10)
		{
			xhi=K*alpha;
		}
		else
		{
			xhi=(f-K)/x;
		}
	}
	kappa=(((-2. + beta)*beta*alpha*alpha*pow(K,2*beta) +
		6.*alpha*beta*nu*rho*pow(K,1. + beta) + K*K*nu*nu*(2. - 3.*rho*rho)))/(8.*K*K);
	
	return CC_Max(f-K,0.0)+xhi/(2.*ARM_NumericConstants::ARM_SQRT_2_PI)*exp(theta)*SABR_expintegral(x*x/2.,kappa,tex,nbsteps);
	
}



double SABR_NumericalIntegration_ImplicitVol(double f,
       double K,
       double tex, 
       double alpha,
       double beta, 
       double rho, 
       double nu, 
       int nbsteps)
{
	double opt=SABR_NumericalIntegration_Call(f,K,tex,alpha,beta,rho,nu,nbsteps);
	return ARM_CF_BS_Formula::callimplicit_totalvolatility(f, 1., K, 1,opt,1e-8)/sqrt(tex);
}



double SABR_NumericalIntegration_FromATMsigma_to_alpha(double f,
       double K,
       double tex, 
       double sigma, 
       double beta, 
       double rho, 
       double nu,
       int nbsteps)
{
	class SABR_vol: public DoubleToDoubleFunc 
	{
	public: 
		double f0;
		double K0;
		double tex0;
		double beta0;
		double rho0;
		double nu0;
		int nbsteps0;
		SABR_vol(double f1,
            double K1,
            double tex1,
            double beta1, 
            double rho1, 
            double nu1,
            int nbsteps1)
      :
		f0(f1),
            K0(K1),
            tex0(tex1),
            beta0(beta1),
            rho0(rho1),
            nu0(nu1),
            nbsteps0(nbsteps1) 
        {}

		virtual double operator() (double alpha0)  const
		{
			return SABR_NumericalIntegration_ImplicitVol(f0,K0,tex0,alpha0,beta0,rho0,nu0,nbsteps0);
		}
	};

	SABR_vol x(f,K,tex,beta,rho,nu,nbsteps);
	double gess=sigma*pow(f,1.-beta);
	return Inverse(x,Inverse::ALWAYSPOSITIVE)(sigma,gess,gess/10,1e-12); 


}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///				Numerical integration method; zp =Z((f+K)/2,K)
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/// 
double SABR_NumericalIntegration_Call2(double f, 
         double K, 
         double tex, 
         double alpha,
         double beta, 
         double rho, 
         double nu, 
         int nbsteps)
{
	double z,zp,fp,x,b1,theta,xhi,kappa,theta0,theta2,theta3,xhi0,xhi1,xhi2,xhi3;
		z=pow(alpha,-1)*pow(1 - beta,-1)*
			(pow(f,1 - beta) - pow(K,1 - beta));	
		
		x=log(pow(1 - rho,-1)*(-rho + nu*z + 
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			0.5)))*pow(nu,-1);
		fp=(f+K)/2.;				// choix heuristque
		zp=pow(alpha,-1)*pow(1. - beta,-1)*
			(pow(fp,1. - beta) - pow(K,1. - beta));
		b1=beta*pow(alpha*(1 - beta)*zp + pow(K,1 - beta),-1);
		if (fabs(f-K)<0.002)
		{
			theta0=log(pow(K,-beta)*pow(pow(pow(K,1. - beta),-2*beta*pow(-1. + beta,-1)),0.5));
			theta2=(pow(alpha,-2)*pow(K,-2*(1. + beta))*((-2 + beta)*beta*alpha*alpha*
				pow(K,2*beta) + 6.*alpha*beta*nu*rho*pow(K,1. + beta) + 
				pow(K,2)*nu*nu*(2. - 3.*rho*rho)))/24.;
			theta3=(pow(alpha,-3)*pow(K,-3*(1 + beta))*(-((-2 + beta)*
				beta*pow(alpha,3)*pow(K,3*beta)) - 6*nu*rho*alpha*alpha*pow(beta,2)*pow(K,1. + 2*beta) + 
				rho*pow(K,3)*pow(nu,3)*(5. - 6.*rho*rho) + alpha*beta*pow(K,2 + beta)*nu*nu*(-2. + 3.*rho*rho)))/24.;
			theta=theta0+theta2*(f-K)*(f-K)+theta3*(f-K)*(f-K)*(f-K);
			xhi0=alpha*pow(K,beta);
			xhi1=(-(nu*rho) + alpha*beta*pow(K,-1. + beta))/2.;
			xhi2=(pow(alpha,-1)*pow(K,-2 - beta)*((-2. + beta)*beta*alpha*alpha*
				pow(K,2*beta) + pow(K,2)*nu*nu*(2. - 3.*rho*rho)))/12.;
			xhi3=(pow(alpha,-2)*pow(K,-3 - 2*beta)*(-((-2. + beta)*beta*pow(alpha,3)*
				pow(K,3*beta)) + rho*pow(K,3)*pow(nu,3)*(5. - 6.*rho*rho) + 
				alpha*beta*pow(K,2 + beta)*nu*nu*(-2. + 3.*rho*rho)))/24.;
			xhi=xhi0+(f-K)*xhi1+(f-K)*(f-K)*xhi2+(f-K)*(f-K)*(f-K)*xhi3;
		}
		else
		{
		
		theta=log(alpha*z*pow(f - K,-1)*pow(f*K,beta/2.)) + 
			log(x*pow(z,-1)*pow(1 - 2*nu*rho*z + 
			nu*nu*pow(z,2),0.25)) + 
			(alpha*b1*nu*rho*pow(z,2))/4.;
		xhi=(f-K)/x;
		}
		
		kappa=((-2. + beta)*beta*alpha*alpha*pow(K,2.*beta)*
			pow(K - alpha*(-1. + beta)*zp*pow(K,beta),-2) + 
			6*alpha*beta*nu*rho*pow(K,beta)*
			pow(K - alpha*(-1 + beta)*zp*pow(K,beta),-1) + 
			nu*nu*(2 + 2*nu*rho*zp - 3*rho*rho - 
			nu*nu*pow(zp,2))*
			pow(1. - 2.*nu*rho*zp + nu*nu*pow(zp,2),-1))/8.;

	return CC_Max(f-K,0.0)+xhi/(2.*ARM_NumericConstants::ARM_SQRT_2_PI)*exp(theta)*SABR_expintegral(x*x/2.,kappa,tex,nbsteps);

}



double SABR_NumericalIntegration_ImplicitVol2(double f, double K, double tex, double alpha, double beta, double rho, double nu, int nbsteps)
{
	double opt=SABR_NumericalIntegration_Call2(f,K,tex,alpha,beta,rho,nu,nbsteps);
	return ARM_CF_BS_Formula::callimplicit_totalvolatility(f, 1., K, 1,opt,1e-8)/sqrt(tex);
}



double SABR_NumericalIntegration_FromATMsigma_to_alpha2(double f,double K,double tex, double sigma, double beta, double rho, double nu,int nbsteps)
{
	class SABR_vol: public DoubleToDoubleFunc 
	{
	public: 
		double f0;
		double K0;
		double tex0;
		double beta0;
		double rho0;
		double nu0;
		int nbsteps0;
		SABR_vol(double f1,double K1,double tex1,  double beta1, double rho1, double nu1,int nbsteps1):
		f0(f1),K0(K1),tex0(tex1),beta0(beta1),rho0(rho1),nu0(nu1),nbsteps0(nbsteps1) {}

		virtual double operator() (double alpha0)  const
		{
			return SABR_NumericalIntegration_ImplicitVol2(f0,K0,tex0,alpha0,beta0,rho0,nu0,nbsteps0);
		}
	};

	SABR_vol x(f,K,tex,beta,rho,nu,nbsteps);
	double gess=sigma*pow(f,1.-beta);
	return Inverse(x,Inverse::ALWAYSPOSITIVE)(sigma,gess,gess/10,1e-12); 


}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///						zs beta=1 (SABR_A in ARM)
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////

double SABR_BetEqOne_CompatibleKernel_ImplVol(double f,double K,double tex, double alpha,double rho, double nu)
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



double SABR_BetaEqOne_CompatibleKernel_FromATMsigma_to_alpha(double f,
															 double K,
															 double tex,
															 double sigma,  
															 double rho, 
															 double nu)
{
	class SABR_vol: public DoubleToDoubleFunc 
	{
	public: 
		double f0;
		double K0;
		double tex0;
		double beta0;
		double rho0;
		double nu0;
		SABR_vol(double f1,double K1,double tex1, double rho1, double nu1):
		f0(f1),K0(K1),tex0(tex1),rho0(rho1),nu0(nu1) {}

		virtual double operator() (double alpha0)  const
		{
			return SABR_BetEqOne_CompatibleKernel_ImplVol(f0,K0,tex0,alpha0,rho0,nu0);
		}
	};

	SABR_vol x(f,K,tex,rho,nu);
	double gess=sigma*f;
	return Inverse(x,Inverse::ALWAYSPOSITIVE)(sigma,gess,gess/10,1e-12); 


}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///						Derivatives;
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Nota Bene : All the following calculations have been computed by Mathematica,
///  So despite their ugly error prone appearence, 
/// they are are certified exact. The associated Mathematica notebook is under /ARMDev
/// 	they are valid only for the flag ANALYTIC)
double SABR_ImpVol_DerMaturity(double f, double K, double tex, double alpha, double beta, double rho, double nu)

{
	if (fabs(f-K)<(f*ARM_GP_CF_SABR_SERIE_USAGE_LIMIT)) 
	{
		return (alpha*pow(f,-3 + beta)*
			(alpha*alpha*pow(-1 + beta,2)*pow(f,2*beta) + 
			6*alpha*beta*nu*rho*pow(f,1 + beta) + 
			pow(f,2)*nu*nu*(2 - 3*rho*rho)))/24.;
	}
	else
	{
		
		double z,x,b1,theta;
		z=pow(alpha,-1)*pow(1 - beta,-1)*
			(pow(f,1 - beta) - pow(K,1 - beta));
		
		
		x=log(pow(1 - rho,-1)*(-rho + nu*z + 
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			0.5)))*pow(nu,-1);
		
		b1=beta*pow(2,1 - beta)*pow(f + K,-1 + beta);
		
		theta=log(alpha*z*pow(f - K,-1)*pow(f*K,beta/2.)) + 
			log(x*pow(z,-1)*pow(1 - 2*nu*rho*z + 
			nu*nu*pow(z,2),0.25)) + 
			(alpha*b1*nu*rho*pow(z,2))/4.;
		
		return log(f*pow(K,-1))*(theta - 
			log(log(f*pow(K,-1))*pow(f - K,-1)*pow(f*K,0.5)))*
			pow(x,-3)*pow(pow(x,-2)*
			(-2*tex*theta + 2*tex*
			log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5)) + pow(x,2)),-1.5);
	}
}

double SABR_ImpVol_DerForward(double f, 
		  double K,
		  double tex,
		  double alpha,
		  double beta, 
		  double rho,
		  double nu)

{
	if (fabs(f-K)<(f*ARM_GP_CF_SABR_SERIE_USAGE_LIMIT)) 
	{
		return (SABR_ImpVol_DerForward(f, f+0.0001, tex, alpha, beta, rho, nu)+
			SABR_ImpVol_DerForward(f, f-0.0001, tex, alpha, beta, rho, nu))/2.;
	}
	else
	{
		double z,x,b1,theta;
		z=pow(alpha,-1)*pow(1 - beta,-1)*
			(pow(f,1 - beta) - pow(K,1 - beta));
		
		
		x=log(pow(1 - rho,-1)*(-rho + nu*z + 
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			0.5)))*pow(nu,-1);
		
		b1=beta*pow(2,1 - beta)*pow(f + K,-1 + beta);
		
		theta=log(alpha*z*pow(f - K,-1)*pow(f*K,beta/2.)) + 
			log(x*pow(z,-1)*pow(1 - 2*nu*rho*z + 
			nu*nu*pow(z,2),0.25)) + 
			(alpha*b1*nu*rho*pow(z,2))/4.;
		
		
		
		double DerImpVoltheta=log(f*pow(K,-1))*pow(x,-3)*
			pow(1 - 2*theta*pow(x,-2) + 
			2*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-2),-1.5);
		
		
		double DerImpVolx= -(log(f*pow(K,-1))*(4*theta*pow(x,-3) - 
			4*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-3))*pow(x,-1)*
			pow(1 - 2*theta*pow(x,-2) + 
			2*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-2),-1.5))/2. - 
			log(f*pow(K,-1))*pow(x,-2)*
			pow(1 - 2*theta*pow(x,-2) + 
			2*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-2),-0.5);
		
		double Derxz=pow(nu,-1)*(nu + ((-2*nu*rho + 2*z*nu*nu)*
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			-0.5))/2.)*pow(-rho + nu*z + 
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),0.5),
			-1);
		double DerImpVolf=(pow(f,-1)*pow(f - K,-1)*pow(x,-3)*
			((f + K)*tex*log(f*pow(K,-1)) - 
			2*(f - K)*(tex + 2*tex*theta - 
			2*tex*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5)) - pow(x,2)))*
			pow(pow(x,-2)*(-2*tex*theta + 
			2*tex*log(log(f*pow(K,-1))*pow(f - K,-1)*
            pow(f*K,0.5)) + pow(x,2)),-1.5))/2.;
		
		double Derthetab1=(alpha*nu*rho*pow(z,2))/4.;

		double Derb1f=(-1 + beta)*beta*pow(2,1 - beta)*pow(f + K,-2 + beta);
		
		double Derthetax=1/x;
		
		double Derthetaalpha=pow(alpha,-1) + (b1*nu*rho*pow(z,2))/4.;
		
		double Derzalpha=-(pow(alpha,-2)*pow(1 - beta,-1)*
			(pow(f,1 - beta) - pow(K,1 - beta)));
		
		
		double Derthetaz=(alpha*b1*nu*rho*z)/2. + pow(z,-1) + 
			z*pow(x,-1)*pow(1 - 2*nu*rho*z + 
			nu*nu*pow(z,2),-0.25)*
			((x*(-2*nu*rho + 2*z*nu*nu)*pow(z,-1)*
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			-0.75))/4. - 
			x*pow(z,-2)*pow(1 - 2*nu*rho*z + 
			nu*nu*pow(z,2),0.25));

		double Derthetaf=(2*f - beta*f + beta*K)*pow(2*f*K - 2*pow(f,2),-1);

		double Derzf=pow(alpha,-1)*pow(f,-beta);
		
		return DerImpVoltheta*( 
      Derthetaf + Derthetaz*Derzf + Derthetax*Derxz*Derzf + 
        Derthetab1*Derb1f) + DerImpVolx*Derxz*Derzf + DerImpVolf;
	}
}

double SABR_ImpVol_DerStrike(double f, 
		double K, 
		double tex,
		double alpha,
		double beta, 
		double rho,
		double nu)

{
	if (fabs(f-K)<(f*ARM_GP_CF_SABR_SERIE_USAGE_LIMIT)) 
	{
		return (SABR_ImpVol_DerStrike(f, f+0.0001, tex, alpha, beta, rho, nu)+
			SABR_ImpVol_DerStrike(f, f-0.0001, tex, alpha, beta, rho, nu))/2.;
	}
	else
	{
		double z,x,b1,theta;
		z=pow(alpha,-1)*pow(1 - beta,-1)*
			(pow(f,1 - beta) - pow(K,1 - beta));
		
		
		x=log(pow(1 - rho,-1)*(-rho + nu*z + 
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			0.5)))*pow(nu,-1);
		
		b1=beta*pow(2,1 - beta)*pow(f + K,-1 + beta);
		
		theta=log(alpha*z*pow(f - K,-1)*pow(f*K,beta/2.)) + 
			log(x*pow(z,-1)*pow(1 - 2*nu*rho*z + 
			nu*nu*pow(z,2),0.25)) + 
			(alpha*b1*nu*rho*pow(z,2))/4.;
		
		
		
		double DerImpVoltheta=log(f*pow(K,-1))*pow(x,-3)*
			pow(1 - 2*theta*pow(x,-2) + 
			2*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-2),-1.5);
		
		
		double DerImpVolx= -(log(f*pow(K,-1))*(4*theta*pow(x,-3) - 
			4*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-3))*pow(x,-1)*
			pow(1 - 2*theta*pow(x,-2) + 
			2*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-2),-1.5))/2. - 
			log(f*pow(K,-1))*pow(x,-2)*
			pow(1 - 2*theta*pow(x,-2) + 
			2*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-2),-0.5);
		
		double Derxz=pow(nu,-1)*(nu + ((-2*nu*rho + 2*z*nu*nu)*
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			-0.5))/2.)*pow(-rho + nu*z + 
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),0.5),
			-1);
		double DerImpVolK=(pow(K,-1)*pow(-f + K,-1)*pow(x,-3)*
			((f + K)*tex*log(f*pow(K,-1)) - 
			2*(f - K)*(tex + 2*tex*theta - 
			2*tex*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5)) - pow(x,2)))*
			pow(pow(x,-2)*(-2*tex*theta + 
			2*tex*log(log(f*pow(K,-1))*pow(f - K,-1)*
            pow(f*K,0.5)) + pow(x,2)),-1.5))/2.;
		
		double Derthetab1=(alpha*nu*rho*pow(z,2))/4.;
		
		double Derb1K=(-1 + beta)*beta*pow(2,1 - beta)*pow(f + K,-2 + beta);
		
		double Derthetax=1/x;
		
		double Derthetaalpha=pow(alpha,-1) + (b1*nu*rho*pow(z,2))/4.;
		
		double Derzalpha=-(pow(alpha,-2)*pow(1 - beta,-1)*
			(pow(f,1 - beta) - pow(K,1 - beta)));
		
		
		double Derthetaz=(alpha*b1*nu*rho*z)/2. + pow(z,-1) + 
			z*pow(x,-1)*pow(1 - 2*nu*rho*z + 
			nu*nu*pow(z,2),-0.25)*
			((x*(-2*nu*rho + 2*z*nu*nu)*pow(z,-1)*
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			-0.75))/4. - 
			x*pow(z,-2)*pow(1 - 2*nu*rho*z + 
			nu*nu*pow(z,2),0.25));
		
		double DerthetaK=-((beta*f + 2*K - beta*K)*pow(2*f*K - 2*pow(K,2),-1));
		
		double DerzK=-(pow(alpha,-1)*pow(K,-beta));
		
		return DerImpVoltheta*( 
			DerthetaK + Derthetaz*DerzK + Derthetax*Derxz*DerzK + 
			Derthetab1*Derb1K) + DerImpVolx*Derxz*DerzK + DerImpVolK;
	}
}



double SABR_ImpVol_DerNu(double f,
		double K, 
		double tex,
		double alpha,
		double beta,
		double rho,
		double nu)

{
	if (fabs(f-K)<(f*ARM_GP_CF_SABR_SERIE_USAGE_LIMIT)) 
	{
		return (alpha*tex*pow(f,-2 + beta)*
			(3*alpha*beta*rho*pow(f,beta) + 
			f*nu*(2 - 3*rho*rho)))/12.;
	}
	else
	{
		
		double z,x,b1,theta;
		z=pow(alpha,-1)*pow(1 - beta,-1)*
			(pow(f,1 - beta) - pow(K,1 - beta));
		
		
		x=log(pow(1 - rho,-1)*(-rho + nu*z + 
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			0.5)))*pow(nu,-1);
		
		b1=beta*pow(2,1 - beta)*pow(f + K,-1 + beta);
		
		theta=log(alpha*z*pow(f - K,-1)*pow(f*K,beta/2.)) + 
			log(x*pow(z,-1)*pow(1 - 2*nu*rho*z + 
			nu*nu*pow(z,2),0.25)) + 
			(alpha*b1*nu*rho*pow(z,2))/4.;
		
		
		
		double DerImpVoltheta=log(f*pow(K,-1))*pow(x,-3)*
			pow(1 - 2*theta*pow(x,-2) + 
			2*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-2),-1.5);
		
		
		double DerImpVolx= -(log(f*pow(K,-1))*(4*theta*pow(x,-3) - 
			4*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-3))*pow(x,-1)*
			pow(1 - 2*theta*pow(x,-2) + 
			2*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-2),-1.5))/2. - 
			log(f*pow(K,-1))*pow(x,-2)*
			pow(1 - 2*theta*pow(x,-2) + 
			2*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-2),-0.5);
		
		double Derxz=pow(nu,-1)*(nu + ((-2*nu*rho + 2*z*nu*nu)*
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			-0.5))/2.)*pow(-rho + nu*z + 
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),0.5),
			-1);
		
		double Derthetax=1/x;
		
		double Derthetanu=(z*(alpha*b1*rho*z + 2*(-rho + nu*z)*
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),-1)))/
			4.;
		
		double Derxnu=pow(nu,-2)*(-log(pow(-1 + rho,-1)*
			(rho - nu*z - 
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),0.5)
			)) + nu*z*pow(1 - 2*nu*rho*z + 
			nu*nu*pow(z,2),-0.5));
		
		return DerImpVoltheta*(Derthetanu + Derthetax*Derxnu) + DerImpVolx*Derxnu;
	}
}



double SABR_ImpVol_DerRho(double f, 
	  double K, 
	  double tex, 
	  double alpha, 
	  double beta, 
	  double rho, 
	  double nu)

{
	if (fabs(f-K)<(f*ARM_GP_CF_SABR_SERIE_USAGE_LIMIT)) 
	{
		return (alpha*nu*tex*pow(f,-2 + beta)*
			(-(f*nu*rho) + alpha*beta*pow(f,beta)))/4.;
	}
	else
	{
		
		double z,x,b1,theta;
		z=pow(alpha,-1)*pow(1 - beta,-1)*
			(pow(f,1 - beta) - pow(K,1 - beta));
		
		
		x=log(pow(1 - rho,-1)*(-rho + nu*z + 
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			0.5)))*pow(nu,-1);
		
		b1=beta*pow(2,1 - beta)*pow(f + K,-1 + beta);
		
		theta=log(alpha*z*pow(f - K,-1)*pow(f*K,beta/2.)) + 
			log(x*pow(z,-1)*pow(1 - 2*nu*rho*z + 
			nu*nu*pow(z,2),0.25)) + 
			(alpha*b1*nu*rho*pow(z,2))/4.;
		
		
		
		double DerImpVoltheta=log(f*pow(K,-1))*pow(x,-3)*
			pow(1 - 2*theta*pow(x,-2) + 
			2*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-2),-1.5);
		
		
		double DerImpVolx= -(log(f*pow(K,-1))*(4*theta*pow(x,-3) - 
			4*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-3))*pow(x,-1)*
			pow(1 - 2*theta*pow(x,-2) + 
			2*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-2),-1.5))/2. - 
			log(f*pow(K,-1))*pow(x,-2)*
			pow(1 - 2*theta*pow(x,-2) + 
			2*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-2),-0.5);
		
		double Derxz=pow(nu,-1)*(nu + ((-2*nu*rho + 2*z*nu*nu)*
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			-0.5))/2.)*pow(-rho + nu*z + 
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),0.5),
			-1);
		
		double Derthetax=1/x;
		
		double Derthetarho=(alpha*b1*nu*pow(z,2))/4. - 
			(nu*z*pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),-1))/
			2.;
		
		double Derxrho=(1 - rho)*pow(nu,-1)*pow(-1 + rho,-2)*
			(-rho + nu*z + (1 - rho)*
			(-1 - nu*z*pow(1 - 2*nu*rho*z + 
			nu*nu*pow(z,2),-0.5)) + 
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),0.5))*
			pow(-rho + nu*z + pow(1 - 2*nu*rho*z + 
			nu*nu*pow(z,2),0.5),-1);
		
		return DerImpVoltheta*(Derthetarho + Derthetax*Derxrho) + DerImpVolx*Derxrho;
	}
}



double SABR_ImpVol_DerBeta(double f, 
		double K, 
		double tex, 
		double alpha, 
		double beta, 
		double rho, 
		double nu)

{
	if (fabs(f-K)<(f*ARM_GP_CF_SABR_SERIE_USAGE_LIMIT)) 
	{
		return (alpha*pow(f,-3 + beta)*
			(2*alpha*tex*pow(f,beta)*
			(3*f*nu*rho + alpha*(-1 + beta)*pow(f,beta))\
			+ log(f)*(3*tex*alpha*alpha*pow(-1 + beta,2)*
			pow(f,2*beta) + 
			12*alpha*beta*nu*rho*tex*pow(f,1 + beta) + 
			pow(f,2)*(24 + 
			tex*nu*nu*(2 - 3*rho*rho)))))/24.;
	}
	else
	{
		
		double z,x,b1,theta;
		z=pow(alpha,-1)*pow(1 - beta,-1)*
			(pow(f,1 - beta) - pow(K,1 - beta));
		
		
		x=log(pow(1 - rho,-1)*(-rho + nu*z + 
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			0.5)))*pow(nu,-1);
		
		b1=beta*pow(2,1 - beta)*pow(f + K,-1 + beta);
		
		theta=log(alpha*z*pow(f - K,-1)*pow(f*K,beta/2.)) + 
			log(x*pow(z,-1)*pow(1 - 2*nu*rho*z + 
			nu*nu*pow(z,2),0.25)) + 
			(alpha*b1*nu*rho*pow(z,2))/4.;
		
		
		
		double DerImpVoltheta=log(f*pow(K,-1))*pow(x,-3)*
			pow(1 - 2*theta*pow(x,-2) + 
			2*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-2),-1.5);
		
		
		double DerImpVolx= -(log(f*pow(K,-1))*(4*theta*pow(x,-3) - 
			4*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-3))*pow(x,-1)*
			pow(1 - 2*theta*pow(x,-2) + 
			2*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-2),-1.5))/2. - 
			log(f*pow(K,-1))*pow(x,-2)*
			pow(1 - 2*theta*pow(x,-2) + 
			2*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-2),-0.5);
		
		double Derxz=pow(nu,-1)*(nu + ((-2*nu*rho + 2*z*nu*nu)*
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			-0.5))/2.)*pow(-rho + nu*z + 
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),0.5),
			-1);
		
		double Derthetab1=(alpha*nu*rho*pow(z,2))/4.;
		
		double Derthetax=1/x;
		
		double Derthetabeta=log(f*K)/2.;
		
		double Derzbeta=pow(alpha,-1)*pow(1 - beta,-2)*
			(pow(f,1 - beta) - pow(K,1 - beta)) + 
			pow(alpha,-1)*pow(1 - beta,-1)*
			(-(log(f)*pow(f,1 - beta)) + 
			log(K)*pow(K,1 - beta));
		// FIXMEFRED: mig.vc8 (23/05/2007 11:36:04):cast
		double Derb1beta=(1 - beta*log(2.) + beta*log(f + K))*pow(2,1 - beta)*
			pow(f + K,-1 + beta);
		
		double Derthetaz=(alpha*b1*nu*rho*z)/2. + pow(z,-1) + 
			z*pow(x,-1)*pow(1 - 2*nu*rho*z + 
			nu*nu*pow(z,2),-0.25)*
			((x*(-2*nu*rho + 2*z*nu*nu)*pow(z,-1)*
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			-0.75))/4. - 
			x*pow(z,-2)*pow(1 - 2*nu*rho*z + 
			nu*nu*pow(z,2),0.25));
		
		return DerImpVoltheta*(Derthetabeta + Derthetaz*Derzbeta + Derthetax*Derxz*Derzbeta + 
			Derthetab1*Derb1beta) + DerImpVolx*Derxz*Derzbeta;
	}
}
  

double SABR_ImpVol_DerAlpha(double f,
		double K,
		double tex,
		double alpha,
		double beta,
		double rho,
		double nu)

{
	if (fabs(f-K)<(f*ARM_GP_CF_SABR_SERIE_USAGE_LIMIT)) 
	{
		return (pow(f,-3 + beta)*(3*tex*alpha*alpha*
			pow(-1 + beta,2)*pow(f,2*beta) + 
			12*alpha*beta*nu*rho*tex*pow(f,1 + beta) + 
			pow(f,2)*(24 + 
			tex*nu*nu*(2 - 3*rho*rho))))/24.;
	}
	else
	{
		double z,x,b1,theta;
		z=pow(alpha,-1)*pow(1 - beta,-1)*
			(pow(f,1 - beta) - pow(K,1 - beta));
		
		
		x=log(pow(1 - rho,-1)*(-rho + nu*z + 
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			0.5)))*pow(nu,-1);
		
		b1=beta*pow(2,1 - beta)*pow(f + K,-1 + beta);
		
		theta=log(alpha*z*pow(f - K,-1)*pow(f*K,beta/2.)) + 
			log(x*pow(z,-1)*pow(1 - 2*nu*rho*z + 
			nu*nu*pow(z,2),0.25)) + 
			(alpha*b1*nu*rho*pow(z,2))/4.;
		
		
		
		double DerImpVoltheta=log(f*pow(K,-1))*pow(x,-3)*
			pow(1 - 2*theta*pow(x,-2) + 
			2*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-2),-1.5);
		
		
		double DerImpVolx= -(log(f*pow(K,-1))*(4*theta*pow(x,-3) - 
			4*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-3))*pow(x,-1)*
			pow(1 - 2*theta*pow(x,-2) + 
			2*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-2),-1.5))/2. - 
			log(f*pow(K,-1))*pow(x,-2)*
			pow(1 - 2*theta*pow(x,-2) + 
			2*log(log(f*pow(K,-1))*pow(f - K,-1)*
			pow(f*K,0.5))*pow(x,-2),-0.5);
		
		double Derxz=pow(nu,-1)*(nu + ((-2*nu*rho + 2*z*nu*nu)*
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			-0.5))/2.)*pow(-rho + nu*z + 
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),0.5),
			-1);
		
		double Derthetab1=(alpha*nu*rho*pow(z,2))/4.;
		
		double Derthetax=1/x;
		
		double Derthetaalpha=pow(alpha,-1) + (b1*nu*rho*pow(z,2))/4.;
		
		double Derzalpha=-(pow(alpha,-2)*pow(1 - beta,-1)*
			(pow(f,1 - beta) - pow(K,1 - beta)));
		
		
		double Derthetaz=(alpha*b1*nu*rho*z)/2. + pow(z,-1) + 
			z*pow(x,-1)*pow(1 - 2*nu*rho*z + 
			nu*nu*pow(z,2),-0.25)*
			((x*(-2*nu*rho + 2*z*nu*nu)*pow(z,-1)*
			pow(1 - 2*nu*rho*z + nu*nu*pow(z,2),
			-0.75))/4. - 
			x*pow(z,-2)*pow(1 - 2*nu*rho*z + 
			nu*nu*pow(z,2),0.25));
		
		return DerImpVoltheta*(Derthetaalpha + Derthetaz*Derzalpha + 
			Derthetax*Derxz*Derzalpha) + DerImpVolx*Derxz*Derzalpha;
	}
}

double SABR_implicit_vol_direct_DerAlpha(double f,
		                                 double K, 
		                                 double tex, 
		                                 double alpha, 
		                                 double beta, 
		                                 double rho, 
		                                 double nu,
		                                 int flag)

{
    double res =  CptSABR_implicit_vol_direct_DerAlpha(f, K, tex,alpha, beta,rho, nu,flag);

    return res;
}



double SABR_implicit_vol_normal_DerAlpha(double f,
		                                 double K, 
		                                 double tex, 
		                                 double alpha,
		                                 double beta,
		                                 double rho, 
		                                 double nu,
		                                 int flag)

{
     double res = CptSABR_implicit_vol_normal_DerAlpha(f, K, tex, alpha, beta, rho, nu,flag);

     return res;
}


////////////////////////////////////////////////////////////////////////////////////////////
///
///				computation of sum{0,t} of exp(-a/x+b*x)/sqrt(x)
///
////////////////////////////////////////////////////////////////////////////////////////////

double SABR_expintegral(double a, double b, double t,int n)
{
	complex<double> I(0,1);
	complex<double> ca_sqrt=sqrt(a);
	complex<double> cb_sqrt=b;cb_sqrt=std::sqrt(cb_sqrt);
	complex<double> ct_sqrt=sqrt(t);
	complex<double> expab=std::exp(2.*I*ca_sqrt*cb_sqrt);
	complex<double> erf1=cerf(-ca_sqrt/ct_sqrt+I*cb_sqrt*ct_sqrt,n);
	complex<double> erf2=cerf(ca_sqrt/ct_sqrt+I*cb_sqrt*ct_sqrt,n);
	double res = (((erf1+1.)/expab+(erf2-1.)*expab)/cb_sqrt).imag()*ARM_NumericConstants::ARM_SQRT_PI/2.;

	return res;

}

double sabr2b_h1(double x, double beta1, double beta2, double lambda)
{
	double aux = 1.+pow(x,beta1-beta2)*(exp(lambda*x)-1);
	return -1./aux/aux;
}

double sabr2b_h2(double x, double beta1, double beta2, double lambda)
{
	return pow(x,beta1-beta2)*((beta1-beta2)/x*(exp(lambda*x)-1.)+lambda*exp(lambda*x));
}

double sabr2b_dh1(double x, double beta1, double beta2, double lambda)
{
	double res;
	res=2./pow(1.+pow(x,beta1-beta2)*(exp(lambda*x)-1.),3.);
	res*=sabr2b_h2(x,beta1,beta2,lambda);
	return res;
}

double sabr2b_dh2(double x, double beta1, double beta2, double lambda)
{
	double res,aux;
	res=(beta1-beta2)/x*sabr2b_h2(x,beta1,beta2,lambda);
	aux=(beta1-beta2)/x/x*(1-exp(lambda*x));
	aux+=lambda*(beta1-beta2)/x*exp(lambda*x);
	aux+=lambda*lambda*exp(lambda*x);
	res+=pow(x,beta1-beta2)*aux;
	return res;
}


double sabr2b_f(double x, double beta1, double beta2, double lambda)
{
	if (x==0) return 1.;
	double x1=pow(x,beta1),x2=pow(x,beta2);
	double res=1./(x1-x2)*(1./(1./x1*exp(-lambda*x)+1./x2*(1-exp(-lambda*x)))-x2);
	return res;
}

double sabr2b_df(double x, double beta1, double beta2, double lambda)
{
	return sabr2b_h1(x,beta1,beta2,lambda)*sabr2b_h2(x,beta1,beta2,lambda);
}

double sabr2b_d2f(double x, double beta1, double beta2, double lambda)
{
	return sabr2b_dh1(x,beta1,beta2,lambda)*sabr2b_h2(x,beta1,beta2,lambda)+sabr2b_h1(x,beta1,beta2,lambda)*sabr2b_dh2(x,beta1,beta2,lambda);
}

double sabr2b_C(double x, double beta1, double beta2, double lambda)
{
	double fx=sabr2b_f(x,beta1,beta2,lambda);
	double res=pow(x,beta1)*fx+pow(x,beta2)*(1-fx);
	return res;
}

double sabr2b_dC(double x, double beta1, double beta2, double lambda)
{
	double x1=pow(x,beta1),x2=pow(x,beta2);
	return (beta1*x1/x-beta2*x2/x)*sabr2b_f(x,beta1,beta2,lambda)+(x1-x2)*sabr2b_df(x,beta1,beta2,lambda)+beta2*x2/x;
}

double sabr2b_d2C(double x, double beta1, double beta2, double lambda)
{
	double x1=pow(x,beta1),x2=pow(x,beta2);
	double res;
	res=(beta1*(beta1-1)*x1/x/x-beta2*(beta2-1)*x2/x/x)*sabr2b_f(x,beta1,beta2,lambda);
	res+=2*(beta1*x1/x-beta2*x2/x)*sabr2b_df(x,beta1,beta2,lambda)+(x1-x2)*sabr2b_d2f(x,beta1,beta2,lambda);
	res+=beta2*(beta2-1)*x2/x/x;
	return res;
}

double sabr2b_I(double f,double k, double beta1, double beta2, double lambda)
{
	double g1=exp(gammalog(1-beta1)),g2=exp(gammalog(1-beta2));
	double g_1f=gammp(1-beta1,lambda*f),g_1k=gammp(1-beta1,lambda*k);
	double g_2f=gammp(1-beta2,lambda*f),g_2k=gammp(1-beta2,lambda*k);
	double res;
	res=pow(lambda,beta1-1)*g1*(g_1f-g_1k)-pow(lambda,beta2-1)*g2*(g_2f-g_2k);
	res+=1./(1-beta2)*(pow(f,1-beta2)-pow(k,1-beta2));
	return res;
}

double sabr2b_atmvol(	double fwd, double maturity,
						double alpha, double beta1, double beta2, double rho,
						double nu, double zero, double lambda)
{
	double fav,g1,g2,cfav;
	double sabr3;
	double forward=fwd-zero;

	fav=forward;
	cfav=sabr2b_C(fav,beta1,beta2,lambda);
	g1=sabr2b_dC(fav,beta1,beta2,lambda)/cfav;g2=sabr2b_d2C(fav,beta1,beta2,lambda)/cfav;
	sabr3=(2*g2-g1*g1+1./fav/fav)/24*alpha*cfav*alpha*cfav+0.25*rho*nu*alpha*g1+(2-3*rho*rho)/24*nu*nu;

	double res=alpha*cfav/fav*(1+maturity*sabr3);
	if (res>0.)
		return res;
	return res;
}

double sabr2b_implicit_vol(	double fwd, double k, double maturity,
							double alpha, double beta1, double beta2, double rho,
							double nu, double zero, double lambda)
{
	if (fabs(fwd-k)<K_NEW_DOUBLE_TOL)
		return sabr2b_atmvol(fwd,maturity,alpha,beta1,beta2,rho,nu,zero,lambda);

	double fav,g1,g2,zeta,xzeta,ifk,cfav;
	double sabr1,sabr2,sabr3;
	double forward=fwd-zero,strike=k-zero;

	fav=sqrt(forward*strike);
	cfav=sabr2b_C(fav,beta1,beta2,lambda);
	g1=sabr2b_dC(fav,beta1,beta2,lambda)/cfav;g2=sabr2b_d2C(fav,beta1,beta2,lambda)/cfav;
	zeta=nu/alpha*(forward-strike)/cfav;
	xzeta=log((sqrt(1-2.*rho*zeta+zeta*zeta)-rho+zeta)/(1-rho));
	ifk=sabr2b_I(forward,strike,beta1,beta2,lambda);
	sabr1=alpha*log(forward/strike)/ifk;
	sabr2=zeta/xzeta;
	sabr3=(2*g2-g1*g1+1./fav/fav)/24*alpha*cfav*alpha*cfav+0.25*rho*nu*alpha*g1*cfav+(2-3*rho*rho)/24*nu*nu;
	double res=sabr1*sabr2*(1+maturity*sabr3);
	if (res>0.)
		return res;
	return res;
}


CC_END_NAMESPACE()
 

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/