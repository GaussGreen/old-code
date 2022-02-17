/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file sabr_calibration.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include "firsttoinc.h"
#include "gpbase/port.h"
#include "gpbase/gpvector.h"


#include <cmath>
#include <complex>
#include "gpbase/gpmatrix.h"

#include "gpbase/removenagwarning.h"
#include "nag.h"
#include "nage04.h"

#include "gpnumlib/gaussiananalytics.h"
#include "gpbase/numericconstant.h"

#include "gpclosedforms/heston.h"

#include "gpclosedforms/heston_calibration.h"
#include "gpclosedforms/optimization1.h"
#include "gpclosedforms/heston_formula.h"


#include "expt.h"   // for the exceptions

using namespace std;

CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_TOLERANCE 1.0e-13
#define ARM_CF_MAXIT 150




/// Calibration of the 8 parameters
GeneralizedHeston_ParameterSet*  GeneralizedHeston_CalibrateToSmile(std::vector<double>* K_Vec,std::vector<double>* ImpVol_Vec,double F,double tex,int nbsteps,int algorithm,
												double V0_0,double omega_0,double theta_0,double ksi_0,double rho_0,double muJ_0,double sigmaJ_0,double lambda_0)
{
	int i;
	int n=8;
	int m=K_Vec->size();
	if(ImpVol_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"GeneralizedHeston_CalibrateToSmile: K_Vec and ImpVol_Vec do not have the same size!" );
	}
	std::vector<double>* strikelist= new std::vector<double>(m);
	std::vector<double>* pricelist= new std::vector<double>(m);
	std::vector<double>* weightlist= new std::vector<double>(m);
	std::vector<double>* initialparamlist= new std::vector<double>(n);
	std::vector<double>* lowerboundaryparamlist= new std::vector<double>(n);
	std::vector<double>* upperboundaryparamlist= new std::vector<double>(n);
	for(i=0;i<m;i++)
	{
		(*strikelist)[i]=(*K_Vec)[i];
		(*pricelist)[i]=(*ImpVol_Vec)[i];
		(*weightlist)[i]=1.;
	}
	(*initialparamlist)[0]=	V0_0			;
	(*initialparamlist)[1]=	omega_0		;
	(*initialparamlist)[2]=	theta_0		;
	(*initialparamlist)[3]=	ksi_0			;
	(*initialparamlist)[4]=	rho_0			;
	(*initialparamlist)[5]=	muJ_0			;
	(*initialparamlist)[6]=	sigmaJ_0		;
	(*initialparamlist)[7]=	lambda_0		;
	(*lowerboundaryparamlist)[0]=0.0001;
	(*lowerboundaryparamlist)[1]=0.0001;
	(*lowerboundaryparamlist)[2]=0.0001;
	(*lowerboundaryparamlist)[3]=0.00001;
	(*lowerboundaryparamlist)[4]=0.00001;
	(*lowerboundaryparamlist)[5]=-10000000.0;
	(*lowerboundaryparamlist)[6]=0.00001;
	(*lowerboundaryparamlist)[7]=0.00001;
	(*upperboundaryparamlist)[0]=10000000.0;;
	(*upperboundaryparamlist)[1]=10000000.0;
	(*upperboundaryparamlist)[2]=10000000.0;
	(*upperboundaryparamlist)[3]=10000000.0;
	(*upperboundaryparamlist)[4]=0.9999;
	(*upperboundaryparamlist)[5]=10000000.0;
	(*upperboundaryparamlist)[6]=10000000.0;
	(*upperboundaryparamlist)[7]=10000000.0;
	class objectiveFuntion : public Optimization_ObjectiveFuntion
	{
	private:
		std::vector<double>* strike_list;
		double forward;
		double maturity;
		int nbsteps;
	public:
		void NAG_CALL operator() (Integer m, Integer n, 
			double x[], /// input
			double f[],	/// output (f(x))
			double fjac[],  /// output  (Df(x,i))
			Integer tdfjac, Nag_Comm *comm)
			
		{
			int i;
			int m_max=(m<strike_list->size())?m:strike_list->size();
			double V0,K,omega,theta,ksi,rho,muJ,sigmaJ,lambda,implvol,der_V0,der_omega,der_theta,der_ksi,der_rho,der_muJ,der_sigmaJ,der_lambda;
			for(i=0;i<m_max;i++)
			{
				K=(*strike_list)[i];
				V0			=x[0];
				omega		=x[1];
				theta		=x[2];
				ksi			=x[3];
				rho			=x[4];
				muJ			=x[5];
				sigmaJ		=x[6];
				lambda		=x[7];
				GeneralizedHeston_DetermineDerivatives(forward, K, V0,  maturity,
							 omega, theta, ksi, rho,  muJ,  sigmaJ,  lambda,
							 &implvol,  &der_V0,&der_omega, &der_theta, &der_ksi, &der_rho,  &der_muJ,  &der_sigmaJ, &der_lambda,
					1,1,1,1,1,1,1,1, nbsteps);
				f[i]		=implvol;
				fjac[i*8]	=der_V0			;
				fjac[i*8+1]	=der_omega		;	
				fjac[i*8+2]	=der_theta		;
				fjac[i*8+3]	=der_ksi		;
				fjac[i*8+4]	=der_rho		;
				fjac[i*8+5]	=der_muJ		;
				fjac[i*8+6]	=der_sigmaJ		;
				fjac[i*8+7]	=der_lambda		;
			}
			
		}
		objectiveFuntion(std::vector<double>* strike_list0,double f,double tex0,int nbsteps0):
		strike_list(strike_list0),forward(f),maturity(tex0),nbsteps(nbsteps0)
		{}
		
	};
	objectiveFuntion func(strikelist,F,tex, nbsteps);
	string tracefile("C:\\Nag_Heston_Calibrate.txt");
	Optimization_Result_Set* result=OptimizeWithDerivatives(
		strikelist,
		pricelist,
		weightlist,
		&func,
		initialparamlist,
		lowerboundaryparamlist,
		upperboundaryparamlist,
		algorithm,
///		FALSE,
		TRUE,
		tracefile
		);
	GeneralizedHeston_ParameterSet* setptr= new GeneralizedHeston_ParameterSet(*result);
	return setptr;

	delete strikelist;
	delete pricelist;
	delete weightlist;
	delete initialparamlist;
	delete lowerboundaryparamlist;
	delete upperboundaryparamlist;
	

}

/// Calibration reduite (Pure heston)
GeneralizedHeston_ParameterSet*  GeneralizedHeston_CalibrateToSmile(std::vector<double>* K_Vec,std::vector<double>* ImpVol_Vec,double F,double tex,double muJ,double sigmaJ,double lambda,int nbsteps,int algorithm,
												double V0_0,double omega_0,double theta_0,double ksi_0,double rho_0)
{
	int i;
	int n=5;
	int m=K_Vec->size();
	if(ImpVol_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"GeneralizedHeston_CalibrateToSmile: K_Vec and ImpVol_Vec do not have the same size!" );
	}
	std::vector<double>* strikelist= new std::vector<double>(m);
	std::vector<double>* pricelist= new std::vector<double>(m);
	std::vector<double>* weightlist= new std::vector<double>(m);
	std::vector<double>* initialparamlist= new std::vector<double>(n);
	std::vector<double>* lowerboundaryparamlist= new std::vector<double>(n);
	std::vector<double>* upperboundaryparamlist= new std::vector<double>(n);
	for(i=0;i<m;i++)
	{
		(*strikelist)[i]=(*K_Vec)[i];
		(*pricelist)[i]=(*ImpVol_Vec)[i];
		(*weightlist)[i]=1.;
	}
	(*initialparamlist)[0]=	V0_0			;
	(*initialparamlist)[1]=	omega_0		;
	(*initialparamlist)[2]=	theta_0		;
	(*initialparamlist)[3]=	ksi_0			;
	(*initialparamlist)[4]=	rho_0			;
	(*lowerboundaryparamlist)[0]=0.0001;
	(*lowerboundaryparamlist)[1]=0.0001;
	(*lowerboundaryparamlist)[2]=0.0001;
	(*lowerboundaryparamlist)[3]=0.00001;
	(*lowerboundaryparamlist)[4]=0.00001;
	
	(*upperboundaryparamlist)[0]=10000000.0;;
	(*upperboundaryparamlist)[1]=10000000.0;
	(*upperboundaryparamlist)[2]=10000000.0;
	(*upperboundaryparamlist)[3]=10000000.0;
	(*upperboundaryparamlist)[4]=0.9999;
	
	class objectiveFuntion : public Optimization_ObjectiveFuntion
	{
	private:
		std::vector<double>* strike_list;
		double forward;
		double maturity;
		double muJ;
		double sigmaJ;
		double lambda;
		int nbsteps;
	public:
		void NAG_CALL operator() (Integer m, Integer n, 
			double x[], /// input
			double f[],	/// output (f(x))
			double fjac[],  /// output  (Df(x,i))
			Integer tdfjac, Nag_Comm *comm)
			
		{
			int i;
			int m_max=(m<strike_list->size())?m:strike_list->size();
			double V0,K,omega,theta,ksi,rho,implvol,der_V0,der_omega,der_theta,der_ksi,der_rho,der_muJ,der_sigmaJ,der_lambda;
			for(i=0;i<m_max;i++)
			{
				K=(*strike_list)[i];
				V0			=x[0];
				omega		=x[1];
				theta		=x[2];
				ksi			=x[3];
				rho			=x[4];
				
				GeneralizedHeston_DetermineDerivatives(forward, K, V0,  maturity,
							 omega, theta, ksi, rho,  muJ,  sigmaJ,  lambda,
							 &implvol,  &der_V0,&der_omega, &der_theta, &der_ksi, &der_rho,  &der_muJ,  &der_sigmaJ, &der_lambda,
					1,1,1,1,1,0,0,0, nbsteps);
				f[i]		=implvol;
				fjac[i*5]	=der_V0			;
				fjac[i*5+1]	=der_omega		;	
				fjac[i*5+2]	=der_theta		;
				fjac[i*5+3]	=der_ksi		;
				fjac[i*5+4]	=der_rho		;
				
			}
			
		}
		objectiveFuntion(std::vector<double>* strike_list0,double f,double tex0,double muJ0,double sigmaJ0,double lambda0,int nbsteps0):
		strike_list(strike_list0),forward(f),maturity(tex0),muJ(muJ0),sigmaJ(sigmaJ0),lambda(lambda0),nbsteps(nbsteps0)
		{}
		
	};
	objectiveFuntion func(strikelist,F,tex,muJ,sigmaJ,lambda, nbsteps);
	string tracefile("C:\\Nag_Heston_Calibrate.txt");
	Optimization_Result_Set* result=OptimizeWithDerivatives(
		strikelist,
		pricelist,
		weightlist,
		&func,
		initialparamlist,
		lowerboundaryparamlist,
		upperboundaryparamlist,
		algorithm,
///		FALSE,
		TRUE,
		tracefile
		);
	GeneralizedHeston_ParameterSet* setptr= new GeneralizedHeston_ParameterSet(*result);
	return setptr;

	delete strikelist;
	delete pricelist;
	delete weightlist;
	delete initialparamlist;
	delete lowerboundaryparamlist;
	delete upperboundaryparamlist;
	

}


void GeneralizedHeston_DetermineDerivatives( double F,double K,double V0, double t,
							double omega,double theta,double ksi,double rho, double muJ, double sigmaJ, double lambda,
							double* impvol, double* der_V0,double* der_omega,double* der_theta,double* der_ksi,double* der_rho, double* der_muJ, double* der_sigmaJ, double* der_lambda,
							bool V0flag,bool omegaflag,bool thetaflag,bool ksiflag,bool rhoflag, bool muJflag, bool sigmaJflag, bool lambdaflag,int nbsteps)

{

	ArgumentList a(F,K,V0,t,omega,theta,ksi,rho,lambda,muJ,sigmaJ,K_CALL,nbsteps);
	
	Power_Expression<ARM_CF_Heston_JumpDiffusion_Formula> y;
(*impvol)=y(a);
if(V0flag) (*der_V0)=y(2,a);
if(omegaflag)(*der_omega)=y(4,a);
if(thetaflag)(*der_theta)=y(5,a);
if(ksiflag)(*der_ksi)=y(6,a);
if(rhoflag)(*der_rho)=y(7,a);
if(muJflag)(*der_muJ)=y(9,a);
if(sigmaJflag)(*der_sigmaJ)=y(10,a);
if(lambdaflag)(*der_lambda)=y(8,a);

}


CC_END_NAMESPACE()


#undef ARM_CF_EPS
#undef ARM_CF_MAXIT

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
