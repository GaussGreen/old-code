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

#include "gpclosedforms/sabr_calibration.h"
#include "gpclosedforms/extendedsabrformula.h"
#include "gpclosedforms/optimization1.h"
#include "gpclosedforms/sabrbdiff1.h"

#include <glob/expt.h>   // for the exceptions

using namespace std;

CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_TOLERANCE 1.0e-13
#define ARM_CF_MAXIT 2000
#define ARM_GP_CF_SABR_SERIE_USAGE_LIMIT 0.05

///   Calibration of the 4 parameters
SABR_ParameterSet* SABR_CalibrateToSmile(ARM_GP_Vector* K_Vec,
	  ARM_GP_Vector* ImpVol_Vec,ARM_GP_Vector* Weigth_Vec,
	  double f,
	  double tex,
	  int flag,
	  int nbsteps, 
	  int algorithm,
	  double alpha0,
	  double beta0,
	  double rho0,
	  double nu0)
{
	int n=4;
	int m=K_Vec->size();
	if(ImpVol_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"SABR_CalibrateToSmile: K_Vec and ImpVol_Vec do not have the same size!" );
	}
	ARM_GP_Vector weightlist(m,1);
	ARM_GP_Vector initialparamlist(n);
	ARM_GP_Vector lowerboundaryparamlist(n);
	ARM_GP_Vector upperboundaryparamlist(n);

	initialparamlist[0]=alpha0;
	initialparamlist[1]=beta0;
	initialparamlist[2]=rho0;
	initialparamlist[3]=nu0;
	lowerboundaryparamlist[0]=0.001;
	lowerboundaryparamlist[1]=0.5001;
	lowerboundaryparamlist[2]=-0.9999;
	lowerboundaryparamlist[3]=0.001;
	upperboundaryparamlist[0]=1000;;
	upperboundaryparamlist[1]=0.9999;
	upperboundaryparamlist[2]=0.9999;
	upperboundaryparamlist[3]=1000;
	class objectiveFuntion : public Optimization_ObjectiveFuntion
	{
	private:
		ARM_GP_Vector* strike_list;
		double forward;
		double maturity;
		int flag;
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
			double K,alpha,beta,rho,nu,implvol,der_alpha,der_beta,der_rho,der_nu,der_f;
			for(i=0;i<m_max;i++)
			{
				K=(*strike_list)[i];
				alpha=x[0];
				beta=x[1];
				rho=x[2];
				nu=x[3];
				double beta_eff;
				if (beta>0.9999)
				{
					beta_eff=0.9999;
				}
				else
				{
					beta_eff=beta;
				}
				SABR_DetermineAllDerivatives(forward,
					K,
					maturity,
					alpha,
					beta_eff,
					rho,
					nu,
					flag,
					&implvol,
					&der_alpha,
					&der_beta,
					&der_rho,
					&der_nu,
					&der_f,
					1,
					1,
					1,
					1,
					0);
				f[i]		=implvol;
				fjac[i*4]	=der_alpha;
				fjac[i*4+1]	=der_beta;
				fjac[i*4+2]	=der_rho;
				fjac[i*4+3]	=der_nu;
			}
			
		}
		objectiveFuntion(ARM_GP_Vector* strike_list0,
			double f0,
			double tex0,
			int flag0,
			int nbsteps0)
		:
			strike_list(strike_list0),
			forward(f0),
			maturity(tex0),
			nbsteps(nbsteps0),
			flag(flag0)
		{}
		
	};
	objectiveFuntion func(K_Vec,f,tex, flag, nbsteps);
	string tracefile("C:\\Nag_SABR_Calibrate.txt");
	Optimization_Result_Set* result=OptimizeWithDerivatives(
		K_Vec,
		ImpVol_Vec,
		Weigth_Vec,
		&func,
		&initialparamlist,
		&lowerboundaryparamlist,
		&upperboundaryparamlist,
		algorithm,
		FALSE,
		tracefile);
	SABR_ParameterSet* setptr= new SABR_ParameterSet(*result,f);
	return setptr;


}

/// Calibration of the 4+1 parameters (f0 is considered as a parameter)
SABR_ParameterSet* SABR_CalibrateToSmile(ARM_GP_Vector* K_Vec,
		  ARM_GP_Vector* ImpVol_Vec,ARM_GP_Vector* Weigth_Vec,
		  double tex,
		  int flag,
		  int nbsteps,
		  int algorithm,
		  double alpha0,
		  double beta0,
		  double rho0,
		  double nu0,
		  double f0)
{
	int n=5;
	int m=K_Vec->size();
	if(ImpVol_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"SABR_CalibrateToSmile: K_Vec and ImpVol_Vec do not have the same size!" );
	}
	ARM_GP_Vector weightlist(m,1);
	ARM_GP_Vector initialparamlist(n);
	ARM_GP_Vector lowerboundaryparamlist(n);
	ARM_GP_Vector upperboundaryparamlist(n);

	initialparamlist[0]=alpha0;
	initialparamlist[1]=beta0;
	initialparamlist[2]=rho0;
	initialparamlist[3]=nu0;
	initialparamlist[4]=f0;
	lowerboundaryparamlist[0]=0.001;
	lowerboundaryparamlist[1]=0.5001;
	lowerboundaryparamlist[2]=-0.9999;
	lowerboundaryparamlist[3]=0.00001;
	lowerboundaryparamlist[4]=0.00001;
	upperboundaryparamlist[0]=100000000;;
	upperboundaryparamlist[1]=0.9999;
	upperboundaryparamlist[2]=0.9999;
	upperboundaryparamlist[3]=10000000;
	upperboundaryparamlist[4]=10000000;
	class objectiveFuntion : public Optimization_ObjectiveFuntion
	{
	private:
		ARM_GP_Vector* strike_list;
		double maturity;
		int flag;
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
			double K,alpha,beta,rho,nu,forward,implvol,der_alpha,der_beta,der_rho,der_nu,der_f;
			for(i=0;i<m_max;i++)
			{
				K=(*strike_list)[i];
				alpha	=x[0];
				beta	=x[1];
				rho		=x[2];
				nu		=x[3];
				forward	=x[4];
				double beta_eff;
				if (beta>0.9999)
				{
					beta_eff=0.9999;
				}
				else
				{
					beta_eff=beta;
				}
				SABR_DetermineAllDerivatives(forward,
					K,
					maturity,
					alpha,
					beta_eff,
					rho,
					nu,
					flag,
					&implvol,
					&der_alpha,
					&der_beta,
					&der_rho,
					&der_nu,
					&der_f,
					1,
					1,
					1,
					1,
					1);
				f[i]		=implvol;
				fjac[i*5]	=der_alpha;
				fjac[i*5+1]	=der_beta;
				fjac[i*5+2]	=der_rho;
				fjac[i*5+3]	=der_nu;
				fjac[i*5+4]	=der_f;
			}
			
		}
		objectiveFuntion(ARM_GP_Vector* strike_list0,
			double tex0,
			int flag0,
			int nbsteps0)
		:
		strike_list(strike_list0),
		maturity(tex0),
		nbsteps(nbsteps0),
		flag(flag0)
		{}
		
	};
	objectiveFuntion func(K_Vec,tex, flag, nbsteps);
	string tracefile("C:\\Nag_SABR_Calibrate.txt");
	Optimization_Result_Set* result=OptimizeWithDerivatives(
		K_Vec,
		ImpVol_Vec,
		Weigth_Vec,
		&func,
		&initialparamlist,
		&lowerboundaryparamlist,
		&upperboundaryparamlist,
		algorithm,
		FALSE,
		tracefile
		);
	SABR_ParameterSet* setptr= new SABR_ParameterSet(*result);
	return setptr;

}


/// Calibration of 3 parameters : Beta Fixed
SABR_ParameterSet*  SABR_CalibrateToSmile(ARM_GP_Vector* K_Vec,ARM_GP_Vector* ImpVol_Vec,
										  ARM_GP_Vector* Weigth_Vec,
										  double f,double beta,double tex,int flag,int nbsteps,int algorithm,
												double alpha0,double rho0,double nu0)
{
	int i;
	int n=3;
	int m=K_Vec->size();
	if(ImpVol_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"SABR_CalibrateToSmile: K_Vec and ImpVol_Vec do not have the same size!" );
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
	(*lowerboundaryparamlist)[0]=0.00001;
	(*lowerboundaryparamlist)[1]=-0.9999;
	(*lowerboundaryparamlist)[2]=0.00001;
	(*upperboundaryparamlist)[0]=0.5;
	(*upperboundaryparamlist)[1]=0.9999;
	(*upperboundaryparamlist)[2]=1.;

	class objectiveFuntion : public Optimization_ObjectiveFuntion
	{
	private:
		ARM_GP_Vector* strike_list;
		double maturity;
		int flag;
		int nbsteps;
		double forward;
		double beta;
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
				SABR_DetermineAllDerivatives(forward,
                    K,
                    maturity,
                    alpha,
                    beta_eff,
                    rho,
                    nu,
                    flag,
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
		objectiveFuntion(ARM_GP_Vector* strike_list0,double f0,double beta0,double tex0,int flag0,int nbsteps0):
		strike_list(strike_list0),maturity(tex0),nbsteps(nbsteps0),flag(flag0),forward(f0),beta(beta0)
		{}
		
	};
	objectiveFuntion func(strikelist,f,beta,tex, flag, nbsteps);
	string tracefile("C:\\Nag_SABR_Calibrate.txt");

	Optimization_Result_Set* result=OptimizeWithDerivatives(
		strikelist,
		pricelist,
		Weigth_Vec,
		&func,
		initialparamlist,
		lowerboundaryparamlist,
		upperboundaryparamlist,
		algorithm,
//		TRUE,
		FALSE,
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


/// Calibration of 3 parameters (beta fixed) + weigths for the alpha, rhos and nus separatly
SABR_ParameterSet*  SABR_CalibrateToSmile(ARM_GP_Vector* K_Vec,ARM_GP_Vector* ImpVol_Vec,ARM_GP_Vector* Weigth_Vec,
										  double f,double beta,double tex,int flag,int nbsteps,int algorithm,
												double alpha0,double rho0,double nu0,double alphap1,double rhop1,double nup1,
												double rweight_alpha1,double rweight_rho1,double rweight_nu1)
{
	int i;
	int n=3;
	int m=K_Vec->size();
	if(ImpVol_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"SABR_CalibrateToSmile: K_Vec and ImpVol_Vec do not have the same size!" );
	}
	ARM_GP_Vector* strikelist= new ARM_GP_Vector(m+n);
	ARM_GP_Vector* pricelist= new ARM_GP_Vector(m+n);
	ARM_GP_Vector* weightlist= new ARM_GP_Vector(m+n);
	ARM_GP_Vector* initialparamlist= new ARM_GP_Vector(n);
	ARM_GP_Vector* lowerboundaryparamlist= new ARM_GP_Vector(n);
	ARM_GP_Vector* upperboundaryparamlist= new ARM_GP_Vector(n);
	for(i=0;i<m;i++)
	{
		(*strikelist)[i]=(*K_Vec)[i];
		(*pricelist)[i]=(*ImpVol_Vec)[i];
		(*weightlist)[i]=1.;
	}
	for(i=0;i<n;i++)
	{
		(*strikelist)[m+i]=0;
	}
	(*pricelist)[m]=alphap1;
	(*pricelist)[m+1]=rhop1;
	(*pricelist)[m+2]=nup1;

	(*weightlist)[m]=rweight_alpha1;
	(*weightlist)[m+1]=rweight_rho1;
	(*weightlist)[m+2]=rweight_nu1;

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
		double alphap;
		double rhop;
		double nup;
		double rweight_alpha;
		double rweight_rho;
		double rweight_nu;
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
			for(i=0;i<m_max-3;i++)
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
				SABR_DetermineAllDerivatives(forward,
                    K,
                    maturity,
                    alpha,
                    beta_eff,
                    rho,
                    nu,
                    flag,
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
			f[m-3]=x[0];
			f[m-2]=x[1];
			f[m-1]=x[2];
			fjac[(m-3)*3]	=1;
			fjac[(m-3)*3+1]	=0;
			fjac[(m-3)*3+2]	=0;
			fjac[(m-2)*3]	=0;
			fjac[(m-2)*3+1]	=1;
			fjac[(m-2)*3+2]	=0;
			fjac[(m-1)*3]	=0;
			fjac[(m-1)*3+1]	=0;
			fjac[(m-1)*3+2]	=1;
			
		}
		objectiveFuntion(ARM_GP_Vector* strike_list0,double f0,double beta0,double tex0,int flag0,int nbsteps0,
			double alphap0,double rhop0,double nup0,double rweight_alpha0,double rweight_rho0,double rweight_nu0):
		strike_list(strike_list0),maturity(tex0),nbsteps(nbsteps0),flag(flag0),forward(f0),beta(beta0),alphap(alphap0),
			rhop(rhop0),nup(nup0),rweight_alpha(rweight_alpha0),rweight_rho(rweight_rho0),rweight_nu(rweight_nu0)
		{}
		
	};
	objectiveFuntion func(strikelist,f,beta,tex, flag, nbsteps,alphap1,rhop1,nup1,rweight_alpha1,rweight_rho1,rweight_nu1);

	Optimization_Result_Set* result=OptimizeWithDerivatives(
		strikelist,
		pricelist,
		/*Weigth_Vec,*/weightlist,
		&func,
		initialparamlist,
		lowerboundaryparamlist,
		upperboundaryparamlist,
		algorithm,
		FALSE,
		""
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

SABR_ParameterSet * SABR_CalibrateToSmileBetaFixedToOne(ARM_GP_Vector* K_Vec,ARM_GP_Vector* ImpVol_Vec,
														ARM_GP_Vector* Weight_Vec,double f,double t,double atmvol)
{
	class ARM_SABR_CalibrationBetaOne : public ARM_LEVMARQFunc
	{
	private:
		bool calibatm;
		ARM_GP_Vector *	KVec;
		ARM_GP_Vector *	ImpVolVec;
		ARM_GP_Vector *	WeightVec;
		double			f;
		double			t;
		double			atmvol;

		double *		alpha;
		double *		rho;
		double *		nu;

	public:
		bool	calibATM() const {return calibatm;};
		double	Alpha() const {return *alpha;};
		double	Rho() const {return *rho;};
		double	Nu() const {return *nu;};

		ARM_SABR_CalibrationBetaOne(ARM_GP_Vector * K_Vec, ARM_GP_Vector * ImpVol_Vec,
					ARM_GP_Vector * Weight_Vec, double fwd, double ttm, double atm)
		{
			KVec		= K_Vec;
			ImpVolVec	= ImpVol_Vec;
			WeightVec	= Weight_Vec;
			f			= fwd;
			t			= ttm;
			atmvol		= atm;
			calibatm	= atmvol > K_DOUBLE_TOL;
			alpha		= new double;
			rho			= new double;
			nu			= new double;
		}

		~ARM_SABR_CalibrationBetaOne()
		{
			delete alpha;
			delete rho;
			delete nu;
		}

		void operator()(double p[], double hx[], int m, int n, void * adata = NULL) const
		{
			(*rho) = p[0] > 0.999 ? 0.999 : p[0] < -0.999 ? -0.999 : p[0];
			(*nu) = p[1] < 1e-8 ? 1e-8 : p[1];
			(*alpha) = calibatm ? GetAlpha() : p[2] < 1e-8 ? 1e-8 : p[2];
			
			for(int k = 0; k < n; k++)
			{
				hx[k] = (*ImpVolVec)[k] - vol((*KVec)[k]);
				
				if(WeightVec != NULL)
					if(WeightVec->size() == KVec->size()) hx[k] *= (*WeightVec)[k];
			}
		}

	private:

		double	GetAlpha() const
		{
			double a, b, c, d;

			a	= 0.25 * t * (*rho) * (*nu);

			b	= 1. + t * (2. - 3. * (*rho) * (*rho)) * (*nu) * (*nu) / 24.;

			c	= - atmvol;

			d	= b * b - 4. * a * c;

			if(d < 0. || fabs(a) < 1e-12)
			{
				return atmvol;
			}
			else
			{
				double rac1 = (- b + sqrt(d)) / (2. * a);
				double rac2 = (- b - sqrt(d)) / (2. * a);

				return (rac1 > 0. ? rac1 : rac2 > 0. ? rac2 : atmvol);
			}
		}

		double vol(double strike) const
		{
			return CptSABR_BetEqOne_ImplVol(f,strike,t,*alpha,*rho,*nu);
		}
	};

	ARM_SABR_CalibrationBetaOne func(K_Vec, ImpVol_Vec, Weight_Vec, f, t, atmvol);

	ARM_GP_Vector x(2 + (func.calibATM() == false));
	x[0] = -0.2;
	x[1] = 0.3;
	if(func.calibATM() == false) x[2] = 0.15;
	ARM_GP_Vector fx(K_Vec->size(), 0.);

	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];

	opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
	opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference jacobian version is used 

	int status = LEVMARQMinization_WithNumDerivatives(func, x, fx, info, 500, opts);

	SABR_ParameterSet * set = new SABR_ParameterSet(f,func.Alpha(),1.,func.Rho(),func.Nu(),info[2]);

	return set;
}

CC_END_NAMESPACE()


#undef ARM_CF_EPS
#undef ARM_CF_MAXIT

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
