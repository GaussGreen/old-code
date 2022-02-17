/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file glambda_calibration.cpp
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
#include "gpbase/gpvector.h"
#include "gpbase/removenagwarning.h"
#include "nag.h"
#include "nage04.h"

#include "gpnumlib/gaussiananalytics.h"
#include "gpbase/numericconstant.h"

#include "gpclosedforms/sabrvanilla.h"
#include "gpclosedforms/sabrimpliedvol.h"
#include "gpclosedforms/smile_sabr.h"
#include "gpclosedforms/smile_glambda.h"

#include "gpclosedforms/extended_sabr_interface.h"

#include "gpclosedforms/glambda_calibration.h"
#include "gpclosedforms/extendedsabrformula.h"
#include "gpclosedforms/optimization1.h"


#include "expt.h"   // for the exceptions

using namespace std;

CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_TOLERANCE 1.0e-13
#define ARM_CF_MAXIT 2000
#define ARM_GP_CF_SABR_SERIE_USAGE_LIMIT 0.05




GLambda_ParameterSet* GLambda_CalibrateFromSABR(double f0,
												double alpha0,
												double beta0,
												double rho0,
												double nu0, 
												double t0,
												int flag0,
												int nbsteps0,
												double scope,	/// scope=0.2 means we use strikes between f*0.8 and f*1.2
												double initial_l1,
												double initial_l2,
												double initial_l3,
												double initial_l4,
												double initial_l5,
												double initial_l6,
												int algorithm)
{
	int n=6;
	int m=7;
	double xp,xq;
	std::vector<double> Proba_list_(m);
	std::vector<double> Quantile_list_(m);
	int i;

	for(i=0;i<m;i++)
	{
		Quantile_list_[i]=f0*(1.+(2.*i+1.-m)/(m-1.)*scope);
		Proba_list_[i]=1.+Export_SABR_VanillaOption(ARM_CF_SABR_VanillaOption_Formula::STRIKE,
							    f0,Quantile_list_[i],
								t0,alpha0,beta0,rho0,nu0,K_CALL,flag0,nbsteps0);
		xq=Quantile_list_[i];xp=Proba_list_[i];
	}
	std::vector<double> weightlist(m,1);
	std::vector<double> initialparamlist(n);
	std::vector<double> lowerboundaryparamlist(n);
	std::vector<double> upperboundaryparamlist(n);

	initialparamlist[0]=initial_l1;
	initialparamlist[1]=initial_l2;
	initialparamlist[2]=initial_l3;
	initialparamlist[3]=initial_l4;
	initialparamlist[4]=initial_l5;
	initialparamlist[5]=initial_l6;
	lowerboundaryparamlist[0]=0.;
	lowerboundaryparamlist[1]=0;
	lowerboundaryparamlist[2]=0;
	lowerboundaryparamlist[3]=0.09;
	lowerboundaryparamlist[4]=0.01;
	lowerboundaryparamlist[5]=0.05;
	upperboundaryparamlist[0]=100000.;;
	upperboundaryparamlist[1]=100000.;
	upperboundaryparamlist[2]=100000.;
	upperboundaryparamlist[3]=100000.;
	upperboundaryparamlist[4]=100000.;
	upperboundaryparamlist[5]=100000.;
	class objectiveFuntion : public Optimization_ObjectiveFuntion
	{
	private:
		std::vector<double>* Proba_list;
	public:
		void NAG_CALL operator() (Integer m, Integer n, 
			double x[],				/// input
			double f[],				/// output (f(x))
			double fjac[],			/// output  (Df(x,i))
			Integer tdfjac, Nag_Comm *comm)
			
		{
			int i;
			double proba,l1,l2,l3,l4,l5,l6,quantile,der_l1,der_l2,der_l3,der_l4,der_l5,der_l6;
			for(i=0;i<m;i++)
			{
				proba=(*Proba_list)[i];
				l1	=x[0];
				l2	=x[1];
				l3	=x[2];
				l4	=x[3];
				l5	=x[4];
				l6	=x[5];
				GLambda_Smile::QuantileAndAllDerivatives(
					proba,
					l1,
					l2,
					l3,
					l4,
					l5,
					l6,
					&quantile,
					&der_l1,
					&der_l2,
					&der_l3,
					&der_l4,
					&der_l5,
					&der_l6);
				f[i]		=quantile;
				fjac[i*6]	=der_l1;
				fjac[i*6+1]	=der_l2;
				fjac[i*6+2]	=der_l3;
				fjac[i*6+3]	=der_l4;
				fjac[i*6+4]	=der_l5;
				fjac[i*6+5]	=der_l6;
			}
			
		}
		objectiveFuntion(
			std::vector<double>* Proba_list0
			)
		:
		Proba_list(Proba_list0)
		{}
		
	};
	objectiveFuntion func(&Proba_list_);
//	string tracefile("C:\\Documents and Settings\\ocroissant\\My Documents\\Nag_Trace.txt");
	string tracefile("");
	double tolerance=0.00001;
	int max_iter=500;
	Optimization_Result_Set* result=OptimizeWithDerivatives(
		&Proba_list_,
		&Quantile_list_,
		&weightlist,
		&func,
		&initialparamlist,
		&lowerboundaryparamlist,
		&upperboundaryparamlist,
		algorithm,
		FALSE,
		tracefile,
		tolerance,
		max_iter
		);
	GLambda_ParameterSet* setptr= new GLambda_ParameterSet(*result);
	return setptr;

}



CC_END_NAMESPACE()


#undef ARM_CF_EPS
#undef ARM_CF_MAXIT

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
