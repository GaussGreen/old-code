/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 03/15/2006
 *
 *  basic functions for the closed form framework 
 *
 *	\file bisabr_calibration.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date March 2006
 */
#include <glob/firsttoinc.h>
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

#include "gpclosedforms/bisabr_spreadoption.h"
#include "gpclosedforms/bisabr_spreadoption_formula.h"

#include "gpclosedforms/bisabr_calibration.h"
#include "gpclosedforms/optimization1.h"
#include "gpclosedforms/eigenvalues.h"


#include <glob/expt.h>   // for the exceptions

using namespace std;

CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_TOLERANCE 1.0e-13
#define ARM_CF_MAXIT 150


///        Creation of a set of correlations acceptables (positve definite correlation matrix )
///        from rho1, rho2 and rhos imposed, other correlations just modified at the margin
/// We start from rho1,rho2,rhos,rhov,rhoc12,rhoc21 which are acceptable by hypothesis
/// we return newrhov,newrhoc12, newrhoc21

void BiSABR_CorrelationEvolution(double rho1,double rho2,double rhos,double rhov,double rhoc12,double rhoc21,
								 double newrho1, double newrho2,
								 double* newrhov,double* newrhoc12,double* newrhoc21)

{
	double r=sqrt(1.-rhoc12*rhoc12);
	double l3=(rho2-rhoc12*rhos)*r;
	double f=sqrt(1.-rhos*rhos-l3*l3);
	double l5=(rhov-rhoc12*rho1)/r;
	double l6=(rhoc21-rhos*rho1-l3*l5)/f;
	double d1=sqrt((1.-newrho1*newrho1)/(1.-rho1*rho1));
	double d2=sqrt((1.-newrho2*newrho2)/(1.-rho2*rho2));
	*newrhov = l5*d1*(l3*(newrho2-rho2*d2)+d2*r)+newrho1*(d2*rhoc12+(newrho2-rho2*d2)*rhos)+l6*d1*(newrho2-rho2*d2)*f;
	*newrhoc12 = d2*rhoc12+(newrho2-rho2*d2)*rhos;
	*newrhoc21 = l3*l5*d1+newrho1*rhos+l6*d1*f;
	return;
}


BiSABR_ParameterSet*  BiSABR_CalibrateToSmile(
			ARM_GP_Vector* F1_Vec,ARM_GP_Vector* Alpha1_Vec,ARM_GP_Vector* Beta1_Vec,ARM_GP_Vector* Rho1_Vec,ARM_GP_Vector* Nu1_Vec,
			ARM_GP_Vector* F2_Vec,ARM_GP_Vector* Alpha2_Vec,ARM_GP_Vector* Beta2_Vec,ARM_GP_Vector* Rho2_Vec,ARM_GP_Vector* Nu2_Vec,
			ARM_GP_Vector* K_Vec,ARM_GP_Vector* T_Vec,ARM_GP_Vector* Price_Vec,ARM_GP_Vector* Weight_Vec,
						double rhos_0,double rhov_0,double rhoc12_0,double rhoc21_0,
						double converg_prec, int nbIter_max, double first_step_max,int flag)
{
	int i,j;
	int n=4; /// number of calibrated parameters
	int m;
	m=F1_Vec->size();
	if(Alpha1_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: alpha1_Vec and F1_Vec do not have the same size!" );
	}
	if(Beta1_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: beta1_Vec and F1_Vec do not have the same size!" );
	}
	if(Rho1_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: rho1_Vec and F1_Vec do not have the same size!" );
	}
	if(Nu1_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: nu1_Vec and F1_Vec do not have the same size!" );
	}
	if(Alpha2_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: alpha2_Vec and F1_Vec do not have the same size!" );
	}
	if(Beta2_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: beta2_Vec and F1_Vec do not have the same size!" );
	}
	if(Rho2_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: rho2_Vec and F1_Vec do not have the same size!" );
	}
	if(Nu2_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: nu2_Vec and F1_Vec do not have the same size!" );
	}
	if(K_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: K_Vec and F1_Vec do not have the same size!" );
	}
	if(T_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: T_Vec and F1_Vec do not have the same size!" );
	}
	if(Price_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: T_Vec and F1_Vec do not have the same size!" );
	}
		if(Weight_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: Weight_Vec and F1_Vec do not have the same size!" );
	}


	ARM_GP_Vector* F1list= new ARM_GP_Vector(m);
	ARM_GP_Vector* Alpha1list= new ARM_GP_Vector(m);
	ARM_GP_Vector* Beta1list= new ARM_GP_Vector(m);
	ARM_GP_Vector* Rho1list= new ARM_GP_Vector(m);
	ARM_GP_Vector* Nu1list= new ARM_GP_Vector(m);
	ARM_GP_Vector* F2list= new ARM_GP_Vector(m);
	ARM_GP_Vector* Alpha2list= new ARM_GP_Vector(m);
	ARM_GP_Vector* Beta2list= new ARM_GP_Vector(m);
	ARM_GP_Vector* Rho2list= new ARM_GP_Vector(m);
	ARM_GP_Vector* Nu2list= new ARM_GP_Vector(m);
	ARM_GP_Vector* Klist= new ARM_GP_Vector(m);
	ARM_GP_Vector* Tlist= new ARM_GP_Vector(m);
	ARM_GP_Vector* pricelist= new ARM_GP_Vector(m);
	ARM_GP_Vector* weightlist= new ARM_GP_Vector(m);
	ARM_GP_Vector* initialparamlist= new ARM_GP_Vector(n);
	ARM_GP_Vector* lowerboundaryparamlist= new ARM_GP_Vector(n);
	ARM_GP_Vector* upperboundaryparamlist= new ARM_GP_Vector(n);
	ARM_GP_Vector* constraintlowerboundaryparamlist= new ARM_GP_Vector(4*m);
	ARM_GP_Vector* constraintupperboundaryparamlist= new ARM_GP_Vector(4*m);
	for(i=0;i<m;i++)
	{
		(*F1list)[i]=(*F1_Vec)[i];
		(*Alpha1list)[i]=(*Alpha1_Vec)[i];
		(*Beta1list)[i]=(*Beta1_Vec)[i];
		(*Rho1list)[i]=(*Rho1_Vec)[i];
		(*Nu1list)[i]=(*Nu1_Vec)[i];
		(*F2list)[i]=(*F2_Vec)[i];
		(*Alpha2list)[i]=(*Alpha2_Vec)[i];
		(*Beta2list)[i]=(*Beta2_Vec)[i];
		(*Rho2list)[i]=(*Rho2_Vec)[i];
		(*Nu2list)[i]=(*Nu2_Vec)[i];
		(*Klist)[i]=(*K_Vec)[i];
		(*Tlist)[i]=(*T_Vec)[i];
		(*pricelist)[i]=(*Price_Vec)[i];
		(*weightlist)[i]=(*Weight_Vec)[i];
	}
	for(i=0;i<4;i++)
	{
		(*lowerboundaryparamlist)[i]=-0.9999;
		(*upperboundaryparamlist)[i]=0.9999;
		for(j=0;j<m;j++)
		{
			(*constraintlowerboundaryparamlist)[j*4+i]=0.00000000001;			/// eigenvalue lower limit 
			(*constraintupperboundaryparamlist)[j*4+i]=1000000000000;			/// eigenvalue upper limit
		}
	}
	(*initialparamlist)[0]=	rhos_0	;
	(*initialparamlist)[1]=	rhov_0	;
	(*initialparamlist)[2]=	rhoc12_0;
	(*initialparamlist)[3]=	rhoc21_0;


	class objectiveFuntion : public Optimization_ObjectiveFuntion
	{
	private:
		ARM_GP_Vector* F1_list;
		ARM_GP_Vector* Alpha1_list;
		ARM_GP_Vector* Beta1_list;
		ARM_GP_Vector* Rho1_list;
		ARM_GP_Vector* Nu1_list;
		ARM_GP_Vector* F2_list;
		ARM_GP_Vector* Alpha2_list;
		ARM_GP_Vector* Beta2_list;
		ARM_GP_Vector* Rho2_list;
		ARM_GP_Vector* Nu2_list;
		ARM_GP_Vector* K_list;
		ARM_GP_Vector* T_list;
		int Flag;
	
	public:
		void NAG_CALL operator() (Integer m, Integer n, 
			double x[], /// input
			double f[],	/// output (f(x))
			double fjac[],  /// output  (Df(x,i))
			Integer tdfjac, Nag_Comm *comm)
			
		{
			int i;
			int m_max=(m<K_list->size())?m:K_list->size();
			double price,F1,Alpha1,Beta1,Rho1,Nu1,F2,Alpha2,Beta2,Rho2,Nu2,K,T,Rhos,Rhov,Rhoc12,Rhoc21,der_Rhos,der_Rhov,der_Rhoc12,der_Rhoc21;
			for(i=0;i<m_max;i++)
			{
				F1=(*F1_list)[i];
				Alpha1=(*Alpha1_list)[i];
				Beta1=(*Beta1_list)[i];
				Rho1=(*Rho1_list)[i];
				Nu1=(*Nu1_list)[i];
				F2=(*F2_list)[i];
				Alpha2=(*Alpha2_list)[i];
				Beta2=(*Beta2_list)[i];
				Rho2=(*Rho2_list)[i];
				Nu2=(*Nu2_list)[i];
				K=(*K_list)[i];
				T=(*T_list)[i];
				Rhos		=x[0];
				Rhov		=x[1];
				Rhoc12		=x[2];
				Rhoc21		=x[3];
			
				BiSABR_DetermineDerivatives(F1, Alpha1,Beta1,Rho1,Nu1,F2, Alpha2,Beta2,Rho2,Nu2,K,T,Rhos,Rhov,Rhoc12,Rhoc21,
							 &price,  &der_Rhos,&der_Rhov, &der_Rhoc12, &der_Rhoc21,1,1,1,1,Flag);
				f[i]		=price;
				fjac[i*4]	=der_Rhos			;
				fjac[i*4+1]	=der_Rhov		;	
				fjac[i*4+2]	=der_Rhoc12		;
				fjac[i*4+3]	=der_Rhoc21		;
			}
			
		}
		objectiveFuntion(
			ARM_GP_Vector* F1_list0,ARM_GP_Vector* Alpha1_list0,ARM_GP_Vector* Beta1_list0,ARM_GP_Vector* Rho1_list0,ARM_GP_Vector* Nu1_list0,
			ARM_GP_Vector* F2_list0,ARM_GP_Vector* Alpha2_list0,ARM_GP_Vector* Beta2_list0,ARM_GP_Vector* Rho2_list0,ARM_GP_Vector* Nu2_list0,
			ARM_GP_Vector* K_list0,ARM_GP_Vector* T_list0,int flag0):
		F1_list(F1_list0),Alpha1_list(Alpha1_list0),Beta1_list(Beta1_list0),Rho1_list(Rho1_list0),Nu1_list(Nu1_list0),
		F2_list(F2_list0),Alpha2_list(Alpha2_list0),Beta2_list(Beta2_list0),Rho2_list(Rho2_list0),Nu2_list(Nu2_list0),
		K_list(K_list0),T_list(T_list0),Flag(flag0)
		{}
		
	};
	class constraintFuntion : public Optimization_ConstraintFuntion
	{
	private:
		
		ARM_GP_Vector* Rho1_list;
		ARM_GP_Vector* Rho2_list;
	
	
	public:
		void NAG_CALL operator() (Integer ncnl, Integer n, 
			double x[], /// input
			double conf[],	/// output (f(x))
			double conjac[],  /// output  (Df(x,i))
		 Nag_Comm *comm)
			
		{
			int i;
// FIXMEFRED: mig.vc8 (22/05/2007 18:26:22):floor (Integer/Integer) doesnt mean anything
			int m=ncnl/n;
			double Rho1,Rho2,Rhos,Rhov,Rhoc12,Rhoc21,
				e1_der_rho12,e2_der_rho12,e3_der_rho12,e4_der_rho12,
				e1_der_rho13,e2_der_rho13,e3_der_rho13,e4_der_rho13,
				e1_der_rho14,e2_der_rho14,e3_der_rho14,e4_der_rho14,
				e1_der_rho23,e2_der_rho23,e3_der_rho23,e4_der_rho23,
				e1_der_rho24,e2_der_rho24,e3_der_rho24,e4_der_rho24,
				e1_der_rho34,e2_der_rho34,e3_der_rho34,e4_der_rho34;
			
				double e1,e2,e3,e4;
			
				for(i=0;i<m;i++)
				{
					Rho1=(*Rho1_list)[i];
					Rho2=(*Rho2_list)[i];
					Rhos		=x[0];
					Rhov		=x[1];
					Rhoc12		=x[2];
					Rhoc21		=x[3];
					SymmetricPolynomOfEigenValues4_AllDerivatives(Rho1,Rhos, Rhoc12,Rhoc21,Rhov,Rho2,
						&e1,&e2,&e3,&e4,
						&e1_der_rho12,&e2_der_rho12,&e3_der_rho12,&e4_der_rho12,
						&e1_der_rho13,&e2_der_rho13,&e3_der_rho13,&e4_der_rho13,
						&e1_der_rho14,&e2_der_rho14,&e3_der_rho14,&e4_der_rho14,
						&e1_der_rho23,&e2_der_rho23,&e3_der_rho23,&e4_der_rho23,
						&e1_der_rho24,&e2_der_rho24,&e3_der_rho24,&e4_der_rho24,
						&e1_der_rho34,&e2_der_rho34,&e3_der_rho34,&e4_der_rho34,
						0,1,1,1,1,0);
					conf[(i*4)]			=e1;
					conf[(i*4+1)]		=e2;
					conf[(i*4+2)]		=e3;
					conf[(i*4+3)]		=e4;
					
					conjac[i*16]	=e1_der_rho13;
					conjac[i*16+1]	=e1_der_rho24;
					conjac[i*16+2]	=e1_der_rho14;
					conjac[i*16+3]	=e1_der_rho23;
					conjac[i*16+4]	=e2_der_rho13;
					conjac[i*16+5]	=e2_der_rho24;
					conjac[i*16+6]	=e2_der_rho14;
					conjac[i*16+7]	=e2_der_rho23;
					conjac[i*16+8]	=e3_der_rho13;
					conjac[i*16+9]	=e3_der_rho24;
					conjac[i*16+10]	=e3_der_rho14;
					conjac[i*16+11]	=e3_der_rho23;
					conjac[i*16+12]	=e4_der_rho13;
					conjac[i*16+13]	=e4_der_rho24;
					conjac[i*16+14]	=e4_der_rho14;
					conjac[i*16+15]	=e4_der_rho23;
				}
				
		}
		constraintFuntion(
			ARM_GP_Vector* Rho1_list0,ARM_GP_Vector* Rho2_list0):
		Rho1_list(Rho1_list0),Rho2_list(Rho2_list0)
		{}
		
	};
	objectiveFuntion func(F1list,Alpha1list,Beta1list,Rho1list,Nu1list,F2list,Alpha2list,Beta2list,Rho2list,Nu2list,Klist,Tlist,flag);

	constraintFuntion confunc(Rho1list,Rho2list);

	Optimization_Result_Set* result=OptimizeWithDerivatives_Constraint(
		Klist,
		pricelist,
		weightlist,
		&func,
		&confunc,
		initialparamlist,
		lowerboundaryparamlist,
		upperboundaryparamlist,
		constraintlowerboundaryparamlist,
		constraintupperboundaryparamlist,
		Optimization_ObjectiveFuntion::NAG_OPT_NLIN_LSQ_CONSTRAINT,
		FALSE,
		"C:\\NagTrace.txt",converg_prec,nbIter_max,first_step_max);

	BiSABR_ParameterSet* setptr= new BiSABR_ParameterSet(*result);
	

	delete F1list;
	delete Alpha1list;
	delete Beta1list;
	delete Rho1list;
	delete Nu1list;
	delete F2list;
	delete Alpha2list;
	delete Beta2list;
	delete Rho2list;
	delete Nu2list;
	delete Klist;
	delete Tlist;
	delete pricelist;
	delete weightlist;
	delete initialparamlist;
	delete lowerboundaryparamlist;
	delete upperboundaryparamlist;

	return setptr;
	

}

BiSABR_ParameterSet*  BiSABR_CalibrateToSmile_ConstantRhos(
			ARM_GP_Vector* F1_Vec,ARM_GP_Vector* Alpha1_Vec,ARM_GP_Vector* Beta1_Vec,ARM_GP_Vector* Rho1_Vec,ARM_GP_Vector* Nu1_Vec,
			ARM_GP_Vector* F2_Vec,ARM_GP_Vector* Alpha2_Vec,ARM_GP_Vector* Beta2_Vec,ARM_GP_Vector* Rho2_Vec,ARM_GP_Vector* Nu2_Vec,
			ARM_GP_Vector* K_Vec,ARM_GP_Vector* T_Vec,ARM_GP_Vector* Price_Vec,ARM_GP_Vector* Weight_Vec,
						double rhos,double rhov_0,double rhoc12_0,double rhoc21_0,
						double converg_prec, int nbIter_max, double first_step_max,int flag)
{
	int i,j;
	int n=3; /// number of calibrated parameters
	int m;
	m=F1_Vec->size();
	if(Alpha1_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: alpha1_Vec and F1_Vec do not have the same size!" );
	}
	if(Beta1_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: beta1_Vec and F1_Vec do not have the same size!" );
	}
	if(Rho1_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: rho1_Vec and F1_Vec do not have the same size!" );
	}
	if(Nu1_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: nu1_Vec and F1_Vec do not have the same size!" );
	}
	if(Alpha2_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: alpha2_Vec and F1_Vec do not have the same size!" );
	}
	if(Beta2_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: beta2_Vec and F1_Vec do not have the same size!" );
	}
	if(Rho2_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: rho2_Vec and F1_Vec do not have the same size!" );
	}
	if(Nu2_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: nu2_Vec and F1_Vec do not have the same size!" );
	}
	if(K_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: K_Vec and F1_Vec do not have the same size!" );
	}
	if(T_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: T_Vec and F1_Vec do not have the same size!" );
	}
	if(Price_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: T_Vec and F1_Vec do not have the same size!" );
	}
		if(Weight_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"BiSABR_CalibrateToSmile: Weight_Vec and F1_Vec do not have the same size!" );
	}


	ARM_GP_Vector* F1list= new ARM_GP_Vector(m);
	ARM_GP_Vector* Alpha1list= new ARM_GP_Vector(m);
	ARM_GP_Vector* Beta1list= new ARM_GP_Vector(m);
	ARM_GP_Vector* Rho1list= new ARM_GP_Vector(m);
	ARM_GP_Vector* Nu1list= new ARM_GP_Vector(m);
	ARM_GP_Vector* F2list= new ARM_GP_Vector(m);
	ARM_GP_Vector* Alpha2list= new ARM_GP_Vector(m);
	ARM_GP_Vector* Beta2list= new ARM_GP_Vector(m);
	ARM_GP_Vector* Rho2list= new ARM_GP_Vector(m);
	ARM_GP_Vector* Nu2list= new ARM_GP_Vector(m);
	ARM_GP_Vector* Klist= new ARM_GP_Vector(m);
	ARM_GP_Vector* Tlist= new ARM_GP_Vector(m);
	ARM_GP_Vector* pricelist= new ARM_GP_Vector(m);
	ARM_GP_Vector* weightlist= new ARM_GP_Vector(m);
	ARM_GP_Vector* initialparamlist= new ARM_GP_Vector(n);
	ARM_GP_Vector* lowerboundaryparamlist= new ARM_GP_Vector(n);
	ARM_GP_Vector* upperboundaryparamlist= new ARM_GP_Vector(n);
	ARM_GP_Vector* constraintlowerboundaryparamlist= new ARM_GP_Vector(4*m);
	ARM_GP_Vector* constraintupperboundaryparamlist= new ARM_GP_Vector(4*m);
	for(i=0;i<m;i++)
	{
		(*F1list)[i]=(*F1_Vec)[i];
		(*Alpha1list)[i]=(*Alpha1_Vec)[i];
		(*Beta1list)[i]=(*Beta1_Vec)[i];
		(*Rho1list)[i]=(*Rho1_Vec)[i];
		(*Nu1list)[i]=(*Nu1_Vec)[i];
		(*F2list)[i]=(*F2_Vec)[i];
		(*Alpha2list)[i]=(*Alpha2_Vec)[i];
		(*Beta2list)[i]=(*Beta2_Vec)[i];
		(*Rho2list)[i]=(*Rho2_Vec)[i];
		(*Nu2list)[i]=(*Nu2_Vec)[i];
		(*Klist)[i]=(*K_Vec)[i];
		(*Tlist)[i]=(*T_Vec)[i];
		(*pricelist)[i]=(*Price_Vec)[i];
		(*weightlist)[i]=(*Weight_Vec)[i];
	}
	for(i=0;i<3;i++)
	{
		(*lowerboundaryparamlist)[i]=-0.9999;
		(*upperboundaryparamlist)[i]=0.9999;
	}
	for(j=0;j<m;j++)
	{
		for(i=0;i<4;i++)
		{
			(*constraintlowerboundaryparamlist)[j*4+i]=0.00000000001;			/// eigenvalue lower limit 
			(*constraintupperboundaryparamlist)[j*4+i]=1000000000000;			/// eigenvalue upper limit
		}
	}
	(*initialparamlist)[0]=	rhov_0	;
	(*initialparamlist)[1]=	rhoc12_0;
	(*initialparamlist)[2]=	rhoc21_0;
	
	
	class objectiveFuntion : public Optimization_ObjectiveFuntion
	{
	private:
		ARM_GP_Vector* F1_list;
		ARM_GP_Vector* Alpha1_list;
		ARM_GP_Vector* Beta1_list;
		ARM_GP_Vector* Rho1_list;
		ARM_GP_Vector* Nu1_list;
		ARM_GP_Vector* F2_list;
		ARM_GP_Vector* Alpha2_list;
		ARM_GP_Vector* Beta2_list;
		ARM_GP_Vector* Rho2_list;
		ARM_GP_Vector* Nu2_list;
		ARM_GP_Vector* K_list;
		ARM_GP_Vector* T_list;
		double rhos;
		int Flag;
	
	public:
		void NAG_CALL operator() (Integer m, Integer n, 
			double x[], /// input
			double f[],	/// output (f(x))
			double fjac[],  /// output  (Df(x,i))
			Integer tdfjac, Nag_Comm *comm)
			
		{
			int i;
			int m_max=(m<K_list->size())?m:K_list->size();
			double price,F1,Alpha1,Beta1,Rho1,Nu1,F2,Alpha2,Beta2,Rho2,Nu2,K,T,Rhos,Rhov,Rhoc12,Rhoc21,der_Rhos,der_Rhov,der_Rhoc12,der_Rhoc21;
			for(i=0;i<m_max;i++)
			{
				F1=(*F1_list)[i];
				Alpha1=(*Alpha1_list)[i];
				Beta1=(*Beta1_list)[i];
				Rho1=(*Rho1_list)[i];
				Nu1=(*Nu1_list)[i];
				F2=(*F2_list)[i];
				Alpha2=(*Alpha2_list)[i];
				Beta2=(*Beta2_list)[i];
				Rho2=(*Rho2_list)[i];
				Nu2=(*Nu2_list)[i];
				K=(*K_list)[i];
				T=(*T_list)[i];
				Rhos		=rhos;
				Rhov		=x[0];
				Rhoc12		=x[1];
				Rhoc21		=x[2];
			
				BiSABR_DetermineDerivatives(F1, Alpha1,Beta1,Rho1,Nu1,F2, Alpha2,Beta2,Rho2,Nu2,K,T,Rhos,Rhov,Rhoc12,Rhoc21,
							 &price,  &der_Rhos,&der_Rhov, &der_Rhoc12, &der_Rhoc21,1,1,1,1,Flag);
				f[i]		=price;
				fjac[i*3]	=der_Rhov			;
				fjac[i*3+1]	=der_Rhoc12		;	
				fjac[i*3+2]	=der_Rhoc21		;
					;
			}
			
		}
		objectiveFuntion(
			ARM_GP_Vector* F1_list0,ARM_GP_Vector* Alpha1_list0,ARM_GP_Vector* Beta1_list0,ARM_GP_Vector* Rho1_list0,ARM_GP_Vector* Nu1_list0,
			ARM_GP_Vector* F2_list0,ARM_GP_Vector* Alpha2_list0,ARM_GP_Vector* Beta2_list0,ARM_GP_Vector* Rho2_list0,ARM_GP_Vector* Nu2_list0,
			double rhos0,ARM_GP_Vector* K_list0,ARM_GP_Vector* T_list0,int flag0):
		F1_list(F1_list0),Alpha1_list(Alpha1_list0),Beta1_list(Beta1_list0),Rho1_list(Rho1_list0),Nu1_list(Nu1_list0),
		F2_list(F2_list0),Alpha2_list(Alpha2_list0),Beta2_list(Beta2_list0),Rho2_list(Rho2_list0),Nu2_list(Nu2_list0),
		K_list(K_list0),T_list(T_list0),rhos(rhos0),Flag(flag0)
		{}
		
	};
	class constraintFuntion : public Optimization_ConstraintFuntion
	{
	private:
		
		ARM_GP_Vector* Rho1_list;
		ARM_GP_Vector* Rho2_list;
		double rhos;
	
	
	public:
		void NAG_CALL operator() (Integer ncnl, Integer n, 
			double x[], /// input
			double conf[],	/// output (f(x))
			double conjac[],  /// output  (Df(x,i))
		 Nag_Comm *comm)
			
		{
			int i;
// FIXMEFRED: mig.vc8 (22/05/2007 18:45:57):floor(int/int)
			int m=ncnl/4;
			double Rho1,Rho2,Rhos,Rhov,Rhoc12,Rhoc21,
				e1_der_rho12,e2_der_rho12,e3_der_rho12,e4_der_rho12,
				e1_der_rho13,e2_der_rho13,e3_der_rho13,e4_der_rho13,
				e1_der_rho14,e2_der_rho14,e3_der_rho14,e4_der_rho14,
				e1_der_rho23,e2_der_rho23,e3_der_rho23,e4_der_rho23,
				e1_der_rho24,e2_der_rho24,e3_der_rho24,e4_der_rho24,
				e1_der_rho34,e2_der_rho34,e3_der_rho34,e4_der_rho34;
			
				double e1,e2,e3,e4;
			
				for(i=0;i<m;i++)
				{
					Rho1=(*Rho1_list)[i];
					Rho2=(*Rho2_list)[i];
					Rhos		=rhos;
					Rhov		=x[0];
					Rhoc12		=x[1];
					Rhoc21		=x[2];
					SymmetricPolynomOfEigenValues4_AllDerivatives(Rho1,Rhos, Rhoc12,Rhoc21,Rhov,Rho2,
						&e1,&e2,&e3,&e4,
						&e1_der_rho12,&e2_der_rho12,&e3_der_rho12,&e4_der_rho12,
						&e1_der_rho13,&e2_der_rho13,&e3_der_rho13,&e4_der_rho13,
						&e1_der_rho14,&e2_der_rho14,&e3_der_rho14,&e4_der_rho14,
						&e1_der_rho23,&e2_der_rho23,&e3_der_rho23,&e4_der_rho23,
						&e1_der_rho24,&e2_der_rho24,&e3_der_rho24,&e4_der_rho24,
						&e1_der_rho34,&e2_der_rho34,&e3_der_rho34,&e4_der_rho34,
						0,0,1,1,1,0);
					conf[(i*4)]			=e1;
					conf[(i*4+1)]		=e2;
					conf[(i*4+2)]		=e3;
					conf[(i*4+3)]		=e4;
					
					conjac[i*12]	=e1_der_rho24;			/// contrainte 1 : e1
					conjac[i*12+1]	=e1_der_rho14;
					conjac[i*12+2]	=e1_der_rho23;
					conjac[i*12+3]	=e2_der_rho24;
					conjac[i*12+4]	=e2_der_rho14;
					conjac[i*12+5]	=e2_der_rho23;
					conjac[i*12+6]	=e3_der_rho24;
					conjac[i*12+7]	=e3_der_rho14;
					conjac[i*12+8]	=e3_der_rho23;
					conjac[i*12+9]	=e4_der_rho24;
					conjac[i*12+10]	=e4_der_rho14;
					conjac[i*12+11]	=e4_der_rho23;
				}
				
		}
		constraintFuntion(
			ARM_GP_Vector* Rho1_list0,ARM_GP_Vector* Rho2_list0,double rhos0):
		Rho1_list(Rho1_list0),Rho2_list(Rho2_list0),rhos(rhos0)
		{}
		
	};
	objectiveFuntion func(F1list,Alpha1list,Beta1list,Rho1list,Nu1list,F2list,Alpha2list,Beta2list,Rho2list,Nu2list,rhos,Klist,Tlist,flag);

	constraintFuntion confunc(Rho1list,Rho2list,rhos);

	Optimization_Result_Set* result=OptimizeWithDerivatives_Constraint(
		Klist,
		pricelist,
		weightlist,
		&func,
		&confunc,
		initialparamlist,
		lowerboundaryparamlist,
		upperboundaryparamlist,
		constraintlowerboundaryparamlist,
		constraintupperboundaryparamlist,
		Optimization_ObjectiveFuntion::NAG_OPT_NLIN_LSQ_CONSTRAINT,
		FALSE,
		"C:\\NagTrace.txt",converg_prec,nbIter_max,first_step_max);

	BiSABR_ParameterSet* setptr= new BiSABR_ParameterSet(*result,rhos);
	

	delete F1list;
	delete Alpha1list;
	delete Beta1list;
	delete Rho1list;
	delete Nu1list;
	delete F2list;
	delete Alpha2list;
	delete Beta2list;
	delete Rho2list;
	delete Nu2list;
	delete Klist;
	delete Tlist;
	delete pricelist;
	delete weightlist;
	delete initialparamlist;
	delete lowerboundaryparamlist;
	delete upperboundaryparamlist;

	return setptr;
	

}



void BiSABR_DetermineDerivatives( double F1,double alpha1,double beta1,double rho1,double nu1,
								 double F2,double alpha2,double beta2,double rho2,double nu2,
								double K,double T,
								double rhos,double rhov,double rhoc12,double rhoc21,
							double* price, double* der_rhos,double* der_rhov,double* der_rhoc12,double* der_rhoc21,
							bool rhos_flag,bool rhov_flag,bool rhoc12_flag,bool rhoc21_flag,int flag)
{

	ArgumentList a(F1,alpha1,beta1,rho1,nu1,
		F2,alpha2,beta2,rho2,nu2,
		rhos,rhov,rhoc12,rhoc21,
		K,T,K_CALL,flag);
	
	Power_Expression<ARM_CF_BiSABR_SpreadOption_Formula> y;
	double e1,e2,e3,e4;
	eigenvalues4(rho1,rhos, rhoc12,rhoc21,rhov,rho2,&e1,&e2,&e3,&e4);
	if((e1<=0) || (e2 <=0) || (e3<=0) || (e4<=0))
	{
		(*price)=-1.;
		if(rhos_flag) (*der_rhos)=0;
		if(rhov_flag) (*der_rhov)=0;
		if(rhoc12_flag) (*der_rhoc12)=0;
		if(rhoc21_flag) (*der_rhoc21)=0;
		
		
	}
	else
	{
		
		(*price)=y(a);
		if(rhos_flag) (*der_rhos)=y(10,a);
		if(rhov_flag) (*der_rhov)=y(11,a);
		if(rhoc12_flag) (*der_rhoc12)=y(12,a);
		if(rhoc21_flag) (*der_rhoc21)=y(13,a);
	}
	
}


CC_END_NAMESPACE()


#undef ARM_CF_EPS
#undef ARM_CF_MAXIT

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
