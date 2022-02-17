/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_glambda_interface.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date sep 2005
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpclosedforms/spreadoption_glambda_student_formula.h"
#include "gpclosedforms/spreadoption_glambda_interface.h"
#include "gpbase/gpvector.h"
#include "gpbase/gplinalgtypedef.h"
#include <cmath>
#include <complex>
#include "gpbase/gpmatrix.h"
#include <glob/expt.h>


CC_BEGIN_NAMESPACE(ARM)

double Export_Student_GLambda_Power_SpreadOption(
								   double l1a, double l2a, double l3a, double l4a, double l5a,double l6a,
								   double l1b, double l2b, double l3b, double l4b, double l5b,double l6b,
								   double copula_corr,double copula_degre,
								   double a10,double b10,double k10,double a20,double b20,double k20,int n)
{
	ArgumentList a( l1a,l2a,l3a,l4a,l5a,l6a,
					l1b,l2b,l3b,l4b,l5b,l6b,
					copula_corr,copula_degre,
					a10, b10, k10, a20, b20, k20,n);

	Power_Expression<ARM_CF_GLambda_Student_PowerSpreadOption_Formula> y;
	return y(a);
}





double Export_Student_GLambda_Power_Digital_SpreadOption(
								   double l1a, double l2a, double l3a, double l4a, double l5a,double l6a,
								   double l1b, double l2b, double l3b, double l4b, double l5b,double l6b,
								   double copula_corr,double copula_degre,
								   double a20,double b20,double k20,int n)
{
	ArgumentList a( l1a,l2a,l3a,l4a,l5a,l6a,
					l1b,l2b,l3b,l4b,l5b,l6b,
					copula_corr,copula_degre,
					a20, b20, k20,n);

	Power_Expression<ARM_CF_GLambda_Student_DigitalPowerSpreadOption_Formula> y;
	return y(a);
}

double Export_Student_GLambda_Power_Index2Digital_SpreadOption(
								   double l1a, double l2a, double l3a, double l4a, double l5a,double l6a,
								   double l1b, double l2b, double l3b, double l4b, double l5b,double l6b,
								   double copula_corr,double copula_degre,
								   double a20,double b20,double k20,int n)
{
	ArgumentList a( l1a,l2a,l3a,l4a,l5a,l6a,
					l1b,l2b,l3b,l4b,l5b,l6b,
					copula_corr,copula_degre,
					a20, b20, k20,n);

	Power_Expression<ARM_CF_GLambda_Student_Index2DigitalPowerSpreadOption_Formula> y;
	return y(a);
}


double Export_Student_GLambda_Power_Index1Digital_SpreadOption(
								   double l1a, double l2a, double l3a, double l4a, double l5a,double l6a,
								   double l1b, double l2b, double l3b, double l4b, double l5b,double l6b,
								   double copula_corr,double copula_degre,
								   double a20,double b20,double k20,int n)
{
	ArgumentList a( l1a,l2a,l3a,l4a,l5a,l6a,
					l1b,l2b,l3b,l4b,l5b,l6b,
					copula_corr,copula_degre,
					a20, b20, k20,n);

	Power_Expression<ARM_CF_GLambda_Student_Index1DigitalPowerSpreadOption_Formula> y;
	return y(a);
}


double Export_Student_GLambda_SpreadOption(
								   double l1a, double l2a, double l3a, double l4a, double l5a,double l6a,
								   double l1b, double l2b, double l3b, double l4b, double l5b,double l6b,
								   double copula_corr,double copula_degre,
								   double k20,int n)
{
	ArgumentList a( l1a,l2a,l3a,l4a,l5a,l6a,
					l1b,l2b,l3b,l4b,l5b,l6b,
					copula_corr,copula_degre,
					k20,n);

	Power_Expression<ARM_CF_GLambda_Student_SpreadOption_Formula> y;
	return y(a);
}


double Export_GLambda_Distribution(  double l1, double l2, double l3, double l4, double l5,double l6,double x)
{
	ArgumentList a(l1,l2,l3,l4,l5,l6);
	return GLambda_Smile::probability_distribution(a,x,1.);
}


double Export_GLambda_Quantile(  double l1, double l2, double l3, double l4, double l5,double l6,double x)
{
	ArgumentList a(l1,l2,l3,l4,l5,l6);
	return GLambda_Smile::quantile(a,x,1.);
}


double Export_GLambda_CompleteSpreadoption(
		ARM_GP_Vector* l1a_Vec,
		ARM_GP_Vector* l2a_Vec,
		ARM_GP_Vector* l3a_Vec,
		ARM_GP_Vector* l4a_Vec,
		ARM_GP_Vector* l5a_Vec,
		ARM_GP_Vector* l6a_Vec,
		ARM_GP_Vector* l1b_Vec,
		ARM_GP_Vector* l2b_Vec,
		ARM_GP_Vector* l3b_Vec,
		ARM_GP_Vector* l4b_Vec,
		ARM_GP_Vector* l5b_Vec,
		ARM_GP_Vector* l6b_Vec,
		ARM_GP_Vector* Dis_Vec,
	 
	  double copula_corr,
	  double copula_degre,
	  double k,
	  double n
	 )
{
	int m=l1a_Vec->size();
	if(l2a_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"Export_GLambda_CompleteSpreadoption: l2a_Vec do not have the same size!" );
	}
	if(l3a_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"Export_GLambda_CompleteSpreadoption: l3a_Vec do not have the same size!" );
	}
	if(l4a_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"Export_GLambda_CompleteSpreadoption: l4a_Vec do not have the same size!" );
	}
	if(l5a_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"Export_GLambda_CompleteSpreadoption: l5a_Vec do not have the same size!" );
	}
	if(l6a_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"Export_GLambda_CompleteSpreadoption: l6a_Vec do not have the same size!" );
	}
	if(l1b_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"Export_GLambda_CompleteSpreadoption: l1b_Vec do not have the same size!" );
	}
	if(l2b_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"Export_GLambda_CompleteSpreadoption: l2b_Vec do not have the same size!" );
	}
	if(l3b_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"Export_GLambda_CompleteSpreadoption: l3b_Vec do not have the same size!" );
	}
	if(l4b_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"Export_GLambda_CompleteSpreadoption: l4b_Vec do not have the same size!" );
	}
	if(l5b_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"Export_GLambda_CompleteSpreadoption: l5b_Vec do not have the same size!" );
	}
	if(l6b_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"Export_GLambda_CompleteSpreadoption: l6b_Vec do not have the same size!" );
	}
	if(Dis_Vec->size() != m)
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +"Export_GLambda_CompleteSpreadoption: Dis_Vec  do not have the same size!" );
	}
	double sum=0;
	int i;
	for(i=0;i<m;i++)
	{
		sum+=Export_Student_GLambda_SpreadOption(
								    (*l1a_Vec)[i],  (*l2a_Vec)[i],  (*l3a_Vec)[i],  (*l4a_Vec)[i],  (*l5a_Vec)[i], (*l6a_Vec)[i],
								    (*l1b_Vec)[i],  (*l2b_Vec)[i],  (*l3b_Vec)[i],  (*l4b_Vec)[i],  (*l5b_Vec)[i], (*l6b_Vec)[i],
								    copula_corr, copula_degre,k, n)*(*Dis_Vec)[i];
	}
	return sum;


}

CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/