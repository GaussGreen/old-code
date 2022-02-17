/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_glambda_interface.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date Sep 2005
 */
 
#ifndef _GP_CF_SPREADOPTION_GLAMBDA_INTERFACE_H
#define _GP_CF_SPREADOPTION_GLAMBDA_INTERFACE_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include "gpclosedforms/basic_distributions.h" /// for ArgumentList_Checking_Result
#include "gpclosedforms/powerspreadoption.h"		/// for instanciation of the templates
#include "gpclosedforms/smile_glambda.h"
#include "gpbase/gpvector.h"

using namespace std;

CC_BEGIN_NAMESPACE(ARM)

double Export_Student_GLambda_Power_SpreadOption(
								   double l1a, double l2a, double l3a, double l4a, double l5a, double l6a,
								   double l1b, double l2b, double l3b, double l4b, double l5b, double l6b,
								   double copula_corr,double copula_degre,
								   double a10,double b10,double k10,double a20,double b20,double k20,int n);


double Export_Student_GLambda_Power_Digital_SpreadOption(double l1a,double l2a,double l3a,double l4a,double l5a,double l6a,
												double l1b,double l2b,double l3b,double l4b,double l5b,double l6b,
												double copula_corr,double copula_degre,
												double a20,double b20,double k20,int n);


double Export_Student_GLambda_Power_Index2Digital_SpreadOption(
								   double l1a, double l2a, double l3a, double l4a, double l5a,double l6a,
								   double l1b, double l2b, double l3b, double l4b, double l5b,double l6b,
								   double copula_corr,double copula_degre,
								   double a20,double b20,double k20,int n);


double Export_Student_GLambda_Power_Index1Digital_SpreadOption(
								   double l1a, double l2a, double l3a, double l4a, double l5a,double l6a,
								   double l1b, double l2b, double l3b, double l4b, double l5b,double l6b,
								   double copula_corr,double copula_degre,
								   double a20,double b20,double k20,int n);

double Export_Student_GLambda_SpreadOption(
								   double l1a, double l2a, double l3a, double l4a, double l5a,double l6a,
								   double l1b, double l2b, double l3b, double l4b, double l5b,double l6b,
								   double copula_corr,double copula_degre,
								   double k20,int n);


double Export_GLambda_Distribution(double l1, double l2, double l3, double l4, double l5,double l6,double x);


double Export_GLambda_Quantile(  double l1, double l2, double l3, double l4, double l5,double l6,double x);

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
	 );


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

