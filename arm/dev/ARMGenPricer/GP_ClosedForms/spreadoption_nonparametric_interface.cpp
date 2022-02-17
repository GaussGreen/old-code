/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 03/15/2007
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_nonparametric_interface.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date March 2007
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpclosedforms/spreadoption_nonparametric_formula.h"
#include "gpclosedforms/spreadoption_nonparametric_interface.h"
#include "gpbase/gpvector.h"
#include "gpbase/gplinalgtypedef.h"
#include <cmath>
#include <complex>
#include "gpbase/gpmatrix.h"
#include <glob/expt.h>


CC_BEGIN_NAMESPACE(ARM)


double Export_Nonparametric_CompleteSpreadoption(
		ARM_GP_Vector* strike_Vec1,
		ARM_GP_Vector* vol_Vec1,
		ARM_GP_Vector* strike_Vec2,
		ARM_GP_Vector* vol_Vec2,
		double S1,double S2,
		double index_begin1,double index_end1,int flag_begin1,int flag_end1,
		double index_begin2,double index_end2,int flag_begin2,int flag_end2,
		double correlation,double maturity,double a1,double b1,double k1,double a2,double b2,double k2,
		int nbsteps,int algorithm,int smiletype
	 )
{
	ArgumentList a(
		strike_Vec1,vol_Vec1,strike_Vec2, vol_Vec2,
		S1,S2,
		 index_begin1, index_end1, flag_begin1,flag_end1,
		 index_begin2, index_end2, flag_begin2,flag_end2,
		 correlation, maturity, a1, b1, k1, a2, b2, k2,
		 nbsteps, algorithm);
	if(smiletype==0)
	{
		Power_Expression<ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula> y;
		return y(a);
	}
	else if(smiletype==1)
	{
		Power_Expression<ARM_CF_NonParametric_Gaussian_PowerSpreadOption2_Formula> y;
		return y(a);
	}
	else if(smiletype==2)
	{
		Power_Expression<ARM_CF_NonParametric_Gaussian_PowerSpreadOption3_Formula> y;
		return y(a);
	}
	else
	{
		Power_Expression<ARM_CF_NonParametric_Gaussian_PowerSpreadOption4_Formula> y;
		return y(a);
	}
}





CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/