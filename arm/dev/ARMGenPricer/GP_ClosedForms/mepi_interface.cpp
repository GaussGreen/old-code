/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file mepi_interface.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2006
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpclosedforms/Mepi_Option.h"
#include "gpclosedforms/mepi_interface.h"

#include <cmath>

CC_BEGIN_NAMESPACE(ARM)



 ///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///  
///			Begining Exportable  Pricing Functions 
///
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

 

vector<double>* Export_Mepi_EDP_VanillaOption(
							  double P0,
							  double f0,
							  double T,
							  double K,
							  double R,
							  double Emin,
							  double Lmax,
							  double gamma0,
							  double gamma1,  
							  double sig,
							  double lambda,
							  double sigJ,
							  double r,
							  double s,
							  double mu,
							  double fees,
							  double volDrift,
							  double volVol,
							  int CallPut,
							  ARM_GP_Vector* A_params
							  )
{	
	const vector<double> a_Vec=A_params->GetValues();
	int schema=a_Vec[0];
	double TimeStepsNb=a_Vec[1];
	double alpha=a_Vec[2];
	int nbSpaceSteps=a_Vec[3];
	double limit_divisor=a_Vec[4];
	double nbSpaceStdDev=a_Vec[5];
	int nbHermiteComponent=a_Vec[6];
	double relativeShift=a_Vec[7];

	vector<double>* timegrid = mepi_timegrid_create( TimeStepsNb, alpha,T);

	vector<double>* result= Mepi_EDP_VanillaOption(
		P0,
		f0,
		T,
		K,
		R,
		Emin,
		Lmax,
		gamma0,
		gamma1,  
		sig,
		lambda,
		sigJ,
		r,
		s,
		mu,
		fees,
		nbSpaceSteps,
		nbSpaceStdDev,
		volDrift,
		volVol,
		limit_divisor,
		relativeShift,
		nbHermiteComponent,
		CallPut,
		timegrid,
		schema
		);
	delete timegrid;
	return result;
}

CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
