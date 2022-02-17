/*!
 *
 * Copyright (c) CDC IXIS CM July 2007 Paris
 *
 * Version initiale 05/02/2007
 *
 *  basic functions for the closed form framework 
 *
 *	\file sabr_newversion.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date May 2007
 */
 
#ifndef _GP_CF_SABR_NEWVERSION_H
#define _GP_CF_SABR_NEWVERSION_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include "gpclosedforms/sabr_calibration.h"

CC_BEGIN_NAMESPACE(ARM)

double SABR_StrikeCuterExp_ImplicitVol(
									   double f,
									   double K,
									   double T,
									   double alpha,
									   double beta,
									   double rho,
									   double nu,
									   double alpha_exp);


double SABR_StrikeCuterTanh_ImplicitVol(
									   double f,
									   double K,
									   double T,
									   double alpha,
									   double beta,
									   double rho,
									   double nu,
									   double alpha_tanh,
									   double kb_tanh);



double SABR_StrikeCuterExp_VanillaOption(
									   double f,
									   double K,
									   double T,
									   double alpha,
									   double beta,
									   double rho,
									   double nu,
									   double alpha_exp,
									   int callput);

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
									   int callput);


/// Calibration of 3 parameters : Beta Fixed
SABR_ParameterSet*  SABR_CalibrateToSmile_withstrikecuter(ARM_GP_Vector* K_Vec,ARM_GP_Vector* ImpVol_Vec,
										  ARM_GP_Vector* Weigth_Vec,
										  double f,double beta,double tex,int flag,
										  double alpha_exp,double alpha_tanh, double kb_tanh,
									      int nbsteps,int algorithm,	double alpha0,double rho0,double nu0);

CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


