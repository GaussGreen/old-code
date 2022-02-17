/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 03/15/2006
 *
 *  basic functions for the closed form framework 
 *
 *	\file trisabr_interface.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date March 2007
 */
 
#ifndef _GP_CF_EXTENDED_TRISABR_INTERFACE_H
#define _GP_CF_EXTENDED_TRISABR_INTERFACE_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpbase/gpvector.h"
#include "gpbase/typedef.h"



CC_BEGIN_NAMESPACE(ARM)


double Export_TriSABR_VanillaOption(
			 double F1,
			 double alpha1,
			 double beta1,
			 double rho1,
			 double nu1,
								 
			 double F2,
			 double alpha2,
			 double beta2,
			 double rho2,
			 double nu2,
								 
			 double F3,
			 double alpha3,
			 double beta3,
			 double rho3,
			 double nu3,
			 
			 double C_rhos12,
			 double C_rhos23,
			 double C_rhos13,
			 double C_rhov12,
			 double C_rhov23,
			 double C_rhov13,
			 double C_rhoc12,
			 double C_rhoc21,
			 double C_rhoc23,
			 double C_rhoc32,
			 double C_rhoc13,
			 double C_rhoc31,
			 double C_K,
			 double C_T,
			 int C_callput,
			 int C_flag,
			 double C_nbsteps);

void Export_TriSABR_Eigenvalues(
			double rho1,double rho2,double rho3,
			double rhos12, double rhos23, double rhos13,
			double rhov12, double rhov23, double rhov13,
			double rhoc12, double rhoc13,
			double rhoc21, double rhoc23,
			double rhoc31, double rhoc32,
			double* e1,double* e2,double* e3,double* e4,double* e5,double* e6);

CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

