/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file tri_sabr.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date Feb 2006
 */
 
#ifndef _GP_CF_TRI_SABR_H
#define _GP_CF_TRI_SABR_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include <complex>
using namespace std; 

CC_BEGIN_NAMESPACE(ARM)


 double TriSABR_VanillaOption(
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


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


