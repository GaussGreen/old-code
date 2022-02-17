/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 03/15/2006
 *
 *  basic functions for the closed form framework 
 *
 *	\file trisabr_interface.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date March 2007
 */

#include "firsttoinc.h"
#include "gpbase/port.h"
#include "gpclosedforms/tri_sabr.h"
#include "gpclosedforms/trisabr_interface.h"
#include "gpclosedforms/eigenvalues.h"



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
			 double C_nbsteps)
{
return 	 TriSABR_VanillaOption(
		 F1,
		 alpha1,
		 beta1,
		 rho1,
		 nu1,	
		 F2,
		 alpha2,
		 beta2,
		 rho2,
		 nu2,		
		 F3,
		 alpha3,
		 beta3,
		 rho3,
		 nu3,		
		 C_rhos12,
		 C_rhos23,
		 C_rhos13,
		 C_rhov12,
		 C_rhov23,
		 C_rhov13,
		 C_rhoc12,
		 C_rhoc21,
		 C_rhoc23,
		 C_rhoc32,
		 C_rhoc13,
		 C_rhoc31,
		 C_K,
		 C_T,
		 C_callput,
		 C_flag,
		 C_nbsteps);

}


void Export_TriSABR_Eigenvalues(
			double rho1,double rho2,double rho3,
			double rhos12, double rhos23, double rhos13,
			double rhov12, double rhov23, double rhov13,
			double rhoc12, double rhoc13,
			double rhoc21, double rhoc23,
			double rhoc31, double rhoc32,
			double* e1,double* e2,double* e3,double* e4,double* e5,double* e6)
{
	 TriSABR_Eigenvalues(
			 rho1,	  rho2,	   rho3,
			 rhos12,  rhos23,  rhos13,
			 rhov12,  rhov23,  rhov13,
			 rhoc12,  rhoc13,
			 rhoc21,  rhoc23,
			 rhoc31,  rhoc32,
			e1,e2,e3,e4,e5,e6);
	 return;
}


CC_END_NAMESPACE()
 

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/