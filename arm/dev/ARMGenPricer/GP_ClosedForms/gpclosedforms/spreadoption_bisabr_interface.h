/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_bisabr_interface.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2006
 */
 
#ifndef _GP_CF_SPREADOPTION_BISABR_INTERFACE_H
#define _GP_CF_SPREADOPTION_BISABR_INTERFACE_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"



CC_BEGIN_NAMESPACE(ARM)


double Export_BiSABR_SpreadOption(	double F1,double alpha1,double beta1,double rho1,double nu1,
									double F2,double alpha2,double beta2,double rho2,double nu2,
									double rhos,double rhov,double rhoc12,double rhoc21,double K,double T,int callput,int flag=2);

double Export_BiSABR_S3_SpreadOption(
									 double F1,double alpha1,double beta1,double rho1,double nu1,
									 double F2,double alpha2,double beta2,double rho2,double nu2,
									 double F3,double alpha3,double beta3,double rho3,double nu3,
									 double rhos,double rhov,double rhoc12,double rhoc21,double copularho,
									 double T,double a1,double b1,double k1,double a2,double b2,double k2,int flag,int nbsteps);


double Export_Complexified_BiSABR_SpreadOption(double F1,double alpha1,double beta1,double rho1,double nu1,
								  double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int callput,double rhos,double rhov,double rhoc12,double rhoc21,int flag = 2);

double Export_SABR_BetaEqualZero_Option (double f,double K,double T,double mu,double alpha,double rho,double nu,int callput);


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

