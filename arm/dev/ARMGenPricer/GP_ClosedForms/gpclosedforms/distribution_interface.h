/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file distribution_interface.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_DISTRIBUTION_INTERFACE_H
#define _GP_CF_DISTRIBUTION_INTERFACE_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include <complex>
CC_USING_NS(std,complex)

CC_BEGIN_NAMESPACE(ARM)

double Export_ShiftedLogNormal_Quantile(double f,double K,double tex,double sigma, double alpha);

double Export_SABR_Quantile(double f,double k,double T,double alpha,double beta, double rho, double nu, double Sabr_Type, double nbsteps);

double Export_SABR_Quantile(double f,double k,double T,double alpha,double beta, double rho, double nu, double Sabr_Type, double nbsteps,
							double a,double b,double c);


double Export_ShiftedLogNormal_Distribution(double f,double K,double tex,double sigma, double alpha);

double Export_SABR_Distribution(double f,double k,double T,double alpha,double beta, double rho, double nu, double Sabr_Type, double nbsteps);

double Export_SABR_Distribution(double f,double k,double T,double alpha,double beta, double rho, double nu, double Sabr_Type, double nbsteps,
								double a,double b,double c);


double Export_SABR_Density(double f,double x,double t,double alpha,double beta, double rho, double nu, double Sabr_Type, double nbsteps);

double Export_SABR_Density(double f,double x,double t,double alpha,double beta, double rho, double nu, double Sabr_Type, double nbsteps,
						   double a,double b,double c);


double Export_BiSABR_Quantile(double f1,double alpha1,double beta1, double rho1, double nu1,
		double f2,double alpha2,double beta2, double rho2, double nu2,
		double rhos,double rhov,double rhoc12, double rhoc21, double x,double T,int flag );

double Export_BiSABR_Quantile(double f1,double alpha1,double beta1, double rho1, double nu1,
		double f2,double alpha2,double beta2, double rho2, double nu2,
		double rhos,double rhov,double rhoc12, double rhoc21, double x,double T,int flag ,
		double a,double b,double c);


double Export_BiSABR_Distribution(double f1,double alpha1,double beta1, double rho1, double nu1,
		double f2,double alpha2,double beta2, double rho2, double nu2,
		double rhos,double rhov,double rhoc12, double rhoc21, double x,double T,int flag );

double Export_BiSABR_Distribution(double f1,double alpha1,double beta1, double rho1, double nu1,
		double f2,double alpha2,double beta2, double rho2, double nu2,
		double rhos,double rhov,double rhoc12, double rhoc21, double x,double T,int flag ,
		double a,double b,double c);


double Export_BiSABR_Density(double f1,double alpha1,double beta1, double rho1, double nu1,
		double f2,double alpha2,double beta2, double rho2, double nu2,
		double rhos,double rhov,double rhoc12, double rhoc21, double x,double T,int flag );

double Export_BiSABR_Density(double f1,double alpha1,double beta1, double rho1, double nu1,
		double f2,double alpha2,double beta2, double rho2, double nu2,
		double rhos,double rhov,double rhoc12, double rhoc21, double x,double T,int flag,
		double a,double b,double c);

double Export_Shifted2LogNormal_Quantile(double f1,double sigma1,double f2,double sigma2,
										 double alpha,double rho,double x,double T,int n);


double Export_Shifted2LogNormal_Distribution(double f1,double sigma1,double f2,double sigma2,
											 double alpha,double rho,double x,double T,int n);


double Export_Shifted2LogNormal_Density(double f1,double sigma1,double f2,double sigma2,
										double alpha,double rho,double x, double T,int n);


double Export_Student_Quantile(double rank, double x);

double Export_Student_Distribution(double rank, double x);

double Export_NonParametric_Distribution(ARM_GP_Vector* x,ARM_GP_Vector* y,ARM_GP_Vector* y2,
							double index_begin,double index_end,
							int beginflag,int endflag, double S,double T, double strike);


double Export_NonParametric_Quantile(ARM_GP_Vector* x,ARM_GP_Vector* y,ARM_GP_Vector* y2,
							double index_begin,double index_end,
							int beginflag,int endflag,double S,double T, double proba);


double Export_NonParametric_LogVolatility(ARM_GP_Vector* x,ARM_GP_Vector* y,
							double index_begin,double index_end,
							int beginflag,int endflag, double strike);


double Export_NonParametric_NormalVolatility(ARM_GP_Vector* x,ARM_GP_Vector* y,
							double index_begin,double index_end,
							int beginflag,int endflag, double strike);


double Export_NonParametric_N_Distribution(ARM_GP_Vector* x,ARM_GP_Vector* y,
							double index_begin,double index_end,
							int beginflag,int endflag, double S,double T, double strike);

double Export_NonParametric_LN_Distribution(ARM_GP_Vector* x,ARM_GP_Vector* y,
							double index_begin,double index_end,
							int beginflag,int endflag, double S,double T, double strike);

double Export_NonParametric_N_Quantile(ARM_GP_Vector* x,ARM_GP_Vector* y,
							double index_begin,double index_end,
							int beginflag,int endflag,double S,double T, double proba);

double Export_NonParametric_LN_Quantile(ARM_GP_Vector* x,ARM_GP_Vector* y,
							double index_begin,double index_end,
							int beginflag,int endflag,double S,double T, double proba);




CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

