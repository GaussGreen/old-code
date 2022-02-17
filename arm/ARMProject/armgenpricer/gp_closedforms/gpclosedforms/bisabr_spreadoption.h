/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file bisabr_spreadoption.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2006
 */
 
#ifndef _GP_CF_BISABR_SPREADOPTION_H
#define _GP_CF_BISABR_SPREADOPTION_H

#include "firsttoinc.h"
#include "gpbase/port.h"
#include <complex>
using namespace std; 

CC_BEGIN_NAMESPACE(ARM)

double GIntegralAnalytical(double a1, double b1,double x1);

double BiSABRIntegral(double A,double B,double G,double F,double gamma);

double BetaEqualZeroSABR(double f,double K,double T,double mu,double alpha,double rho,double nu,int callput);

double BiSABR_SpreadOption(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,double rhos,double rhov,double rhoc12,double rhoc21);

double BiSABR_SpreadOption2(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,double rhos,double rhov,double rhoc12,double rhoc21);

double BiSABR_SpreadOption3(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,double rhos,double rhov,double rhoc12,double rhoc21);

double BiSABR_SpreadOption4(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,double rhos,double rhov,double rhoc12,double rhoc21);

complex<double> Complexified_BiSABR_SpreadOption(complex<double> F1,complex<double> alpha1,complex<double> beta1,complex<double> rho1,complex<double> nu1,complex<double> F2,complex<double> alpha2,complex<double> beta2,complex<double> rho2,complex<double> nu2,
					complex<double> K,complex<double> T,complex<double> rhos,complex<double> rhov,complex<double> rhoc12,complex<double> rhoc21);


double Complete_BiSABR_SpreadOption(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21);

double Complete_BiSABR_SpreadOption_2(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21);


double Complete_BiSABR_SpreadOption2(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21);

double BiSABR_SpreadOption_bilog_corrected(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,double rhos,double rhov,double rhoc12,double rhoc21);

double Complete_BiSABR_SpreadOption_bilog_corrected(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21);

double Packaged_BiSABR_SpreadOption(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
								  double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21,int flag);

double Packaged_Complexified_BiSABR_SpreadOption(double F1,double alpha1,double beta1,double rho1,double nu1,
								  double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int callput,double rhos,double rhov,double rhoc12,double rhoc21,int flag);





CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


