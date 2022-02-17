/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file bisabr_digital_spreadoption.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date August 2006
 */
 
#ifndef _GP_CF_BISABR_DIGITAL_SPREADOPTION_H
#define _GP_CF_BISABR_DIGITAL_SPREADOPTION_H

#include "firsttoinc.h"
#include "gpbase/port.h"
#include <cmath>
#include <complex>

using std::sqrt;
using std::log;
using std::real;
using std::imag;
using std::exp;
using std::complex;

CC_BEGIN_NAMESPACE(ARM)


complex<double> GIntegralAnalytical(complex<double> a1, complex<double> b1,complex<double> x1);

complex<double> BiSABRIntegral(complex<double> A,complex<double> B,complex<double> G,complex<double> F,complex<double> gamma);

complex<double> Complexified_BiSABR_SpreadOption(complex<double> F1,complex<double> alpha1,complex<double> beta1,complex<double> rho1,complex<double> nu1,
													  complex<double> F2,complex<double> alpha2,complex<double> beta2,complex<double> rho2,complex<double> nu2,
					complex<double> K,complex<double> T,complex<double> rhos,complex<double> rhov,complex<double> rhoc12,complex<double> rhoc21);

double BiSABR_Digital_SpreadOption(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21);


double BiSABR_Digital_SpreadOption_2(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21);

double BiSABR_Digital_SpreadOption_3(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21);

double BiSABR_Digital_SpreadOption_4(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21);

double BiSABR_Digital_SpreadOption_4(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,complex<double> alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21);

/// compute the SABR parameters associated with  S1/S2 and K/S2 for the change of numeraires

void BiSABREqivalents(double F1,double alpha1,double beta1,double rho1,double nu1,		
				double F2,double alpha2,double beta2,double rho2,double nu2,
				double rhos,double rhov,double rhoc12,double rhoc21	,double K,
				/* output */
				double& alphaZ,double& betaZ,double& rhoZ, double& nuZ,
				complex<double>& alphaY,double& betaY,double& rhoY, double& nuY,
				double& rhoYZ,double& rhoalphaYalphaZ,double& rhoalphaYZ,double& rhoYalphaZ
				);

double BiSABR_Digital_SpreadOption(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,complex<double> alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21);

double BiSABR_Digital_SpreadOption_PayS1(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21);

double BiSABR_Digital_SpreadOption_PayS1_4(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21);

double BiSABR_Digital_SpreadOption_PayS2(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21);

double BiSABR_Digital_SpreadOption_PayS2_4(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21);


double Corrected_BiSABR_Digital_SpreadOption_PayS1(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21);

double Corrected_BiSABR_Digital_SpreadOption_PayS2(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21);


double BiSABR_Digital_SpreadOption_PayS3(
									   double F1,double alpha1,double beta1,double rho1,double nu1,
									   double F2,double alpha2,double beta2,double rho2,double nu2,
									   double rhos,double rhov,double rhoc12,double rhoc21,
									   double F3,
									   double sigma3,double rho13,double rho23,
									   double K,double T,int CallPut);

double BiSABR_Digital_SpreadOption_PayS3_4(
									   double F1,double alpha1,double beta1,double rho1,double nu1,
									   double F2,double alpha2,double beta2,double rho2,double nu2,
									   double rhos,double rhov,double rhoc12,double rhoc21,
									   double F3,
									   double sigma3,double rho13,double rho23,
									   double K,double T,int CallPut);


double Corrected_BiSABR_Digital_SpreadOption_PayS3(double F1,double alpha1,double beta1,double rho1,double nu1,
									   double F2,double alpha2,double beta2,double rho2,double nu2,
									   double rhos,double rhov,double rhoc12,double rhoc21,
									   double F3,
									   double sigma3,double rho13,double rho23,
									   double K,double T,int CallPut);

CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


