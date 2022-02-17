/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 03/15/2006
 *
 *  basic functions for the closed form framework 
 *
 *	\file extended_sabr.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date March 2006
 */
 
#ifndef _GP_CF_EXTENDED_BISABR_INTERFACE_H
#define _GP_CF_EXTENDED_BISABR_INTERFACE_H

#include "firsttoinc.h"
#include "gpbase/port.h"
#include "gpbase/gpvector.h"
#include "gpbase/typedef.h"

//#include "gpclosedforms/bisabr_calibration.h"


CC_BEGIN_NAMESPACE(ARM)



////////////////////////////////////////////////////////////////////////////////////////////
///
///				Calibration of BiSABR
///
////////////////////////////////////////////////////////////////////////////////////////////

class BiSABR_ParameterSet;

BiSABR_ParameterSet*  Export_BiSABR_CalibrateToSmile(
			std::vector<double>* F1_Vec,std::vector<double>* Alpha1_Vec,std::vector<double>* Beta1_Vec,std::vector<double>* Rho1_Vec,std::vector<double>* Nu1_Vec,
			std::vector<double>* F2_Vec,std::vector<double>* Alpha2_Vec,std::vector<double>* Beta2_Vec,std::vector<double>* Rho2_Vec,std::vector<double>* Nu2_Vec,
			std::vector<double>* K_Vec,std::vector<double>* T_Vec,std::vector<double>* Price_Vec,std::vector<double>* Weight_Vec,
						double rhos_0,double rhov_0,double rhoc12_0,double rhoc21_0,
						double converg_prec, int nbIter_max, double first_step_max,int rhos_flag,int flag = 2);


void Export_BiSABR_CorrelationEvolution(double rho1,double rho2,double rhos,double rhov,double rhoc12,double rhoc21,
								 double newrho1, double newrho2,
								 double* newrhov,double* newrhoc12,double* newrhoc21);


double Export_BiSABR_Digital_SpreadOption(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21,int flag = 2);


double Export_BiSABR_Digital_SpreadOption_PayS1(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21,int flag = 2);

double Export_BiSABR_Digital_SpreadOption_PayS2(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21,int flag = 2);


double Export_BiSABR_Digital_SpreadOption_PayS3(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
					double rhos,double rhov,double rhoc12,double rhoc21,double S3,double sigma3, double rho13, double rho23, double K,double T,int CallPut,int flag = 2);



CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

