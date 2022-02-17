/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 03/15/2006
 *
 *  basic functions for the closed form framework 
 *
 *	\file bisabr_calibration.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date March 2006
 */
 
#ifndef _GP_CF_BISABR_CALIBRATION_H
#define _GP_CF_BISABR_CALIBRATION_H

#include "firsttoinc.h"
#include "gpbase/port.h"
#include "gpbase/gpvector.h"
#include "gpbase/typedef.h"
#include "gpbase/removenagwarning.h"
#include "nag.h"
#include "nage04.h"
#include "gpclosedforms/basic_distributions.h"

#include "gpclosedforms/optimization1.h"
#include <complex>
CC_USING_NS(std,complex)

CC_BEGIN_NAMESPACE(ARM)

void BiSABR_CorrelationEvolution(double rho1,double rho2,double rhos,double rhov,double rhoc12,double rhoc21,
								 double newrho1, double newrho2,
								 double* newrhov,double* newrhoc12,double* newrhoc21);


class BiSABR_ParameterSet
{
	private:
			
		double rhos;
		double rhov;
		double rhoc12;
		double rhoc21;

		double objective;
	public:
		
		double get_rhos() {return		rhos		;}
		double get_rhov() {return		rhov		;}
		double get_rhoc12() {return		rhoc12		;}
		double get_rhoc21() {return		rhoc21		;}
		
		double get_objective() {return	objective;}

		BiSABR_ParameterSet():rhos(0),rhov(0),rhoc12(0),rhoc21(0),objective(0)
		{}
		BiSABR_ParameterSet(double rhos_a,double rhov_a,double rhoc12_a,double rhoc21_a):
		rhos(rhos_a),rhov(rhov_a),rhoc12(rhoc12_a),rhoc21(rhoc21_a),objective(0)
		{}
		BiSABR_ParameterSet(Optimization_Result_Set r)
		{
			rhos			=(*(r.OptimalParamSet))[0];
			rhov			=(*(r.OptimalParamSet))[1];
			rhoc12			=(*(r.OptimalParamSet))[2];
			rhoc21			=(*(r.OptimalParamSet))[3];
			objective		=   r.OptimalObjective;
		}
		BiSABR_ParameterSet(Optimization_Result_Set r,double rhos0)
		{
			rhos			=rhos0;
			rhov			=(*(r.OptimalParamSet))[0];
			rhoc12			=(*(r.OptimalParamSet))[1];
			rhoc21			=(*(r.OptimalParamSet))[2];
			objective		=   r.OptimalObjective;
		}
};

void BiSABR_DetermineDerivatives( double F1,double alpha1,double beta1,double rho1,double nu1,
								 double F2,double alpha2,double beta2,double rho2,double nu2,
								double K,double T,
								double rhos,double rhov,double rhoc12,double rhoc21,
							double* price, double* der_rhos,double* der_rhov,double* der_rhoc12,double* der_rhoc21,
							bool rhos_flag,bool rhov_flag,bool rhoc12_flag,bool rhoc21_flag,int flag);

BiSABR_ParameterSet*  BiSABR_CalibrateToSmile(
			std::vector<double>* F1_Vec,std::vector<double>* Alpha1_Vec,std::vector<double>* Beta1_Vec,std::vector<double>* Rho1_Vec,std::vector<double>* Nu1_Vec,
			std::vector<double>* F2_Vec,std::vector<double>* Alpha2_Vec,std::vector<double>* Beta2_Vec,std::vector<double>* Rho2_Vec,std::vector<double>* Nu2_Vec,
			std::vector<double>* K_Vec,std::vector<double>* T_Vec,std::vector<double>* Price_Vec,std::vector<double>* Weight_Vec,
						double rhos_0,double rhov_0,double rhoc12_0,double rhoc21_0,
						double converg_prec, int nbIter_max, double first_step_max,int flag);

BiSABR_ParameterSet*  BiSABR_CalibrateToSmile_ConstantRhos(
			std::vector<double>* F1_Vec,std::vector<double>* Alpha1_Vec,std::vector<double>* Beta1_Vec,std::vector<double>* Rho1_Vec,std::vector<double>* Nu1_Vec,
			std::vector<double>* F2_Vec,std::vector<double>* Alpha2_Vec,std::vector<double>* Beta2_Vec,std::vector<double>* Rho2_Vec,std::vector<double>* Nu2_Vec,
			std::vector<double>* K_Vec,std::vector<double>* T_Vec,std::vector<double>* Price_Vec,std::vector<double>* Weight_Vec,
						double rhos,double rhov_0,double rhoc12_0,double rhoc21_0,
						double converg_prec, int nbIter_max, double first_step_max,int flag);


CC_END_NAMESPACE()


#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

