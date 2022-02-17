/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file heston_calibration.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_HESTON_CALIBRATION_H
#define _GP_CF_HESTON_CALIBRATION_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpbase/gpvector.h"
#include "gpbase/typedef.h"
#include "gpbase/removenagwarning.h"
#include "nag.h"
#include "nage04.h"
#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/normal_heston.h"
#include "gpnumlib/levmarq.h"
#include "gpnumlib/solver.h"
#include "gpclosedforms/vanilla_normal.h"
#include "gpclosedforms/inverse.h"

#include "gpclosedforms/optimization1.h"
#include <complex>
CC_USING_NS(std,complex)

CC_BEGIN_NAMESPACE(ARM)


class GeneralizedHeston_ParameterSet
{
	private:
			
		double V0;
		double omega;
		double theta;
		double ksi;
		double rho;
		double muJ;
		double sigmaJ;
		double lambda;

		double objective;
	public:
		
		double get_V0() {return			V0		;}
		double get_omega() {return		omega	;}
		double get_theta() {return		theta	;}
		double get_ksi() {return		ksi		;}
		double get_rho() {return		rho		;}
		double get_muJ() {return		muJ		;}
		double get_sigmaJ() {return		sigmaJ	;}
		double get_lambda() {return		lambda	;}
		double get_objective() {return	objective;}

		GeneralizedHeston_ParameterSet():V0(0),omega(0),theta(0),ksi(0),rho(0),muJ(0),sigmaJ(0),lambda(0),objective(0)
		{}
		GeneralizedHeston_ParameterSet(double V0_a,double omega_a,double theta_a,double ksi_a,double rho_a,double muJ_a,double sigmaJ_a,double lambda_a):
		V0(V0_a),omega(omega_a),theta(theta_a),ksi(ksi_a),rho(rho_a),muJ(muJ_a),sigmaJ(sigmaJ_a),lambda(lambda_a),objective(0)
		{}
		GeneralizedHeston_ParameterSet(Optimization_Result_Set r)
		{
			V0			=(*(r.OptimalParamSet))[0];
			omega		=(*(r.OptimalParamSet))[1];
			theta		=(*(r.OptimalParamSet))[2];
			ksi			=(*(r.OptimalParamSet))[3];
			rho			=(*(r.OptimalParamSet))[4];
			muJ			=(*(r.OptimalParamSet))[5];
			sigmaJ		=(*(r.OptimalParamSet))[6];
			lambda		=(*(r.OptimalParamSet))[7];
			objective=r.OptimalObjective;
		}
};

void GeneralizedHeston_DetermineDerivatives( double F,double K,double V0, double t,
							double omega,double theta,double ksi,double rho, double muJ, double sigmaJ, double lambda,
							double* impvol, double* der_V0,double* der_omega,double* der_theta,double* der_ksi,double* der_rho, double* der_muJ, double* der_sigmaJ, double* der_lambda,
							bool V0flag,bool omegaflag,bool thetaflag,bool ksiflag,bool rhoflag, bool muJflag, bool sigmaJflag, bool lambdaflag,int nbsteps);

GeneralizedHeston_ParameterSet*  GeneralizedHeston_CalibrateToSmile(ARM_GP_Vector* K_Vec,ARM_GP_Vector* ImpVol_Vec,double F, double t,int nb,int algorithm,
												double V0_0,double omega_0,double theta_0,double ksi_0,double rho_0,double muJ_0,double sigmaJ_0,double lambda_0);



CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

