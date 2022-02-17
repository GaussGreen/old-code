/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file sabr_calibration.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_SABR_CALIBRATION_H
#define _GP_CF_SABR_CALIBRATION_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpbase/gpvector.h"
#include "gpbase/typedef.h"
#include "gpbase/removenagwarning.h"
#include "nag.h"
#include "nage04.h"
#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/sabrvanilla.h"
#include "gpclosedforms/optimization1.h"
#include <complex>
#include "gpnumlib/levmarq.h"

CC_USING_NS(std,complex)

CC_BEGIN_NAMESPACE(ARM)


class SABR_ParameterSet
{
	private:
	
		double alpha;
		double beta;
		double rho;
		double nu;
		double forward;
		double objective;
	public:
		
		double get_forward() {return forward;}
		double get_alpha() {return alpha;}
		double get_beta() {return beta;}
		double get_rho() {return rho;}
		double get_nu() {return nu;}
		double get_objective() {return objective;}

		SABR_ParameterSet():forward(0),alpha(0),beta(0),rho(0),nu(0),objective(0)
		{}
		SABR_ParameterSet(double f,double a,double b,double r,double n,double o):forward(f),alpha(a),beta(b),rho(r),nu(n),objective(o)
		{}
		SABR_ParameterSet(Optimization_Result_Set r)
		{
			alpha=(*(r.OptimalParamSet))[0];
			beta=(*(r.OptimalParamSet))[1];
			rho=(*(r.OptimalParamSet))[2];
			nu=(*(r.OptimalParamSet))[3];
			forward=(*(r.OptimalParamSet))[4];
			objective=r.OptimalObjective;
		}

	

		SABR_ParameterSet(Optimization_Result_Set r,double f, double beta0)
		{
			alpha=(*(r.OptimalParamSet))[0];
			beta=beta0;
			rho=(*(r.OptimalParamSet))[1];
			nu=(*(r.OptimalParamSet))[2];
			forward=f;
			objective=r.OptimalObjective;
		}

			SABR_ParameterSet(Optimization_Result_Set r,double f)
		{
			alpha=(*(r.OptimalParamSet))[0];
			beta=(*(r.OptimalParamSet))[1];
			rho=(*(r.OptimalParamSet))[2];
			nu=(*(r.OptimalParamSet))[3];
			forward=f;
			objective=r.OptimalObjective;
		}
};

SABR_ParameterSet*  SABR_CalibrateToSmile(ARM_GP_Vector* K_Vec,ARM_GP_Vector* ImpVol_Vec,ARM_GP_Vector* Weigth_Vec,double f,double t,int flag,int nbsteps,int algorithm,
												double alpha0,double beta0,double rho0,double nu0);
SABR_ParameterSet*  SABR_CalibrateToSmile(ARM_GP_Vector* K_Vec,ARM_GP_Vector* ImpVol_Vec,ARM_GP_Vector* Weigth_Vec,double t,int flag,int nbsteps,int algorithm,
												double alpha0,double beta0,double rho0,double nu0,double f0);
SABR_ParameterSet*  SABR_CalibrateToSmile(ARM_GP_Vector* K_Vec,ARM_GP_Vector* ImpVol_Vec,ARM_GP_Vector* Weigth_Vec,double f,double beta,double t,int flag,int nbsteps,int algorithm,
												double alpha0,double rho0,double nu0);
SABR_ParameterSet*  SABR_CalibrateToSmile(ARM_GP_Vector* K_Vec,ARM_GP_Vector* ImpVol_Vec,ARM_GP_Vector* Weigth_Vec,double f,double beta,double t,int flag,int nbsteps,int algorithm,
												double alpha0,double rho0,double nu0,double alphap,double rho0p,double nup,double rweight_alpha,double rweight_rho,double rweight_nu);

SABR_ParameterSet*	SABR_CalibrateToSmileBetaFixedToOne(ARM_GP_Vector* K_Vec,ARM_GP_Vector* ImpVol_Vec,ARM_GP_Vector* Weight_Vec,double f,double t,double atmvol);




CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

