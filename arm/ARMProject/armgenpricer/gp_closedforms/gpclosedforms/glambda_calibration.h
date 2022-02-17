/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file GLambda_calibration.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date Jal5ary 2004
 */
 
#ifndef _GP_CF_GLambda_CALIBRATION_H
#define _GP_CF_GLambda_CALIBRATION_H

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


class GLambda_ParameterSet
{
	private:
	
		double l1;
		double l2;
		double l3;
		double l4;
		double l5;
		double l6;
		double objective;
	public:
		
		double get_l1() {return l1;}
		double get_l2() {return l2;}
		double get_l3() {return l3;}
		double get_l4() {return l4;}
		double get_l5() {return l5;}
		double get_l6() {return l6;}
		double get_objective() {return objective;}

		GLambda_ParameterSet():l1(0),l2(0),l3(0),l4(0),l5(0),l6(0),objective(0)
		{}
		GLambda_ParameterSet(double l10,double l20,double l30,double l40,double l50,double l60,double o):l1(l10),l2(l20),l3(l30),l4(l40),l5(l50),l6(l60),objective(o)
		{}
		GLambda_ParameterSet(Optimization_Result_Set r)
		{
			l1=(*(r.OptimalParamSet))[0];
			l2=(*(r.OptimalParamSet))[1];
			l3=(*(r.OptimalParamSet))[2];
			l4=(*(r.OptimalParamSet))[3];
			l5=(*(r.OptimalParamSet))[4];
			l6=(*(r.OptimalParamSet))[5];
			objective=r.OptimalObjective;
		}
};

GLambda_ParameterSet*  GLambda_CalibrateFromSABR(double f,double alpha,double beta,double rho,double nu, double t,int flag,int nbsteps,double scope,
													double initial_l1,double initial_l2,double initial_l3,double initial_l4,double initial_l5,double initial_l6,
													int algorithm);


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

