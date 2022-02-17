/*
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 * $Log: calibdirection.h,v $
 * Revision 1.1  2005/06/16 14:51:19  emezzine
 * Initial revision
 *
 *
 */

/*! \file enumcalibration.h
 *
 *  \brief 
 *
 *	\author  E.M Ezzine 
 *	\version 1.0
 *	\date June 2005
 */


#ifndef _INGPCALIB_ENUMCALIBRATION_H
#define _INGPCALIB_ENUMCALIBRATION_H

#include "gpbase/port.h"

CC_BEGIN_NAMESPACE( ARM )

struct ARM_CalibMethodType
{
	enum MethodType
		{
			Bootstrap1D,
			BootstrapND,
			Optimize,
			Optimize1D,
			Numerical,
			HW2FOnly,
			Unknown
		};
};

struct ARM_ModelFitterSolverType
{
    enum SolverType
	{
		NewtonRaphson = 0,
		NewtonRaphsonWithRetrial,
		SmoothNewthonRhaphson,
		NewthonRhaphsonNoThrow,
		NewtonRaphsonWithDichotomy,
		Dichotomy,
		NagSolver,
		Brent,
		NoSolverType
	};
};

struct ARM_ModelFitterOptimizerType
{
   enum OptimizerType
	{
		bounds_no_deriv = 0,
		nlin_lsq,
		lsq_deriv,
		lsq_check_deriv,
		OptimiseWithBrent,
		NoOptimizerType
	};
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/