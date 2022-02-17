/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file solver.h
 *
 *  \brief 
 *	\author  A. Triki
 *	\version 1.0
 *	\date October 2005
 */

#ifndef _INGPNUMLIB_RUNGEKUTTA_H
#define _INGPNUMLIB_RUNGEKUTTA_H

///GP Base
#include "gpbase/ostringstream.h"
#include "gpbase/warningkeeper.h"

#include "odefunctions.h"

#include <glob/expt.h>
#include <cmath>
#include <iomanip> /// for setprecision()
using namespace std;

CC_BEGIN_NAMESPACE( ARM )

extern const int RK4_NBTMP_VAR;
extern const int RK5_NBTMP_VAR;

void Rk4ToNextStep(ARM_GP_Vector* yInOut, ARM_GP_Vector* dydx, int n, double x, double h,
			const ARM_ODEFunc& function, vector<ARM_GP_Vector>& solverVars);

void RkFunction(ARM_GP_Vector* yInOut, int nvar, double x1,double x2, int nstep, ARM_GP_Vector* xp,
			const ARM_ODEFunc& function, vector<ARM_GP_Vector>& solverVars);

void rkckToNextStep(ARM_GP_Vector* yInOut, ARM_GP_Vector* dydx, int n, double x, double h,
			ARM_GP_Vector* yerr, const ARM_ODEFunc& function, vector<ARM_GP_Vector>& solverVars);

void RkckIntermediateFunction(ARM_GP_Vector* yInOut, ARM_GP_Vector* dydx, int n, double* x, double htry, double eps,
							ARM_GP_Vector* yscal, double* hdid, double* hnext,
							const ARM_ODEFunc& function, vector<ARM_GP_Vector>& solverVars);

void odeint(ARM_GP_Vector* yInOut, int nvar, double x1, double x2, double eps, double tiny, double h1,
			double hmin, int kmax, double dxsav,
			const ARM_ODEFunc& function, vector<ARM_GP_Vector>& solverVars,
			int *nok, int *nbad, int *kount, int& nstp,
			ARM_GP_Vector* xp, ARM_GP_Vector* yp);

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
