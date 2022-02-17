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

#include "expt.h"
#include <cmath>
#include <iomanip> /// for setprecision()
using namespace std;

CC_BEGIN_NAMESPACE( ARM )

extern const int RK4_NBTMP_VAR;
extern const int RK5_NBTMP_VAR;

void Rk4ToNextStep(std::vector<double>* yInOut, std::vector<double>* dydx, int n, double x, double h,
			const ARM_ODEFunc& function, vector<std::vector<double>>& solverVars);

void RkFunction(std::vector<double>* yInOut, int nvar, double x1,double x2, int nstep, std::vector<double>* xp,
			const ARM_ODEFunc& function, vector<std::vector<double>>& solverVars);

void rkckToNextStep(std::vector<double>* yInOut, std::vector<double>* dydx, int n, double x, double h,
			std::vector<double>* yerr, const ARM_ODEFunc& function, vector<std::vector<double>>& solverVars);

void RkckIntermediateFunction(std::vector<double>* yInOut, std::vector<double>* dydx, int n, double* x, double htry, double eps,
							std::vector<double>* yscal, double* hdid, double* hnext,
							const ARM_ODEFunc& function, vector<std::vector<double>>& solverVars);

void odeint(std::vector<double>* yInOut, int nvar, double x1, double x2, double eps, double tiny, double h1,
			double hmin, int kmax, double dxsav,
			const ARM_ODEFunc& function, vector<std::vector<double>>& solverVars,
			int *nok, int *nbad, int *kount, int& nstp,
			std::vector<double>* xp, std::vector<double>* yp);

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
