/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 * \file rungekutta.cpp
 *
 *  \brief
 *
 *	\author  A. Triki
 *	\version 1.0
 *	\date October 2005
 */

#include "gpnumlib/rungekutta.h"
#include "gpbase/gpvector.h"
#include "gpbase/utilityport.h"
#include <math.h>

CC_BEGIN_NAMESPACE( ARM )

inline double SIGN(const double &a, const double &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}


const double SAFETY = 0.9;
const double PGROW  = -0.2;
const double PSHRNK = -0.25;
const double ERRCON = 1.89e-4;
const double MAXSTP = 10000;
const double TINY	= 1.0e-30;

/// Internal variables aliases for RK4
const int RK4_DYM		= 0;
const int RK4_DYT		= 1;
const int RK4_YT		= 2;
const int RK4_DXDY		= 3;
const int RK4_NBTMP_VAR	= 4;

/// Internal variables aliases for RK5
const int RK5_AK2		= 0;
const int RK5_AK3		= 1;
const int RK5_AK4		= 2;
const int RK5_AK5		= 3;
const int RK5_AK6		= 4;
const int RK5_YTEMP		= 5;
const int RK5_YOUT		= 6;
const int RK5_YERR		= 7;
const int RK5_YSCAL		= 8;
const int RK5_DXDY		= 9;
const int RK5_NBTMP_VAR	= 10;

///The value ERRCON equals (5/SAFETY) raised to the power (1/PGROW), see use below.

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
///	Function    : Rk4ToNextStep
/// Description : Standard 4 Step Runge Kutta
///					
///	Action		: Given values for the variables y[1..n] and their derivatives dydx[1..n] known at x, use the
///				fourth-order Runge-Kutta method to advance the solution over an interval h and return the
///				incremented variables as yout[1..n], which need not be a distinct array from y. The user
///				supplies the routine derivs(x,y,dydx), which returns derivatives dydx at x.
///				The next step is stored in yout
////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Rk4ToNextStep(std::vector<double>* yInOut, std::vector<double>* dydx, int n, double x, double h,
			const ARM_ODEFunc& function, vector<std::vector<double>>& solverVars)
{
	int i;
	double xh,hh,h6;
	std::vector<double>* dym = & solverVars[RK4_DYM];
	std::vector<double>* dyt = & solverVars[RK4_DYT];
	std::vector<double>* yt  = & solverVars[RK4_YT];
				
	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;

	// First step
	for (i=0;i<n;i++)
		(*yt)[i]=(*yInOut)[i]+hh*((*dydx)[i]);

	// Second step
	function.derivs(xh,yt,dyt);
	for (i=0;i<n;i++) 
		(*yt)[i]=(*yInOut)[i]+hh*(*dyt)[i];

	// Third step
	function.derivs(xh,yt,dym);
	for (i=0;i<n;i++) 
	{
		(*yt)[i]	=(*yInOut)[i]+h*(*dym)[i];
		(*dym)[i]  += (*dyt)[i];
	}
	
	// Fourth step : accumulate increments with proper weights
	function.derivs(x+h,yt,dyt);
	for (i=0;i<n;i++)
		(*yInOut)[i] += h6*((*dydx)[i]+(*dyt)[i]+2.0*(*dym)[i]);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
///	Function    : RkFunction
///	Description : Simple function without Step Size Control
///					
///	Action		: Starting from initial values vstart[1..nvar] known at x1 use fourth-order Runge-Kutta
///				to advance nstep equal increments to x2. The user-supplied routine derivs(x,v,dvdx)
///				evaluates derivatives. Results are stored in the global variables y[1..nvar][1..nstep+1]
///				and xx[1..nstep+1].
///				
///				WARNING :: If we store intermediate values y size = nvar * (nbStep+1) else y size = nvar
///							y1t1,y2t1,...ynt1,y1t2,y2t2.....
///				Vstart =y0 ; nvar = number of variables yi; x1 = xstart, x2 = xend,nbstep, y (output), xx(output) 
////////////////////////////////////////////////////////////////////////////////////////////////////////////

void RkFunction(std::vector<double>* yInOut, int nvar, double x1,double x2, int nstep, std::vector<double>* xp,
			const ARM_ODEFunc& function, vector<std::vector<double>>& solverVars)

{
	int k;
	double x,h;
	std::vector<double>* dxdy	= & solverVars[RK4_DXDY];;

	(*xp)[0]=x1;

	x=x1;
	h=(x2-x1)/nstep;
	for (k=0;k<nstep;k++) 
	{ 
		function.derivs(x,yInOut,dxdy);
		Rk4ToNextStep(yInOut,dxdy,nvar,x,h,function,solverVars);
		if ((float)(x+h) == x) 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"Step size too small in routine RkFunction");
		x += h;

		(*xp)[k+1]=x; //Store intermediate steps.
	}
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////
///	Function    : rkckToNextStep
/// Description : fth-order Cash-Karp Runge-Kutta method
///					
///	Action		: Given values for n variables y[1..n] and their derivatives dydx[1..n] known at x, use
///				the fth-order Cash-Karp Runge-Kutta method to advance the solution over an interval h
///				and return the incremented variables as yout[1..n]. Also return an estimate of the local
///				truncation error in yout using the embedded fourth-order method. The user supplies the routine
///				derivs(x,y,dydx), which returns derivatives dydx at x.
////////////////////////////////////////////////////////////////////////////////////////////////////////////
const double a2	=	0.2;
const double a3	=	0.3;
const double a4	=	0.6;
const double a5	=	1.0;
const double a6	=	0.875;
const double b21=	0.2;
const double b31=	3.0/40.0;
const double b32=	9.0/40.0;
const double b41=	0.3;
const double b42=	-0.9;
const double b43=	1.2;
const double b51=	-11.0/54.0;
const double b52=	2.5;
const double b53=	-70.0/27.0;
const double b54=	35.0/27.0;
const double b61=	1631.0/55296.0;
const double b62=	175.0/512.0;
const double b63=	575.0/13824.0;
const double b64=	44275.0/110592.0;
const double b65=	253.0/4096.0;
const double c1	=	37.0/378.0;
const double c3	=	250.0/621.0;
const double c4	=	125.0/594.0;
const double c6	=	512.0/1771.0;
const double dc1=	c1-2825.0/27648.0;
const double dc3=	c3-18575.0/48384.0;
const double dc4=	c4-13525.0/55296.0;
const double dc5=	-277.00/14336.0;
const double dc6=	c6-0.25;

/// Internal variables aliases for RK5

void rkckToNextStep(std::vector<double>* yIn, std::vector<double>* dydx, int n, double x, double h,
			std::vector<double>* yOut, std::vector<double>* yerr, const ARM_ODEFunc& function, vector<std::vector<double>>& solverVars)
{
	int i;
		
	std::vector<double>* ak2		= & solverVars[RK5_AK2];
	std::vector<double>* ak3		= & solverVars[RK5_AK3];
	std::vector<double>* ak4		= & solverVars[RK5_AK4];
	std::vector<double>* ak5		= & solverVars[RK5_AK5];
	std::vector<double>* ak6		= & solverVars[RK5_AK6];
	std::vector<double>* ytemp	= & solverVars[RK5_YTEMP];

	// First step
	for (i=0;i<n;i++)
		(*ytemp)[i]=(*yIn)[i]+b21*h*(*dydx)[i];
	
	// Second step
	function.derivs(x+a2*h,ytemp,ak2);
	for (i=0;i<n;i++)
		(*ytemp)[i]=(*yIn)[i]+h*(b31*(*dydx)[i]+b32*(*ak2)[i]);
	
	// Third step
	function.derivs(x+a3*h,ytemp,ak3); 
	for (i=0;i<n;i++)	
		(*ytemp)[i]=(*yIn)[i]+h*(b41*(*dydx)[i]+b42*(*ak2)[i]+b43*(*ak3)[i]);
	
	// Fourth step
	function.derivs(x+a4*h,ytemp,ak4);
	for (i=0;i<n;i++)
		(*ytemp)[i]=(*yIn)[i]+h*(b51*(*dydx)[i]+b52*(*ak2)[i]+b53*(*ak3)[i]+b54*(*ak4)[i]);
	
	// Fifth step
	function.derivs(x+a5*h,ytemp,ak5);
	for (i=0;i<n;i++)
		(*ytemp)[i]=(*yIn)[i]+h*(b61*(*dydx)[i]+b62*(*ak2)[i]+b63*(*ak3)[i]+b64*(*ak4)[i]+b65*(*ak5)[i]);
	
	// Sixth step : accumulate increments with proper weights
	function.derivs(x+a6*h,ytemp,ak6);
	for (i=0;i<n;i++)
		(*yOut)[i] = (*yIn)[i]+h*(c1*(*dydx)[i]+c3*(*ak3)[i]+c4*(*ak4)[i]+c6*(*ak6)[i]);
	
	// Estimate error as diference between 4th and 5th order methods
	for (i=0;i<n;i++)
		(*yerr)[i]=h*(dc1*(*dydx)[i]+dc3*(*ak3)[i]+dc4*(*ak4)[i]+dc5*(*ak5)[i]+dc6*(*ak6)[i]);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
///	Function    : RkckIntermediateFunction
///	Description : fth-order Cash-Karp Runge-Kutta method with error control
///					
///	Action		: Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure accuracya nd
///				adjust stepsize. Input are the dependent variable vector y[1..n] and its derivative dydx[1..n]
///				at the starting value of the independent variable x. Also input are the stepsize to be attempted
///				htry, the required accuracy eps, and the vector yscal[1..n] against which the error is
///				scaled. On output, y and x are replaced bythei r new values, hdid is the stepsize that was
///				actuallyac complished, and hnext is the estimated next stepsize. derivs is the user-supplied
///				routine that computes the right-hand side derivatives.
////////////////////////////////////////////////////////////////////////////////////////////////////////////

void RkckIntermediateFunction(std::vector<double>* yInOut, std::vector<double>* dydx, int n, double* x, double htry, double eps,
							std::vector<double>* yscal, double* hdid, double* hnext,
							const ARM_ODEFunc& function, vector<std::vector<double>>& solverVars)
{
	int i;
	double  errmax,h,htemp;
	std::vector<double>* yout  = & solverVars[RK5_YOUT];
	std::vector<double>* yerr  = & solverVars[RK5_YERR];

	h=htry; //Set stepsize to the initial trial value.
	for (;;) 
	{
		rkckToNextStep(yInOut,dydx,n,*x,h,yout,yerr,function,solverVars); //Take a step.
		errmax=0.0; //Evaluate accuracy.
		for (i=0;i<n;i++) 
			errmax=CC_Max(errmax,fabs((*yerr)[i]/(*yscal)[i]));
		errmax /= eps; //Scale relative to required tolerance.
		if (errmax <= 1.0) break; //Step succeeded. Compute size of next step.
		
		htemp=SAFETY*h*pow(errmax,PSHRNK);
		//Truncation error too large, reduce stepsize.
		h=(h >= 0.0 ? CC_Max(htemp,0.1*h) : CC_Min(htemp,0.1*h));
		//No more than a factor of 10.
		if (h==0.0) 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"stepsize underflow in rkqs");
	}
	if (errmax > ERRCON) 
		*hnext=SAFETY*h*pow(errmax,PGROW);
	else 
		*hnext=5.0*h; //No more than a factor of 5 increase.
	*x += (*hdid=h);

	for (i=0;i<n;i++)
		(*yInOut)[i] = (*yout)[i];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
///	Function    : odeint
///	Description : fth-order Cash-Karp Runge-Kutta method with error control
///					
///	Action		: Runge-Kutta driver with adaptive stepsize control. Integrate starting values ystart[1..nvar]
///				from x1 to x2 with accuracy eps, storing intermediate results in global variables. h1 should
///				be set as a guessed .rst stepsize, hmin as the minimum allowed stepsize (can be zero). On
///				output nok and nbad are the number of good and bad (but retried and .xed) steps taken, and
///				ystart is replaced byv alues at the end of the integration interval. derivs is the user-supplied
///				routine for calculating the right-hand side derivative, while rkqs is the name of the stepper
///				routine to be used.
////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*User storage for intermediate results. Preset kmax and dxsav in the calling program. If kmax .=
0 results are stored at approximate intervals dxsav in the arrays xp[1..kount], yp[1..nvar]
[1..kount], where kount is output by odeint. De.ning declarations for these variables, with
memoryallo cations xp[1..kmax] and yp[1..nvar][1..kmax] for the arrays, should be in
the calling program.*/

//x1 = xstart, x2 = xend, eps = error, h1 = first guess step, hmin = minimum allawed stepsize (0),
//  nok(output), nbad(output),kmax (maximum storage), kount (output), yp (output), dxsav (approximate intervals if no storage)
//  nstp(nbre of steps used except internal tries)
void odeint(std::vector<double>* yInOut, int nvar, double x1, double x2, double eps, double tiny, double h1,
			double hmin, int kmax, double dxsav,
			const ARM_ODEFunc& function, vector<std::vector<double>>& solverVars,
			int *nok, int *nbad, int *kount, int& nstp,
			std::vector<double>* xp, std::vector<double>* yp)
{
	int i;
	double xsav,x,hnext,hdid,h;

	std::vector<double>* yscal	= & solverVars[RK5_YSCAL];
	std::vector<double>* dydx		= & solverVars[RK5_DXDY];

	x=x1;
	h=SIGN(h1,x2-x1);
	*nok = (*nbad) = (*kount) = -1;
	
	if (kmax >= 0) 
		xsav=x-dxsav*2.0; //Assures storage of .rst step.
	
	for (nstp=0;nstp<MAXSTP;nstp++) 
	{ 
		//Take at most MAXSTP steps.
		function.derivs(x,yInOut,dydx);
		for (i=0;i<nvar;i++)
			//Scaling used to monitor accuracy. This general-purpose choice can be modi.ed if need be.
			(*yscal)[i]=fabs((*yInOut)[i])+fabs((*dydx)[i]*h)+tiny;
		if (kmax > 0 && (*kount) < kmax-1 && fabs(x-xsav) > fabs(dxsav)) 
		{
			(*xp)[++(*kount)]=x; //Store intermediate results.
			if(yp)
			{
				for (i=0;i<nvar;i++) 
					(*yp)[i + nvar * (*kount)] = (*yInOut)[i];
			}
			xsav=x;
		}
		if ((x+h-x2)*(x+h-x1) > 0.0) 
			h=x2-x; //If stepsize can overshoot, decrease.
	
		RkckIntermediateFunction (yInOut,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,function,solverVars);
		if (hdid == h) 
			++(*nok);
		else 
			++(*nbad);
		
		if ((x-x2)*(x2-x1) >= 0.0) 
		{
			//Are we done?
			if (kmax) 
			{
				(*xp)[++(*kount)]=x; //Save .nal step.
				if(yp)
				{
					for (i=0;i<nvar;i++)
						(*yp)[i + nvar * (*kount)]=(*yInOut)[i];
				}
			}
			return; //Normal exit.
		}
		if (fabs(hnext) <= hmin) 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					"Step size too small in odeint");
		h=hnext;
	}
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					"Too many steps in routine odeint");

	
}

CC_END_NAMESPACE()
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

