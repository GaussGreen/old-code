/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: QModelAnalytics.cpp,v $
 * Revision 1.1  2004/09/22 10:15:09  ebenhamou
 * Initial revision
 *
 *
 *
 */


/*! \file QModelAnalytics.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */

#include "gpmodels/QModelAnalytics.h"

/// gpbase
#include "gpbase/numericconstant.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/vanilla_normal.h"
#include "gpclosedforms/normal.h"
#include <cmath>

/// gpnumlib
#include "gpnumlib/gaussiananalytics.h"
#include "gpnumlib/newtonraphson.h"
#include "gpnumlib/normalinvcum.h"


CC_BEGIN_NAMESPACE( ARM )




////////////////////////////////////////////////////
///	Static routine to price with a Qmodel 1F
////////////////////////////////////////////////////
double QModelAnalytics::BSFormulaForQModel(double forward,
							double totalvolatility,
							double bondprice,
							double strike,
							double CallPut)
{
    /// Same function as BlackSholes_Formula() but call to NormalCDF2
    /// far better near 0 than cdfnormal()
	if(strike<=K_DOUBLE_TOL)
		return CallPut == K_CALL? bondprice*(forward-strike): 0.0;

	if(totalvolatility<=K_DOUBLE_TOL)
	{
		double intrinsicValue = CallPut*bondprice*(forward-strike);
		return (intrinsicValue>K_DOUBLE_TOL) ? intrinsicValue: 0.0;
	}

	/// Standard Black Scholes case
    double limit=0.001;
    double nd1,d1 = (log(forward/strike))/totalvolatility+0.5*totalvolatility ;
    double nd2,d2 = d1-totalvolatility;

    if(d1>limit || d1<-limit)
        nd1 = ARM_GaussianAnalytics::cdfNormal(CallPut*d1);
    else
        nd1 = RegularNormalCDF2(CallPut*d1);
    if(d2>limit || d2<-limit)
        nd2 = ARM_GaussianAnalytics::cdfNormal(CallPut*d2);
    else
        nd2 = RegularNormalCDF2(CallPut*d2);

    double value = CallPut*bondprice*(forward*nd1-strike*nd2);

    return value;
}


double QModelAnalytics::BSQFunction( 
		double forward, 
		double strike, 
		double qSigma,
		double maturity,
		double qParameter,
		double df,
		int callOrPut,
        double qForward,
        double* vega)
{
    /// The underlying is shifted log with an absolute shift m = qForward*(1.0/qParameter-1)

    double price;
    double var,k;
    int cp;

    double m = qForward*(1.0/qParameter-1);
    double fwd = forward + m;
    bool isNotNullQ = (fabs(qParameter)>1.0e-10);

	if( isNotNullQ && fwd > 0.0 )
    {
        k   = strike + m;
        var = qSigma*sqrt(maturity)*fabs(qParameter);
		price = BSFormulaForQModel(fwd,var,df,k,callOrPut);
        if(vega)
            *vega = BlackSholes_Derivative_2(fwd,var,df,k,callOrPut);
    }
	else if( isNotNullQ && fwd < 0.0)
    {
        fwd = -fwd;
        k   = -(strike+m);
        var = qSigma*sqrt(maturity)*fabs(qParameter);
        cp  = callOrPut==K_CALL ? K_PUT : K_CALL;
        price = BSFormulaForQModel(fwd,var,df,k,cp);
        if(vega)
            *vega = BlackSholes_Derivative_2(fwd,var,df,k,cp);
    }
	else if(!isNotNullQ)
    {
        /// Q parameter is very closed to 0 => normal pricing formula
        double stdDev = qSigma*forward;
		price = df*VanillaOption_N(forward,stdDev,strike,maturity,callOrPut);
        if(vega)
            *vega = df*VegaVanillaOption_N(forward,stdDev,strike,maturity,callOrPut);
    }
    else
    {
        /// fwd=0.0 !
        price = (strike+m)*callOrPut < 0 ? callOrPut*df*(forward-strike) : 0.0;
        if(vega)
            *vega=0.0;
    }
    return price;
}



////////////////////////////////////////////////////
///	Static routine to price with a Qmodel 1F
////////////////////////////////////////////////////
double QModelAnalytics::DigitalBSQFunction( 
		double forward, 
		double strike, 
		double qSigma,
		double maturity,
		double qParameter,
		double df,
		int callOrPut,
        double qForward)
{
    /// The underlying is shifted log with an absolute shift m = qForward*(1.0/qParameter-1)

    double m = qForward*(1.0/qParameter-1);
    double fwd = forward + m;

	if( fwd > 0.0 )
		return DigitalBlackSholes_Formula(fwd, qSigma*sqrt(maturity)*fabs(qParameter), df, strike+m, callOrPut);

	else if( fwd < 0.0 )
        return DigitalBlackSholes_Formula(-fwd, qSigma*sqrt(maturity)*fabs(qParameter), df, -(strike+m), callOrPut==K_CALL ? K_PUT : K_CALL);

	else
        /// Q parameter is very closed to 0 => normal pricing formula
		return df*VanillaDigitalOption_N(forward,qSigma*forward,strike,maturity,callOrPut);
}


double QModelAnalytics::NormalVolToQVol(
		double forward, 
        double norSigma,
        double maturity,
        double qParameter)
{
    /// Compute ATM normal price
    double sqrtMat = sqrt(maturity);

	///Price ATM
    double target=norSigma*sqrtMat*ARM_NumericConstants::ARM_INVSQRT2PI;

    /// Build & call a NR solver
    QPriceFunct qPriceFct(forward,forward,maturity,qParameter,forward);
    UnaryFuncWithNumDerivative<double,double> fctToSolve( qPriceFct );

    double qVol0;
	if( fabs(qParameter) > K_NEW_DOUBLE_TOL )
		qVol0 = ARM_NormalInvCum::Inverse_erf_Moro(0.5*(target*qParameter/forward+1.0))
		*2.0/(qParameter*sqrtMat);
	else
		qVol0 = target*ARM_NumericConstants::ARM_SQRT_2_PI/sqrtMat/forward;

    /// Create a NR solver then find the solution
	const double ftol = 1.0e-8;
    T_NewtonRaphsonSolver< UnaryFuncWithNumDerivative<double,double> > solver(fctToSolve,target,ftol);

	///To initialize departure point
	solver.setInitialGuess(qVol0);

    double qVol=solver.Solve();

    return qVol;
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

