/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  Header for the   BlackSholes templates and related
 *
 *	\file vanilla_bs.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

#ifndef _GP_CF_VANILLA_BS_H
#define _GP_CF_VANILLA_BS_H

#include "gpbase/port.h"
#include <gpbase/argconv.h>
#include "basic_distributions.h"
#include "gpclosedforms/inverse.h"
#include <cmath>

CC_BEGIN_NAMESPACE( ARM )

double BS(double f, 
		  double k, 
		  double t, 
		  double v, 
		  int callput = K_CALL,
		  double* delta = NULL,
		  double* gamma	= NULL,
		  double* vega  = NULL,
		  double* theta = NULL);

double VanillaImpliedVol_BS (double f, double k, double t, double target, int callput, double* guess = NULL, bool* success = NULL);

class ARM_ImpliedVolBS
{
private:
	double	F;
	double	K;
	double	T;
	int		CallPut;
	
public:
	ARM_ImpliedVolBS(double f, double k, double t, int callPut) : 
	  F(f), K(k), T(t), CallPut(callPut)
	  {}

public:

	double	vol(double target, bool * success = NULL)
	{
		double best;
		RootFunc func(this);
		double res = brentSolve(func, target, 1e-8,2.,1e-10,100,0,&best);
		if(success != NULL) *success = fabs(res-best) < K_DOUBLE_TOL ? false : true;
		return res;
	}

	double vol(double target, double k, bool * success = NULL)
	{
		K = k;
		return vol(target, success);
	}

private:
	
	class RootFunc : public DoubleToDoubleFunc
	{
	private:
		ARM_ImpliedVolBS * _this;

	public:
		RootFunc(ARM_ImpliedVolBS * myClass) : _this(myClass)
		{}

	public:
		double operator()(double x) const
		{
			return BS(_this->F, _this->K, _this->T, x, _this->CallPut);
		}
	};

	friend class RootFunc;
};

double BSImpliedVol(double f, double k, double t, double target, int callPut, bool * success = NULL);

double VanillaOptionFromCall(double forward,double k,double call,int CallOrPutFlag);

double BlackSholes_Formula(	double forward,
							double totalvolatility,
							double bondprice,
							double strike,
							double CallPut);


double BlackSholes_Formula(	double forward,
							double volatility,
							double bondprice,
							double strike,
							double maturity,
							double CallPut);

/// version for the digital
double DigitalBlackSholes_Formula(	double forward,
							double totalvolatility,
							double bondprice,
							double strike,
							double CallPut);

double DigitalBlackSholesSmooth_Formula(	double forward,
											double totalvolatility_plus,
											double totalvolatility_moins,
											double bondprice,
											double strike,
											double CallPut,
											double spread);

double BlackSholes_Derivative_1(	double forward,			/// delta
									double totalvolatility,
									double bondprice,
									double strike,
									double CallPut);

double BlackSholes_Derivative_2(	double forward,			/// vega total
									double totalvolatility,
									double bondprice,
									double strike,
									double CallPut);

double BlackSholes_Derivative_3(	double forward,			/// bond price total
									double totalvolatility,
									double bondprice,
									double strike,
									double CallPut);

double BlackSholes_Derivative_4(	double forward,			// strike
									double totalvolatility,
									double bondprice,
									double strike,
									double CallPut);

double BlackSholes_Derivative_1_1(	double forward,			/// gamma
									double totalvolatility,
									double bondprice,
									double strike,
									double CallPut);

double Asymptotic_BlackSholesTimeValue(	double forward,
										double strike,
										double totalvolatility);

double Asymptotic_BlackSholesTimeValue_ImplicitVol(
										double forward,
										double strike,
										double optprice);

double 	BS_RatioOption(	double S1,double Mu1,double Sigma1,
					   double S2,double Mu2,double Sigma2,double Rho,double K,double T,int CallPut);


double 	BS_ProductOption(	double S1,double Mu1,double Sigma1,
						 double S2,double Mu2,double Sigma2,double Rho,double K,double T,int CallPut);

double 	BS_PayFloat_DigitalOption(	double S1,double Mu1,double Sigma1,
		double S2,double Mu2,double Sigma2,double Rho,double K,double T,int CallPut);

double BS_PayFloat_DigitalOptionSmooth_Formula(	double forward1,
												double totalvolatility1,
												double forward2,
												double totalvolatility2,
												double correl,
												double bondprice,
												double strike,
												double CallPut,
												double spread);
						 
						 CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
