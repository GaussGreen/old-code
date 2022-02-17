/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file bondanalytic.cpp
 *
 *  \brief file for all the bond analytics (inspired from bondmath.cpp)
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date May 2004
 */


#include "gpclosedforms/bondanalytic.h"

/// gpbase
#include "gpbase/functor.h"

/// gpnumlib
#include "gpnumlib/solver.h"
#include "gpnumlib/newtonraphson.h"
#include "gpnumlib/numfunction.h"

/// kernel
#include <glob/expt.h>   /// for constants

/// standard libraries
#include <cmath>


CC_BEGIN_NAMESPACE( ARM )

const double DEF_PRECISION = 1e-8;

const double PriceToYieldInit = 0.03;

class BondPriceToInverse : public ARM_GP::UnaryFunc<double,double>
{
public: 
		BondPriceToInverse(
			double coupon,
			double redemptionValue,
			int frequency,
			double settlToNext,
			double settlToMaturity,
			double accruingTime) :
		itsCoupon(coupon),
		itsRedemptionValue(redemptionValue),
		itsFrequency(frequency),
		itsSettlToNext(settlToNext),
		itsSettlToMaturity(settlToMaturity),
		itsAccruingTime(accruingTime)
		{
		};

		virtual double operator() (double yield) const 
		{
			return ARM_BondAnalytics::YieldToPrice(
				itsCoupon,
				itsRedemptionValue,
				yield,
				itsFrequency,
				itsSettlToNext,
				itsSettlToMaturity,
				itsAccruingTime);
		}
private:
	double itsCoupon;
	double itsRedemptionValue;
	int itsFrequency;
	double itsSettlToNext;
	double itsSettlToMaturity;
	double itsAccruingTime;
};

////////////////////////////////////////////////////
///	Class  : ARM_BondAnalytics
///	Routine: PriceToYield
///	Returns: double
///	Action : compute the yield of a bond based on its clean price
////////////////////////////////////////////////////
double ARM_BondAnalytics::PriceToYield( 
	double coupon,
    double redemptionValue,
    double cleanPrice,
    int frequency,
    double settlToNext,
    double settlToMaturity,
	double accruingTime )
{
	/// standard arm basis
	const double ARM_STD_BASIS = 100.0;

	double pricetoyield = 0.0;

	if ((cleanPrice > K_NEW_DOUBLE_TOL) && (cleanPrice < YieldToPrice(coupon,redemptionValue,0.0,frequency,settlToNext,settlToMaturity,accruingTime)-K_NEW_DOUBLE_TOL))
	{	
		BondPriceToInverse func(coupon,redemptionValue,frequency,settlToNext,settlToMaturity,accruingTime);

		UnaryFuncWithNumDerivative<double> funcWithDeriv(func);

		T_NewtonRaphsonSolver<UnaryFuncWithNumDerivative<double> > solver(funcWithDeriv,cleanPrice,DEF_PRECISION,DEF_PRECISION);

		///To initialize departure point
		solver.setInitialGuess(PriceToYieldInit);

		pricetoyield = solver.Solve();
	}
		return pricetoyield;
}


////////////////////////////////////////////////////
///	Class  : ARM_BondAnalytics
///	Routine: YieldToPrice
///	Returns: double
///	Action : compute the clean price of a bond
////////////////////////////////////////////////////

double ARM_BondAnalytics::YieldToPrice(
	double coupon,
    double redemptionValue,
    double yield,
    int frequency,
    double settlToNext,
    double settlToMaturity,
	double accruingTime )
{
	/// are we done?
	if ( settlToMaturity < -K_NEW_DOUBLE_TOL )
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, ARM_USERNAME + ": maturity is in the past!" );

	double couponPeriod = 1.0/(double)(frequency);
    settlToNext      *= (double)frequency;
    settlToMaturity  *= (double)frequency;
    coupon           *= couponPeriod;
    yield            *= couponPeriod;

    double factor = 1.0 + yield;
	if ( factor < K_NEW_DOUBLE_TOL  )
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, ARM_USERNAME + ": yield below -1.0 makes no sense!" );

	/// to compute the number of coupons, add 0197 (7 days) to ge exact period rounding
    int leftCouponNb	= (int)( settlToMaturity - settlToNext + 0.0197 ) + 1;
    double numerator, idenominator, dirtyPrice;

	///  Various cases:
	///		1) last coupon period (use simple interest discounting)
	///		2) very small yield
	///		3) large discounting factor
	///		4) normal case
	///		5) huge negative yield and long term.

	const double ONE_PERIOD = 1.0;

	///	1) last coupon period (use simple interest discounting)
    if ( settlToMaturity <= ONE_PERIOD )
	{
		if( settlToMaturity < K_NEW_DOUBLE_TOL )
		{
			dirtyPrice = redemptionValue;
		}
		else
		{
			dirtyPrice = (redemptionValue + coupon) / (1.0 + yield * settlToMaturity);
		}
	}
	///	2) very small yield
    else if( fabs(yield) < K_NEW_DOUBLE_TOL )
    {
        dirtyPrice = redemptionValue + coupon * (double) leftCouponNb;
    }

	///	3) large discounting factor
    else if ( (numerator= pow(factor, (double) leftCouponNb)) > K_HUGE_DOUBLE )
    {
        dirtyPrice = coupon * pow(factor,((double)leftCouponNb-settlToMaturity)) / yield;
    }
	
	///	4) normal case
    else if ( (idenominator = pow( factor, -settlToMaturity )) < K_HUGE_DOUBLE )
    {
        dirtyPrice = idenominator * (redemptionValue + coupon * (numerator-1.0) / yield);
    }
	
	///	5) huge negative yield and long term.
    else
    {
        dirtyPrice = K_HUGE_DOUBLE;
    }

	double accrued = coupon * accruingTime / couponPeriod;
	
	return dirtyPrice-accrued;
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/



