//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : ImpliedLossModel.hpp
//
//   Description : Gives index prices of off-market tranches by interpolating set of CDOQuotes
//
//   Author      : Matthias Arnsdorf
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_IMPLIED_LOSS_MODEL_HPP
#define QLIB_IMPLIED_LOSS_MODEL_HPP

#include "edginc/CDOQuotes.hpp"
#include "edginc/ITrancheQuoteInterpolator.hpp"
#include "edginc/IExpectedLossInterpolator.hpp"


DRLIB_BEGIN_NAMESPACE


class PRODUCTS_DLL ImpliedLossModel: public CObject
{

public:
    static CClassConstSP const TYPE;


	/** Destructor */
    virtual ~ImpliedLossModel();

	

	// possible output types
	static const string SPREAD;		// par spread  
	static const string ANNUITY;	// anuuity / duration  
	static const string CONT_LEG;	// contingent leg  
	static const string UPFRONT;	//  upfront  
	static const string LOSS;		// loss	 
	static const string RISKY_ZERO;	// risky zero	
	
	// default running coupon
	static const double DEFAULT_COUPON;


	/** return array of outputType for set of instruments defined by 4-tuple {lowStrike, highStrike, expiry, outputType} */
	// if outputType = UPRONT then default coupon of 500bps is assumed
	DoubleArraySP getValues(
		const CDOQuotes & marketQuotes,		// set of market tranche qutoes  
		const ICDSParSpreads & index,		// index swap spreads			 
		const ExpiryArray & expiries,		// set of expiries for output   
		const DoubleArray & lowStrikes,		// set of low strikes for output  
		const DoubleArray & highStrikes,	// set of high strikes for output   
		const StringArray &	outputType,     // output type  
        IForwardRatePricerSP model          /** for calculating fees */
		) const;

	// Overload to ake in coupons
	/** return array of outputType for set of instruments defined by 5-tuple {lowStrike, highStrike, expiry, outputType, coupon} */
	DoubleArraySP getValues(
		const CDOQuotes & marketQuotes,		// set of market tranche qutoes  
		const ICDSParSpreads & index,		// index swap spreads			 
		const ExpiryArray & expiries,		// set of expiries for output   
		const DoubleArray & lowStrikes,		// set of low strikes for output  
		const DoubleArray & highStrikes,	// set of high strikes for output   
		const StringArray &	outputType,		// output type  
		const DoubleArray & coupons,        // coupons needed if outputType = upfront  	
        IForwardRatePricerSP model          /** for calculating fees */
		) const;


protected:
	
	/** private constructor */
	ImpliedLossModel();

private:
    // For reflection
    static void load (CClassSP& clazz);

	static IObject* defaultImpliedLossModel();


	// FIELDS ---------------------------------------

	/** tranche interpolation model */
	ITrancheQuoteInterpolatorSP trancheQuoteInterpolator;

	/** expected loss interpolation model [optional] */
	IExpectedLossInterpolatorSP expLossInterpolator;

};

DECLARE(ImpliedLossModel);

DRLIB_END_NAMESPACE

#endif
