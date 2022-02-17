//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : ExpLossPrior.hpp
//
//   Description : Custom expected loss prior. Derives of ITrancheQuoteInterpolator
//
//   Author      : Matthias Arnsdorf
//
//   Date        : June 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_EXP_LOSS_PRIOR_HPP
#define QLIB_EXP_LOSS_PRIOR_HPP

#include "edginc/ITrancheQuoteInterpolator.hpp"


DRLIB_BEGIN_NAMESPACE

/** Class to produce an expected loss surface consisting of a custom expected loss surface */
class PRODUCTS_DLL ExpLossPrior: public CObject,
	public virtual ITrancheQuoteInterpolator 
{

public:
	static CClassConstSP const TYPE;


	/** Destructor */
	virtual ~ExpLossPrior();

	/** Called immediately after object constructed */
	virtual void validatePop2Object();    

	/** Returns expected loss surface for set of strikes and dates */
	virtual ExpectedLossSurfaceSP getELSurface(
		const CDOQuotes & marketQuotes,			// market tranche quotes. Ignored for trivial prior 
		const ICDSParSpreads & indexSwapSpreads	// index swap spreads 
		) const;


private:
	/** private constructor */
	ExpLossPrior();

	// For reflection
	static void load (CClassSP& clazz);

	static IObject* defaultExpLossPrior();

	

	// FIELDS ---------------------------------------

	/** dates */
	ExpiryArraySP expiries;

	/** strikes */
	DoubleArraySP strikes;

	/** losses */
	DoubleArrayArraySP baseLosses;

};

DECLARE(ExpLossPrior);

DRLIB_END_NAMESPACE

#endif


