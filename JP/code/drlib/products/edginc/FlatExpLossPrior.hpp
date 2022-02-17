//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : FlatExpLossPrior.hpp
//
//   Description : Flat expected loss prior. Derives of ITrancheQuoteInterpolator
//
//   Author      : Matthias Arnsdorf
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_FLAT_EXP_LOSS_PRIOR_HPP
#define QLIB_FLAT_EXP_LOSS_PRIOR_HPP

#include "edginc/ITrancheQuoteInterpolator.hpp"


DRLIB_BEGIN_NAMESPACE

/** Class to produce an expected loss surface consisting of a single, constant expected loss */
class PRODUCTS_DLL FlatExpLossPrior: public CObject,
						public virtual ITrancheQuoteInterpolator 
{

public:
    static CClassConstSP const TYPE;

	/** public constructor */
	FlatExpLossPrior(double expLoss);

	/** Destructor */
    virtual ~FlatExpLossPrior();

	/** Called immediately after object constructed */
    virtual void validatePop2Object();    

	/** Returns expected loss surface for set of strikes and dates */
    virtual ExpectedLossSurfaceSP getELSurface(
		const CDOQuotes & marketQuotes,			// market tranche quotes. Ignored for trivial prior 
		const ICDSParSpreads & indexSwapSpreads	// index swap spreads 
			) const;
		

private:
	/** private constructor */
	FlatExpLossPrior();

    // For reflection
    static void load (CClassSP& clazz);

	static IObject* defaultFlatExpLossPrior();

	/** default expected loss as percentage of tranche size*/
	static const double DEFAULT_LOSS;

	// FIELDS ---------------------------------------

	/** constant expected loss value */
	double expLoss;

};

DECLARE(FlatExpLossPrior);

DRLIB_END_NAMESPACE

#endif


