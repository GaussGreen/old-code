//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : IExpectedLossPrior.hpp
//
//   Description : Base class for expected loss priors
//
//   Author      : Matthias Arnsdorf
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_IEXPECTED_LOSS_PRIOR_HPP
#define QLIB_IEXPECTED_LOSS_PRIOR_HPP

#include "edginc/Object.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/ExpectedLossSurface.hpp"
#include "edginc/CDOQuotes.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include "edginc/DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE


class PRODUCTS_DLL IExpectedLossPrior: public virtual IObject {

public:
    static CClassConstSP const TYPE;
    
	virtual ~IExpectedLossPrior();


	
	/** Returns expected loss surface for set of strikes and dates */
    virtual ExpectedLossSurfaceSP getELSurface(
		const CDOQuotes & marketQuotes,			/** market tranche quotes */
		const ICDSParSpreads & indexSwapSpreads,/** index swap spreads */
		const DoubleArray & strikes,			/** strikes for EL surface (needs to include 0 and 1) */
		const DateTimeArray & dates				/** dates for EL surface */	
		) const = 0;
    

private:
    // For reflection
    static void load (CClassSP& clazz);
};

DECLARE(IExpectedLossPrior);

DRLIB_END_NAMESPACE

#endif
