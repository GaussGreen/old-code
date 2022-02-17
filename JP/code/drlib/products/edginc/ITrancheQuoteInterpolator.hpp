//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : ITrancheQuoteInterpolator.hpp
//
//   Description : Interface for classes that can interpolate tranche quotes and produce an ExpectedLossSurface
//
//   Author      : Matthias Arnsdorf
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_ITRANCHE_QUOTE_INTERPOLATOR_HPP
#define QLIB_ITRANCHE_QUOTE_INTERPOLATOR_HPP

#include "edginc/ExpectedLossSurface.hpp"
#include "edginc/CDOQuotes.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include "edginc/DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE


class PRODUCTS_DLL ITrancheQuoteInterpolator: public virtual IObject {

public:
    static CClassConstSP const TYPE;
    
	virtual ~ITrancheQuoteInterpolator();


	/** Returns expected loss surface for set of strikes and dates
		This surface should reprice the input trancheQuotes */
    virtual ExpectedLossSurfaceSP getELSurface(
		const CDOQuotes & quotes,					/** market tranche quotes */
		const ICDSParSpreads & indexSwapSpreads		/** index swap spreds */
		) const  = 0;
    

private:
    // For reflection
    static void load (CClassSP& clazz);
};

DECLARE(ITrancheQuoteInterpolator);

DRLIB_END_NAMESPACE

#endif

