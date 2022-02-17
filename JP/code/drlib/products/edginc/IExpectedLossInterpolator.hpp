//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : IExpectedLossInterpolator.hpp
//
//   Description : Interface for expected loss interpolators
//
//   Author      : Matthias Arnsdorf
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_IEXPECTED_LOSS_INTERPOLATOR_HPP
#define QLIB_IEXPECTED_LOSS_INTERPOLATOR_HPP

#include "edginc/Object.hpp"
#include "edginc/ExpectedLossSurface.hpp"
#include "edginc/DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE


class PRODUCTS_DLL IExpectedLossInterpolator: public virtual IObject {

public:
    static CClassConstSP const TYPE;
    
	/** Destructor */
	virtual ~IExpectedLossInterpolator();

	/** Returns an expected loss surface that interpolates/extrapolates the targetExpLosses */
    virtual ExpectedLossSurfaceSP getELSurface(
		const ExpectedLossSurface & targetExpLosses	/** Expected loss points to interpolate */
		) const = 0;
    

private:
    // For reflection
    static void load (CClassSP& clazz);
};

DECLARE(IExpectedLossInterpolator);

DRLIB_END_NAMESPACE

#endif


