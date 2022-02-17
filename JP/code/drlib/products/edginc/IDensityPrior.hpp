//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : IDensityPrior.hpp
//
//   Description : Base class for probability density priors
//
//   Author      : Matthias Arnsdorf
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_IDENSITY_PRIOR_HPP
#define QLIB_IDENSITY_PRIOR_HPP

#include "edginc/Object.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE


class PRODUCTS_DLL IDensityPrior: public virtual IObject {

public:
    static CClassConstSP const TYPE;
    
	virtual ~IDensityPrior();


	/** Returns gridPoints (fractional losses) in [start,end] at which probability density is defined.
		Density is assumed piecewise linear between these points 
		grid points can depend on the date*/
	virtual DoubleArraySP getGridPoints(
		double start, 
		double end,
		const DateTime & date
		) = 0;

	/** Returns values of probability density at gridPoints and at given date. */
    virtual DoubleArraySP getDensity(
		const DoubleArray & gridPoints,	/** points at which to return density */
		const DateTime & date			/** date at which want density */
			) const = 0;
    
	/** set range on which density is defined */
	virtual void setRange(double lowBound, double highBound) = 0;

private:
    // For reflection
    static void load (CClassSP& clazz);
};

DECLARE(IDensityPrior);

DRLIB_END_NAMESPACE

#endif


