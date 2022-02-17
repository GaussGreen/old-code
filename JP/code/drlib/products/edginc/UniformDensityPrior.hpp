//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : UniformDensityPrior.hpp
//
//   Description : Uniform probability density prior
//
//   Author      : Matthias Arnsdorf
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_UNIFORM_DENSITY_PRIOR_HPP
#define QLIB_UNIFORM_DENSITY_PRIOR_HPP

#include "edginc/Object.hpp"
#include "edginc/IDensityPrior.hpp"
#include "edginc/DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE


class PRODUCTS_DLL UniformDensityPrior:  public CObject,
							public virtual IDensityPrior 
{

public:
    static CClassConstSP const TYPE;

	/** public constructor */
	UniformDensityPrior(double lowBound, double highBound);

	/** Destructor */
    virtual ~UniformDensityPrior();

	/** Called immediately after object constructed */
    virtual void validatePop2Object();    


	/** Returns gridPoints (fractional losses) in [start,end] at which probability density is defined.
		Density is assumed piecewise linear between these points */
	virtual DoubleArraySP getGridPoints(
		double start, 
		double end,
		const DateTime & date);

	/** Returns values of probability density at gridPoints. */
    virtual DoubleArraySP getDensity(
		const DoubleArray & gridPoints,	/** points at which to return density */
		const DateTime & date			/** date at which want density */
		) const;
    
	
	/** set range on which density is defined */
	virtual void setRange(double lowBound, double highBound);



private:

	/** private constructor */
	UniformDensityPrior();

    // For reflection
    static void load (CClassSP& clazz);

	static IObject* defaultUniformDensityPrior();


	// FIELDS ---------------------------------------

	/** low bound of range on which density is defined */
	double lowBound;
	/** high bound of range on which density is defined */
	double highBound;

};

DECLARE(UniformDensityPrior);

DRLIB_END_NAMESPACE

#endif
