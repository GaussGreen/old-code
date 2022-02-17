//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : LHPADensityPrior.hpp
//
//   Description : Probability density prior using large homogeneous pool approximation
//					Currently supports Gaussian Copula
//   Author      : Matthias Arnsdorf
//
//   Date        : June 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_LHPA_DENSITY_PRIOR_HPP
#define QLIB_LHPA_DENSITY_PRIOR_HPP

#include "edginc/Object.hpp"
#include "edginc/IDensityPrior.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include "edginc/DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE


class LHPADensityPrior:  public CObject,
						public virtual IDensityPrior 
{

public:
    static CClassConstSP const TYPE;

	/** public constructor */
	LHPADensityPrior(
		double lowBound, 
		double highBound, 
		int numIntervals,
		double sigma,
		double floor);

	/** Destructor */
    virtual ~LHPADensityPrior();

	/** Called immediately after object constructed */
    void validatePop2Object();    


	/** Returns gridPoints (fractional losses) in [start,end] at which probability density is defined.
		Density is assumed piecewise linear between these points */
	DoubleArraySP getGridPoints(
		double start, 
		double end,
		const DateTime & date);

	/** calculate the grid points at date */
	void calcGridPoints(const DateTime & date);

	/** Returns values of probability density at gridPoints. */
    DoubleArraySP getDensity(
		const DoubleArray & gridPoints,	/** points at which to return density */
		const DateTime & date			/** date at which want density */
		) const;
    
	
	/** set range on which density is defined */
	void setRange(double lowBound, double highBound);



private:

	/** private constructor */
	LHPADensityPrior();

    // For reflection
    static void load (CClassSP& clazz);

	static IObject* defaultLHPADensityPrior();

	/** maximum value for any exponent */
	static const double maxExp;

	/** high bound for density */
	static const double huge;

	/** default numIntervals */
	static const int def_numIntervals;

	/** default sigma */
	static const double def_sigma;

	// FIELDS ---------------------------------------

	/** low bound of range on which density is defined */
	// this needs to be 0 for now
	double lowBound;
	/** high bound of range on which density is defined */
	double highBound;

	/** number of intervals between low and high bound */
	int numIntervals;

	/* standard deviation for grid determination */
	double sigma;

	/** expiries for beta term structure */
	DateTimeArraySP betaDates;

	/** betas */
	DoubleArraySP betas;

	/** density floor */
	double floor;


	/**index swap spreads */
	ICDSParSpreadsWrapper index;

	/** market data */
	MarketDataSP market;

	// Transient fields

	// At any given time we store one integration grid per date
	// the idea is that the gird is calculated when getGridPoints is called
	// and is then available for a getDensity call at the same date

	/** discretisation grid cache*/
	DoubleArraySP grid;
	/** date at which grid is defined */
	DateTimeSP gridDate;

	/** homogeneous grid of cumulative density points */
	DoubleArraySP cummDensity;
};

DECLARE(LHPADensityPrior);

DRLIB_END_NAMESPACE

#endif
