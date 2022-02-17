//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : EntropyExpLossInterpolator.hpp
//
//   Description : Interpolates expected losses using entropy minimisation
//
//   Author      : Matthias Arnsdorf
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_ENTROPY_EXP_LOSS_INTERPOLATOR_HPP
#define QLIB_ENTROPY_EXP_LOSS_INTERPOLATOR_HPP

#include "edginc/Object.hpp"
#include "edginc/IExpectedLossInterpolator.hpp"
#include "edginc/IDensityPrior.hpp"
#include "edginc/DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE


class PRODUCTS_DLL EntropyExpLossInterpolator:  public CObject,
							public virtual IExpectedLossInterpolator {

public:
    static CClassConstSP const TYPE;

	/** public constructor */
	EntropyExpLossInterpolator(
		const DoubleArray & strikes,	/** strikes at which to intepolate */
		IDensityPriorSP prior			/** density prior	*/
		);		

	/** Destructor */
    virtual ~EntropyExpLossInterpolator();

	/** Returns an expected loss surface.
		The surface interpolates/extrapolates the targetExpLosses */
    virtual ExpectedLossSurfaceSP getELSurface(
		const ExpectedLossSurface & targetExpLosses	/** Expected loss points to interpolate */
		) const;
    
	


private:

	/** minimum separation between tranches */
	static const double minSep; 

	/** separtion of strike in homogeneous part of default storage grid */
	static const double elStep;

	/** Private constructor */
	EntropyExpLossInterpolator();

    // For reflection
    static void load (CClassSP& clazz);

	static IObject* defaultEntropyExpLossInterpolator();

	// METHODS --------------------------------------

	/** create strike at which expected losses for final ExpectedLossSurface are defined  
		This includs 0, 1, baseStrikes and strikes  */
	DoubleArraySP createStorageStrikes(
		const DoubleArray& baseStrikes, // base strikes of target tranches
		double RR						// recovery
		) const;

	/** return constraint coefficient matrix needed for strike interpolation */
	DoubleArrayArraySP constraintCoeffMatrix(
		const DoubleArray& integrationGrid,
		const DoubleArray& baseStrikes
		) const;

	/** Routine to interpolate expected losses accross strikes using entropy minimisation */
	void interpolate(
		const DoubleArray& baseStrikes,			// (I) strikes for target base expected losses
		const DoubleArray& baseEL,				// (I) target base expected loses. 
		const DoubleArray& integrationGrid,		// (I) integration grid. Assume h_coeffs and prior are piecewise linear between points 
		const DoubleArray& prior,				// (I) prior probability density defined at integration grid points. Assumed piecewise linear inbetweem
		const DoubleArrayArray& h_coeffs,		// (I) constraint coefficient matrix returned by constraintCoeffMatrix()
		const DoubleArray& calculationGrid,		// (I) points at which we want output losses	
		DoubleArray& baseELout					// (O) output expected losses
		) const;

	// FIELDS -------------------------
	/** strikes at which to inteprolate */
	DoubleArraySP strikes;
	
	/** prior for density */
	IDensityPriorSP prior;

	/** scale factor to multiply all input expected losses 
	This can be useful for generating wide spread portfolios */
	double scale;

	/** first optimizer */
	OptimizerNDSP optimizer1; 
	/** second optimizer */
	OptimizerNDSP optimizer2;

	/** True if only want to use unmodified prior density */
	bool showPrior;

};

DECLARE(EntropyExpLossInterpolator);

DRLIB_END_NAMESPACE

#endif
