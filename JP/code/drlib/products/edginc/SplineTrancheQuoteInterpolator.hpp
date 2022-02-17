//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : SplineTrancheQuoteInterpolator.hpp
//
//   Description : Interpolation of CDOQuotes using mulit-spline approach
//
//   Author      : Matthias Arnsdorf
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_SPLINE_TRANCHE_QUOTE_INTERPOLATOR_HPP
#define QLIB_SPLINE_TRANCHE_QUOTE_INTERPOLATOR_HPP


#include "edginc/ITrancheQuoteInterpolator.hpp"
#include "edginc/CreditFeeLegWithPV.hpp"
#include "edginc/ConstrainedCubicSpline.hpp"


DRLIB_BEGIN_NAMESPACE

/** Class to interpolate CDO tranche quotes. 

	Inputs are a CDO quotes object, an ICDSparSpreads object
	representing the index swap par spreads and a set of model parameters.

	The ouput is an expected loss surface defined at the chosen time interpolation dates and a set of base 
	strikes which are constructed as follows:
	- include all strikes EXCEPT 0 AND <= 1-RR which define the tranche quotes. 
	- include 1-RR and 1.

	A constant recovery rate, RR, is assumed, which means that there can be no loss above 1-RR 
	(Given a unit portfolio notional).

	The expected loss surface is constructed such that the input tranche quotes are repriced. This assumes
	LINEAR interpolation of the loss surface.

	The algorithm used for the interpolation is a multi-spline using quadratic programming. This means that we
	do a cubic spline interpolation of the expected losses. The splines are solved via a constrained minimisation
	algortihm (quadratic programming). The spline curves represent the expected losses of the normalised or scaled 
	tranches (i.e. tranches with unit size). 

	Multi-spline means that we solve a single minimisation problem where the variable is the sequence of all tranche curves.
	This allows us to impose the consistency constraints between the tranches.

	Note that although we are solving for a cubic spline our final curves are constructed by taking the spline values at
	the node points (time interpolation dates) and assuming piecewise linear interpolation in between. This is important
	because it ensures that the tranche par spread is linear in the expected loss curve, which is necessary since quadratic
	programming only allows linear constraints.

*/
class PRODUCTS_DLL SplineTrancheQuoteInterpolator:  public CObject,
							public virtual ITrancheQuoteInterpolator
{

public:
    static CClassConstSP const TYPE;

	/** public constructor */
	SplineTrancheQuoteInterpolator(
		const DateTimeArray & interpDates,
		ITrancheQuoteInterpolatorSP  prior,
		double smoothingWeight,
		bool calibrateToIndex,
		bool calibrateToMid,
		double minSepTime,
		double minSepStrike
		);

	/** Destructor */
    virtual ~SplineTrancheQuoteInterpolator();

	/** Called immediately after object constructed */
    virtual void validatePop2Object();    

	/** Returns expected loss surface for set of strikes and dates
		This surface should reprice the input trancheQuotes */
    virtual ExpectedLossSurfaceSP getELSurface(
		const CDOQuotes & quotes,					/** market tranche quotes */
		const ICDSParSpreads & indexSwapSpreads		/** index swap spreds */
		) const;
	
  

private:

	// DEFAULTS ----------------

	/** default smoothing weight */
	static const double DEFAULT_SMOOTHING;

	/** private constructor */
	SplineTrancheQuoteInterpolator();

    // For reflection
    static void load (CClassSP& clazz);

	static IObject* defaultSplineTrancheQuoteInterpolator();

	// METHODS ------------------------------------------------------

	/** create interpDates from interpInterval */
	DateTimeArraySP createSplineDates(
		const CDOQuotes & quotes,
		const ICDSParSpreads & index
		) const;

	/** return object containing spline curves */ 
	ConstrainedCubicSplineSP getSpline( 
		  const CDOQuotes & quotes,				// CDO quotes												 
		  const ICDSParSpreads & index,			// index swap spreads										 
		  const DoubleArray & lowS,				// low strikes for tranches which we want to interpolate	 
		  const DoubleArray & highS,			// corresponding high strikes								 
		  DateTimeArraySP	splineDates,		// splineDates												 
		  const DoubleArray & timeGrid			// yearFracs (ACT/365) corresponding to interpDates		 
		  ) const;	

	/** construct and return EL surface from solved spline  */
	ExpectedLossSurfaceSP populateELSurface(
		const DoubleArray & lowS,				// low strikes						 
		const DoubleArray & highS,				// high strikes					 
		ConstrainedCubicSplineSP  spline,		// spline that has been solved		 
		DateTimeArraySP	splineDates,			// dates at which spline is defined	 
		double RR								// recovery rate					 
		) const;


	/** get prior loss points array from prior class  */
	void setPrior( 
		  const CDOQuotes & quotes,				// tranche quotes  
		  const ICDSParSpreads & index,			// index swap spreads  
		  const DoubleArray & lowS,				// low strikes  
		  const DoubleArray & highS,			// high strikes  
		  DateTimeArraySP splineDates,			// splineDates  
		  DoubleArray & prior,					// (O) prior array  
		  DoubleArray &	priorWeights			// (O) prior weights  
		  ) const;


	/** populate the linear spread constraint and constraint values  */
	void getSpreadConstraints(
		const CDOQuotes & quotes,				// CDO quotes	 
		const ICDSParSpreads & index,			// index swap spreads	 
		const DoubleArray & lowS,				// low strikes		 
		const DoubleArray & highS,				// high strikes	 
		DateTimeArraySP	splineDates,			// splineDates  
		const DoubleArray & timeGrid,			// interpolation time grid	 
		const string & quoteType,				// MID, BID or ASK	 
		DoubleArray & spreadConstraintVals,		// (O) constraint values	 [v]	 
		DoubleArray & spreadConstraintFuncs		// (O) constraint functionals [c(i)]	 
		) const;


	/** populate concavity constraints  */
	void getConcavityConstraints(
		const DoubleArray & lowS,				// low strikes	 
		const DoubleArray & highS,				// high strikes	 
		double minSep,							// absolute minimum separation between tranche EL's (unit notional tranches) 
		double RR,								// recovery rate  
		int numSplineDates,						// size of splineDates  
		DoubleArray & constraintVals,			// (O) constraint values	 
		DoubleArray & constraintFuncs			// (O) constraint functionals	 
		) const;	

	/** populate bound constraint, i.e. all normalised EL's <= 1  */
	void getBoundConstraint(
		 int numStrikes,						// number of tranches	 
		 int numSplineDates,					// size of splineDates  
		 DoubleArray & constraintVals,			// constraint values	 
		 DoubleArray & constraintFuncs			// constraint functionals	 
		 ) const;


	/** returns constraint value and populates constraint functional 
		for linear spread constraint  */
	double linearSpreadConstraint(
		double spread,							// target spread	 
		double upfront,							// target upfront	 	
		const DateTime & quoteDate,				// date corersponding to quotes  
		double lowStrikeTarget,					// low strike corresponding to quotes. 0 for index swap quote  
		double highStrikeTarget,				// high strike corresponding to quotes. 1 for index swap quote  
		const DoubleArray & lowS,				// low strikes	 
		const DoubleArray & highS,				// high strikes	 
		DateTimeArraySP splineDates,			// splineDates. Need to start at valueDate  
		const DoubleArray & timeGrid,			// time grid for interpolation	 
		 YieldCurveConstSP & yieldCurveSP,		// yield curve	 
		 const DoubleArray & yearFracsDf,		// Act/365F year fraction corresponding to zeroDates	 
		 const DoubleArray & dfCurve,			// array of discount factors defined at yearFracsDf times	 
		 CreditFeeLegWithPVSP indexFeeLegSP,	// index fee leg	 
		 double RR,								// recovery  
		const string & quoteType,				// MID, BID or ASK	 
		int startIndex,							// first index in constraintArray to populate. This should be start of constraintFunc row	 
		DoubleArray & constraintArray			// (O) constraint array to populate	 
		) const;


	/**	function that produces contingent and fee leg effective curve values 
		of a chosen tranche with strikes: lowStrikeTarget, highStrikeTarget */  
	void effCurveValue(
		double lowStrikeTarget,		// low strike of tranche for which we want effective curve  
		double highStrikeTarget,	// high strike of tranche for which we want effective curve  
		double lowStrikeTest,		// low strike of tranche which has 0 effective curve val  
		double highStrikeTest,		// high strike of tranche that has 0 effective curve val  
		double RR,					// recovery rate  
		double & contLegValue,		// (O) contingent leg eff curve value  
		double & feeLegValue		// (O) fee leg eff curve value  
		) const;

	
	// FIELDS --------------------------------
	
	/** prior loss surface */
	ITrancheQuoteInterpolatorSP prior;

	/** smoothing weight. Lies in [0,1[ */
	double smoothingWeight;
	
	/** time interpolation discretisation interval. */
	string interpInterval;

	/** Override dates at which to interpolate [default: quote dates] */
	DateTimeArraySP interpDates;


	/** calibrate to index swap spread (TRUE) or super senior quote (FALSE) */
	bool calibrateToIndex;

	/** calibrate to Mid (TRUE) or to bid/ask quotes (FALSE) */
	bool calibrateToMid;

	/** minimum separation between losses of all normalised tranches at consecutive interpDate points (bps) */
	double minSepTime;
	
	/** minimum separation between losses of consecutive normalised tranches at every time point (bps) */
	double minSepStrike;

	/** weights for prior per tranche at last date in time line (introduced for zero stress)*/
	DoubleArraySP lastWeights;

	

};

DECLARE(SplineTrancheQuoteInterpolator);

DRLIB_END_NAMESPACE

#endif


