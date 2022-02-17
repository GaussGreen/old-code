//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research (Credit Hybrids)
//
//   Filename    : ConstrainedCubicSpline.hpp
//
//   Description : finds cubic spline under linear equality and inequality constraints using quadratic programming
//
//   Author      : Matthias Arnsdorf
//
//   Date        : 5 October 2005
//
//
//----------------------------------------------------------------------------
#include "edginc/AtomicArray.hpp"
#include "edginc/QuadraticProg.hpp"

#ifndef EDR_CONSTRAINED_CUBIC_SPLINE_HPP
#define EDR_CONSTRAINED_CUBIC_SPLINE_HPP

DRLIB_BEGIN_NAMESPACE


/*************************************************************************** 

	Implementation of algortihm described in:
	
	B.A. Turlach (1997), "Constrained Smoothing Splines Revisited", 
	Australian National University Preprint
	
	Uses quadratic programming to solve for cubic spline under linear equality
	and inequality constraints.

****************************************************************************/


class UTIL_DLL ConstrainedCubicSpline
{
public:

	/** Destructor */
	virtual ~ConstrainedCubicSpline();

	/** Constructor */
	ConstrainedCubicSpline(
		const DoubleArray & xGrid,			/** node points of spline */
		const DoubleArray & prior,			/** prior values at node points */
		const DoubleArray &	priorWeights,	/** weight for each prior value */					
		double smoothingWeight,				/** smoothing weight */
		int numCurves						/** number of splines in multi-spline */
		);
		
		

	/** get array of all values. If have more than one spline
	values are returned consecutively */
	void getValueArray(DoubleArray & valueArray) const;

	/** overloaded version returns only values for spline indexed by splineNum (zero offset)*/
	void getValueArray(int splineNum, DoubleArray & valueArray) const;


	/** Returns values of spline curveNum at given time points 
		Note that curveNum has unit offset, i.e. 1<= curveNum <=numCurves */
	DoubleArraySP values(int curveNum, DoubleArray& times);

	
	
	/** find the spline with no addtional constraints 
		Throws IMSLException if not succesful*/
	void solve();

	/** solve for spline with additional constraints 
			Throws IMSLException if not succesful */
	void solve(
		const DoubleArray& equalConstraintFuncs,	// flattened matrix of dim =  (rows x cols) = (numEqConstraints x numVariables)
		const DoubleArray& equalConstraintVals,		// Array dim = numEqConstraints
		const DoubleArray& inequalConstraintFuncs,	// flattened matrix dim =(rows x cols) = (numInEqConstraints x numVariables)
		const DoubleArray& inequalConstraintVals	// Array dim = numInEqConstraints
		);

	/** set prior values */
	void setPrior(const DoubleArray & prior);

	/** set target values that spline should hit */
	void setTargets(
		const DoubleArray & targetValues,		/** values spline should go through */	
		const IntArray & targetIndexInGrid,	    /** index of position in xGrid corresponding to targetValue */
		const IntArray & splineIndex			/** for each target value this gives the spline number (zero offset) the target refers to */
		);

	/** overload of setTarget for single spline */
	void setTargets(
		const DoubleArray & targetValues,			/** values spline should go through */
		const IntArray & targetIndexInGrid		/** index of position in xGrid corresponding to targetValue */
		);

	/** set splines to be increasing with minimum separation */
	void setIncreasing(double minSep);

	/** set spline to be concave (curvature < 0) or convex (curvature >=0) */
	void setCurvature(double curvature);

	/** set discrete curvature (i.e. only consider spline values at xGrid points not 2nd derivatives) */
	void setDiscreteCurvature(double curvature, double minSep);


	/** set upper bound for values */
	void setHighBound(const DoubleArrayArray & bounds); // bounds[splineNum][gridPoint]
	
	/** overload of setHighBound for single spline */
	void setHighBound(const DoubleArray & bounds);		// bounds[gridPoint]

	/** set lower bound for values */
	void setLowBound(const DoubleArrayArray & bounds);	// bounds[splineNum][gridPoint]

	/** overload of setLow Bound for single spline */
	void setLowBound(const DoubleArray & bounds);		// bounds[gridPoint]

	

	

private:
	
	// Methods /////////////////////////////////////////////
	void calcSplineConstraintFuncs();
	void calcHessian();
	void calcLinCoeffs(const DoubleArray & prior);

	void calcSplineParams();

	/** returns value of spline @curveIndex given xGridIndex I and time t such that xGrid[I] <= t <= xGrid[I+1] 
	0 < I < N-1 */
	// 0 <= curveIndex < numCurves
	double val(int curveIndex, int xGridIndex, double time) const;

	
	
	// Fields ////////////////////////////////////////////////////////////////
	
	/** Node points for spline */
	DoubleArray xGrid;
	/** weights for prior */
	DoubleArray  priorWeights;
	/** smoothig weight */
	double  smoothingWeight;
	/** nnumber of splines in multi-spline */
	int numCurves;
	
	/** array of differeces in xGrid: dx[i] = xGrid[i+1] - xGrid[i]*/
	DoubleArray dx;

	/** flattend matrix that contains spline consistency constraints only */
	DoubleArray splineConstraintFuncs;
	DoubleArray splineConstraintVals;
	DoubleArray hessian;
	DoubleArray linCoeffs;

	/** number of xGrid points */
	int N;
	/** number of spline points */
	int numVariables;
	/** total number of (linear) constraints for spline (not additional constraints) */
	int numSplineConstraints;
	/** number of spline equality constraints (linear) */
	int numSplineEqualityConstraints;

	/** target constraints functions*/
	DoubleArray targetFuncs;
	/** Additional equality constraint values */
	DoubleArray targetVals;

	/** increasing constraints functions*/
	DoubleArray increaseFuncs;
	/** increasing constraint values */
	DoubleArray increaseVals;

	/** curvature constraints functions*/
	DoubleArray curvatureFuncs;
	/** curvature constraint values */
	DoubleArray curvatureVals;

	/** highBound constraints functions*/
	DoubleArray highBoundFuncs;
	/** highBound constraint values */
	DoubleArray highBoundVals;

	/** highBound constraints functions*/
	DoubleArray lowBoundFuncs;
	/** highBound constraint values */
	DoubleArray lowBoundVals;

	
	

	// ----------------- Results --------------------------------------------------
	/** Array which holds N value points and (N-2) 2nd derivatives consecutively */
	DoubleArraySP result;

	/** Flag to indicate if result has been calculated */
	bool resultIsCalculated;

	
	/** Arrays of length N to hold spline coefficients */
	// spline in segment (xGrid[i], xGrid[i+1]) has form:
	// S_i(t) = value[i] + b[i](t-t[i]) + c[i](t-t[i])^2 + d[i](t-t[i])^3

	DoubleArray b;
	DoubleArray c;
	DoubleArray d;

	/** Flag to indeicate that spline params have been calculated */
	bool splineParamsAreCalculated;

	

};

typedef refCountPtr<ConstrainedCubicSpline> ConstrainedCubicSplineSP;

DRLIB_END_NAMESPACE

#endif

