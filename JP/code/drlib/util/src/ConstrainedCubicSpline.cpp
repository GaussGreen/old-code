//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research (Credit Hybrids)
//
//   Filename    : ConstrainedCubicSpline.cpp
//
//   Description : finds cubic spline under linear equality and inequality constraints using quadratic programming
//
//   Author      : Matthias Arnsdorf
//
//   Date        : 5 October 2005
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/ConstrainedCubicSpline.hpp"
#include "edginc/Format.hpp"




DRLIB_BEGIN_NAMESPACE

/** Destructor */
ConstrainedCubicSpline::~ConstrainedCubicSpline()
{}

/** Constructor */
ConstrainedCubicSpline::ConstrainedCubicSpline(
		const DoubleArray & _xGrid,
		const DoubleArray & _prior,
		const DoubleArray &	_priorWeights,					
		double _smoothingWeight,
		int _numCurves
		) 
		: 
		xGrid(_xGrid),
		priorWeights(_priorWeights),
		smoothingWeight(_smoothingWeight),
		numCurves(_numCurves),
		dx(0),
		splineConstraintFuncs(0),
		splineConstraintVals(0),
		hessian(0),
		linCoeffs(0),
		N(0),
		numVariables(0),
		numSplineConstraints(0),
		numSplineEqualityConstraints(0),
		targetFuncs(0),
		targetVals(0),
		increaseFuncs(0),
		increaseVals(0),
		curvatureFuncs(0),
		curvatureVals(0),
		highBoundFuncs(0),
		highBoundVals(0),
		lowBoundFuncs(0),
		lowBoundVals(0),
		result(0),
		resultIsCalculated(false),
		b(0),
		c(0),
		d(0),
		splineParamsAreCalculated(false)
		
		
	
{
	static const string method = "ConstrainedCubicSpline::ConstrainedCubicSpline";
	try
	{

		int i; //counters

		N = xGrid.size();
		// -------------- validation
		if(N<2) throw ModelException("Size of xGrid is less than 2");
		if(priorWeights.size() != numCurves*N) throw ModelException("priorWeights size different to xGrid size*numCurves");
		if(_prior.size() != numCurves*N) throw ModelException("prior size different to xGrid size");
		if(smoothingWeight <0) throw ModelException("smoothingWeight needs to be positive");

			
		// Define dx
		dx.resize(N-1);
		for(i = 0;i<N-1;i++)
		{
			dx[i] = xGrid[i+1] - xGrid[i];
			
			if(dx[i] <= 0)
			{
				throw ModelException("dx["+Format::toString(i)+"] = "+Format::toString(dx[i])+" <= 0");
			}
		}
		
		/** number of variables needed to describe a multi spline */
		// This is N values and N-2 2nd derivatives per curve
		// A spline is defined by N values and (N-2) 2nd derivatives per curve.
		// The 2nd derivative at the end points of the spline is assumed 0
		numVariables = numCurves*(2*N-2);							
		
		

		/** total number of constraints needed to define spline */
		numSplineConstraints = numCurves*(N - 2);					// "num rows" in splineConstraintFuncs
		
		// A spline is defined using quadratic programming only using equality constraints
		// Hence the number of equality constraints = total number of constraints
		// Additioanal inequality constraints that the splines should satisdy can be set when calling slove(..)

		/** number of equality constraints */
		numSplineEqualityConstraints = numSplineConstraints;		 
		
	
		// calculate the 'design matrices' 
		calcHessian();
		calcSplineConstraintFuncs(); 
		calcLinCoeffs(_prior);
		
		
		// define spline constraint values
		// Have \sum_i splineConstraintFuncs[i][j] x[i] = splineConstraintVals[j]
		// where x are the spline values
		// j runs over the number of constraints
		
		splineConstraintVals.resize(numSplineConstraints, 0.0); // spline constraints all equate to 0
		
		

	} catch (exception & e)
	{
		throw ModelException(e,method);
	}
}
		

//////////////////////////////////////////////////////////////////////////////////////////////////
/** calcSplineConstraintFuncs

	Denote splineConstraintFunc matrix as A(i,j)  single spine dim = (n-2)x(2n-2)
	Defined by
	
	\sum_j A(i,j)x[j] = constraintVal[i] = 0;
	   
	Form is A = (Q^T, -R)
	 
	 R as defined as:
			 R(i,i) = (dx(i) + dx(i+1))/3
			 R(i,i+1) = R(i+1,i) = dx(i+1)/6
			 otherwise R = 0
	
	  Q^T is defined as:	
			 Q^T(i,i) = dx(i)^{-1} 
			 Q^T(i,i+1) = -dx(i)^{-1}-dx(i+1)^{-1}
			 Q^T(i,i+2) = dx(i+1)^{-1}
			 else 0

	in case of multiple spline the constraint matrix consists of multiple (block diagonal) copies of the 
	spline constraint matrix

  Note that matrix is represented by flattening to DoubleArray
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void ConstrainedCubicSpline::calcSplineConstraintFuncs()
{
	static const string method = "ConstrainedCubicSpline::calcSplineConstraintFuncs";
	try
	{
		int i,j,k; // counters
		// for individual spline constraint matrices
		int numRows = N-2;
		int numCols = 2*N-2;
		// for total constraint matrix
		int totNumRows = numCurves*numRows;
		int totNumCols = numCurves*numCols;

		DoubleArray singleSplineConstraint(numRows*numCols);
		splineConstraintFuncs.resize(totNumRows*totNumCols);

		//first N cols given by Q^T
		for(j=0; j< N;j++)
		{
			for(i = 0;i<numRows;i++)
			{
				if(i==j)
				{
					singleSplineConstraint[numCols*i + j] = 1/dx[i];
				}
				else if(j == i+1)
				{
					singleSplineConstraint[numCols*i + j] = -1/dx[i] - 1/dx[i+1];
				}
				else if(j == i+2)
				{
					singleSplineConstraint[numCols*i + j] = 1/dx[i+1];
				}
				else singleSplineConstraint[numCols*i + j] = 0.0;
			}
		}
		// next n-2 cols given by -R
		int col; // column number in R matrix
		for(j = N; j< numCols;j++)
		{
			col = j-N;
			for(i = 0; i < numRows; i++)
			{
				if(col == i)
				{
					singleSplineConstraint[numCols*i + j] = -(dx[i]+dx[i+1])/3;
				}
				else if(col == i+1)
				{
					singleSplineConstraint[numCols*i + j] = -dx[col]/6;
				}
				else if(col == i-1)
				{
					singleSplineConstraint[numCols*i + j] = -dx[i]/6;
				}
				else singleSplineConstraint[numCols*i + j] = 0.0;
			}
		}

		// copy singleSplineConstraints into splineConstraintFuncs
		int row;
		for(k = 0;k<numCurves;k++)
		{
			for(i= k*numRows;i<(k+1)*numRows;i++) // loop through rows
			{
				for(j = 0;j<totNumCols;j++) // loop through columns
				{
					if(j<k*numCols || j >=(k+1)*numCols) splineConstraintFuncs[i*totNumCols + j] = 0.0;
					
					else
					{
						row = i - k*numRows;
						col = j - k*numCols;
						splineConstraintFuncs[i*totNumCols + j] = singleSplineConstraint[row*numCols+col];
					}
				}
			}
		}
					

	} catch (exception & e)
	{
		throw ModelException(e, method);
	}
}


/////////////////////////////////////////////////////////////////////////////////////////
/** calcSplineHessian 

  Denote Hessian by H(i,j).  Single spline Hessian has dim (2n-2)x(2n-2)

  Defined by:					W	0
				H(i,j)  =	(	0	s*R	)

  where W = diag(priorWeights) and s = smoothingWeight. R is as defined above.

  For multi-splines H is a contains block diagonal conpies of the above - 1 per curve.

  Note that matirx is represented by flattening to DoubleArray
*/
////////////////////////////////////////////////////////////////////////////////////////////////
void ConstrainedCubicSpline::calcHessian()
{
	static const string method = "ConstrainedCubicSpline::calcHessian";
	try
	{
		int i,j,k,l;
		int numRows = 2*N-2;
		int numCols = 2*N-2;
		int totNumRows = numCurves*numRows;
		int totNumCols = numCurves*numCols;
		
		DoubleArray singleHessian(numRows*numCols);
		hessian.resize(totNumRows*totNumCols);
		
		//first N rows given by diagonal weights
		for(i = 0; i < N;i++)
		{
			for(j = 0;j < numCols;j++)
			{
				if(i==j) singleHessian[numCols*i + j] = priorWeights[i];
				else singleHessian[numCols*i + j] =0.0;
			}
		}
		// next N row given by smoothingWeight*R, where R is a design matrix
		// R(i,i) = (dx(i) + dx(i+1))/3
		// R(i,i+1) = R(i+1,i) = dx(i+1)/6
		// otherwise R = 0
		for(i = N; i < numRows;i++)
		{
			for(j = 0;j < numRows;j++)
			{
				if(j==i)
				{
					singleHessian[numCols*i + j] = smoothingWeight*(dx[i-N]+dx[i-N+1])/3;
				}
				else if(j==i+1)
				{
					singleHessian[numCols*i + j] = smoothingWeight*dx[j-N]/6;
				}
				else if(j==i-1 && j >= N)
				{
					singleHessian[numCols*i + j] = smoothingWeight*dx[i-N]/6;
				}
				else singleHessian[numCols*i + j] =0.0;
			}
		}

		// Do multiple diagonal copies into hessian
		int row;
		int col;
		for(k = 0;k<numCurves;k++)
		{
			//copy prior weights into top diagonal of single Hessian if have more than one curve
			if(k>0)
			{
				for(l=0;l<N;l++)
				{
					singleHessian[numCols*l+l] = priorWeights[k*N+l];
				}
			}

			for(i= k*numRows;i<(k+1)*numRows;i++) // loop through rows
			{
				for(j = 0;j<totNumCols;j++) // loop through columns
				{
					if(j<k*numCols || j >=(k+1)*numCols) hessian[i*totNumCols + j] = 0.0;
					
					else
					{
						row = i - k*numRows;
						col = j - k*numCols;
						hessian[i*totNumCols + j] = singleHessian[row*numCols+col];
					}
				}
			}
		}
		

	} catch (exception & e)
	{
		throw ModelException(e, method);
	}
}

///////////////////////////////////////////////////////////////////////////////////
/** calcLinCoeffs 
 denote vector of linear coeffs by f(i). Dim (2n-2) per Curve

  f(i)  = - priorWeights[i]*prior[i] for i<N 
  f(i) = 0 otherwise

  for multi-spline f(i) contains consecutive copies of the above
*/
//////////////////////////////////////////////////////////////////////////////////
void ConstrainedCubicSpline::calcLinCoeffs(const DoubleArray & prior)
{
	static const string method = "ConstrainedCubicSpline::calcLinCoeffs";
	try
	{
		int i,k;
		
		int numRows = 2*N-2;
		int totNumRows = numCurves*numRows;

		linCoeffs.resize(totNumRows);

		for(k = 0;k<numCurves;k++)
		{
			for (i = k*numRows ;i < (k+1)*numRows;i++)
			{
				if(i < k*numRows+N)
				{
					linCoeffs[i] = -prior[k*N + i-k*numRows]*priorWeights[k*N + i-k*numRows];
				}
				else linCoeffs[i] = 0.0;
			}
		}

	} catch (exception & e)
	{
		throw ModelException(e, method);
	}
}

////////////////////////////////////////////////////////////////////////////////
/** solve 
 solve the QP problem without additional constraints
*/
////////////////////////////////////////////////////////////////////////////////
void ConstrainedCubicSpline::solve()
{				
	solve(DoubleArray(0),DoubleArray(0),DoubleArray(0),DoubleArray(0));	 
}

//////////////////////////////////////////////////////////////////////////////////
/** solve with additonal constraints */
////////////////////////////////////////////////////////////////////////////////////
void ConstrainedCubicSpline::solve(
		const DoubleArray& equalConstraintFuncs,  // flattened matrix dim rows x cols = numEqConstraints x numVariables
		const DoubleArray& equalConstraintVals,	 // Array dim numEqConstraints
		const DoubleArray& inequalConstraintFuncs, //flattened matrix dim rows x cols = numInEqConstraints x numVariables
		const DoubleArray& inequalConstraintVals // Array dim numInEqConstraints
		)
{
	static const string method = "ConstrainedCubicSpline::solve(.,.,.,.)";
	try
	{
		// validation
		int numEq = equalConstraintVals.size();
		int numInEq = inequalConstraintVals.size();

		if(equalConstraintFuncs.size() != numEq*numVariables)
		{
			throw ModelException("equalConstraintFuncs not of size numEq*numVariables");
		}
		if(inequalConstraintFuncs.size() != numInEq*numVariables)
		{
			throw ModelException("inqualConstraintFuncs not of size numInEq*numVariables");
		}

		// create full constraint matrix by prepending. Note equality constraints have to come before ineq constraints
		DoubleArray constraints(0);
		DoubleArray constraintVals(0);

		// inequality constraints -------------------------------------------
		constraints.insert(constraints.begin(),inequalConstraintFuncs.begin(),inequalConstraintFuncs.end());
		constraintVals.insert(constraintVals.begin(),inequalConstraintVals.begin(),inequalConstraintVals.end());
		
		// increasing constraint
		constraints.insert(constraints.begin(),increaseFuncs.begin(),increaseFuncs.end());
		constraintVals.insert(constraintVals.begin(),increaseVals.begin(),increaseVals.end());
		
		// curvature constraint
		constraints.insert(constraints.begin(),curvatureFuncs.begin(),curvatureFuncs.end());
		constraintVals.insert(constraintVals.begin(),curvatureVals.begin(),curvatureVals.end());
		
		/// high  bound
		constraints.insert(constraints.begin(),highBoundFuncs.begin(),highBoundFuncs.end());
		constraintVals.insert(constraintVals.begin(),highBoundVals.begin(),highBoundVals.end());

		// low bound
		constraints.insert(constraints.begin(),lowBoundFuncs.begin(),lowBoundFuncs.end());
		constraintVals.insert(constraintVals.begin(),lowBoundVals.begin(),lowBoundVals.end());
		
		
		// equality constraints ----------------------------------------------------------
		constraints.insert(constraints.begin(),equalConstraintFuncs.begin(),equalConstraintFuncs.end());
		constraintVals.insert(constraintVals.begin(),equalConstraintVals.begin(),equalConstraintVals.end());
		
		// add target constraints
		constraints.insert(constraints.begin(),targetFuncs.begin(),targetFuncs.end());
		constraintVals.insert(constraintVals.begin(),targetVals.begin(),targetVals.end());
		
		// add spline constraints
		constraints.insert(constraints.begin(),splineConstraintFuncs.begin(),splineConstraintFuncs.end());
		constraintVals.insert(constraintVals.begin(),splineConstraintVals.begin(),splineConstraintVals.end());
		
	
		/** total number of equality constraints */
		int totalNumEqConstraints = numSplineEqualityConstraints;
		totalNumEqConstraints += numEq + targetVals.size();

		// calculate result array
		result = QuadraticProg::minimize(
			constraintVals.size(),			// Total number of linear constraints (m)
			numVariables,					// number of variables (n)
			totalNumEqConstraints,			// number of linear equality constraints (meq)
			constraints,					// flattened array of size m x n containing equality constraints in first meq rows followed by inequality constraints
			constraintVals,					// array with m components containing RHS of linear constraints (b)
			linCoeffs,						// array with n components containing the coeffs of linear term of objective function (g)
			hessian						// array (m x n) containing Hessian matix of objective function (h)
		);


		// set calc flags
		resultIsCalculated = true;
		splineParamsAreCalculated = false; 
				

	} 
	catch(IMSLException& eIMSL)
	{
		throw IMSLException(eIMSL);
	}
	catch (exception & e)
	{
		throw ModelException(e,method);
	}
}


//////////////////////////////////////////////////////////////////////
/** setPrior */
void ConstrainedCubicSpline::setPrior(const DoubleArray & prior)
{
	static const string method = "constrainedCubicSpline::setPrior";
	try
	{
		if(prior.size() != priorWeights.size())
		{
			throw ModelException("Prior size is not equal to priorWeights size");
		}

		calcLinCoeffs(prior);

	} catch (exception & e)
	{
		throw ModelException(e,method);
	}
}

///////////////////////////////////////////////////////////////////////////////
/**
  add target constraints for mulit-spline
*/
void ConstrainedCubicSpline::setTargets(
		const DoubleArray & targetValues,
		const IntArray & targetIndexInGrid,		//
		const IntArray & splineNumber			// for each target value this gives the spline number the target refers to
		)
{
	static const string method = "ConstrainedCubicSpline::setTargets()";
	try
	{
		int i;
		
		/** number of rows in (flattened constraint matrix) = number of constraints */
		int numRows = targetValues.size();
		/** number of varaibles for single spline */
		int singleSplineLength = 2*N-2;
		/** number of columns in flattened constraint matrix */
		int numCols = numVariables;
		

		if(targetIndexInGrid.size() != numRows)
		{
			throw ModelException("targetIndexInGrid size (" + Format::toString((int)targetIndexInGrid.size()) +
				") different from targetValues size (" + Format::toString(numRows) +")");
		}
		if(splineNumber.size() != numRows) throw ModelException("splineNumber has incorrect size");


		/** Array (flattened matrix) of target constraint coefficients */
		targetFuncs.resize(numRows*numCols, 0.0);
	

		// set coefficients corresponding to targetValues
		for(i = 0;i< numRows;i++)
		{
			if(targetIndexInGrid[i] <0 || targetIndexInGrid[i] > N-1)
			{
				throw ModelException("targetIndexInGrid[i] ("
					+ Format::toString(targetIndexInGrid[i]) +
					") at position: " + Format::toString(i) +
					" is out of bounds [0, " +
					Format::toString(N-1) + "]"
					);
			}
			
			if(splineNumber[i] < 0 || splineNumber[i] > numCurves-1)
			{
				throw ModelException("splineNumber[i] ("
					+ Format::toString(splineNumber[i]) +
					") at position: " + Format::toString(i) +
					" is out of bounds [0, " +
					Format::toString(numCurves-1) + "]"
					);
			}
				
			targetFuncs[i*numCols+ singleSplineLength*splineNumber[i] + targetIndexInGrid[i]] = 1.0;
		}

		
		// update vals. Note that vals are just given by targetValues
		targetVals = targetValues;
		


	} catch (exception & e)
	{
		throw ModelException(e, method);
	}

}

////////////////////////////////////////////////////////////////////////////////////
/** overload setTargets for single spline */
///////////////////////////////////////////////////////////////////////////////////
void ConstrainedCubicSpline::setTargets(
		const DoubleArray & targetValues,
		const IntArray & targetIndexInGrid
		)
{
	static const string method = "ConstrainedCubicSpline::setTargets(singleSpline)";
	try
	{
		int numTargets = targetValues.size();
		
		/** spline number array needed for setTarget */
		IntArray splineNum(numTargets, 0); // zero offset
		

		setTargets(targetValues, targetIndexInGrid, splineNum);
		

	} catch (exception & e)
	{
		throw ModelException(e, method);
	}

}

////////////////////////////////////////////////////////////////////////////////////////////
/**
  set splines to be increasing with step = minSep
*/
void ConstrainedCubicSpline::setIncreasing(double minSep)
{
	static const string method = "ConstrainedCubicSpline::setIncreasing";
	try
	{
		int i,j, row;
		
		/** number of rows in (flattened constraint matrix) = number of constraints */
		// have one constraint for each value except for frist
		int numRows = (N-1)*numCurves;
		/** number of varaibles for single spline */
		int singleSplineLength = 2*N-2;
		/** number of columns in flattened constraint matrix */
		int numCols = numVariables;
		


		/** Array (flattened matrix) of target constraint coefficients */
		increaseFuncs.resize(numRows*numCols, 0.0);
		increaseVals.resize(numRows);


		// set coefficients in func
		// loop through splines
		for(i = 0;i<numCurves;i++)
		{
			//loop through values starting at 1
			for(j = 1; j< N;j++)
			{
				// row for constraint
				row = i*(N-1) + (j-1);
				
				//set constraint vals
				increaseVals[row] = minSep;
				
				// constraint x_{i} - x_{i-1} >= 0
				increaseFuncs[row*numCols + i*singleSplineLength + (j-1)] = -1.0;
				increaseFuncs[row*numCols +i*singleSplineLength + j] = 1.0;
			}
		}


		
		

	} catch (exception & e)
	{
		throw ModelException(e, method);
	}

}

////////////////////////////////////////////////////////////////////////////////////////////
/**
  set splines to be concave or convex
*/
///////////////////////////////////////////////////////////////////////////////////////////
void ConstrainedCubicSpline::setCurvature(double curvature)
{
	static const string method = "ConstrainedCubicSpline::setCurvature";
	try
	{
		int i,j, row;
		
		/** number of rows in (flattened constraint matrix) = number of constraints */
		// have one constraint for each 2nd derivative
		int numRows = (N-2)*numCurves;
		/** number of variables for single spline */
		int singleSplineLength = 2*N-2;
		/** number of columns in flattened constraint matrix */
		int numCols = numVariables;
		


		/** Array (flattened matrix) of target constraint coefficients */
		curvatureFuncs.resize(numRows*numCols);
		curvatureVals.resize(numRows);

		for(i = 0;i<curvatureFuncs.size();i++) curvatureFuncs[i] = 0.0; // initialise

		// set coefficients in func
		// loop through splines
		for(i = 0;i<numCurves;i++)
		{
			//loop through 2nd derivs starting
			for(j = 0; j< N-2;j++)
			{
				// row for constraint
				row = i*(N-2) + j;
				
				//set constraint vals
				curvatureVals[row] = 0.0 ;
				
				// constraint: gamma_{i}  >= 0 for convex or -gamma_{i} >= 0 for concave, (gamma is 2nd deriv)
				curvatureFuncs[row*numCols + i*singleSplineLength + N + j] = curvature<0 ? -1.0 : 1.0;
				
			}
		}
		

	} catch (exception & e)
	{
		throw ModelException(e, method);
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////
/**
  set discrete curvatrure -i.e. compare numerical 2nd derivative base on values at node points only.
  cannot set discrete curvature and curvature at the same time
*/
/////////////////////////////////////////////////////////////////////////////////////////////
void ConstrainedCubicSpline::setDiscreteCurvature(double curvature, double minSep)
{
	static const string method = "ConstrainedCubicSpline::setDiscreteCurvature";
	try
	{
		int i,j, row;
		
		/** number of rows in (flattened constraint matrix) = number of constraints */
		// have for each differnce in (discrete) first derivs i.e. N-2
		int numRows = (N-2)*numCurves;
		/** number of variables for single spline */
		int singleSplineLength = 2*N-2;
		/** number of columns in flattened constraint matrix */
		int numCols = numVariables;
		
		// set curvature to -1 or 1
		curvature = curvature<0 ? -1.0 : 1.0;
		
		/** Array (flattened matrix) of target constraint coefficients */
		curvatureFuncs.resize(numRows*numCols, 0.0);
		curvatureVals.resize(numRows);


		// set coefficients in func
		// loop through splines
		for(i = 0;i<numCurves;i++)
		{
			//loop through values starting at 0
			for(j = 0; j < N-2;j++)
			{
				// row for constraint
				row = i*(N-2) + j;
				
				//set constraint vals
				curvatureVals[row] = minSep;
				
				
				// constraint: first derivative is increasing or decreasing monotonically 
				// (x_{i+1}-x_i)/dx  >= (<=) (x_i-x_{i-1})/dx + (-) minSep  for convex (concave)
				curvatureFuncs[row*numCols + i*singleSplineLength + j+2] = curvature/dx[j+1];
				curvatureFuncs[row*numCols + i*singleSplineLength  + j+1] = -curvature*(1/dx[j]+1/dx[j+1]);
				curvatureFuncs[row*numCols + i*singleSplineLength  + j] = curvature/dx[j];
				
			}
		}
		

	} catch (exception & e)
	{
		throw ModelException(e, method);
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////
/**
  set upper bound for splines
*/
//////////////////////////////////////////////////////////////////////////////////////////////
void ConstrainedCubicSpline::setHighBound(const DoubleArrayArray & bounds)
{
	static const string method = "ConstrainedCubicSpline::setHighBound";
	try
	{
		int i,j, row;
		
		/** number of rows in (flattened constraint matrix) = number of constraints */
		// have one constraint for each value in each spline
		int numRows = N*numCurves;
		/** number of variables for single spline */
		int singleSplineLength = 2*N-2;
		/** number of columns in flattened constraint matrix */
		int numCols = numVariables;
		
		if(bounds.size() != numCurves) throw ModelException("Need bound array for each spline");
		for(i=0;i<numCurves;i++)
		{
			if(bounds[i].size() != N) throw ModelException("bound array has incorrect gridPoint dimension");
		}


		/** Array (flattened matrix) of target constraint coefficients */
		highBoundFuncs.resize(numRows*numCols, 0.0);
		highBoundVals.resize(numRows);


		// set coefficients in func
		// loop through splines
		for(i = 0;i<numCurves;i++)
		{
			//loop through values
			for(j = 0; j< N;j++)
			{
				// row for constraint
				row = i*N + j;
				
				//set constraint vals
				highBoundVals[row] = -bounds[i][j];
				
				// constraint: - x_i >= -bound
				highBoundFuncs[row*numCols + i*singleSplineLength  + j] = - 1.0;
				
			}
		}
		

	} catch (exception & e)
	{
		throw ModelException(e, method);
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////
/** overload for single spline */
//////////////////////////////////////////////////////////////////////////////////////////////
void ConstrainedCubicSpline::setHighBound(const DoubleArray & bounds)
{
	static const string method = "ConstrainedCubicSpline::setHighBound(DoubleArray)";
	try
	{
		DoubleArrayArray multiBound(1);
		multiBound[0] = bounds;
		setHighBound(multiBound);

	} catch (exception & e)
	{
		throw ModelException(e, method);
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////
/**
  set lower bound for splines
*/
//////////////////////////////////////////////////////////////////////////////////////////////
void ConstrainedCubicSpline::setLowBound(const DoubleArrayArray & bounds)
{
	static const string method = "ConstrainedCubicSpline::setLowBound";
	try
	{
		int i,j, row;
		
		/** number of rows in (flattened constraint matrix) = number of constraints */
		// have one constraint for each value in each spline
		int numRows = N*numCurves;
		/** number of variables for single spline */
		int singleSplineLength = 2*N-2;
		/** number of columns in flattened constraint matrix */
		int numCols = numVariables;
		
		if(bounds.size() != numCurves) throw ModelException("Need bound array for each spline");
		for(i=0;i<numCurves;i++)
		{
			if(bounds[i].size() != N) throw ModelException("bound array has incorrect gridPoint dimension");
		}


		/** Array (flattened matrix) of target constraint coefficients */
		lowBoundFuncs.resize(numRows*numCols, 0.0);
		lowBoundVals.resize(numRows);

		// set coefficients in func
		// loop through splines
		for(i = 0;i<numCurves;i++)
		{
			//loop through values
			for(j = 0; j< N;j++)
			{
				// row for constraint
				row = i*N + j;
				
				//set constraint vals
				lowBoundVals[row] = bounds[i][j];
				
				// constraint:  x_i >= bound
				lowBoundFuncs[row*numCols + i*singleSplineLength  + j] = 1.0;
				
			}
		}
		

	} catch (exception & e)
	{
		throw ModelException(e, method);
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////
/** overload for single spline */
//////////////////////////////////////////////////////////////////////////////////////////////
void ConstrainedCubicSpline::setLowBound(const DoubleArray & bounds)
{
	static const string method = "ConstrainedCubicSpline::setLowBound(DoubleArray)";
	try
	{
		DoubleArrayArray multiBound(1);
		multiBound[0] = bounds;
		setLowBound(multiBound);

	} catch (exception & e)
	{
		throw ModelException(e, method);
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////
/** Spline Params  a,b and c*/
//////////////////////////////////////////////////////////////////////////////////////////////
void ConstrainedCubicSpline::calcSplineParams()
{
	static const string method = "ConstrainedCubicSpline::calcSplineParams()";
	try
	{
		int i,k;

		// number of variables for each curve in result, i.e. N values and N-2 2nd derivs
		int length = 2*N-2;

		b.resize(numCurves*N);
		c.resize(numCurves*N);
		d.resize(numCurves*N);

		for(k=0;k<numCurves;k++)
		{
			for(i = 0;i < N;i++)
			{
				
				// c -------------------
				if(i==0 || i==N-1)
				{
					c[k*N + i] = 0;
				}
				else c[k*N + i] = (*result)[k*length + N + (i-1)]; // Note that first and last 2nd deriv are not stored in result
			}
			for(i = 0;i < N;i++)
			{
				if(i<N-1)
				{
					// b -----------------------
					b[k*N + i] = ((*result)[k*length + (i+1)] - (*result)[k*length + i]) / dx[i] - (2*c[k*N + i]+ c[k*N + (i+1)])*dx[i]/3;

					//d -----------------------------
					d[k*N + i] = (c[k*N + i+1]-c[k*N + i])/(3*dx[i]);
				
				}
				else   // i == N-1
				{
					//b ----------------------------
					b[k*N + i] = ((*result)[k*length+i]-(*result)[k*length+(i-1)])/dx[i-1] + c[k*N+(i-1)]*dx[i-1]/3;

					//d -----------------------------
					d[k*N+i] = 0;
				}
			}
		}

		splineParamsAreCalculated = true;

	} catch (exception & e)
	{
		throw ModelException(e,method);
	}
	
}



//////////////////////////////////////////////////////////////////////////////////////////////
/**
 populates valueArray with values of splines consecutively 
*/
//////////////////////////////////////////////////////////////////////////////////////////////
void ConstrainedCubicSpline::getValueArray(DoubleArray & valueArray) const
{
	static const string method = "ConstrainedCubicSpline::getValueArray";
	try
		{
			if(valueArray.size() != numCurves*N) throw ModelException("valueArray has incorrect size");
				
			if(!resultIsCalculated) throw ModelException("Spline not solved"); // check if result has been calculated
			
			int i,k;
			for(k = 0;k<numCurves;k++)
			{
				for(i = 0;i < N;i++)
				{
					valueArray[k*N + i] = (*result)[k*(2*N-2) + i];
				}
			}

		} catch (exception & e)
		{
			throw ModelException(e,method);
		}

}

// overloaded version that give ouput just for single spline denoted by splineNum
void ConstrainedCubicSpline::getValueArray(int splineNum, DoubleArray & valueArray) const
{
	static const string method = "ConstrainedCubicSpline::getValueArray";
	try
		{
			if(valueArray.size() != N) throw ModelException("valueArray has incorrect size");
				
			if(!resultIsCalculated) throw ModelException("Spline not solved"); // check if result has been calculated
			
			int i;
			for(i = 0;i < N;i++)
			{
				valueArray[i] = (*result)[splineNum*(2*N-2) + i];
			}
			

		} catch (exception & e)
		{
			throw ModelException(e,method);
		}

}




//////////////////////////////////////////////////////////////////////////////////////////////
/** val
 inline for performance , no validation
 spline is of form a + b(t-ti) + c(t-ti)^2 + d(t-ti)^3

 need 0 <= I < N !!!!
*/
//////////////////////////////////////////////////////////////////////////////////////////////
inline double ConstrainedCubicSpline::val(int curveIndex, int I, double t) const
{
	double Dt = (t - xGrid[I]);

	return (*result)[curveIndex*(2*N-2) + I] + b[curveIndex*N + I]*Dt + c[curveIndex*N + I]*(Dt*Dt) + d[curveIndex*N + I]*(Dt*Dt*Dt);

}

//////////////////////////////////////////////////////////////////////////////////////////////
/** values 

 have linear extrapolation outside xGrid
 assume that t is ordered
*/
//////////////////////////////////////////////////////////////////////////////////////////////

DoubleArraySP ConstrainedCubicSpline::values(int curveNum, DoubleArray& t) 
{
	static const string method = "ConstrainedCubicSpline::values";
	try
		{
			// note curveNum has unit offset	
			if(curveNum > numCurves) throw ModelException("curveNum is too high");			
			if(curveNum < 1) throw ModelException("curveNum is too low");
			int curveIndex = curveNum-1;

			if(!resultIsCalculated) throw ModelException("Spline not solved"); // check if result has been calculated
			if(!splineParamsAreCalculated) calcSplineParams();

			int i;
			
			int numTimes = t.size();
			
			if(numTimes<1) throw ModelException("times array is empty");
			
			DoubleArraySP out = DoubleArraySP(new DoubleArray(numTimes));
			
			// find first xGrid point strictly after t[0];
			int iHigh = 0;
			while(iHigh < N && xGrid[iHigh] <= t[0]) iHigh++;


			for(i = 0; i < numTimes; i++)
			{	
				if(iHigh > N) throw ModelException("iHigh > N"); // This can't happen
				if(iHigh == 0) // t < xGrid[0] ; linear extrapolation a + b(t - xGrid[0])
				{
					(*out)[i] = val(curveIndex,0, t[i]);
				}
				else
				{
					(*out)[i] = val(curveIndex, iHigh - 1, t[i]);
				}

				// check if need to move iHigh
				if(i+1 < numTimes)
				{
					if(iHigh < N && xGrid[iHigh] <= t[i+1]) iHigh++;
				}
			}

			return out;

		} catch (exception & e)
		{
			throw ModelException(e,method);
		}

}


 
DRLIB_END_NAMESPACE


