/*
 *
 * Copyright (c) IXIS CIB November 2004 Paris
 *
 *	\file amc_andersen.cpp
 *  \brief
 *
 *	\author  A Schauly
 *	\version 1.0
 *	\date November 2004
 */


/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

#include "gpnummethods/amc_andersen.h"

/// gpbase
#include "gpbase/specialsort.h"
#include "gpbase/gpmatrixlinalg.h"

/// gpnumlib
#include "gpnumlib/optimizer.h"

/// gpnummethods
#include "gpnummethods/amcmethod.h"
#include "gpnummethods/amc_exercboundcalc.h"

CC_BEGIN_NAMESPACE( ARM )

const double MaxValue = 1e15;
const double MinValue = -1e15;

////////////////////////////////////////////////////
///	Class  : ARM_AMCAndersen
///	Routine: ARM_AMCAndersen
///	Returns:
///	Action : Default constructor
////////////////////////////////////////////////////
ARM_AMCAndersen::ARM_AMCAndersen( size_t ItersNb, bool sortedMaximization  )
:	ARM_ExerciseBoundaryCalc(), itsUseSortedMaximisation(sortedMaximization)
{
	itsItersNb = ItersNb;
	SetName(ARM_AMCANDERSEN);
}

////////////////////////////////////////////////////
///	Class  : ARM_AMCAndersen
///	Routine: Clone,View, toString
///	Returns:
///	Action : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_AMCAndersen::Clone() const
{
	return new ARM_AMCAndersen(*this);
}


string ARM_AMCAndersen::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << " =======> Andersen METHOD <====== " << CC_NS(std,endl);
	os << indent << " Andersen Algortihm using strategy 1" << CC_NS(std,endl);
	os << indent << itsItersNb << "Iterations for exercise boundary calculation" << CC_NS(std,endl);
	os << indent << " Using ";
	if(itsUseSortedMaximisation) 
		os << " Sorted Maximisation algorithm";
	else 
		os << " Standard Brent algorithm";
	os << CC_NS(std,endl);
	return os.str();
}

////////////////////////////////////////////////////
///	Class  : ARM_AndersenCalculate
///	Routine: ComputeExerciseBoundary
///	Returns: void 
///	Action : to be given to Brent or another optimization
/// Algorithm. 
////////////////////////////////////////////////////

class ARM_AndersenCalculate
{
private: 
	ARM_VectorPtr itsTrigValue, itsPayoff, itsContOpt;
	size_t size;
public: 
	ARM_AndersenCalculate( const ARM_VectorPtr& trigValue, const ARM_VectorPtr& payoff, const ARM_VectorPtr& contOpt );
	double operator () ( double x ) const;
};


ARM_AndersenCalculate::ARM_AndersenCalculate( const ARM_VectorPtr& trigValue, const ARM_VectorPtr& payoff, const ARM_VectorPtr& contOpt )
: itsTrigValue(trigValue),
itsPayoff(payoff),
itsContOpt(contOpt),
size(payoff->size())
{
	
}

double ARM_AndersenCalculate::operator () ( double x ) const
{
	double result;
	unsigned int i;

	result = 0;

	for ( i = 0 ; i < size ; i ++ )
		result += ( itsTrigValue->Elt(i) >= x ) ? itsPayoff->Elt( i ) : itsContOpt->Elt( i );

	return - result;
}

////////////////////////////////////////////////////
///	Class  : ARM_AndersenReverseCalculate
///	Routine: ComputeExerciseBoundary
///	Returns: void 
///	Action : to be given to Brent or another optimization
/// Algorithm. 
////////////////////////////////////////////////////

class ARM_AndersenReverseCalculate
{
private: 
	ARM_VectorPtr itsTrigValue, itsPayoff, itsContOpt;
	size_t size;
public: 
	ARM_AndersenReverseCalculate( const ARM_VectorPtr& trigValue, const ARM_VectorPtr& payoff, const ARM_VectorPtr& contOpt );
	double operator () ( double x ) const;
};


ARM_AndersenReverseCalculate::ARM_AndersenReverseCalculate( const ARM_VectorPtr& trigValue, const ARM_VectorPtr& payoff, const ARM_VectorPtr& contOpt )
: itsTrigValue(trigValue),
itsPayoff(payoff),
itsContOpt(contOpt),
size(payoff->size())
{
}

double ARM_AndersenReverseCalculate::operator () ( double x ) const
{
	double result;
	unsigned int i;

	result = 0;

	for ( i = 0 ; i < size ; i ++ )
		result += ( itsTrigValue->Elt(i) <= x ) ? itsPayoff->Elt( i ) : itsContOpt->Elt( i );

	return - result;
}

////////////////////////////////////////////////////
///	Class  : ARM_AMCAndersen
///	Routine: ComputeStdExerciseBoundary
///	Returns: double
///	Action : Computes the exercise boundary using std Brent bisection maximization
////////////////////////////////////////////////////

pair<double,bool> ARM_AMCAndersen::ComputeStdExerciseBoundary( const ARM_VectorPtr& payoff, const ARM_VectorPtr& contOpt, const ARM_VectorPtr& trigValue )
{
	double maxopt, minopt, exerciseboundary;
	double plainExerciseBoundary, reverseExerciseBoundary;
	bool boundaryType;

	size_t size,i;

	// If the trig value argument is null we use the payoff
	ARM_AndersenCalculate aac((trigValue.IsNull()?payoff:trigValue), payoff, contOpt );
	ARM_AndersenReverseCalculate aarc((trigValue.IsNull()?payoff:trigValue), payoff, contOpt );

	/// finds upper and lower bounds for the optmization algorithm
	maxopt = (trigValue.IsNull()?payoff:trigValue)->Elt( 0 );
	minopt = (trigValue.IsNull()?payoff:trigValue)->Elt( 0 );
	size = (trigValue.IsNull()?payoff:trigValue)->size();

	for ( i = 0 ; i < size ; i ++ )
	{
		maxopt = MAX( maxopt, (trigValue.IsNull()?payoff:trigValue)->Elt( i ) );
		minopt = MIN( minopt, (trigValue.IsNull()?payoff:trigValue)->Elt( i ) );
	}

	if( (maxopt == minopt) && (minopt > MinValue) && (maxopt < MaxValue))
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Intrinsinc value seems to be constant");

	// To handle very large fees
	if( (maxopt == minopt) && ((minopt <= MinValue) || (maxopt >= MaxValue)) )
	{
		if (maxopt >=  MaxValue)
		{
			exerciseboundary = maxopt*1,0001;
			boundaryType = true;
		}

		if (minopt <=  MinValue)
		{
			exerciseboundary = minopt*1.0001;
			boundaryType = false;
		}
	}
	else
	{
		/// Does the optimization using Brent
		/// ExerciseBoundary of the form payoff >= a
		T_GoldenSectionOptimizer<ARM_AndersenCalculate> Optimizer( aac, minopt, 0.5*(minopt+maxopt), maxopt, 0.00001, 100 );
		/// ExerciseBoundary of the from payoff <= a
		T_GoldenSectionOptimizer<ARM_AndersenReverseCalculate> ReverseOptimizer( aarc, minopt, 0.5*(minopt+maxopt), maxopt, 0.00001, 100 );
		plainExerciseBoundary = Optimizer.Optimize();
		reverseExerciseBoundary = ReverseOptimizer.Optimize();

		boundaryType = (aac( plainExerciseBoundary ) <= aarc(reverseExerciseBoundary)) ? true : false;
		exerciseboundary = boundaryType ? plainExerciseBoundary : reverseExerciseBoundary;

		/// check if exerciseboundary == minopt or == maxopt.
		/// In such a case, the option is always/never exercised. 
		if( exerciseboundary <= minopt ) exerciseboundary = -1E20;
		if( exerciseboundary >= maxopt ) exerciseboundary = 1E20;
	}

	return pair<double,bool> (exerciseboundary,boundaryType);
}



////////////////////////////////////////////////////
///	Class  : ARM_AMCAndersen
///	Routine: ComputeExerciseBoundary
///	Returns: ARM_ExerciseBoundaryPtr
///	Action : Computes the exercise boundary using sorted maximization
////////////////////////////////////////////////////

pair<double,bool> ARM_AMCAndersen::ComputeSortedExerciseBoundary( const ARM_VectorPtr& payoff, const ARM_VectorPtr& contOpt, const ARM_VectorPtr& trigValue)
{	
	double	maxopt = (trigValue.IsNull()?payoff:trigValue)->Elt( 0 ),
			minopt = (trigValue.IsNull()?payoff:trigValue)->Elt( 0 ),
			sum	   = 0.0,
			sumCont= 0.0;

	/// use int because we use --i!
	int i, j, size = payoff->size();

	/// computes the sum of continuation value
	for ( i = 0 ; i < size ; i ++ )
	{
		maxopt = MAX( maxopt, (trigValue.IsNull()?payoff:trigValue)->Elt( i ) );
		minopt = MIN( minopt, (trigValue.IsNull()?payoff:trigValue)->Elt( i ) );
		sum   += (*payoff)[i];
		sumCont += (*contOpt)[i];
	}

	double exerciseboundary;
	bool exerciseType;

	/// validation for asburd minimization!
	if( maxopt == minopt  && (minopt > MinValue) && (maxopt < MaxValue))
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Intrinsinc value seems to be constant");


	// To handle very large fees
	if( (maxopt == minopt) && ((minopt <= MinValue) || (maxopt >= MaxValue)) )
	{
		if (maxopt >=  MaxValue)
		{
			exerciseboundary = maxopt*1,0001;
			exerciseType = true;
		}

		if (minopt <=  MinValue)
		{
			exerciseboundary = minopt*1.0001;
			exerciseType = false;
		}

	}
	else
	{
		/// sort in increasing order based on the trigger value!
		ARM_GP_VectorPtr trigValueSorted(static_cast<ARM_GP_Vector*>((trigValue.IsNull()?payoff:trigValue)->Clone()));
		ARM_GP_VectorPtr payoffSorted(static_cast<ARM_GP_Vector*>(payoff->Clone()));
		ARM_GP_VectorPtr contOptSorted(static_cast<ARM_GP_Vector*>(contOpt->Clone()));
		ARM_GP_VectorPtr indexSorted(new ARM_GP_Vector(size));

		for (i = 0; i < size; ++i)
			(*indexSorted)[i] = i;

		ARM_T_Sort<double,double>::sortTwoVectorsWithSameSize( *trigValueSorted, *indexSorted );

		for (i = 0; i < size; ++i)
		{
			(*payoffSorted)[i] = (*payoff)[(*indexSorted)[i]];
			(*contOptSorted)[i] = (*contOpt)[(*indexSorted)[i]];
		}

		exerciseboundary=minopt;
		double valueMax=sum, currentValue=sum, valueMin=sumCont, currentMinValue=sumCont;
		double doubleSum=sumCont+sum;


		ARM_GP_VectorPtr valuesMin( new ARM_GP_Vector(size ) );
		ARM_GP_VectorPtr valuesMax( new ARM_GP_Vector(size ) );

		int maxIndex = 0; // Index for Exercise if TrigValue >= a
		int minIndex = 0; // Index for Exercise if TrigValue <= a

		for( i=0; i<size; ++i )
		{
			currentValue += (*contOptSorted)[i]-(*payoffSorted)[i];
			currentMinValue= doubleSum-currentValue;
			(*valuesMax)[i] = currentValue;
			(*valuesMin)[i] = currentMinValue;

			if( currentValue > valueMax)
			{	
				valueMax		= currentValue; 
				maxIndex		= i;
			}
			if( currentMinValue > valueMin )
			{
				valueMin		= currentMinValue;
				minIndex		= i;
			}
		}

		exerciseType = (valueMax >  valueMin) ? true : false;
		
		exerciseboundary = exerciseType ? (*trigValueSorted)[maxIndex] : (*trigValueSorted)[minIndex];
		ARM_GP_VectorPtr values = exerciseType ? valuesMax : valuesMin;

		const int nbPoints = size/10;
		const int degre = 2;

		if ((maxIndex >= nbPoints) && (maxIndex <= size-nbPoints-1))
		{
			int Rows = 2*nbPoints+1;
			int Cols = degre+1;

			ARM_GP_Matrix Mat( Rows, Cols );
			ARM_GP_Vector Y(Rows);

			for ( i = 0 ; i <Rows ; i++ )
			{
				Y[i] = (*values)[i+maxIndex-nbPoints];
				for ( j = 0 ; j< Cols ; j++ )
					Mat(i,j) = pow((*trigValueSorted)[i+maxIndex-nbPoints],j);
			}

			ARM_GP_Vector* X = LeastSquareRegression( Mat, Y );

			double a = (*X)[2];
			double b = (*X)[1];

			exerciseboundary = -b/2/a;
		}

		/// to handle degenerated cases 
		if( exerciseboundary <= minopt ) exerciseboundary = -1E20;
		if( exerciseboundary >= maxopt ) exerciseboundary =  1E20;
	}

	return pair<double,bool> (exerciseboundary,exerciseType); 
}


////////////////////////////////////////////////////
///	Class  : ARM_AMCAndersen
///	Routine: ComputeExerciseBoundary
///	Returns: ARM_ExerciseBoundaryPtr
///	Action : Computes the exercise Boundary using Brent
////////////////////////////////////////////////////

ARM_ExerciseBoundary * ARM_AMCAndersen::ComputeExerciseBoundary( const ARM_VectorPtr& payoff, 
	const ARM_VectorPtr& contOpt,
	const ARM_GP_MatrixPtr& StatesVector )
{
	ARM_ExerciseBoundary * exercbound;

	if (!payoff.IsNull()&&!contOpt.IsNull()&&!StatesVector.IsNull())
	{
		// We extract the first row as a trigger value for andersen
		ARM_GP_VectorPtr firstCol(NULL);
		if (StatesVector->cols())
			firstCol = ARM_GP_VectorPtr(StatesVector->GetColumn(0));

		pair<double,bool> exerciseboundary = exerciseboundary = itsUseSortedMaximisation? 
			ComputeSortedExerciseBoundary( payoff, contOpt, firstCol)
		:	ComputeStdExerciseBoundary(payoff, contOpt, firstCol);

		/// Create an ExerciseBoundary Ptr which will be returned
		exercbound = new ARM_AndersenExerciseBoundary( exerciseboundary.first, exerciseboundary.second );
	}
	// Create an empy exercise boundary
	else
	{
		exercbound = new ARM_AndersenExerciseBoundary( 0.0, true );
	}

	return exercbound;
}


CC_END_NAMESPACE()
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
