/*
 *
 * Copyright (c) IXIS CIB November 2004 Paris
 *
 *	\file exerciseboundary.cpp
 *	\author  A Schauly
 *	\version 1.0
 *	\date November 2004
 */


#include "gpbase/ostringstream.h"
#include "gpnumlib/typedef.h"
#include "gpbase/ostringstream.h"
#include "gpinfra/exerciseboundary.h"
#include "gpbase/checkarg.h"
#include "gpbase/gpmatrix.h"

#include "gpnumlib/regression.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_AndersenExerciseBoundary
///	Routine: Clone, toString
///	Returns:
///	Action : Standard ARM object support
////////////////////////////////////////////////////

string ARM_AndersenExerciseBoundary::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << " =======> Andersen Exercise Boundary <====== " << CC_NS(std,endl);
	os << indent << " No Comments yet" << CC_NS(std,endl);

	return os.str();
}

////////////////////////////////////////////////////
///	Class  : ARM_AndersenExerciseBoundary
///	Routine: GetExerciseBoundary
///	Returns: ARM_GP_VectorPtr&
///	Action : Exercise Boundary (in for of a vector)
////////////////////////////////////////////////////

ARM_GP_VectorPtr ARM_AndersenExerciseBoundary::GetExerciseBoundary( void ) const
{
	ARM_GP_Vector * vec = new ARM_GP_Vector(1);
	ARM_GP_VectorPtr * result = new ARM_GP_VectorPtr( vec );
	(*vec)[0] = itsValue;
	return *result;
}



////////////////////////////////////////////////////
///	Class  : ARM_LSExerciseBoundary
///	Routine: Clone, toString
///	Returns:
///	Action : Standard ARM object support
////////////////////////////////////////////////////


string ARM_LSExerciseBoundary::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << " =======> LS Exercise Boundary <====== " << CC_NS(std,endl);
	os << indent << " No Comments yet" << CC_NS(std,endl);

	return os.str();
}

////////////////////////////////////////////////////
///	Class  : ARM_LSExerciseBoundary
///	Routine: ComputeContinationValues
///	Returns: ARM_GP_VectorPtr & which contains pseudo 
///  continuation values computed according to the LS Algorithm
///	Action : computes the ARM_GP_VectorPtr& returned;
////////////////////////////////////////////////////

ARM_GP_VectorPtr ARM_LSExerciseBoundary::ComputePseudoContinuationValues( 
		const ARM_GP_MatrixPtr& StatesVector )
{
	return itsRegression->ComputeValues(StatesVector);
}


////////////////////////////////////////////////////
///	Class  : ARM_LSExerciseBoundary
///	Routine: EvalInPlace
///	Returns: void 
///	Action : 
////////////////////////////////////////////////////

void ARM_LSExerciseBoundary::EvalInPlace( const ARM_VectorPtr& payoff, const ARM_VectorPtr& contOpt, 
										 const ARM_GP_MatrixPtr& StatesVector )
{
	size_t size = payoff->size();
	size_t i;

	#ifdef __GP_STRICT_VALIDATION
		if ( StatesVector == ARM_GP_MatrixPtr( NULL ) )
			throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
            "SatesVector is null " );

		if ( StatesVector->rows() != size || StatesVector->cols() != itsRegression->GetFactorNb() 
			       || contOpt->size() != size )
			throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
            "SatesVector, itsContinuationValueCoefficients and payoff are supposed to be the same size. " );
	#endif

	ARM_VectorPtr PseudoContinuationValues( ComputePseudoContinuationValues( StatesVector ) );

	for ( i=0 ; i < size ; i++ )
		(*payoff)[i ] = (*payoff)[i ] >= ( *PseudoContinuationValues )[i] ? (*payoff)[i ] : (*contOpt)[i ];

	return;
}


////////////////////////////////////////////////////
///	Class  : ARM_LSExerciseBoundary
///	Routine: GetExerciseBoundary
///	Returns: ARM_GP_VectorPtr 
///	Action : Return the result of the regression
////////////////////////////////////////////////////

ARM_GP_VectorPtr ARM_LSExerciseBoundary::GetExerciseBoundary( void ) const
{
	ARM_LSRegression* reg = dynamic_cast<ARM_LSRegression*>(&*itsRegression);

	if (reg)
		return reg->GetRegressedCoeffs();
	else
		return ARM_GP_VectorPtr(NULL);
}

////////////////////////////////////////////////////
///	Class  : ARM_AndersenExerciseBoundary
///	Routine: EvalInPlace
///	Returns: void 
///	Action : 
////////////////////////////////////////////////////

void ARM_AndersenExerciseBoundary::EvalInPlace( const ARM_VectorPtr& payoff, const ARM_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector )
{
	// We extract the first row as a trigger value for andersen
	ARM_GP_VectorPtr trigValue(NULL);
	if (StatesVector->cols())
		trigValue = ARM_GP_VectorPtr(StatesVector->GetColumn(0));
	else
	// By default it uses the payoff
		trigValue = ARM_GP_VectorPtr( new ARM_GP_Vector(*payoff));

	size_t size, i;

	size = payoff->size();

	#ifdef __GP_STRICT_VALIDATION
		if ( size != contOpt->size() )  // compares payoff->size() with contOpt->size()
			throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
            "vectors payoff and contOpt are supposed to be the same size! " );
	#endif

	/// Replaces payoff values with new continuation values	
	if( itsExerciseIfPayoffIsGreaterThanValue )
	{
		for ( i = 0 ; i < size ; i ++ )
			(*payoff)[i ] = trigValue->Elt( i ) >= itsValue ? (*payoff)[i ] : (*contOpt)[i ];
	}
	else
	{
		for ( i = 0 ; i < size ; i ++ )
			(*payoff)[i ] = trigValue->Elt( i ) <= itsValue ? (*payoff)[i ] : (*contOpt)[i ];
	}

	return;
}

////////////////////////////////////////////////////
///	Class  : ARM_SmoothTreeExerciseBoundary
///	Routine: EvalInPlace
///	Returns: void 
///	Action : 
////////////////////////////////////////////////////

void ARM_SmoothTreeExerciseBoundary::EvalInPlace( const ARM_VectorPtr& payoff, const ARM_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector )
{
	size_t size, i;

	size = payoff->size();

	#ifdef __GP_STRICT_VALIDATION
		if ( size != contOpt->size() )  // compares payoff->size() with contOpt->size()
			throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
            "vectors payoff and contOpt are supposed to be the same size! " );
	#endif

	/// Do the smoothed resulting payoff	
	for ( i = 0 ; i < size ; i ++ )
		(*payoff)[i] = (*itsSmoothValues)[i] + ((*payoff)[i] >= (*contOpt)[i] ? (*payoff)[i] : (*contOpt)[i]);

	return;
}


CC_END_NAMESPACE()
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
