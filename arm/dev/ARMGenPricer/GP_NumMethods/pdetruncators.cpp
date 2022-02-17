/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pdetrunactors.cpp
 *  \brief
 *
 *	\author  A Schauly
 *	\version 1.0
 *	\date July 2005
 */

#include "gpnummethods/pdetruncators.h"
#include "gpinfra/pricingstates.h"
#include <math.h>

CC_USING_NS(std,vector)
#include <vector>

CC_BEGIN_NAMESPACE( ARM )

/////////////////////////////////////////////////////////////////////////////////
//// Global for all truncators
/////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : PDE_Truncator
///	Routine: TruncateTranstionMatrixes
///	Returns: void
///	Action : troncate tridiagonal transition Matrix
///    and adds new boundary conditions (giben by the two last parameters)
////////////////////////////////////////////////////
void PDE_Truncator::TruncateTranstionMatrixes( size_t toTimeIdx, const ARM_GP_VectorPtr& UpperTerms, const ARM_GP_VectorPtr& LowerTerms, const ARM_GP_VectorPtr& DiagTerms, 
	const ARM_GP_VectorPtr& newUpperLeftLimitConditions, const ARM_GP_VectorPtr& newLowerRightLimitConditions )
{
	size_t newSize = getToTimeSize( toTimeIdx ), 
		oldSize = DiagTerms->size();

	if( newSize != oldSize )
	{
		size_t i=(oldSize-newSize)/2;

		ARM_GP_Vector::iterator 
			iter11 = DiagTerms->begin()+1,
			iter12 = iter11 + i,
			iter21 = UpperTerms->begin(), 
			iter22 = iter21 + i,
			iter31 = LowerTerms->begin(), 
			iter32 = iter31 + i;

		i=1;

		/// translation of the values
		for( ; i < newSize-1 ; ++i, ++iter11, ++iter21, ++iter31,++iter12, ++iter22, ++iter32 )
		{
			(*iter11) = (*iter12);
			(*iter21) = (*iter22);
			(*iter31) = (*iter32);
		}

		/// New boundary conditions
		UpperTerms->Elt(0) = newUpperLeftLimitConditions->Elt(0);
		DiagTerms->Elt(0) = newUpperLeftLimitConditions->Elt(1);
		LowerTerms->Elt(0) = newUpperLeftLimitConditions->Elt(2);
		LowerTerms->Elt(newSize-2) = newLowerRightLimitConditions->Elt(2);
		DiagTerms->Elt(newSize-1) = newLowerRightLimitConditions->Elt(1);
		UpperTerms->Elt(newSize-2) = newLowerRightLimitConditions->Elt(0);

		/// resizing matrix
		DiagTerms->resize( newSize );
		LowerTerms->resize( newSize-1 );
		UpperTerms->resize( newSize-1 );
	}

}

/////////////////////////////////////////////////////////////////////////////////
//// For Dummy truncator
//// This truncator basically does not truncate
/////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : PDE_DummyTruncator
///	Routine: toString
///	Returns: string
///	Action : returns description of the object
////////////////////////////////////////////////////
string PDE_DummyTruncator::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << " Dummy Truncator" << CC_NS(std,endl);
	os << indent << "     Does nothing " << CC_NS(std,endl);

	return os.str();
}

////////////////////////////////////////////////////////////////////////////////
///// End Dummy troncator
////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////
//// Truncation for PDE1F
////  Truncation for one factor models
/////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////
///	Class  : PDE1F_Truncator
///	Routine: TruncatePricingStates
///	Returns: void
///	Action : truncate pricingstates (nummethod states, payoffsstates, payoffs)
////////////////////////////////////////////////////
void PDE1F_Truncator::TruncatePricingStates( size_t toTimeIdx, const ARM_PricingStatesPtr& states )
{
	size_t PayoffsSize = states->GetPayoffsSize();
	size_t PayoffStatesSize = states->GetPayoffStatesSize();
	size_t i=0;

	/// NumMethod states
	TruncateMatrix(toTimeIdx, *(states->GetNumMethodStates()) );

	/// PayoffStates
	for( i=0 ; i<PayoffStatesSize ; ++i)
	{
		ARM_PayoffStates& payoffStates = states->GetPayoffStates(i);

		ARM_GP_Matrix& intermediatePayoffs = const_cast<ARM_GP_Matrix&> (payoffStates.GetIntermediatePayoffs() );
		TruncateMatrix( toTimeIdx, intermediatePayoffs );

		ARM_GP_Vector vec = payoffStates.GetPayoffs();
		TruncateVector( toTimeIdx, vec );
		payoffStates.SetPayoffs( vec );
	}

	/// Payoffs
	for( i=0 ; i<PayoffsSize ; ++i)
	{
		ARM_VectorPtr vec = states->GetPayoffVec(i);
		TruncateVector( toTimeIdx, *vec );

		for( int j=0; j<vec->size() ; ++j)
			states->SetPayoff( j, i, vec->Elt(j) );

	}
	states->resizePayoffs( states->GetPayoffs().GetRowsNb(), getToTimeSize( toTimeIdx ) );

}


////////////////////////////////////////////////////
///	Class  : PDE1F_Truncator
///	Routine: TruncateVector
///	Returns: void
///	Action : truncate vector (from pricing states)
////////////////////////////////////////////////////
void PDE1F_Truncator::TruncateVector( size_t toTimeIdx, ARM_GP_Vector& vec )
{
	size_t newSize = getToTimeSize( toTimeIdx ), 
		oldSize = vec.size();

	if( newSize != oldSize )
	{
		size_t i=(oldSize-newSize)/2;

		ARM_GP_Vector::iterator 
			iter1 = vec.begin(),
			iter2 = iter1 + i;

		i=0;

		/// Elts translation
		for( ; i < newSize ; ++i, ++iter1, ++iter2 )
			(*iter1) = (*iter2);

		/// Resiznig
		vec.resize( newSize );
	}

}


////////////////////////////////////////////////////
///	Class  : PDE1F_Truncator
///	Routine: TruncateMatrix
///	Returns: void
///	Action : truncate matrix (from pricing states)
////////////////////////////////////////////////////
void PDE1F_Truncator::TruncateMatrix( size_t toTimeIdx, ARM_GP_Matrix& matrix )
{
	size_t newSize = getToTimeSize( toTimeIdx ), 
		oldSize = matrix.cols();

	if( newSize != oldSize )
	{
		size_t reductionSize=(oldSize-newSize)/2;
		size_t rowsNb = matrix.rows();
		size_t i=0,j=0;
		ARM_GP_Vector::iterator iter1, iter2;

		/// For each line of the matrix
		for( ; j < rowsNb ; ++j)
		{
			iter1 = matrix.begin() + j*oldSize;
			iter2 = iter1+reductionSize;

			/// Elts translation
			for( i=0 ; i<newSize;++i, ++iter1, ++iter2)
				(*iter1) = (*iter2);
		}

		/// Resizing
		matrix.resize( rowsNb, newSize );
	}
}

////////////////////////////////////////////////////
///	Class  : PDE1F_Truncator
///	Routine: toString
///	Returns: string
///	Action : returns description of the object
////////////////////////////////////////////////////
string PDE1F_Truncator::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << " Dummy 1FTruncator" << CC_NS(std,endl);
	os << indent << " New Sizes are:" << CC_NS(std,endl);
	os << indent; 

	ARM_GP_T_Vector<size_t>::const_iterator iter = getToTimeSizes().begin(), 
		iterEnd = getToTimeSizes().end();


	for( ; iter != iterEnd ; ++iter )
		os << (*iter) << " "; 

	os << CC_NS(std,endl);

	return os.str();
}

/////////////////////////////////////////////////////////////////////////////////
//// End Truncation for PDE1F
/////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////
//// Begin Truncation for PDE2F
////  Truncation for one factor models
////   WARNING: Does not work well !!!
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : PDE2F_Truncator
///	Routine: TruncatePricingStates
///	Returns: void
///	Action : truncate pricingstates
////////////////////////////////////////////////////
void PDE2F_Truncator::TruncatePricingStates( size_t toTimeIdx, const ARM_PricingStatesPtr& states )
{
	size_t PayoffsSize = states->GetPayoffsSize();
	size_t PayoffStatesSize = states->GetPayoffStatesSize();
	size_t i=0;

	/// NumMethod states
	TruncateMatrix(toTimeIdx, *(states->GetNumMethodStates()) );

	/// PayoffStates
	for( i=0 ; i<PayoffStatesSize ; ++i)
	{
		ARM_PayoffStates& payoffStates = states->GetPayoffStates(i);

		ARM_GP_Matrix& intermediatePayoffs = const_cast<ARM_GP_Matrix&> (payoffStates.GetIntermediatePayoffs() );
		TruncateMatrix( toTimeIdx, intermediatePayoffs );

		ARM_GP_Vector vec = payoffStates.GetPayoffs();
		TruncateVector( toTimeIdx, vec );
		payoffStates.SetPayoffs( vec );
	}

	/// Payoffs
	for( i=0 ; i<PayoffsSize ; ++i)
	{
		ARM_VectorPtr vec = states->GetPayoffVec(i);
		TruncateVector( toTimeIdx, *vec );

		for( int j=0; j<vec->size() ; ++j)
			states->SetPayoff( j, i, vec->Elt(j) );

	}

}


////////////////////////////////////////////////////
///	Class  : PDE2F_Truncator
///	Routine: TruncateVector
///	Returns: void
///	Action : truncate vector (from pricing states)
////////////////////////////////////////////////////
void PDE2F_Truncator::TruncateVector( size_t toTimeIdx, ARM_GP_Vector& vec )
{
	size_t newSize = getToTimeSize( toTimeIdx ), 
		oldSize = sqrt(static_cast<double>(vec.size()));

	if( newSize != oldSize )
	{
		size_t newSizeSquare = newSize*newSize;
		size_t Jump = oldSize-newSize;
		size_t i=Jump/2,j=0;

		ARM_GP_Vector::iterator 
			iter1 = vec.begin(),
			iter2 = iter1 + i*(oldSize+1);

		i=0;

		/// Elts translation. This is more complicated in the present case than in 1F
		for( ; i < newSize ; ++i )
		{
			for( ; j < newSize ; ++j, ++iter1, ++iter2 )
				(*iter1) = (*iter2);

			iter2+=Jump;
		}

		/// vector resizing
		vec.resize( newSizeSquare );
	}
}


////////////////////////////////////////////////////
///	Class  : PDE1F_Truncator
///	Routine: TruncateMatrix
///	Returns: void
///	Action : truncate matrix (from pricing states)
////////////////////////////////////////////////////
void PDE2F_Truncator::TruncateMatrix( size_t toTimeIdx, ARM_GP_Matrix& matrix )
{
	size_t newSize = getToTimeSize( toTimeIdx ), 
		oldSize = sqrt(static_cast<double>(matrix.cols()));

	if( newSize != oldSize )
	{
		size_t oldSizeSquare = matrix.cols();
		size_t newSizeSquare = newSize*newSize;
		size_t Jump = oldSize-newSize;
		size_t i=Jump/2,j=0;
		size_t rowsNb = matrix.rows();
		size_t k=0;
		ARM_GP_Vector::iterator iter1,iter2;
		iter1 = matrix.begin();

		for( ; k<rowsNb ; ++k)
		{
			i=Jump/2;
			iter2 = matrix.begin() + i*(oldSize+1)+k*oldSizeSquare;

			i=0;

		/// Elts translation. This is more complicated in the present case than in 1F
			for( ; i < newSize ; ++i )
			{
				for( ; j < newSize ; ++j, ++iter1, ++iter2 )
					(*iter1) = (*iter2);

				iter2+=Jump;
			}
		}

		/// Matrix resizing
		matrix.resize( rowsNb, newSizeSquare );
	}
}


////////////////////////////////////////////////////
///	Class  : PDE2F_Truncator
///	Routine: toString
///	Returns: string
///	Action : returns description of the object
////////////////////////////////////////////////////
string PDE2F_Truncator::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << " PDE 2FTruncator" << CC_NS(std,endl);
	os << indent << " New Sizes are:" << CC_NS(std,endl);
	os << indent; 

	ARM_GP_T_Vector<size_t>::const_iterator iter = getToTimeSizes().begin(), 
		iterEnd = getToTimeSizes().end();


	for( ; iter != iterEnd ; ++iter )
		os << (*iter) << " "; 

	os << CC_NS(std,endl);

	return os.str();
}
/////////////////////////////////////////////////////////////////////////////////
//// End Truncation for PDE2F
/////////////////////////////////////////////////////////////////////////////////

CC_END_NAMESPACE()