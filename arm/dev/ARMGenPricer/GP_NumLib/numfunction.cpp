/*!
 *
 * Copyright (c) IXIS CIB November 2004 Paris
 *
 *	\file numfunction.cpp
 *
 *  \brief implementation of ARM_VectorPtrDbleFuncLinearCombination
 *	\author  A. Schauly
 *	\version 1.0
 *	\date November 2004
 */

#include "gpnumlib/numfunction.h"
#include "gpbase/gpvector.h"
#include <cmath>

CC_BEGIN_NAMESPACE( ARM )

///////////////////////////////////////////////////////////////////////////
///	Class  : struct ARM_VectorPtrDbleFuncLinearCombination
///	Routine: constructor
///	Returns:
///	Action : 
///////////////////////////////////////////////////////////////////////////

ARM_VectorPtrDbleFuncLinearCombination::ARM_VectorPtrDbleFuncLinearCombination
	( const ARM_VectorPtr& Coefficients, 
		const ARM_VectorPtrDbleFuncPtrVector& FuncVect )
{
	itsCoefficients = Coefficients;
	itsFuncVect = FuncVect;
}

///////////////////////////////////////////////////////////////////////////
///	Class  : struct ARM_VectorPtrDbleFuncLinearCombination
///	Routine: operator()
///	Returns:
///	Action : Returns sum( itscoeffs(i)*itsfuncvec(i)(x) )
///////////////////////////////////////////////////////////////////////////

double ARM_VectorPtrDbleFuncLinearCombination::operator ()( ARM_VectorPtr x ) const
{
	size_t size = itsCoefficients->size();
	size_t i; 
	double result = 0;

	#if defined( __GP_STRICT_VALIDATION )
	/// itsCoefficients and itsFuncVect must have the same size 
		if ( size != itsFuncVect.size() )
		{
			throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + "coefficients and funcvect sizes do not match !" );
		}
	#endif
	
	for ( i = 0 ; i < size ; i++ )
	{
		result += itsCoefficients->Elt(i) * ( * itsFuncVect[i]) ( x );
	}

	return result;
}

///////////////////////////////////////////////////////////////////////////
///	Class  : struct ARM_Monomial
///	Routine: operator()
///	Returns:
///	Action : prod( x_i^j_i )
///////////////////////////////////////////////////////////////////////////

double ARM_Monomial::operator ()( ARM_VectorPtr x ) const
{
	double result;
	size_t size = itsRiseToPower->size();
	size_t i;

	if ( size != x->size() )
	{
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			" size of x and monomial are suppose to be the same" );
	}

	result = 1;

	for ( i = 0 ; i < size ; i++ )
	{
		result *= pow( (*x)[i], (*itsRiseToPower)[i] );
	}

	return result;
}

///////////////////////////////////////////////////////////////////////////
///	Class  : struct ARM_Monomial
///	Routine: constructor
///	Returns:
///	Action : 
///////////////////////////////////////////////////////////////////////////

ARM_Monomial::ARM_Monomial( const ARM_IntVectorPtr& RiseToPower )
:itsRiseToPower( RiseToPower ){ }

CC_END_NAMESPACE()

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
