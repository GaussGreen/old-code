/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file binummethod.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpnummethods/binummethod.h"

/// gpbase
#include "gpbase/ostringstream.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/pricerinfo.h"
#include "gpinfra/exerciseboundary.h"

/// numerical method should evaluate in reverse order
/// in Monte Carlo should go last row to first row
/// in backward trees/PDES should go first row to last row
/// recursive call will make sure to call appropriate methods

/// order of backwardisation/MC is the following
/// initialisation
/// browse ParseTree and execute accordingly
/// tree should have only one cashflow!


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_BINumMethod
///	Routine: Constructor
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////
ARM_BINumMethod::ARM_BINumMethod( ): ARM_NumMethod()
{
	CC_ARM_SETNAME( ARM_BINUMMETHOD );
}


////////////////////////////////////////////////////
///	Class  : ARM_BINumMethod
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_BINumMethod::toString(const string& indent, const string& nextIndent) const
{ 
	return indent + string( "Backward Induction Numerical Method\n" ); 
}


////////////////////////////////////////////////////
///	Class  : ARM_BINumMethod
///	Routine: Clone method
///	Returns: a new copy of this
///	Action : Copy this
////////////////////////////////////////////////////
ARM_Object* ARM_BINumMethod::Clone() const
{
	return new ARM_BINumMethod( *this );
}

////////////////////////////////////////////////////
///	Class  : ARM_BINumMethod
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_BINumMethod::~ARM_BINumMethod()
{}


////////////////////////////////////////////////////
///	Class  : ARM_BINumMethod
///	Routine: Copy Constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_BINumMethod::ARM_BINumMethod( const ARM_BINumMethod& rhs)
: ARM_NumMethod( rhs )
{}


////////////////////////////////////////////////////
///	Class  : ARM_BINumMethod
///	Routine: Assignment operator
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_BINumMethod& ARM_BINumMethod::operator=( const ARM_BINumMethod& rhs )
{
	if( this !=	 &rhs )
	{
	    ARM_NumMethod::operator = ( rhs );
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_BINumMethod
///	Routine: Init
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

ARM_PricingStatesPtr ARM_BINumMethod::Init(ARM_PricingModel&, double firstInductTime)
{
	/// nothing yet
    return ARM_PricingStatesPtr(new ARM_PricingStates(1,1));
}



////////////////////////////////////////////////////
///	Class  : ARM_BINumMethod
///	Routine: induction function
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

ARM_PricingStatesPtr ARM_BINumMethod::Induct(
        const ARM_PricingModel& model,	/// model
		ARM_PricingStatesPtr& states,	/// states,
		double	toTime  )
{
	/// nothing as everything is done in the process paid payoff!
	return states;
}

////////////////////////////////////////////////////
///	Class  : ARM_FINumMethod
///	Routine: ComputeExercise
///	Returns: an exercise boundary
///	Action : compute an exercise boundary of an
/// exercise node.
////////////////////////////////////////////////////
ARM_ExerciseBoundary * ARM_BINumMethod::ComputeExercise(const ARM_GP_VectorPtr& payoff, const ARM_GP_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector)
{
	return NULL;
}


////////////////////////////////////////////////////
///	Class  : ARM_BINumMethod
///	Routine: ReInit
///	Returns: 
///	Action : reinitialisation for loop pricing! not supported
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_BINumMethod::ReInit( const ARM_PricingModel& model)
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"Binum method does not support multiple loops!" );
}

////////////////////////////////////////////////////
///	Class  : ARM_BINumMethod
///	Routine: ReInitLoop
///	Returns: ARM_PricingStatesPtr
///	Action : Reinit for loop change (in case of pricing 
/// direction change)
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_BINumMethod::ReInitLoop(const ARM_PricingModel& model)
{ 
	return *(new ARM_PricingStatesPtr( NULL )); 
}


////////////////////////////////////////////////////
///	Class  : ARM_BINumMethod
///	Routine: GetBuckets,GetBucketIndex
///	Returns: 
///	Action : returns the buckets and the index
////////////////////////////////////////////////////
ARM_VectorPtr ARM_BINumMethod::GetBuckets() const
{
	/// one bucket with size 1
	return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,1)); 
}

size_t ARM_BINumMethod::GetBucketIndex() const
{
	/// first bucket only
	return 0;
}



////////////////////////////////////////////////////
///	Class  : ARM_BINumMethod
///	Routine: CreatePricerInfo
///	Returns: ARM_PricerInfo*
///	Action : creates the corresponding pricer info
////////////////////////////////////////////////////
ARM_PricerInfo* ARM_BINumMethod::CreatePricerInfo( const ARM_PricingModel& model ) const
{
	return new ARM_PInfo_SingleLoop;
}


////////////////////////////////////////////////////
///	Class  : ARM_BINumMethod
///	Routine: GetSpotProbabilities
///	Returns: ARM_GP_MatrixPtr
///	Action : Return probas to reach a future state from spot date
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_BINumMethod::GetSpotProbabilities(const ARM_GP_Vector& eventTimes) const
{
    ARM_GP_Matrix* probas = new ARM_GP_Matrix(1,eventTimes.size());

    for(size_t i=0;i<eventTimes.size();++i)
        (*probas)(0,i)=1.0;

    return ARM_GP_MatrixPtr(probas);
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

