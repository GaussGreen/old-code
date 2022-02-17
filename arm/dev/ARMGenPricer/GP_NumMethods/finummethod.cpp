/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file finummethod.cpp
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2003
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpnummethods/finummethod.h"

/// gpbase
#include "gpbase/ostringstream.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/pricerinfo.h"
#include "gpinfra/exerciseboundary.h"


/// proto type numerical method to test 
/// forward induction

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_FINumMethod
///	Routine: Constructor
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////
ARM_FINumMethod::ARM_FINumMethod( ): ARM_NumMethod()
{
	CC_ARM_SETNAME( ARM_FINUMMETHOD );
}


////////////////////////////////////////////////////
///	Class  : ARM_FINumMethod
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_FINumMethod::toString(const string& indent, const string& nextIndent) const
{ 
	return indent + string( "Forward Induction Numerical Method\n" ); 
}


////////////////////////////////////////////////////
///	Class  : ARM_FINumMethod
///	Routine: Clone method
///	Returns: a new copy of this
///	Action : Copy this
////////////////////////////////////////////////////
ARM_Object* ARM_FINumMethod::Clone() const
{
	return new ARM_FINumMethod( *this );
}


////////////////////////////////////////////////////
///	Class  : ARM_FINumMethod
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_FINumMethod::~ARM_FINumMethod()
{}



////////////////////////////////////////////////////
///	Class  : ARM_FINumMethod
///	Routine: Copy Constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_FINumMethod::ARM_FINumMethod( const ARM_FINumMethod& rhs)
: ARM_NumMethod( rhs )
{
}


////////////////////////////////////////////////////
///	Class  : ARM_FINumMethod
///	Routine: Assignment operator
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_FINumMethod& ARM_FINumMethod::operator=( const ARM_FINumMethod& rhs )
{
	if( this !=	 &rhs )
	    ARM_NumMethod::operator = ( rhs );
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_FINumMethod
///	Routine: Init
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

ARM_PricingStatesPtr ARM_FINumMethod::Init( ARM_PricingModel&, double firstInductTime )
{
    return ARM_PricingStatesPtr(new ARM_PricingStates(1,1));
}

////////////////////////////////////////////////////
///	Class  : ARM_FINumMethod
///	Routine: induction function
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

ARM_PricingStatesPtr ARM_FINumMethod::Induct(
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
ARM_ExerciseBoundary * ARM_FINumMethod::ComputeExercise(const ARM_GP_VectorPtr& payoff, const ARM_GP_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector )
{
	return NULL;
}


////////////////////////////////////////////////////
///	Class  : ARM_FINumMethod
///	Routine: ReInitLoop
///	Returns: ARM_PricingStatesPtr
///	Action : Reinit for loop change (in case of pricing 
/// direction change)
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_FINumMethod::ReInitLoop(const ARM_PricingModel& model)
{ 
	return *(new ARM_PricingStatesPtr( NULL )); 
}

////////////////////////////////////////////////////
///	Class  : ARM_FINumMethod
///	Routine: ReInit
///	Returns: 
///	Action : reinitialisation for loop pricing! not supported
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_FINumMethod::ReInit( const ARM_PricingModel& model)
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"forward induction method does not support multiple loops!" );
}

////////////////////////////////////////////////////
///	Class  : ARM_FINumMethod
///	Routine: GetBuckets,GetBucketIndex
///	Returns: 
///	Action : returns the buckets and the index
////////////////////////////////////////////////////
ARM_VectorPtr ARM_FINumMethod::GetBuckets() const
{
	/// one bucket with size 1
	return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,1)); 
}

size_t ARM_FINumMethod::GetBucketIndex() const
{
	/// first bucket only
	return 0;
}



////////////////////////////////////////////////////
///	Class  : ARM_FINumMethod
///	Routine: CreatePricerInfo
///	Returns: ARM_PricerInfo*
///	Action : creates the corresponding pricer info
////////////////////////////////////////////////////
ARM_PricerInfo* ARM_FINumMethod::CreatePricerInfo( const ARM_PricingModel& model ) const
{
	return new ARM_PInfo_SingleLoop;
}

////////////////////////////////////////////////////
///	Class  : ARM_FINumMethod
///	Routine: GetSpotProbabilities
///	Returns: ARM_GP_MatrixPtr
///	Action : Return probas to reach a future state from spot date
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_FINumMethod::GetSpotProbabilities(const ARM_GP_Vector& eventTimes) const
{
    ARM_GP_Matrix* probas = new ARM_GP_Matrix(1,eventTimes.size());

    for(size_t i=0;i<eventTimes.size();++i)
        (*probas)(0,i)=1.0;

    return ARM_GP_MatrixPtr(probas);
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

