/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: mixtenummethod.cpp,v $
 * Revision 1.1  2003/12/30 16:39:43  ebenhamou
 * Initial revision
 *
 *
 *
 */


/*! \file mixtenummethod.cpp
 *
 *  \brief
 *
 *	\author  R Guillemot
 *	\version 1.0
 *	\date Novemmber 2004
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpnummethods/mixtenummethod.h"

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
///	Class  : ARM_MixteNumMethod
///	Routine: Constructor
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////
ARM_MixteNumMethod::ARM_MixteNumMethod( ): ARM_NumMethod(),
itsPricingDirCurrLoop(ARM_NumMethod::GP_FWDBCKWDLOOKING)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_MixteNumMethod
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_MixteNumMethod::toString(const string& indent,const string& nextIndent) const
{ 
	return "\n" + indent + string( "Mixte Induction Numerical Method\n" ); 
}


////////////////////////////////////////////////////
///	Class  : ARM_MixteNumMethod
///	Routine: Clone method
///	Returns: a new copy of this
///	Action : Copy this
////////////////////////////////////////////////////
ARM_Object* ARM_MixteNumMethod::Clone() const
{
	return new ARM_MixteNumMethod( *this );
}


////////////////////////////////////////////////////
///	Class  : ARM_MixteNumMethod
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_MixteNumMethod::~ARM_MixteNumMethod()
{}



////////////////////////////////////////////////////
///	Class  : ARM_MixteNumMethod
///	Routine: Copy Constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_MixteNumMethod::ARM_MixteNumMethod( const ARM_MixteNumMethod& rhs)
: ARM_NumMethod( rhs ),
itsPricingDirCurrLoop( rhs.itsPricingDirCurrLoop )
{}


////////////////////////////////////////////////////
///	Class  : ARM_MixteNumMethod
///	Routine: Assignment operator
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_MixteNumMethod& ARM_MixteNumMethod::operator=( const ARM_MixteNumMethod& rhs )
{
	if( this !=	 &rhs )
	{
	    ARM_NumMethod::operator = ( rhs );
		
		itsPricingDirCurrLoop = rhs.itsPricingDirCurrLoop;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_MixteNumMethod
///	Routine: Init
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

ARM_PricingStatesPtr ARM_MixteNumMethod::Init( ARM_PricingModel&, double firstInductTime )
{
	/// nothing yet
	itsPricingDirCurrLoop = ARM_NumMethod::GP_FWDLOOKING;
    return ARM_PricingStatesPtr(new ARM_PricingStates(1,1));
}

////////////////////////////////////////////////////
///	Class  : ARM_MixteNumMethod
///	Routine: induction function
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

ARM_PricingStatesPtr ARM_MixteNumMethod::Induct(
        const ARM_PricingModel& model,	/// model
		ARM_PricingStatesPtr& states,	/// states,
		double	toTime  )
{
	/// nothing as everything is done in the process paid payoff!
	return states;
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeMethod
///	Routine: ComputeExercise
///	Returns: An exercise boundary
///	Action : Compute an exercise boundary
///          of an exercise node
////////////////////////////////////////////////////
ARM_ExerciseBoundary* ARM_MixteNumMethod::ComputeExercise(const ARM_GP_VectorPtr& payoff, const ARM_GP_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector )
{
	return NULL;
}

////////////////////////////////////////////////////
///	Class  : ARM_MixteNumMethod
///	Routine: ReInit
///	Returns: 
///	Action : reinitialisation for loop pricing! not supported
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_MixteNumMethod::ReInit( const ARM_PricingModel& model)
{
	itsPricingDirCurrLoop = ARM_NumMethod::GP_BCKWDLOOKING;
	return ARM_PricingStatesPtr(new ARM_PricingStates(1,1));
}

////////////////////////////////////////////////////
///	Class  : ARM_MixteNumMethod
///	Routine: ReInitLoop
///	Returns: ARM_PricingStatesPtr
///	Action : Reinit for loop change (in case of pricing 
/// direction change)
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_MixteNumMethod::ReInitLoop(const ARM_PricingModel& model)
{ 
	itsPricingDirCurrLoop = ARM_NumMethod::GP_BCKWDLOOKING;
	return ARM_PricingStatesPtr(new ARM_PricingStates(1,1));
}

////////////////////////////////////////////////////
///	Class  : ARM_MixteNumMethod
///	Routine: ComputeTimeSteps
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_MixteNumMethod::ComputeTimeSteps(ARM_PricingModel& model)
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "ComputeTimeSteps : unimplemented function !");
}

////////////////////////////////////////////////////
///	Class  : ARM_MixteNumMethod
///	Routine: GetBuckets,GetBucketIndex
///	Returns: 
///	Action : returns the buckets and the index
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MixteNumMethod::GetBuckets() const
{
	/// one bucket with size 1
	return static_cast<ARM_VectorPtr>(new std::vector<double>(2,1)); 
}

size_t ARM_MixteNumMethod::GetBucketIndex() const
{
	/// first bucket only
	return 0;
}



////////////////////////////////////////////////////
///	Class  : ARM_MixteNumMethod
///	Routine: CreatePricerInfo
///	Returns: ARM_PricerInfo*
///	Action : creates the corresponding pricer info
////////////////////////////////////////////////////
ARM_PricerInfo* ARM_MixteNumMethod::CreatePricerInfo( const ARM_PricingModel& model ) const
{
	return new ARM_PInfo_SingleLoop;
}

////////////////////////////////////////////////////
///	Class  : ARM_MixteNumMethod
///	Routine: GetSpotProbabilities
///	Returns: ARM_GP_MatrixPtr
///	Action : Return probas to reach a future state from spot date
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_MixteNumMethod::GetSpotProbabilities(const std::vector<double>& eventTimes) const
{
    ARM_GP_Matrix* probas = new ARM_GP_Matrix(1,eventTimes.size());

    for(size_t i=0;i<eventTimes.size();++i)
        (*probas)(0,i)=1.0;

    return ARM_GP_MatrixPtr(probas);
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

