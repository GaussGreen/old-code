/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pricingstates.cpp
 *
 *  \brief pricing states summarizes the state of the  world
 *
 *	\author  main=E. Benhamou, secondary= JM Prie
 *	\version 1.0
 *	\date October 2003
 */


#include "gpinfra/pricingstates.h"

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_PricingStates
///	Routine: Constructor
///	Returns: 
///	Action : Constructor allocating memory
////////////////////////////////////////////////////
ARM_PricingStates::ARM_PricingStates(
	size_t nbStates,
	size_t nbModelStates,
    size_t nbPayoffs,
    size_t nbNumMethodStates,
    size_t nbIntermediatePayoffs,
    size_t nbSnapshots,
	size_t nbPayoffStates,
    bool   otherPayoffsFlag,
	size_t nbProbaChanges)
:
itsModelStates( new ARM_GP_Matrix( nbModelStates,nbStates) ),
itsPayoffs(nbPayoffs,nbStates), 
itsProbaChanges(nbProbaChanges,nbStates),
itsNumMethodStates( new ARM_GP_Matrix( nbNumMethodStates,nbStates) ),
itsPayoffStatesVector(nbPayoffStates)
{
	itsPricingStatesContextVector = ARM_PricingStatesContextPtrVectorPtr( new ARM_PricingStatesContextPtrVector( 0 ) );
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingStates
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_PricingStates::ARM_PricingStates(const ARM_PricingStates& rhs)
:	ARM_RootObject(rhs), 
	itsModelStates( rhs.itsModelStates ),
	itsNumMethodStates( rhs.itsNumMethodStates ),
    itsPayoffs( rhs.itsPayoffs ),
	itsProbaChanges ( rhs.itsProbaChanges ),
	itsPricingStatesContextVector( rhs.itsPricingStatesContextVector ),
    itsPayoffStatesVector(rhs.itsPayoffStatesVector)
{}


////////////////////////////////////////////////////
///	Class  : ARM_PricingStates
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_PricingStates::~ARM_PricingStates()
{}


////////////////////////////////////////////////////
///	Class  : ARM_PricingStates
///	Routine: Affectation
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_PricingStates& ARM_PricingStates::operator = (const ARM_PricingStates& rhs)
{
	if(this != &rhs)
	{
		ARM_RootObject::operator=(rhs);
		itsModelStates					= rhs.itsModelStates;
		itsNumMethodStates				= rhs.itsNumMethodStates;
		itsPayoffs						= rhs.itsPayoffs;
		itsProbaChanges					= rhs.itsProbaChanges;
		itsPricingStatesContextVector	= rhs.itsPricingStatesContextVector;
        itsPayoffStatesVector			= rhs.itsPayoffStatesVector;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingStates
///	Routines: Clone,View
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_PricingStates::Clone() const
{
	return new ARM_PricingStates(*this);
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingStates
///	Routines: GetCopyOfPricingStatesContextVector
///	Returns : ARM_PricingStatesContextPtrVectorPtr
///	Action  : copies itsPricingStatesContextVector's elements
/// (not the pointers, but the underlying elements). 
////////////////////////////////////////////////////
const ARM_PricingStatesContextPtrVectorPtr ARM_PricingStates::GetCopyOfPricingStatesContextVector() const
{
	ARM_PricingStatesContextPtrVector * NewPricingStatesContextPtrVector = new ARM_PricingStatesContextPtrVector( 0 );
	ARM_PricingStatesContextPtrVector::iterator endIter = itsPricingStatesContextVector->end();

	for( ARM_PricingStatesContextPtrVector::iterator iter = itsPricingStatesContextVector->begin() ; iter != endIter ; iter ++ )
		NewPricingStatesContextPtrVector->push_back( ARM_PricingStatesContextPtr( static_cast<ARM_PricingStatesContext*>((*iter)->Clone()) ) );

	return ARM_PricingStatesContextPtrVectorPtr( NewPricingStatesContextPtrVector );
}

////////////////////////////////////////////////////
///	Class   : ARM_PricingStates
///	Routines: AddPricingStatesContextVector
///	Returns : void 
///	Action  : enqueues every element from pricingStatesContextVector
/// in itsPricingStatesContextVector. 
////////////////////////////////////////////////////
void ARM_PricingStates::AddPricingStatesContextVector( const ARM_PricingStatesContextPtrVectorPtr& pricingStatesContextVector ) const
{
	ARM_PricingStatesContextPtrVector::iterator endIter = pricingStatesContextVector->end();

	for( ARM_PricingStatesContextPtrVector::iterator iter = pricingStatesContextVector->begin() ; iter != endIter ; iter++ )
		itsPricingStatesContextVector->push_back( *iter );
}

////////////////////////////////////////////////////
///	Class   : ARM_PricingStates
///	Routines: resizePayoffStatesVector
///	Returns : void 
///	Action  : Resize payoff. 
////////////////////////////////////////////////////
void ARM_PricingStates::resizePayoffStatesVector(size_t nbPayoffStates) 
{
	itsPayoffStatesVector.resize(nbPayoffStates);
	for (size_t i = 0; i < nbPayoffStates; ++i)
		itsPayoffStatesVector[i].resizePayoff(size());
}

////////////////////////////////////////////////////
///	Class   : ARM_PayoffStates
///	Routines: SetOtherPayoffsFlag
///	Returns : 
///	Action  : Accessor. 
////////////////////////////////////////////////////

void ARM_PricingStates::SetOtherPayoffsFlag(bool otherPayoffsFlag)
{
	for (size_t i = 0; i < itsPayoffStatesVector.size(); ++i)
		itsPayoffStatesVector[i].SetOtherPayoffsFlag(otherPayoffsFlag);
}


////////////////////////////////////////////////////
///	Class   : ARM_PayoffStates
///	Routines: SetIVFlag
///	Returns : 
///	Action  : Accessor. 
////////////////////////////////////////////////////

void ARM_PricingStates::SetIVFlag(bool otherPayoffsFlag)
{
	for (size_t i = 0; i < itsPayoffStatesVector.size(); ++i)
		itsPayoffStatesVector[i].SetIVFlag(otherPayoffsFlag);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class   : ARM_PayoffStates
///	Routines: ARM_PayoffStates
///	Returns : 
///	Action  : Constructor. 
////////////////////////////////////////////////////

ARM_PayoffStates::ARM_PayoffStates(
size_t nbStates,
size_t nbIntermediatePayoffs,
size_t nbSnapshots,
bool otherPayoffsFlag,
bool ivFlag)
: itsIntermediatePayoffs(nbIntermediatePayoffs,nbStates),
itsPayoffSnapshots(nbSnapshots),
itsOtherPayoffsFlag(otherPayoffsFlag),
itsIVFlag(ivFlag)
{
}


////////////////////////////////////////////////////
///	Class   : ARM_PayoffStates
///	Routines: ARM_PayoffStates
///	Returns : 
///	Action  : Copy Constructor. 
////////////////////////////////////////////////////
ARM_PayoffStates::ARM_PayoffStates(const ARM_PayoffStates& rhs)
: ARM_RootObject(rhs)
{
	if (this != &rhs)
	{
		itsPayoffs = rhs.itsPayoffs;
		itsIntermediatePayoffs = rhs.itsIntermediatePayoffs;
		itsPayoffSnapshots = rhs.itsPayoffSnapshots;
		itsOtherPayoffsFlag = rhs.itsOtherPayoffsFlag;
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_PayoffStates
///	Routine: Affectation
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_PayoffStates& ARM_PayoffStates::operator=(const ARM_PayoffStates& rhs)
{
	if (this != &rhs)
	{
		ARM_RootObject::operator=(rhs);
		itsPayoffs = rhs.itsPayoffs;
		itsIntermediatePayoffs = rhs.itsIntermediatePayoffs;
		itsPayoffSnapshots = rhs.itsPayoffSnapshots;
		itsOtherPayoffsFlag = rhs.itsOtherPayoffsFlag;
		itsIVFlag = rhs.itsIVFlag;
	}

	return *this;
}

////////////////////////////////////////////////////
///	Class   : ARM_PayoffStates
///	Routines: ~ARM_PayoffStates
///	Returns : 
///	Action  : Destructor. 
////////////////////////////////////////////////////
ARM_PayoffStates::~ARM_PayoffStates()
{
}

////////////////////////////////////////////////////
///	Class   : ARM_PayoffStates
///	Routines: Clone
///	Returns : 
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_PayoffStates::Clone() const
{
	return new ARM_PayoffStates(*this);
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/

/*---- End of file ----*/



