/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file cfmethod.cpp
 *  \brief
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date May 2005
 */


/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpnummethods/cfmethod.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/curve.h"
#include "gpbase/env.h"

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/discretisationscheme.h"
#include "gpinfra/pricerinfo.h"
#include "gpinfra/exerciseboundary.h"

//// gpnumlib
#include "gpnumlib/random.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_CFMethod
///	Routine: Constructor
///	Returns:
///	Action : Standard ARM object support
////////////////////////////////////////////////////

ARM_CFMethod::ARM_CFMethod(	CFMethodType cf_method,  const ARM_GP_Matrix& parameters)
	:
	itsCFmethod(cf_method),
	itsCFparameters(parameters)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_CFMethod
///	Routine: copy constructor
///	Returns:
///	Action : 
////////////////////////////////////////////////////

ARM_CFMethod::ARM_CFMethod(const ARM_CFMethod& rhs)
:
itsCFmethod(rhs.itsCFmethod),
itsCFparameters(rhs.itsCFparameters)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_CFMethod
///	Routine: Clone,View, toString
///	Returns:
///	Action : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_CFMethod::Clone() const
{
	return new ARM_CFMethod(*this);
}


string ARM_CFMethod::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << " =======> CF METHOD <====== " << CC_NS(std,endl);
	os << indent << itsCFparameters.toString(indent) << "\n";
	/// bucket pricer
	return os.str();
}

////////////////////////////////////////////////////
///	Class  : ARM_CFMethod
///	Routine: Init, ReInit
///	Returns: 
///	Action : Initialiation of the method (model is non const because of the postInit
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_CFMethod::Init( ARM_PricingModel& model, double firstInductTime )
{
	return ARM_CFMethod::ReInit(model);
}



/// ReInit is the function to be used before doing a for loop at each beginning
/// of the loop!
ARM_PricingStatesPtr ARM_CFMethod::ReInit( const ARM_PricingModel& model)
{
	ARM_PricingStatesPtr states( new ARM_PricingStates(1,1));
	return states;
}

////////////////////////////////////////////////////
///	Class  : ARM_CFMethod
///	Routine: ReInitLoop
///	Returns: ARM_PricingStatesPtr
///	Action : ReInitLoop
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_CFMethod::ReInitLoop(const ARM_PricingModel& model)
{ 
	return *(new ARM_PricingStatesPtr( NULL )); 
}


////////////////////////////////////////////////////
///	Class  : ARM_CFMethod
///	Routine: Induct
///	Returns: 
///	Action : induct from one time to another!
////////////////////////////////////////////////////

ARM_PricingStatesPtr ARM_CFMethod::Induct(
    const ARM_PricingModel& model,
	ARM_PricingStatesPtr& states,
	double toTime)
{    
	return states;
}



////////////////////////////////////////////////////
///	Class  : ARM_CFMethod
///	Routine: ComputeExercise
///	Returns: an exercise boundary
///	Action : compute an exercise boundary of an
/// exercise node.
////////////////////////////////////////////////////
ARM_ExerciseBoundary * ARM_CFMethod::ComputeExercise(const ARM_GP_VectorPtr& payoff, const ARM_GP_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector )
{
	return NULL;
}

////////////////////////////////////////////////////
///	Class  : ARM_CFMethod
///	Routine: GetBuckets,GetBucketIndex
///	Returns: 
///	Action : returns the buckets and the index
////////////////////////////////////////////////////
ARM_VectorPtr ARM_CFMethod::GetBuckets() const
{
	return ARM_VectorPtr(new ARM_GP_Vector(1,1.0)); 
}


size_t ARM_CFMethod::GetBucketIndex() const
{
	return 0;
}


////////////////////////////////////////////////////
///	Class  : ARM_CFMethod
///	Routine: CreatePricerInfo
///	Returns: ARM_PricerInfo*
///	Action : creates the corresponding pricer info
////////////////////////////////////////////////////
ARM_PricerInfo* ARM_CFMethod::CreatePricerInfo( const ARM_PricingModel& model ) const
{
	//To be Changed
	return new ARM_PInfo_SingleLoop;
}

////////////////////////////////////////////////////
///	Class  : ARM_CFMethod
///	Routine: GetSpotProbabilities
///	Returns: ARM_GP_MatrixPtr
///	Action : Return probas to reach a future state from spot date
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_CFMethod::GetSpotProbabilities(const ARM_GP_Vector& eventTimes) const
{
    return ARM_GP_MatrixPtr(NULL);
}


////////////////////////////////////////////////////
///	Class  : ARM_CFMethod
///	Routine: GetArrowDebreuPrices 
///	Returns: ARM_GP_VectorPtr
///	Action : returns the Arrow Debreu prices at the slice timeIdx
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_CFMethod::GetArrowDebreuPrices(size_t timeIdx, const ARM_PricingModel& model ) const
{
	return ARM_GP_VectorPtr(NULL);
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

