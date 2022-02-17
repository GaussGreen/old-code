/*
 *
 * Copyright (c) IXIS CIB November 2004 Paris
 *
 *	\file amcmethod.cpp
 *  \brief
 *
 *	\author  A Schauly
 *	\version 1.0
 *	\date November 2004
 */


/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpnummethods/amcmethod.h"

// gpinfra
#include "gpinfra/typedef.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/gramfunctorbase.h"
#include "gpinfra/gramfunctorsimple.h"

//gpnummethod
#include "gpnummethods/impsampler.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_AMCMethod
///	Routine: Default constructor
///	Returns:
///	Action : Note that itsLoopNb = 2
////////////////////////////////////////////////////
ARM_AMCMethod::ARM_AMCMethod( size_t itersNb, const ARM_RandomGeneratorPtrVector& randGenVector, ARM_SamplerBase* sampler, ARM_ExerciseBoundaryCalc* exerBoundCalc, size_t MaxBucketSize, const ARM_ImpSamplerPtr& impSampler, const ARM_PathSchemePtr& pathScheme ) 
: ARM_MCMethod( itersNb, randGenVector, sampler, MaxBucketSize, impSampler, pathScheme  ),
itsExerBoundCalc(NULL)
{
	// Number of iteration after exercise boundary has been computed/for ex computation ; nb of buckets
	size_t itersNbNoExerc, itersNbExerc, bucketsNb, EvenBucketSize; 

	itsExerBoundCalc = ARM_ExerciseBoundaryCalcPtr(static_cast<ARM_ExerciseBoundaryCalc*>(exerBoundCalc->Clone()));
	itsPricingDirCurrLoop = ARM_NumMethod::GP_FWDLOOKING;
	itsLoopNb = 2;


	//////////////////////////////////////////////////////////////////////////////////////
	///// This part is for the design of buckets
	//////////////////////////////////////////////////////////////////////////////////////

	EvenBucketSize = 0xfffffffe & MaxBucketSize; 
	/// (We want to have buckets with an even number of simulations only)

	/// Number of trajectories used for exercise boundary
	itersNbExerc = 0xfffffffe & (exerBoundCalc->getItsItersNb());
	/// Remaining trajectories
	itersNbNoExerc = 0xfffffffe & MAX(0, ( itersNb-itersNbExerc ) );

	/// Number of Buckets (complications due to the fact that itersNbNoExerc is an unsigned int)
	bucketsNb = itersNbNoExerc > 1 ? itersNbNoExerc-1 : 0;
	bucketsNb = 2+ bucketsNb/EvenBucketSize;
	itsBuckets = ARM_VectorPtr( new ARM_GP_Vector( bucketsNb ) );

	/// The first bucket has the number of trajectories specified for exercise boundary computation. 
	*(itsBuckets->begin()) = itersNbExerc;
	
	/// There are EvenBucketSize paths in a bucket, from the second to the last-1 bucket. 
	for( ARM_GP_Vector::iterator bucketIter = itsBuckets->begin()+1 ; bucketIter != itsBuckets->end() ; bucketIter++ )
		(*bucketIter) = EvenBucketSize;

	/// The last bucket has the remaining trajectories. 
	if( itersNbNoExerc % EvenBucketSize != 0 ) 
		*(itsBuckets->end()-1) = 0xfffffffe & (itersNbNoExerc % EvenBucketSize);
	else
		if( itersNbNoExerc == 0 )
			*(itsBuckets->end()-1) = 0;

}
////////////////////////////////////////////////////
///	Class  : ARM_AMCMethod
///	Routine: Default constructor with 1 random gen
///	Returns:
///	Action : 
////////////////////////////////////////////////////
ARM_AMCMethod::ARM_AMCMethod( size_t itersNb, const ARM_RandomGeneratorPtr& randGen, ARM_SamplerBase* sampler, ARM_ExerciseBoundaryCalc* exerBoundCalc, size_t MaxBucketSize, const ARM_ImpSamplerPtr& impSampler, const ARM_PathSchemePtr& pathScheme  )
: ARM_MCMethod( itersNb, randGen, sampler, MaxBucketSize, impSampler, pathScheme ),
itsExerBoundCalc(NULL)
{
	// Number of iteration after exercise boundary has been computed/for ex computation ; nb of buckets
	size_t itersNbNoExerc, itersNbExerc, bucketsNb, EvenBucketSize; 

	itsExerBoundCalc = ARM_ExerciseBoundaryCalcPtr(static_cast<ARM_ExerciseBoundaryCalc*>(exerBoundCalc->Clone()));
	itsPricingDirCurrLoop = ARM_NumMethod::GP_FWDLOOKING;
	itsLoopNb = 2;


	//////////////////////////////////////////////////////////////////////////////////////
	///// This part is for the design of buckets
	//////////////////////////////////////////////////////////////////////////////////////

	EvenBucketSize = 0xfffffffe & MaxBucketSize; 
	/// (We want to have buckets with an even number of simulations only)

	/// Number of trajectories used for exercise boundary
	itersNbExerc = 0xfffffffe & (exerBoundCalc->getItsItersNb());
	/// Remaining trajectories
	itersNbNoExerc = 0xfffffffe & MAX(0, ( itersNb-itersNbExerc ) );

	/// Number of Buckets (complications due to the fact that itersNbNoExerc is an unsigned int)
	bucketsNb = itersNbNoExerc > 1 ? itersNbNoExerc-1 : 0;
	bucketsNb = 2+ bucketsNb/EvenBucketSize;
	itsBuckets = ARM_VectorPtr( new ARM_GP_Vector( bucketsNb ) );

	/// The first bucket has the number of trajectories specified for exercise boundary computation. 
	*(itsBuckets->begin()) = itersNbExerc;
	
	/// There are EvenBucketSize paths in a bucket, from the second to the last-1 bucket. 
	for( ARM_GP_Vector::iterator bucketIter = itsBuckets->begin()+1 ; bucketIter != itsBuckets->end() ; bucketIter++ )
		(*bucketIter) = EvenBucketSize;

	/// The last bucket has the remaining trajectories. 
	if( itersNbNoExerc % EvenBucketSize != 0 ) 
		*(itsBuckets->end()-1) = 0xfffffffe & (itersNbNoExerc % EvenBucketSize);
	else
		if( itersNbNoExerc == 0 )
			*(itsBuckets->end()-1) = 0;

}

////////////////////////////////////////////////////
///	Class  : ARM_AMCMethod
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_AMCMethod::ARM_AMCMethod(const ARM_AMCMethod& rhs) : ARM_MCMethod( rhs ),
itsExerBoundCalc(), itsPricingDirCurrLoop( rhs.itsPricingDirCurrLoop ), itsLoopNb( rhs.itsLoopNb ), 
modelStatesVector( rhs.modelStatesVector ), pricingStatesContextVector( rhs.pricingStatesContextVector ),
itsLastInductTime( rhs.itsLastInductTime ), itsPrevLastInductTime( rhs.itsPrevLastInductTime), 
iter( rhs.iter ), iter2( rhs.iter2 )
{
	itsExerBoundCalc = ARM_ExerciseBoundaryCalcPtr(static_cast<ARM_ExerciseBoundaryCalc*>(rhs.itsExerBoundCalc->Clone()));
}

////////////////////////////////////////////////////
///	Class  : ARM_AMCMethod
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_AMCMethod::~ARM_AMCMethod()
{
}


////////////////////////////////////////////////////
///	Class  : ARM_AMCMethod
///	Routine: Induct
///	Returns: ARM_PricingStatesPtr
///	Action : Induct
////////////////////////////////////////////////////

ARM_PricingStatesPtr ARM_AMCMethod::Induct(const ARM_PricingModel& model,ARM_PricingStatesPtr& states,double toTime )
{
	if ( itsPricingDirCurrLoop == ARM_NumMethod::GP_FWDLOOKING )
	{ // in the fwd part, we do like in MCMethod and we save all the model states. 
		ARM_PricingStatesPtr pricingStates( ARM_MCMethod::Induct( model, states, toTime ) );
		//ARM_GP_MatrixPtr mat( static_cast<ARM_GP_Matrix*>( pricingStates->GetModelStates()->Clone() ) );
		ARM_GP_MatrixPtr mat;
		modelStatesVector.push_back( mat ); // modelStates are stored
		ARM_PricingStatesContextPtrVectorPtr context( pricingStates->GetCopyOfPricingStatesContextVector()  );
		pricingStatesContextVector.push_back( context );

		itsLastInductTime = toTime;

		return pricingStates;
	}
	
	return InductBackward( model, states, toTime ); // FIXME
}

////////////////////////////////////////////////////
///	Class  : ARM_AMCMethod
///	Routine: InductBackward
///	Returns: ARM_PricingStatesPtr
///	Action : InductBackWard
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_AMCMethod::InductBackward(const ARM_PricingModel& model,ARM_PricingStatesPtr& states,double toTime )
{ 
	/// nothing to do if we are making no progress
    if( iter == (modelStatesVector.rend()-1) )
        return states; 

	// we use the previously computed modelStates for the backward loop
	// !! Take care, we start on the before last model states (the last one 
	// is never used). 
	if( iter < modelStatesVector.rend()-1 )
	{
		iter++;
		iter2++;
	}
	states->SetModelStates( *iter );
	states->SetPricingStatesContextVector( *iter2 );

	int lastTimeIdx		= GetLastTimeIdx();


	ARM_GP_Vector* timeSteps = GetTimeSteps();
	while ((lastTimeIdx > 0) && ((*timeSteps)[lastTimeIdx] > toTime))
	{
		model.GetNumeraire()->MoveNumeraireBckwd();
		lastTimeIdx--;
	}
	SetLastTimeIdx(lastTimeIdx);

	
	model.GetNumeraire()->Update(model,states,lastTimeIdx-1);

	return states;
}

////////////////////////////////////////////////////
///	Class  : ARM_AMCMethod
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Init
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_AMCMethod::Init( ARM_PricingModel& model, double firstInductTime )
{ 
	// The first loop is forward
	itsPricingDirCurrLoop = ARM_NumMethod::GP_FWDLOOKING;
	// modelStatesVector is initialized at the beginning of the method
	// so that this nummeth can be used many times
	modelStatesVector.resize(0);
	pricingStatesContextVector.resize(0);
	return ARM_MCMethod::Init(model,firstInductTime);
}

////////////////////////////////////////////////////
///	Class  : ARM_AMCMethod
///	Routine: ReInitLoop
///	Returns: ARM_PricingStatesPtr
///	Action : Reinit for loop change. Changes
/// pricing Direction
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_AMCMethod::ReInitLoop(const ARM_PricingModel& model)
{
	// The second loop is backward
	itsPricingDirCurrLoop = ARM_NumMethod::GP_BCKWDLOOKING;
	iter = modelStatesVector.rbegin();
	iter2 = pricingStatesContextVector.rbegin();
	ARM_PricingStatesPtr states( model.FirstPricingStates( (*itsBuckets)[itsBucketIndex] ) );

	/// resets the numeraire
	ARM_NumerairePtr numeraire = model.GetNumeraire();

#if defined(__GP_STRICT_VALIDATION)
	if( numeraire == ARM_NumerairePtr(NULL) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			": Numeraire should not be NULL");
#endif
	numeraire->ResetLoop();

	return Induct( model, states, itsLastInductTime );
}

////////////////////////////////////////////////////
///	Class  : ARM_AMCMethod
///	Routine: ReInit
///	Returns: ARM_PricingStatesPtr
///	Action : Init
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_AMCMethod::ReInit( const ARM_PricingModel& model)
{
	itsPricingDirCurrLoop = ARM_NumMethod::GP_FWDLOOKING;
	modelStatesVector.resize(0);
	pricingStatesContextVector.resize(0);
	return ARM_MCMethod::ReInit( model );
}

////////////////////////////////////////////////////
///	Class  : ARM_AMCMethod
///	Routine: ComputeExercise
///	Returns: An exercise boundary
///	Action : Compute an exercise boundary
///          of an exercise node
////////////////////////////////////////////////////
ARM_ExerciseBoundary * ARM_AMCMethod::ComputeExercise(const ARM_GP_VectorPtr& payoff, const ARM_GP_VectorPtr& contOpt, 
													   const ARM_GP_MatrixPtr& StatesVector )
{
	return itsExerBoundCalc->ComputeExerciseBoundary( payoff,contOpt,StatesVector );
}

////////////////////////////////////////////////////
///	Class  : ARM_AMCMethod
///	Routine: Finalize
///	Returns: void
///	Action : frees memory at end of pricing
////////////////////////////////////////////////////
void ARM_AMCMethod::Finalize()
{
	pricingStatesContextVector.clear();
	modelStatesVector.clear();
}


////////////////////////////////////////////////////
///	Class  : ARM_AMCMethod
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_AMCMethod::toString(const string& indent, const string& nextIndent) const
{
	string str = ARM_MCMethod::toString();
	CC_Ostringstream os;

	os << indent << " This MCMethod has been changed to an American MCMethod " << CC_NS(std,endl);
	os << indent << " " << CC_NS(std,endl);
	os << indent << " =======> Exerise Boundary Calculator <====== " << CC_NS(std,endl);
	str += os.str();
	str += itsExerBoundCalc->toString();

	return str;
}



////////////////////////////////////////////////////
///	Class  : ARM_AMCMethod
///	Routine: ModifyParseTree
///	Returns: void 
///	Action : modifies parse tree. 
/// Changes MAX(a,PV(b) ) into Exercise(a,b)
////////////////////////////////////////////////////
void ARM_AMCMethod::ModifyParseTree( ARM_PricingAdviser * pricingAdviser ) const
{
	CC_STL_VECTOR( ARM_ExpNodePtr ) * PVNodes = pricingAdviser->getPVNodes();

	if( PVNodes->size() != 0 )
		        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : no PV allowed in deal description if it is priced with AMC. " );

}

////////////////////////////////////////////////////
///	Class  : ARM_AMCMethod
///	Routine: NeedToCreateDefaultArgument
///	Returns: bool
///	Action : tells for the exercise function if requires to creae default argument
////////////////////////////////////////////////////
bool ARM_AMCMethod::NeedToCreateDefaultArgument() const
{
	return itsExerBoundCalc->NeedToCreateDefaultArgument();
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

