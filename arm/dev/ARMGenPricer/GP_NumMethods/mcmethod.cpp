/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file mcmethod.cpp
 *  \brief
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpnummethods/mcmethod.h"

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
#include "gpinfra/pricerinfo.h"
#include "gpinfra/exerciseboundary.h"

//// gpnumlib
#include "gpnumlib/random.h"

/// gpnummethods
//#include "gpnummethods/typedef.h"
#include "gpnummethods/sampler.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/impsampler.h"

#include <utility>
using namespace std;


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_MCMethod
///	Routine: Default constructor
///	Returns:
///	Action : Constructor, default is to use superBucket!
////////////////////////////////////////////////////
ARM_MCMethod::ARM_MCMethod( 
size_t itersNb,
const ARM_RandomGeneratorPtrVector& randGenVector, 
ARM_SamplerBase* sampler,
size_t MaxBucketSize,
ARM_ImpSamplerPtr impSampler,
ARM_PathSchemePtr pathScheme)
:
ARM_NumMethod(sampler),
	itsRandGenVector(NULL),
	itsFirstInductTime(0),
	itsBuckets(),
	itsBucketIndex(-1),
	itsOtherPayoffsFlag(true)
{
	// By default we use the dummy importance sampler
	if(!impSampler.IsNull())
		itsImpSampler = ARM_ImpSamplerPtr(static_cast<ARM_ImpSampler*>(impSampler->Clone()));
	else
		itsImpSampler = ARM_ImpSamplerPtr(new ARM_DummyImpSampler);

	if(!pathScheme.IsNull())
		itsPathScheme = ARM_PathSchemePtr(static_cast<ARM_PathScheme*>(pathScheme->Clone()));
	else
		itsPathScheme = ARM_PathSchemePtr(new ARM_IncrementalPathScheme);

	Initialize(itersNb,randGenVector,MaxBucketSize);
}

////////////////////////////////////////////////////
///	Class  : ARM_MCMethod
///	Routine: Default constructor
///	Returns:
///	Action : Constructor, default is to use superBucket!
////////////////////////////////////////////////////
ARM_MCMethod::ARM_MCMethod( 
size_t itersNb, 
const ARM_RandomGeneratorPtr& randGen, 
ARM_SamplerBase* sampler,
size_t MaxBucketSize,
ARM_ImpSamplerPtr impSampler,
ARM_PathSchemePtr pathScheme )
:	ARM_NumMethod(sampler),
	itsRandGenVector(NULL),
	itsFirstInductTime(0),
	itsBuckets(),
	itsBucketIndex(-1),
	itsOtherPayoffsFlag(true)
{
	// By default we use the dummy importance sampler
	if(!impSampler.IsNull())
		itsImpSampler = ARM_ImpSamplerPtr(static_cast<ARM_ImpSampler*>(impSampler->Clone()));
	else
		itsImpSampler = ARM_ImpSamplerPtr(new ARM_DummyImpSampler);

	if(!pathScheme.IsNull())
		itsPathScheme = ARM_PathSchemePtr(static_cast<ARM_PathScheme*>(pathScheme->Clone()));
	else
		itsPathScheme = ARM_PathSchemePtr(new ARM_IncrementalPathScheme);

	ARM_RandomGeneratorPtrVector randGenVector(1);
	randGenVector[0] = randGen;
	Initialize(itersNb,randGenVector,MaxBucketSize);
}


////////////////////////////////////////////////////
///	Class  : ARM_MCMethod
///	Routine: Initialize
///	Returns:
///	Action : Continue the construction of the MC Method
////////////////////////////////////////////////////
void ARM_MCMethod::Initialize(
size_t itersNb,
const ARM_RandomGeneratorPtrVector& randGenVector,
size_t MaxBucketSize )
{
	CC_Ostringstream os;

	int size = randGenVector.size();
	itsRandGenVector.resize(size);
	for(int i=0;i<size;i++)
		itsRandGenVector[i] = ARM_RandomGeneratorPtr((ARM_RandomGenerator*)(randGenVector[i]->Clone()));
	
	size_t bucketsNb;
#if defined(__GP_STRICT_VALIDATION)
	//if( itersNb & 1 )
	//	ARM_THROW( ERR_INVALID_ARGUMENT, "nb of iteration is odd...can only be even!" );
#endif

	for(int j=0;j<itsRandGenVector.size();j++)
	{
		if( itsRandGenVector[j]->GetDistributionType() != ARM_RandomGenerator::ARM_Normal )
			ARM_THROW( ERR_INVALID_ARGUMENT, "waiting for a normal generator!" );
	}

	os << "MaxBucketSize = " << MaxBucketSize << std::endl;
	os << "IterNb = " << itersNb << std::endl;

	if( MaxBucketSize < 1 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "waiting for a max bucket size > 0" );

	/*if( itersNb >= 100000) 
		itsOtherPayoffsFlag = false;*/

	size_t EvenBucketSize = MaxBucketSize & 0xfffffffe;
	/// & 0xfffffffe : make sure that the number of path is even by doing and bitwise

	os << "EvenBucketSize = " << EvenBucketSize << std::endl;

	/// Computes the number of buckets needed (complications since itersNb is an unsigned int)
	bucketsNb = (itersNb > 1) ? itersNb-1 : 0;
	bucketsNb = bucketsNb/EvenBucketSize + 1;
	itsBuckets = ARM_VectorPtr( new ARM_GP_Vector(bucketsNb) );

	for ( ARM_GP_Vector::iterator iter = itsBuckets->begin() ; iter != itsBuckets->end() ; iter++ )
		*iter = EvenBucketSize;

	if( (itersNb % EvenBucketSize) != 0 ) *(itsBuckets->end()-1) = ( itersNb % EvenBucketSize ) & 0xfffffffe;
	
	SetName(ARM_MCMETHOD);
}
////////////////////////////////////////////////////
///	Class  : ARM_MCMethod
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_MCMethod::ARM_MCMethod(const ARM_MCMethod& rhs)
:	ARM_NumMethod(rhs),
itsRandGenVector(rhs.itsRandGenVector), 
itsFirstInductTime(rhs.itsFirstInductTime),
itsBuckets(rhs.itsBuckets),
itsBucketIndex(rhs.itsBucketIndex),
itsOtherPayoffsFlag(rhs.itsOtherPayoffsFlag),
itsImpSampler(ARM_ImpSamplerPtr(static_cast<ARM_ImpSampler*>(rhs.itsImpSampler->Clone()))),
itsPathScheme(ARM_PathSchemePtr(static_cast<ARM_PathScheme*>(rhs.itsPathScheme->Clone())))
{
}


////////////////////////////////////////////////////
///	Class  : ARM_MCMethod
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_MCMethod::~ARM_MCMethod()
{
}



////////////////////////////////////////////////////
///	Class  : ARM_MCMethod
///	Routine: operator =
///	Returns: itself
///	Action : Affectation
////////////////////////////////////////////////////
ARM_MCMethod& ARM_MCMethod::operator=(const ARM_MCMethod& rhs)
{
	if( this != & rhs )
	{
		this->~ARM_MCMethod();
        new (this) ARM_MCMethod(rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_MCMethod
///	Routine: Clone,View, toString
///	Returns:
///	Action : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_MCMethod::Clone() const
{
	return new ARM_MCMethod(*this);
}


string ARM_MCMethod::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << " =======> MC METHOD <====== " << CC_NS(std,endl);
	os << indent << " Pricing direction      : " << GP_PricingDirectionTxt[(int)GetPricingDirection()] << CC_NS(std,endl);
	os << indent << " Nb of break points     : " << GetNbSteps() << CC_NS(std,endl);
	os << indent << " Time Steps			 : ";
	
	if( GetTimeSteps() )
		os << GetTimeSteps()->toString();
	else
		os << "Not initialized!\n";

	for(int i=0;i<itsRandGenVector.size();i++)
	{
		os << indent << " Generator "<<i<<"       : " 
			<< itsRandGenVector[i]->toString() << CC_NS(std,endl);
	}
	os << indent << " First Induct time : " << itsFirstInductTime << CC_NS(std,endl);
	os << indent << " Bucket description: " << itsBuckets->toString( indent, nextIndent );
	os << indent << " Sampler: "			<< GetSampler()->toString() << CC_NS(std,endl);
	os << indent << " Scheduler: "			<< GetSampler()->GetScheduler()->toString() << CC_NS(std,endl);
	if (&*GetImpSampler())
		os << indent << " Importance Sampler: "	<< GetImpSampler()->toString(indent,nextIndent) << CC_NS(std,endl);
	if (&*GetPathScheme())
		os << indent << " Path Scheme: "	<< GetPathScheme()->toString(indent,nextIndent) << CC_NS(std,endl);

	/// bucket pricer
	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_MCMethod
///	Routine: GetRandGen
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_RandomGeneratorPtr ARM_MCMethod::GetRandGen(int i) const
{
	if(i>=itsRandGenVector.size())
		return itsRandGenVector[itsRandGenVector.size()-1];
	else
		return itsRandGenVector[i];
}

////////////////////////////////////////////////////
///	Class  : ARM_MCMethod
///	Routine: Init, ReInit
///	Returns: 
///	Action : Initialiation of the tree (model is non const because of the postInit
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_MCMethod::Init( ARM_PricingModel& model, double firstInductTime )
{
	// The first induct time is going from the generic security
	itsFirstInductTime = firstInductTime;

	/// 1) compute time steps: can it jump from one time to the next one?
	ComputeAndSetTimeSteps(model);
	
	/// 2) Compute local variances and std Dev (does the Cholesky decomposition within the model)
	model.ModelStateLocalVariancesAndStdDev(*(GetTimeSteps()));
	model.NumMethodStateLocalVariancesAndStdDev(*(GetTimeSteps()));

	/// 3) Initialise the sampler
	// We need to do that for communiation reason
	ARM_GP_VectorPtr timeStepCopy(static_cast<ARM_GP_Vector*>(GetTimeSteps()->Clone()));
	ARM_TimeStepsAndSlices* result = GetSampler()->Init(model,timeStepCopy,false);

	// Just used by the tree !!!
	delete result;
	
	/// 4) gives back the hand to the model for post init operation
	model.PostInit();

	/// 5) reset the vector of random nb generator
	ARM_SamplerNDBase* samplerNDBase = GetSampler()->ToSamplerNDBase();

	int totaldim = samplerNDBase->TotalDimension();

	size_t dim=samplerNDBase->dim();
	size_t timeSteps=GetTimeSteps()->size();
	int sizeRandGen = itsRandGenVector.size();
	int i=0;
	while (sizeRandGen>1)
	{
		ARM_GP_T_Vector<size_t> nbOfPathList(1);
		nbOfPathList[0] += (*itsBuckets)[i];
		itsRandGenVector[i]->reset(totaldim, nbOfPathList, dim );
		i++;
		sizeRandGen--;
	}
	int sizeFinalRandom = itsBuckets->size()-i;
	ARM_GP_T_Vector<size_t> nbOfPathList(sizeFinalRandom);
	for(size_t h=0; h<sizeFinalRandom; ++h )
		nbOfPathList[h] += (*itsBuckets)[i+h];
	
	// itsRandGenVector[i]->reset( totaldim, nbOfPathList, dim );
	itsRandGenVector[i]->reset( nbOfPathList, samplerNDBase->nbFactors() );
	// 6) initialize the path scheme;
	itsPathScheme->Init(timeSteps,GetSampler());

    /// 7) finally, reinit the model = get the initial step, set the time index to 0 and induct to the first time
	/// first initialize states to initStates
	itsBucketIndex = -1;
	return ARM_MCMethod::ReInit(model);
}

////////////////////////////////////////////////////
///	Class  : ARM_MCMethod
///	Routine: BuildProcessStates
///	Returns:
///	Action : Build the process states which will be 
/// set in the num method states
////////////////////////////////////////////////////
void ARM_MCMethod::BuildProcessStates( const ARM_PricingModel& model )
{
	size_t factorsNb=model.FactorCount();
	size_t nbTimeSteps=GetTimeSteps()->size();

	size_t bucketSize	= (*itsBuckets)[itsBucketIndex];

	// MC needs absolutely a N dimensional sampler
	ARM_SamplerNDBase* samplerNDBase = GetSampler()->ToSamplerNDBase();
	// this dimension is the number of factors really simulated by the Monte Carlo
	size_t dim = samplerNDBase->dim();

//	ARM_GP_MatrixPtr gaussian( new ARM_GP_Matrix(dim,bucketSize) );

//	ARM_GP_MatrixPtr zStates(new ARM_GP_Matrix(dim, bucketSize));
//	ARM_GP_MatrixPtr integXStates(new ARM_GP_Matrix(factorsNb,bucketSize,0.0));
//	ARM_GP_MatrixPtr incZStates(new ARM_GP_Matrix(dim,bucketSize,0.0));
//	ARM_GP_MatrixPtr incXStates(new ARM_GP_Matrix(factorsNb,bucketSize,0.0));

	ARM_GP_MatrixPtr zStates(new ARM_GP_Matrix);
	ARM_GP_MatrixPtr integXStates(new ARM_GP_Matrix(factorsNb,bucketSize,0.0));
	ARM_GP_MatrixPtr incZStates(new ARM_GP_Matrix);
	ARM_GP_MatrixPtr incXStates(new ARM_GP_Matrix);

	ARM_MatrixPtrVector pathSchemeProcessStates;
	itsProcessStates.resize(nbTimeSteps);

//	ARM_GP_Vector localStdDevs(dim);

	ARM_BoolVector integMCFlags = model.NeedMCIntegProcess();

	size_t timeIdx, stateIdx, factorIdx;

	// By default the numerical method
	if (integMCFlags.size() == 0)
	{
		integMCFlags.resize(factorsNb);

		for (factorIdx = 0; factorIdx < factorsNb; ++factorIdx)
			integMCFlags[factorIdx]=false;
	}
	else
	{
		if (integMCFlags.size() != factorsNb)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				ARM_USERNAME + ": integ MC flags should be equal to the number factors!");
	}

	// The path scheme compute the process increment
	itsPathScheme->ComputeProcessStates(
		bucketSize,
		GetRandGen(itsBucketIndex),
		pathSchemeProcessStates);


	ARM_VectorPtrVector drifts = GetImpSampler()->ComputeDrifts(
		dim,
		(*GetTimeSteps()),
		GetSampler());

	GetImpSampler()->ComputeExps(
		bucketSize,
		dim,
		*GetTimeSteps(),
		itsProcessStates,
		GetSampler());


	for (timeIdx = 0; timeIdx < nbTimeSteps-1; ++timeIdx)
	{
		int locFactorNb = samplerNDBase->GetLocalVar(timeIdx).size();

		zStates->resize(dim, bucketSize);
		incZStates->resize(locFactorNb,bucketSize);
		incXStates->resize(locFactorNb,bucketSize);

		itsProcessStates[timeIdx] = ARM_GP_MatrixPtr(new ARM_GP_Matrix(locFactorNb, bucketSize));
		samplerNDBase->ApplyRelDriftToIntXStates(timeIdx,integXStates);

		if (drifts.size() != 0)
		{	
			for (stateIdx = 0; stateIdx < bucketSize; stateIdx++)
			{
				// We multiply the gaussian with the local volatility of each
				// factor and add the importance sampling drift
				for (factorIdx = 0; factorIdx < dim; ++factorIdx)
					(*incZStates)(factorIdx,stateIdx) = (*pathSchemeProcessStates[timeIdx])(factorIdx,stateIdx)+(*drifts[timeIdx])(factorIdx);
			}
		}
		else
		{
			for (stateIdx = 0; stateIdx < bucketSize; stateIdx++)
			{
				// We multiply the gaussian with the local volatility of each
				// factor
				for (factorIdx = 0; factorIdx < locFactorNb; ++factorIdx)
					(*incZStates)(factorIdx,stateIdx) = (*pathSchemeProcessStates[timeIdx])(factorIdx,stateIdx);
			}
		}
		// Then we apply the rotation
		samplerNDBase->ComputeZtoXStates(timeIdx+1,incZStates,incXStates);

		samplerNDBase->ApplyAbsDriftToXStates(timeIdx,incXStates);

		for (stateIdx = 0; stateIdx < bucketSize; stateIdx++)
		{
			for (factorIdx = 0; factorIdx < locFactorNb; ++factorIdx)
			{
				// Based on the flag we stored the integrated or the increment states.
                // They are used by the model to compute model states at timeIdx+1
                // (be careful : by convention in case of integrated states
                // modelStates at timeIdx+1 = processStates(=numStates) at timeIdx = states valid at timeIdx+1)
				if(!integMCFlags[factorIdx])
					(*itsProcessStates[timeIdx])(factorIdx,stateIdx) = (*incXStates)(factorIdx,stateIdx);
				else
				{
					// These states are valide at timeIdx+1
					(*integXStates)(factorIdx,stateIdx) = (*integXStates)(factorIdx,stateIdx)+(*incXStates)(factorIdx,stateIdx);
					(*itsProcessStates[timeIdx])(factorIdx,stateIdx) = (*integXStates)(factorIdx,stateIdx);
				}
			}
		}
	}
}

/// ReInit is the function to be used before doing a for loop at each beginning
/// of the loop!
ARM_PricingStatesPtr ARM_MCMethod::ReInit( const ARM_PricingModel& model)
{
	/// increment by one the current bucket nb
	++itsBucketIndex;

#if defined(__GP_STRICT_VALIDATION)
	if( itsBucketIndex >= itsBuckets->size() )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			ARM_USERNAME + ": superbucket index is out of range!");
#endif

    /// get the initial step, set the time index to 0 and induct to the first time
	/// first initialize states to initStates
	ARM_PricingStatesPtr states( model.FirstPricingStates( (*itsBuckets)[itsBucketIndex] ) );
	states->SetOtherPayoffsFlag(itsOtherPayoffsFlag);

#if defined(__GP_STRICT_VALIDATION)
	if(GetTimeSteps()->size() < 2 )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"less than 2 time Steps! Makes no sense");
#endif

	int indexPos = 0;
    SetLastTimeIdx(indexPos);

	/// resets the numeraire
	ARM_NumerairePtr numeraire = model.GetNumeraire();

#if defined(__GP_STRICT_VALIDATION)
	if( numeraire == ARM_NumerairePtr(NULL) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			": Numeraire should not be NULL");
#endif
	numeraire->Reset(ARM_NumMethod::GP_FWDLOOKING);
	// Initialize the numeraire if if the first induct time is today
	if (itsFirstInductTime == 0)
	{
		numeraire->MoveNumeraireFwd();
		numeraire->Update(model,states,itsFirstInductTime);
	}
	
	BuildProcessStates(model);

	return Induct(model,states,itsFirstInductTime );
}

////////////////////////////////////////////////////
///	Class  : ARM_MCMethod
///	Routine: ReInitLoop
///	Returns: ARM_PricingStatesPtr
///	Action : ReInitLoop
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_MCMethod::ReInitLoop(const ARM_PricingModel& model)
{ 
	return *(new ARM_PricingStatesPtr( NULL )); 
}

////////////////////////////////////////////////////
///	Class  : ARM_MCMethod
///	Routine: ComputeTimeSteps
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_MCMethod::ComputeAndSetTimeSteps( const ARM_PricingModel& model) const
{
	if( !model.SupportAnalyticMarginal())
	{
		ARM_GP_VectorPtr pNewTimeSteps = GetSampler()->GetScheduler()->ComputeTimeSteps(model);
		const_cast< ARM_MCMethod* >(this)->SetTimeSteps(*pNewTimeSteps);
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_MCMethod
///	Routine: Induct
///	Returns: 
///	Action : induct from one time to another!
////////////////////////////////////////////////////

ARM_PricingStatesPtr ARM_MCMethod::Induct(
    const ARM_PricingModel& model,
	ARM_PricingStatesPtr& states,
	double toTime)
{
    double lastTimeStep	= GetLastTimeStep();
	int lastTimeIdx		= GetLastTimeIdx();
	int maxTimeIdx		= GetTimeSteps()->size();

    if( lastTimeStep > toTime + K_NEW_DOUBLE_TOL )
    {
		CC_Ostringstream os;
		os << ARM_USERNAME << " on " << ARM_COMPUTERNAME << ": Inconsistency in scheduling index";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
    }

	/// nothing to do if we are making no progress
    if( lastTimeStep >= toTime - K_NEW_DOUBLE_TOL )
        return states; 

    /// Find ending time index
    int endTimeIdx=lastTimeIdx;
	while(	endTimeIdx<maxTimeIdx 
		&&	GetTimeStep(endTimeIdx) <= toTime - K_NEW_DOUBLE_TOL )
		++endTimeIdx;

    if(	endTimeIdx==maxTimeIdx ||	GetTimeStep(endTimeIdx) > toTime + K_NEW_DOUBLE_TOL)
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
        "Can't find the date in the forward propagation schedule");

    /// Loop over to the to time!
	ARM_NumerairePtr numeraire = model.GetNumeraire();
	
	/// do we need to keep track of the discounting term as this is the cash numeraire?
    if( model.GetNumeraire()->GetType() == ARM_Numeraire::Cash ||
        model.GetNumeraire()->GetType() == ARM_Numeraire::RollingCash )
	{
		/// induction loop
		for(size_t timeIdx=lastTimeIdx; timeIdx<endTimeIdx; ++timeIdx)
		{
			numeraire->MoveNumeraireFwd();
			numeraire->Update(model,states,timeIdx);
			states->SetNumMethodStates(itsProcessStates[timeIdx]);
			model.MCModelStatesFromToNextTime(states,timeIdx);
		}
	}
	else
	{
		/// induction loop
		numeraire->MoveNumeraireFwd();
		numeraire->Update(model,states,lastTimeIdx);
		for(size_t timeIdx=lastTimeIdx; timeIdx<endTimeIdx; ++timeIdx)
		{
			states->SetNumMethodStates(itsProcessStates[timeIdx]);
			model.MCModelStatesFromToNextTime(states,timeIdx);
		}
	}

	numeraire->FinalizeInduction(toTime,states,model);

    /// Set ending index for next call
    SetLastTimeIdx(endTimeIdx);
    return states;
}



////////////////////////////////////////////////////
///	Class  : ARM_MCMethod
///	Routine: ComputeExercise
///	Returns: an exercise boundary
///	Action : compute an exercise boundary of an
/// exercise node.
////////////////////////////////////////////////////
ARM_ExerciseBoundary * ARM_MCMethod::ComputeExercise(const ARM_GP_VectorPtr& payoff, const ARM_GP_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector )
{
	return NULL;
}

////////////////////////////////////////////////////
///	Class  : ARM_MCMethod
///	Routine: GetBuckets,GetBucketIndex
///	Returns: 
///	Action : returns the buckets and the index
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MCMethod::GetBuckets() const
{
	return itsBuckets; 
}


size_t ARM_MCMethod::GetBucketIndex() const
{
#if defined(__GP_STRICT_VALIDATION)
	if( itsBucketIndex >= itsBuckets->size() )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			ARM_USERNAME + ": superbucket index is out of range!");
#endif
	return itsBucketIndex;
}


////////////////////////////////////////////////////
///	Class  : ARM_MCMethod
///	Routine: CreatePricerInfo
///	Returns: ARM_PricerInfo*
///	Action : creates the corresponding pricer info
////////////////////////////////////////////////////
ARM_PricerInfo* ARM_MCMethod::CreatePricerInfo( const ARM_PricingModel& model ) const
{
	//To be Changed
	return new ARM_PInfo_MultipleLoop(model,itsRandGenVector[itsRandGenVector.size()-1]->StdDevComputationFunc());
}

////////////////////////////////////////////////////
///	Class  : ARM_MCMethod
///	Routine: GetNbOfPaths
///	Returns: size_t 
///	Action : Return nb of paths according to the buckets
////////////////////////////////////////////////////
size_t ARM_MCMethod::GetNbOfPaths() const
{
    size_t nbBuckets=itsBuckets->size(),nbPaths=0;
    for(size_t i=0;i<nbBuckets;++i)
        nbPaths+=(*itsBuckets)[i];
	return nbPaths;
}


////////////////////////////////////////////////////
///	Class  : ARM_MCMethod
///	Routine: GetNumMethodStateGlobalVars
///	Returns: size_t 
///	Action : returns for each time step the global VCV matrix from
///          asOfDate to current time step
////////////////////////////////////////////////////
const ARM_MatrixVector& ARM_MCMethod::GetNumMethodStateGlobalVars() const
{
	return GetSampler()->GetGlobalVCV();
}


////////////////////////////////////////////////////
///	Class  : ARM_MCMethod
///	Routine: GetSpotProbabilities
///	Returns: ARM_GP_MatrixPtr
///	Action : Return probas to reach a future state from spot date
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_MCMethod::GetSpotProbabilities(const ARM_GP_Vector& eventTimes) const
{
    size_t i,j,nbPaths=GetNbOfPaths();

    double cstProba = 1.0/nbPaths;

    ARM_GP_Matrix* probas = new ARM_GP_Matrix(nbPaths,eventTimes.size());
    for(i=0;i<nbPaths;++i)
        for(j=0;j<eventTimes.size();++j)
        (*probas)(i,j)=cstProba;

    return ARM_GP_MatrixPtr(probas);
}


////////////////////////////////////////////////////
///	Class  : ARM_MCMethod
///	Routine: GetArrowDebreuPrices 
///	Returns: ARM_GP_VectorPtr
///	Action : returns the Arrow Debreu prices at the slice timeIdx
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_MCMethod::GetArrowDebreuPrices(size_t timeIdx, const ARM_PricingModel& model ) const
{
	ARM_GP_VectorPtr ArrowDebreuPrices;
	ARM_Numeraire::NumeraireType type = model.GetNumeraire()->GetType();
	
	if( ARM_Numeraire::NumeraireType::Cash == type )
	{
		ARM_NumeraireCash* numeraireCash = static_cast<ARM_NumeraireCash*>(&*model.GetNumeraire());

		double nbPaths = GetNbOfPaths();
		double time = (*GetTimeSteps())(timeIdx);
		ARM_GP_VectorPtr discountingTerm = numeraireCash->GetDiscount(time);

		if( GetLastTimeIdx() == timeIdx && discountingTerm != ARM_GP_VectorPtr(NULL) )
		{
			ArrowDebreuPrices	= ARM_GP_VectorPtr( (ARM_GP_Vector*) discountingTerm->Clone() );
			*ArrowDebreuPrices /= nbPaths;
		}
		else if( (*GetTimeSteps())(timeIdx) < K_DOUBLE_TOL )
		{
			ArrowDebreuPrices = ARM_GP_VectorPtr( new ARM_GP_Vector( nbPaths, 1.0/nbPaths ) );
		}
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, "should only asked at the lastTime idx the Arrow Debreu prices!" );
	}
	else
	{
		size_t i,nbPaths=GetNbOfPaths();
		ArrowDebreuPrices = ARM_GP_VectorPtr( new ARM_GP_Vector(nbPaths) );
		ARM_ZeroCurvePtr ZcCurve=model.GetZeroCurve();
		double zcT = ZcCurve->DiscountPrice((*GetTimeSteps())(timeIdx)/K_YEAR_LEN);

		for(i=0;i<nbPaths;++i)
			(*ArrowDebreuPrices)[i] = zcT/nbPaths;
	}

	return ArrowDebreuPrices;
}

////////////////////////////////////////////////////
///	Class  : ARM_MCMethod
///	Routine: ProcessPaidPayoffs 
///	Returns: 
///	Action : This functions applies the importance 
/// sampling exponential
////////////////////////////////////////////////////

void ARM_MCMethod::ProcessPaidPayoffs(ARM_VectorPtr& payoffs, double evalTime) const
{
	itsImpSampler->ProcessPaidPayoffs(payoffs,evalTime);
}

////////////////////////////////////////////////////
///	Class  : ARM_MCMethod
///	Routine: ProcessUnPaidPayoffs
///	Returns: 
///	Action : This function removes the importance 
/// sampling exponential
////////////////////////////////////////////////////

void ARM_MCMethod::ProcessUnPaidPayoffs(ARM_VectorPtr& payoffs, double evalTime) const
{
	itsImpSampler->ProcessUnPaidPayoffs(payoffs,evalTime);
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

