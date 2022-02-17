/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file treebase.cpp
 *
 *  \brief
 *
 *	\author  JM Prie E Benhamou
 *	\version 1.0
 *	\date November 2004
 */

/// to remove identified warning
#include "gpbase/removeidentifiedwarning.h"

#include "gpnummethods/treebase.h"
#include "gpbase/env.h"
#include "gpbase/timer.h"

/// gpnummethods
#include "gpnummethods/sampler.h"
#include "gpnummethods/slice.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/reconnector.h"
#include "gpnummethods/smoother.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricerinfo.h"
#include "gpinfra/exerciseboundary.h"
#include "gpinfra/pricingstates.h"


static double totalBackwardInductTime=0.0;

CC_BEGIN_NAMESPACE( ARM )


const double ARM_TreeBase::LastSliceCalibrationTermInDays = 90.0; // 3M

////////////////////////////////////////////////////
///	Class  : ARM_TreeBase
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_TreeBase::CopyNoCleanUp(const ARM_TreeBase& rhs)
{
    itsTruncator            = rhs.itsTruncator ? (ARM_TruncatorBase*) rhs.itsTruncator->Clone() : NULL;
    itsReconnector          = rhs.itsReconnector ? (ARM_ReconnectorBase*) rhs.itsReconnector->Clone() : NULL;
    DuplicateCloneablePointorAndNullVectorInPlace<ARM_SliceBase>( *rhs.itsSlices, *itsSlices );
    itsComputeSpotProbas    = rhs.itsComputeSpotProbas;
    itsSmoother             = rhs.itsSmoother ? (ARM_SmootherBase*) rhs.itsSmoother->Clone() : NULL;
    itsStates               = rhs.itsStates != ARM_PricingStatesPtr(NULL) ? ARM_PricingStatesPtr(new ARM_PricingStates(*(rhs.itsStates))) : ARM_PricingStatesPtr(NULL);
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeBase
///	Routine: CleanUp
///	Returns: 
///	Action : Arguments destruction
////////////////////////////////////////////////////
void ARM_TreeBase::CleanUp()
{
    delete itsTruncator;
    delete itsReconnector;
    delete itsSmoother;
    DeletePointorVector<ARM_SliceBase>( *itsSlices );
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeBase
///	Routine: Default constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_TreeBase::ARM_TreeBase( const ARM_SamplerBase* sampler, const ARM_TruncatorBase* truncator, const ARM_ReconnectorBase* reconnector, const ARM_SmootherBase* smoother, bool computeSpotProbas)
:	ARM_NumMethod(sampler),
    itsTruncator( truncator ? (ARM_TruncatorBase*) truncator->Clone() : NULL ),
    itsReconnector( reconnector ? (ARM_ReconnectorBase*) reconnector->Clone() : NULL ),
    itsSlices( ARM_SliceVectorPtr(new ARM_SliceVector(0)) ),
    itsComputeSpotProbas(computeSpotProbas),
    itsSmoother(smoother ? (ARM_SmootherBase*) smoother->Clone() : NULL),
    itsStates( ARM_PricingStatesPtr(new ARM_PricingStates) )
{}


////////////////////////////////////////////////////
///	Class  : ARM_TreeBase
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_TreeBase::ARM_TreeBase(const ARM_TreeBase& rhs)
:   ARM_NumMethod(rhs),
    itsTruncator(NULL),
    itsReconnector(NULL),
    itsSmoother(NULL),
    itsSlices( ARM_SliceVectorPtr(new ARM_SliceVector(0)) )
{
    CopyNoCleanUp(rhs);
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeBase
///	Routine: Assignement operator
///	Returns: ARM_TreeBase&
///	Action : 
////////////////////////////////////////////////////
ARM_TreeBase& ARM_TreeBase::operator= (const ARM_TreeBase& rhs )
{
    if( this != &rhs )
    {
        ARM_NumMethod::operator =(rhs);
        CleanUp();
        CopyNoCleanUp(rhs);
    }
    return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeBase
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_TreeBase::~ARM_TreeBase()
{
    CleanUp();
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeBase
///	Routine: SetSlices
///	Returns: void
///	Action : set a slice vector ptr
////////////////////////////////////////////////////
void ARM_TreeBase::SetSlices( const ARM_SliceVectorPtr sliceVector )
{
    DeletePointorVector<ARM_SliceBase>( *itsSlices );
    itsSlices = sliceVector;
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeBase
///	Routine: ComputeTimeSteps
///	Returns: ARM_GP_VectorPtr
///	Action : compute a schedule
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_TreeBase::ComputeTimeSteps(const ARM_PricingModel& model) const
{
	return GetSampler()->GetScheduler()->ComputeTimeSteps(model);
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeBase
///	Routine: Init
///	Returns: 
///	Action : Initialiation of the tree
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_TreeBase::Init( ARM_PricingModel& model, double firstInductTime, const ARM_VectorPtrVector& prevDriftCorrections, size_t driftOffset, const ARM_GP_Vector& targetLocalPayoffs, ARM_GP_VectorPtr& timeSteps)
{
    /// Check dimension consistency
    size_t nbDims = dim();
    if(model.FactorCount() != nbDims)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
            " : number of factors of the model differs from the tree dimension");

    totalBackwardInductTime=0.0;
	ARM_Timer timer;
	timer.ClockStartTime();

    /// Init the sampler to get the time steps and slices
    ARM_TimeStepsAndSlices* result = GetSampler()->Init(model,timeSteps);
 	timer.ClockEndTime();
    double samplerInitTime = timer.GetDuration();

    SetTimeSteps(*(result->itsTimeSteps));

	timer.ClockStartTime();
    SetSlices(result->itsSlices);
 	timer.ClockEndTime();

	delete result;

    double setSliceTime = timer.GetDuration();


	/// gives back the hand to the model for post init operation
	model.PostInit();

    /// Last index setting for backward induction purpose
    size_t nbSteps = GetTimeSteps()->size();
    if(nbSteps < 2 || GetTimeStep(0) > 0.0)
	    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
            "Tree schedule must contain at least one time step and start at spot date !" );

    int lastTimeIdx = nbSteps-1;
	timer.ClockStartTime();

    /// Truncation initialisation
    GetTruncator()->Init(GetSlices(),GetSampler());

	timer.ClockEndTime();
    double truncatorInitTime = timer.GetDuration();

    /// Build slices and nodes
    ARM_SliceVectorPtr slices=GetSlices();
    ARM_TransitionStates transStates(GetTransitionSize());
    size_t totalNbStates=0;
/***
FILE* f=fopen("c:\\temp\\dumpTree2G.txt","a");
fprintf(f,"\nsetSlice=%10.5lf\tsampler=%10.5lf\ttruncator=%10.5lf\n",setSliceTime,samplerInitTime,truncatorInitTime);
***/
	timer.ClockStartTime();

    double nextTime,targetPayoff;
    ARM_PricingStatesPtr nullStates(NULL);
    ARM_GP_VectorPtr nullDriftCorrection(NULL);
    ARM_GP_VectorPtr df;
    for(size_t timeIdx=0; timeIdx<lastTimeIdx; ++timeIdx)
    {
        /// Set current slice index, it could be useful...
        SetLastTimeIdx(timeIdx);

        if(timeIdx < targetLocalPayoffs.size())
            targetPayoff = targetLocalPayoffs[timeIdx];
        else
        {
            /// Default target = df from asOf to next slice time
            nextTime = GetTimeStep(timeIdx+1);
            df=model.DiscountFactor( model.GetRefModel()->GetModelName() ,0.0,nextTime,nullStates);
            targetPayoff = (*df)[0];
        }
        (*slices)[timeIdx]->LinkedWithNextSlice( (*slices)[timeIdx+1], GetSampler(), GetTruncator(), itsComputeSpotProbas, transStates, (timeIdx < prevDriftCorrections.size() ? prevDriftCorrections[timeIdx] : nullDriftCorrection), driftOffset, targetPayoff);
        totalNbStates += (*slices)[timeIdx]->size();
    }

    /// Calibrate last slice of the tree if necessary using previous calibrated drifts if any
    SetLastTimeIdx(lastTimeIdx);
    if(model.NeedArrowDebreuPrices() && model.NeedLocalDiscount())
    {
        if(prevDriftCorrections.size() >= nbSteps && prevDriftCorrections[lastTimeIdx]->size() == nbDims)
            (*slices)[lastTimeIdx]->SetDriftCorrectionVect(prevDriftCorrections[lastTimeIdx]);

        nextTime = GetTimeStep(lastTimeIdx) + LastSliceCalibrationTermInDays;
        if(targetLocalPayoffs.size() >= nbSteps)
            targetPayoff = targetLocalPayoffs[lastTimeIdx];
        else
        {
            /// Default target = df from asOf to slice time + LastSliceCalibrationTermInDays
            df=model.DiscountFactor( model.GetRefModel()->GetModelName() ,0.0,nextTime,nullStates);
            targetPayoff = (*df)[0];
        }
        bool isLocalDfComputed = false;
        ARM_GP_VectorPtr localDf = (*slices)[lastTimeIdx]->ComputeDriftCorrection(&model,GetSampler(),nextTime,driftOffset,targetPayoff,isLocalDfComputed);

/****
FILE* f=fopen("c:\\temp\\dumpTree2G.txt","a");
ARM_VectorPtr driftCor;
fprintf(f,"#%1dD drift corrections\n",nbDims);
for(size_t timeIdx=0; timeIdx<=lastTimeIdx; ++timeIdx)
{
    driftCor = (*slices)[timeIdx]->GetDriftCorrectionVect();
    fprintf(f,"#%3d\tt=%6.2lf\t",timeIdx,GetTimeStep(timeIdx));
    for( size_t i=0; i<nbDims; ++i )
        fprintf(f,"%15.10lf\t",(*driftCor)[i]);
    fprintf(f,"\n");
}
fclose(f);
****/

    }

	timer.ClockEndTime();
    double totalForwardInductTime = timer.GetDuration();
/***
fprintf(f,"totalNbStates=%6d\ttotalForwardInduct=%10.5lf\n",totalNbStates,totalForwardInductTime);
fclose(f);
***/


    /// Create last slice tree states...
    ARM_PricingStatesPtr initStates( model.FirstPricingStates( (*slices)[lastTimeIdx]->size() ) );

    ///... set a copy to the rolling slice states...
    itsStates = ARM_PricingStatesPtr( new ARM_PricingStates(*initStates) );

    ///... then compute tree states and convert them to model states
	ARM_NumerairePtr numeraire = model.GetNumeraire();
    numeraire->Reset(ARM_NumMethod::GP_BCKWDLOOKING);
    numeraire->MoveNumeraireBckwd();
	numeraire->Update(model,initStates,lastTimeIdx);

    ComputeStateVariables(lastTimeIdx, initStates);
    model.TreeStatesToModelStates(initStates, lastTimeIdx);

    return initStates;
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeBase
///	Routine: ReInitLoop
///	Returns: ARM_PricingStatesPtr
///	Action : Reinit for loop change (in case of pricing 
/// direction change)
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_TreeBase::ReInitLoop(const ARM_PricingModel& model)
{ 
	return *(new ARM_PricingStatesPtr( NULL )); 
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeBase
///	Routine: ReInit
///	Returns: 
///	Action : Reinitialisation for loop pricing
///			 Not supported for tree like method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_TreeBase::ReInit( const ARM_PricingModel& model)
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "Tree method does not support multiple loops" );
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeBase
///	Routine: GetBuckets,GetBucketIndex
///	Returns: 
///	Action : Returns the buckets and the index
////////////////////////////////////////////////////
ARM_VectorPtr ARM_TreeBase::GetBuckets() const
{
	/// one bucket with size 1
	return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,1)); 
}

size_t ARM_TreeBase::GetBucketIndex() const
{
	/// first bucket only
	return 0;
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeBase
///	Routine: CreatePricerInfo
///	Returns: ARM_PricerInfo*
///	Action : Creates the corresponding pricer info
////////////////////////////////////////////////////
ARM_PricerInfo* ARM_TreeBase::CreatePricerInfo( const ARM_PricingModel& model ) const
{
	return new ARM_PInfo_SingleLoop;
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeBase
///	Routine: GetSpotProbabilities
///	Returns: ARM_GP_MatrixPtr
///	Action : Select probas to reach a future states
///			 at given dates from spot date
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_TreeBase::GetSpotProbabilities(const ARM_GP_Vector& eventTimes) const
{
    size_t i,j,nbStates=0,eventIdx=0;
    const ARM_GP_Vector& sched = *(GetTimeSteps());

    /// No probabilities available
    if(! &sched || !itsComputeSpotProbas || itsSlices->size()==0)
        return ARM_GP_MatrixPtr(NULL);

    if(eventTimes.size()==1 && eventTimes[0]==0.0)
        return ARM_GP_MatrixPtr(new ARM_GP_Matrix(1,1,1.0));

    /// Time selection in schedule
    vector<int> timeIdx(eventTimes.size(),-1);
    size_t sliceSize;
    for(i=0;i<sched.size() && eventIdx < eventTimes.size();++i)
    {
        if(eventTimes[eventIdx] <= sched[i])
        {
            if(eventTimes[eventIdx] == sched[i] && (sliceSize=(*itsSlices)[i]->GetSpotProbas()->size()) > 0)
            {
                /// Probabilities computed at this time step
                timeIdx[eventIdx]=i;
                nbStates = CC_Max(nbStates,sliceSize);
            }

            ++eventIdx;
        }
    }

    /// Save spot probabilities
    /// 1st line = event time
    ARM_GP_Matrix* probas = new ARM_GP_Matrix(nbStates+1,timeIdx.size());
    for(j=0;j<timeIdx.size();++j)
    {
        if(timeIdx[j]!=-1)
        {
            (*probas)(0,j)=sched[timeIdx[j]];
            for(i=0;i<(*itsSlices)[timeIdx[j]]->GetSpotProbas()->size();++i)
                (*probas)(i+1,j)=(*((*itsSlices)[timeIdx[j]]->GetSpotProbas()))[i];
        }
        else
        {
            /// Slice not found or probabilities not available
            (*probas)(0,j)=eventTimes[j];
            i=0;
        }

        for(;i<nbStates;++i)
            (*probas)(i+1,j)=0.0;
    }

    return ARM_GP_MatrixPtr(probas);
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeBase
///	Routine: GetSpotProbabilities
///	Returns: ARM_GP_VectorPtr
///	Action : Return probas to reach a future states
///			 at a given time step
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_TreeBase::GetSpotProbabilities(size_t timeIdx) const
{
    /// Test if probabilities are available
    const ARM_GP_Vector& sched = *(GetTimeSteps());
    if(! &sched || !itsComputeSpotProbas || itsSlices->size()==0 || timeIdx >= sched.size())
        return ARM_GP_VectorPtr(NULL);

	return (*itsSlices)[timeIdx]->GetSpotProbas();
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeBase
///	Routine: Induct
///	Returns: List of pricing states
///	Action : backward induction from last time
///          to the next one
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_TreeBase::Induct( const ARM_PricingModel& model, ARM_PricingStatesPtr& states,  double toTime)
{
#ifdef __GP_STRICT_VALIDATION
    if(!GetTimeSteps())
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "No schedule initialised" );
#endif

    double lastTimeStep=GetLastTimeStep();

#ifdef __GP_STRICT_VALIDATION
    if(lastTimeStep < toTime - K_NEW_DOUBLE_TOL)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "Inconsistency in backpropagation schedule index");
#endif


	ARM_Timer timer;
	timer.ClockStartTime();

	int lastTimeIdx=GetLastTimeIdx();
	if( lastTimeStep <= toTime + K_NEW_DOUBLE_TOL )
		return states; /// nothing to do!

    /// Find ending time index
	int endTimeIdx;
    for( endTimeIdx=lastTimeIdx; endTimeIdx>0; --endTimeIdx )
        if( GetTimeStep(endTimeIdx) <= toTime + K_NEW_DOUBLE_TOL )
            break;

	size_t PayoffStatesSize = states->GetPayoffStatesSize();
	size_t p;

    /// Loop backward for each time steps
   
    ARM_PricingStatesPtr nextPricingStates=states;

	bool isOtherPayoffs;

	if (itsStates->GetPayoffStatesSize() != PayoffStatesSize)
		itsStates->resizePayoffStatesVector(PayoffStatesSize);

	for (p = 0; p < PayoffStatesSize; ++p)
	{
		isOtherPayoffs = states->GetPayoffStates(p).GetOtherPayoffsFlag();
		itsStates->GetPayoffStates(p).SetOtherPayoffsFlag(isOtherPayoffs);
	}

    ARM_PricingStatesPtr curPricingStates;
    size_t nbSnapshots;
    size_t nbPayoffs                = nextPricingStates->GetPayoffsSize();
	size_t nbProbaChanges			= nextPricingStates->GetProbaChangesSize();
	size_t nbIntermediatePayoffs;
    size_t nbDims                   = nextPricingStates->NumMethodStatesSize();
	size_t nbModelStates            = nextPricingStates->ModelStatesSize();
    double payoff,probaChange,intermediatePayoff,df;

    size_t nbStates,stateIdx,payoffIdx,probaChangeIdx;
    double time,nextTime=lastTimeStep;

	curPricingStates = itsStates;

    ARM_SliceBase* nextSlice=(*GetSlices())[lastTimeIdx];
    ARM_SliceBase* slice;
    ARM_TransitionStates transStates(GetTransitionSize());
    size_t totalNbStates=0;

    for(int timeIdx=lastTimeIdx-1;timeIdx>=endTimeIdx;--timeIdx)
    {
		GetAuxiliaryPayoffFunc auxPayoffFtor(nextPricingStates);
		GetPayoffFunc payoffFtor(nextPricingStates);
		GetProbaChangeFunc probaChangeFtor(nextPricingStates);
		GetIntermediatePayoffFunc interPayoffFtor(nextPricingStates);

        time=GetTimeStep(timeIdx);

        /// Transfer payoff snapshots because no retro-propagation needed
        slice = (*GetSlices())[timeIdx];
        nbStates=slice->size();

		curPricingStates = itsStates;
        if((curPricingStates->GetStatesOfPayoffsSize() != nbStates) || 
		   (curPricingStates->GetPayoffsSize() != nbPayoffs))
            curPricingStates->resizePayoffs(nbPayoffs,nbStates);

		if((curPricingStates->GetStatesOfProbaChangesSize() != nbStates) 
			|| (curPricingStates->GetProbaChangesSize() != nbProbaChanges))
            curPricingStates->resizeProbaChanges(nbProbaChanges,nbStates);

		for (p = 0; p < PayoffStatesSize; ++p)
		{
			curPricingStates->GetPayoffStates(p).resizePayoff(nbStates);
			isOtherPayoffs = itsStates->GetPayoffStates(p).GetOtherPayoffsFlag();
			if(isOtherPayoffs)
			{
				nbSnapshots = nextPricingStates->GetPayoffStates(p).GetPayoffSnapshotsSize();
				nbIntermediatePayoffs = nextPricingStates->GetPayoffStates(p).GetIntermediatePayoffsSize();

				if(curPricingStates->GetPayoffStates(p).GetIntermediatePayoffStatesSize() != nbStates ||
				   curPricingStates->GetPayoffStates(p).GetIntermediatePayoffsSize() != nbIntermediatePayoffs)
					curPricingStates->GetPayoffStates(p).resizeIntermediatePayoffs(nbIntermediatePayoffs,nbStates);

				if(curPricingStates->GetPayoffStates(p).GetPayoffSnapshotsSize() != nbSnapshots)
					curPricingStates->GetPayoffStates(p).resizePayoffSnapshots(nbSnapshots);

				for(size_t i=0;i<nbSnapshots;++i)
					curPricingStates->GetPayoffStates(p).SetPayoffSnapshot(i,nextPricingStates->GetPayoffStates(p).GetPayoffSnapshot(i));
			}
		}

        if(model.NeedLocalDiscount())
        {
            /// Model requires non trivial local discount values :
            /// get current slice states then convert them to model states...
			/// avoid double action on the last time index
			ComputeStateVariables(timeIdx, curPricingStates);
			model.TreeStatesToModelStates(curPricingStates,timeIdx);

            /// ...and compute local df
			ARM_GP_VectorPtr localDiscounts = model.LocalDiscounts(timeIdx,nextTime-time,curPricingStates);

/****
FILE* f=fopen("c:\\temp\\dumpTree2G.txt","a");
fprintf(f,"\n\nInduct timeIdx=%2d\n",timeIdx);
size_t i;
for(i=0;i<nextPricingStates->GetPayoffStatesSize();++i)
    fprintf(f,"  payoffs=%15.10lf\n",nextPricingStates->GetPayoff(i,0));
fprintf(f,"\n");
for(i=0;i<localDiscounts->size();++i)
    fprintf(f,"  localDf=%15.13lf\n",(*localDiscounts)[i]);
fclose(f);
****/

            /// Update payoffs & intermediate payoffs
            for(stateIdx=0; stateIdx<nbStates;++stateIdx )
            {
                df = (*localDiscounts)[stateIdx];
                for(payoffIdx=0;payoffIdx<nbPayoffs;++payoffIdx)
                {
					auxPayoffFtor.SetPayoffIdx(payoffIdx);
                    payoff = slice->ComputeExpectation(stateIdx,nextSlice,auxPayoffFtor,transStates);
                    curPricingStates->SetPayoff(stateIdx,payoffIdx, df * payoff);
                }
				for (probaChangeIdx=0;probaChangeIdx<nbProbaChanges;++probaChangeIdx)
				{
					probaChangeFtor.SetProbaChangeIdx(probaChangeIdx);
					probaChange = slice->ComputeExpectation(stateIdx,nextSlice,probaChangeFtor,transStates);
					curPricingStates->SetProbaChange(stateIdx,probaChangeIdx, df * probaChange);
				}
				for(p = 0; p < PayoffStatesSize; ++p)
				{
					payoffFtor.SetPayoffIdx(p);
					payoff = slice->ComputeExpectation(stateIdx,nextSlice,payoffFtor,transStates);
					curPricingStates->GetPayoffStates(p).SetPayoff(stateIdx, df * payoff);
					isOtherPayoffs = itsStates->GetPayoffStates(p).GetOtherPayoffsFlag();
					if (isOtherPayoffs)
					{
						for(payoffIdx=0;payoffIdx<nbIntermediatePayoffs;++payoffIdx)
						{
							interPayoffFtor.SetPayoffIdx(p);
							interPayoffFtor.SetInterIdx(payoffIdx);
							intermediatePayoff = slice->ComputeExpectation(stateIdx,nextSlice,interPayoffFtor,transStates);
							curPricingStates->GetPayoffStates(p).SetIntermediatePayoff(stateIdx,payoffIdx, df * intermediatePayoff);
						}
					}
				}

                /// Reset transition states from current node to next slice
                transStates.SetUsedSize(0);
            }
        }
        else
        {
            /// Update payoffs & intermediate payoffs
            for(stateIdx=0; stateIdx<nbStates;++stateIdx )
            {
                for(payoffIdx=0;payoffIdx<nbPayoffs;++payoffIdx)
                {
					auxPayoffFtor.SetPayoffIdx(payoffIdx);
                    payoff = slice->ComputeExpectation(stateIdx,nextSlice,auxPayoffFtor,transStates);
                    curPricingStates->SetPayoff(stateIdx,payoffIdx,payoff);
                }
				for (probaChangeIdx=0;probaChangeIdx<nbProbaChanges;++probaChangeIdx)
				{
					probaChangeFtor.SetProbaChangeIdx(probaChangeIdx);
					probaChange = slice->ComputeExpectation(stateIdx,nextSlice,probaChangeFtor,transStates);
					curPricingStates->SetProbaChange(stateIdx,probaChangeIdx, probaChange);
				}
				for(p = 0; p < PayoffStatesSize; ++p)
				{
					payoffFtor.SetPayoffIdx(p);
					payoff = slice->ComputeExpectation(stateIdx,nextSlice,payoffFtor,transStates);
					curPricingStates->GetPayoffStates(p).SetPayoff(stateIdx, payoff);
					isOtherPayoffs = itsStates->GetPayoffStates(p).GetOtherPayoffsFlag();
					if (isOtherPayoffs)
					{
						for(payoffIdx=0;payoffIdx<nbIntermediatePayoffs;++payoffIdx)
						{
							interPayoffFtor.SetPayoffIdx(p);
							interPayoffFtor.SetInterIdx(payoffIdx);
							intermediatePayoff = slice->ComputeExpectation(stateIdx,nextSlice,interPayoffFtor,transStates);
							curPricingStates->GetPayoffStates(p).SetIntermediatePayoff(stateIdx,payoffIdx,intermediatePayoff);
						}
					}
				}

                /// Reset transition states from current node to next slice
                transStates.SetUsedSize(0);
            }
			if (model.NeedStatesEval(timeIdx))
			{
				ComputeStateVariables(timeIdx, curPricingStates);
				model.TreeStatesToModelStates(curPricingStates, timeIdx);
			}
        }
		nbProbaChanges = curPricingStates->GetProbaChangesSize();
        itsStates=nextPricingStates;
		nextPricingStates=curPricingStates;
        nextTime=time;
        nextSlice=slice;
        totalNbStates += nbStates;
    }

	/// get numeraire
	ARM_NumerairePtr numeraire = model.GetNumeraire();
    if(    ARM_Numeraire::RollingPayment == numeraire->GetType() 
        || ARM_Numeraire::RollingEvent == numeraire->GetType() )
    {
	    if ( timeIdx>=0 )
	    {
		    numeraire->MoveNumeraireBckwd();
		    numeraire->Update(model,states,timeIdx);
	    }
    }

	if (!model.NeedModelTime())
	{
		/// Compute current slice states then convert to model states
		ComputeStateVariables(endTimeIdx, curPricingStates);
		model.TreeStatesToModelStates(curPricingStates, endTimeIdx);
	}

    /// Set ending index for next call
    SetLastTimeIdx(endTimeIdx);

	timer.ClockEndTime();
    double backwardInductTime = timer.GetDuration();
    totalBackwardInductTime += backwardInductTime;

/***
FILE* f=fopen("c:\\temp\\dumpTree2G.txt","a");
fprintf(f,"nbTimeStep=%3d\tnbStates=%6d\tbackwardInduct=%10.5lf\ttotalBackwardInduct=%10.5lf\tnbPayoffs=%3d\tnbInterPayoffs=%3d\tnbSnapshots=%3d\n",
        lastTimeIdx-endTimeIdx,totalNbStates,backwardInductTime,totalBackwardInductTime,nbPayoffs,nbIntermediatePayoffs,nbSnapshots);
fclose(f);
***/

    return curPricingStates;
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeBase
///	Routine: ExerciseSmoothing1D
///	Returns: void
///	Action : Compute smoothing values that correct a 1D
///          exercise boundary transition
////////////////////////////////////////////////////
void ARM_TreeBase::ExerciseSmoothing1D(const ARM_GP_Vector& exerFct, double coef, ARM_GP_Vector& smoothValues, ARM_GP_Vector& exerStates) const
{
    size_t stateIdx,nbStates = exerFct.size();

    /// Locate state index exerIdx such that exerFct[]
    /// changes it sign between exerIdx & exerIdx+1
    double val,lastVal=exerFct[0];
    smoothValues[0] = 0.0;
    ARM_IntVector exerIdx;
    for(stateIdx=1;stateIdx<nbStates;++stateIdx)
    {
        smoothValues[stateIdx]=0.0;
        val = exerFct[stateIdx];
        if((val > 0.0 && lastVal < 0.0) || (val < 0.0 && lastVal > 0.0))
            exerIdx.push_back(stateIdx-1);
        lastVal=val;
    }

    size_t i,nbExerArea=exerIdx.size();
    if(nbExerArea==0)
        return; /// always exercised or never exercised then no smoothing

    /// Theoritical state of exercise boundary is not used at the moment
    /// but it could be saved in exerStates vector
    double unusedTheoExerState;
    for(i=0;i<nbExerArea;++i)
        unusedTheoExerState=itsSmoother->Compute(exerFct,exerIdx[i],coef,smoothValues);
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeBase
///	Routine: ComputeExercise
///	Returns: An exercise boundary
///	Action : Compute an exercise boundary
///          of an exercise node
////////////////////////////////////////////////////
ARM_ExerciseBoundary* ARM_TreeBase::ComputeExercise(const ARM_GP_VectorPtr& payoff, const ARM_GP_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector )
{
    size_t stateIdx,nbStates = payoff->size();
    ARM_GP_VectorPtr exerFct( new ARM_GP_Vector(nbStates) );
    for(stateIdx=0;stateIdx<nbStates;++stateIdx)
        (*exerFct)[stateIdx] = (*payoff)[stateIdx] - (*contOpt)[stateIdx];

    ARM_GP_VectorPtr smoothValues( new ARM_GP_Vector(nbStates,0.0) );

    /// exerStates could collect interpolated exercise boundary
    /// but not used at the moment
    ARM_GP_VectorPtr unusedExerStates( NULL );

    ExerciseSmoothing(exerFct,smoothValues,unusedExerStates);

	return new ARM_SmoothTreeExerciseBoundary(smoothValues,unusedExerStates);
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeBase
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_TreeBase::toString(const string& indent, const string& nextIndent ) const
{
    CC_Ostringstream os;

    os << ARM_NumMethod::toString(indent,nextIndent);

    os << "\n";
    if( itsTruncator)
        os << "Truncator  :  " << itsTruncator->toString( indent, nextIndent ) << "\n";
    if( itsReconnector)
        os << "Reconnector:  " << itsReconnector->toString( indent, nextIndent ) << "\n";
    if( itsSmoother)
        os << "Smoother   :  " << itsSmoother->toString( indent, nextIndent ) << "\n";
    os << "Spot Probas Computation = " << (itsComputeSpotProbas ? "On" : "Off") << "\n";

    if( itsSlices != ARM_SliceVectorPtr(NULL) )
    {    
        os << indent << "Slices     :  "  << "\n";
        for(size_t i=0; i<itsSlices->size(); ++i )
        {
            os << indent << "#" << std::dec << std::setw(3) << i << "\t";
            os << "t=" << std::fixed << std::setprecision(2) << std::setw(8) << GetTimeStep(i) << " => " << (*itsSlices)[i]->toString(indent, nextIndent ) << "\n";
        }
    }
    return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeBase
///	Routine: GetArrowDebreuPrices
///	Returns: ARM_GP_VectorPtr
///	Action : return the AD prices vector at an input slice
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_TreeBase::GetArrowDebreuPrices(size_t timeIdx, const ARM_PricingModel& model ) const
{ 
#ifdef __GP_STRICT_VALIDATION
	if( timeIdx >= itsSlices->size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": out of bounds" );
#endif
	return (*itsSlices)[timeIdx]->GetArrowDebreuPrices();
}


////////////////////////////////////////////////////
///	Class  : ARM_NumMethod
///	Routine: GetNumMethodStateGlobalVars 
///	Returns: ARM_MatrixVector
///	Action : returns for each time step the global VCV matrix from
///          asOfDate to current time step
////////////////////////////////////////////////////
const ARM_MatrixVector& ARM_TreeBase::GetNumMethodStateGlobalVars() const
{
    return GetSampler()->GetGlobalVCV();
}


////////////////////////////////////////////////////
///	Class  : ARM_NumMethod
///	Routine: GetDriftCorrections 
///	Returns: void
///	Action : Collect drift correction for all slices
///          and put them in a vector
////////////////////////////////////////////////////
void ARM_TreeBase::GetDriftCorrections(ARM_VectorPtrVector& driftCorrections ) const
{
    size_t i,nbSlices=itsSlices->size();
    driftCorrections.resize(nbSlices);
    for(i=0;i<nbSlices;++i)
        driftCorrections[i] = (*itsSlices)[i]->GetDriftCorrectionVect();
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
