/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pdemethod.cpp
 *  \brief
 *
 *	\author  A Schauly
 *	\version 1.0
 *	\date July 2005
 */


/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpnummethods/pdemethod.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/vectormanip.h"
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
///	Class  : ARM_PDEMethod
///	Routine: Default constructor
///	Returns:
///	Action : Constructor, default is to use superBucket!
////////////////////////////////////////////////////
ARM_PDEMethod::ARM_PDEMethod( ARM_PDENumericalScheme * numericalScheme, size_t timeItersNb, size_t spaceItersNb, double stdDevNb )
:
	ARM_NumMethod(),
	itsOtherPayoffsFlag(true), 
	itsTimeItersNb(timeItersNb),
	itsSpaceItersNb(spaceItersNb),
	itsStdDevNb(stdDevNb),
	itsPDENumericalScheme(numericalScheme),
	itsPrevLastTimeIdx  (-1),
	itsPrevToTime		(-1)
{
	itsSpaceItersNb = ((itsSpaceItersNb%2) ? itsSpaceItersNb : itsSpaceItersNb+1);

	if( itsPDENumericalScheme )
	{
		itsPDENumericalScheme->setSpaceDiscretizationPointsNb(itsSpaceItersNb);
		itsPDENumericalScheme->setStdDevNb(stdDevNb);
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_PDEMethod
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_PDEMethod::ARM_PDEMethod(const ARM_PDEMethod& rhs)
:	ARM_NumMethod(rhs), 
	itsOtherPayoffsFlag(rhs.itsOtherPayoffsFlag),
	itsTimeItersNb(rhs.itsTimeItersNb),
	itsSpaceItersNb(rhs.itsSpaceItersNb),
	itsMethodName(rhs.itsMethodName),
	itsStdDevNb(rhs.itsStdDevNb),
	itsPDENumericalScheme(NULL),
	itsPrevLastTimeIdx  (rhs.itsPrevLastTimeIdx),
	itsPrevToTime		(rhs.itsPrevToTime)
{
	if( rhs.itsPDENumericalScheme )
		itsPDENumericalScheme = static_cast<ARM_PDENumericalScheme *> ( rhs.itsPDENumericalScheme->Clone() );
}


////////////////////////////////////////////////////
///	Class  : ARM_PDEMethod
///	Routine: Copy no clean up 
///	Returns: 
///	Action : standard copy of the object
////////////////////////////////////////////////////
void ARM_PDEMethod::CopyNoCleanUp(const ARM_PDEMethod& rhs)
{
	itsOtherPayoffsFlag			= rhs.itsOtherPayoffsFlag;
	itsTimeItersNb				= rhs.itsTimeItersNb;
	itsSpaceItersNb				= rhs.itsSpaceItersNb;
	itsPDENumericalScheme		= rhs.itsPDENumericalScheme;
	itsMethodName				= rhs.itsMethodName;
	itsStdDevNb					= rhs.itsStdDevNb;
	itsPrevLastTimeIdx			= rhs.itsPrevLastTimeIdx;
	itsPrevToTime				= rhs.itsPrevToTime;
}


////////////////////////////////////////////////////
///	Class  : ARM_PDEMethod
///	Routine: CleanUp
///	Returns: 
///	Action : deleted allocated objects!
////////////////////////////////////////////////////
void ARM_PDEMethod::CleanUp()
{
	delete itsPDENumericalScheme;
	itsPDENumericalScheme = NULL;
}


////////////////////////////////////////////////////
///	Class  : ARM_PDEMethod
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_PDEMethod::~ARM_PDEMethod()
{
	CleanUp();
}



////////////////////////////////////////////////////
///	Class  : ARM_PDEMethod
///	Routine: operator =
///	Returns: itself
///	Action : Affectation
////////////////////////////////////////////////////
ARM_PDEMethod& ARM_PDEMethod::operator=(const ARM_PDEMethod& rhs)
{
	if(this != &rhs)
	{
		ARM_NumMethod::operator=(rhs);
		CleanUp();
		CopyNoCleanUp(rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_PDEMethod
///	Routine: Clone,View, toString
///	Returns:
///	Action : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_PDEMethod::Clone() const
{
	return new ARM_PDEMethod(*this);
}


string ARM_PDEMethod::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << " =======> PDE METHOD <====== " << CC_NS(std,endl);
	os << indent << " Pricing direction      : " << GP_PricingDirectionTxt[(int)GetPricingDirection()] << CC_NS(std,endl);
	os << CC_NS(std,endl);

	if( itsPDENumericalScheme)
		os << itsPDENumericalScheme->toString();
	else
		os << "No Numerical Scheme Set" <<CC_NS(std,endl);

	/*if( GetTimeSteps() )
		os << indent << "Using Time steps:" << CC_NS(std,endl) << indent << GetTimeSteps()->toString() << CC_NS(std,endl);*/

	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_PDEMethod
///	Routine: Init, ReInit
///	Returns: 
///	Action : Initialiation of the tree (model is non const because of the postInit
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_PDEMethod::Init( ARM_PricingModel& model, double firstInductTime )
{
	itsPrevLastTimeIdx			= -1;
	itsPrevToTime				= -1;

	itsPDENumericalScheme->setSpaceDiscretizationPointsNb(itsSpaceItersNb);
	itsPDENumericalScheme->setStdDevNb(itsStdDevNb);

	// The first induct time is going from the generic security
	itsFirstInductTime = firstInductTime;

	/// 1) compute time steps: can it jump from one time to the next one?
	ComputeAndSetTimeSteps(model);

	/// gives back the hand to the model for post init operation
	model.PostInit();

	SetLastTimeIdx(GetTimeSteps()->size()-1);

	return ARM_PDEMethod::ReInit(model);
}



/// ReInit is the function to be used before doing a for loop at each beginning
/// of the loop!
ARM_PricingStatesPtr ARM_PDEMethod::ReInit( const ARM_PricingModel& model)
{
	/// increment by one the current bucket nb
	if(!itsPDENumericalScheme)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"No Numerical Scheme Set!");

	ARM_PricingStatesPtr states( itsPDENumericalScheme->Init( model ) );
	states->SetOtherPayoffsFlag(itsOtherPayoffsFlag);

#if defined(__GP_STRICT_VALIDATION)
	if(GetTimeSteps()->size() < 2 )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"less than 2 time Steps! Makes no sense");
#endif

	/// resets the numeraire
	ARM_NumerairePtr numeraire = model.GetNumeraire();

#if defined(__GP_STRICT_VALIDATION)
	if( numeraire == ARM_NumerairePtr(NULL) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			": Numeraire should not be NULL");
#endif
    numeraire->Reset(ARM_NumMethod::GP_BCKWDLOOKING);
    numeraire->MoveNumeraireBckwd();
	numeraire->Update(model,states,GetLastTimeIdx());
	model.TreeStatesToModelStates(states, GetLastTimeIdx());

	return states;
}

////////////////////////////////////////////////////
///	Class  : ARM_PDEMethod
///	Routine: ReInitLoop
///	Returns: ARM_PricingStatesPtr
///	Action : ReInitLoop
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_PDEMethod::ReInitLoop(const ARM_PricingModel& model)
{ 
	ARM_THROW( ERR_INVALID_ARGUMENT, "ReInitLoop not implemented for PDEMethod" );
}

////////////////////////////////////////////////////
///	Class  : ARM_PDEMethod
///	Routine: ComputeTimeSteps
///	Returns: void
///	Action : Adds voaltilitiesTimeSteps to pricingTimeSteps
////////////////////////////////////////////////////
void ARM_PDEMethod::ComputeAndSetTimeSteps(const ARM_PricingModel& model) const
{
	/// Adds timesteps from the model changes in volatility
	ARM_GP_Vector* timeSteps = const_cast< ARM_PDEMethod* >(this)->GetTimeSteps();
	double lastTimeStep = (*timeSteps)[timeSteps->size()-1];
	
	/// Get model times
	//ARM_GP_VectorPtr VolTimeSteps = model.VolatilitiesAndCorrelationTimesSteps();
	ARM_GP_VectorPtr VolTimeSteps = ARM_GP_VectorPtr(new ARM_GP_Vector());

	/// Merge
	std::vector<double>* pNewTimeSteps = MergeSortedVectorNoDuplicates( &timeSteps->GetValues(), &(VolTimeSteps->GetValues()) );
	timeSteps = NULL;
		
	/// Remove all time steps beyond last event date
	int resizing  = -1;
	for (size_t i(0); i<pNewTimeSteps->size(); i++)
	{
		if ( (*pNewTimeSteps)[i] > lastTimeStep + 1.e-7 )
		{
			resizing = i;
			break;
		}
	}

	if (resizing != -1)
		pNewTimeSteps->resize(resizing);

	/// Set resulting time steps
	const_cast< ARM_PDEMethod* >(this)->SetTimeSteps(ARM_GP_Vector(*pNewTimeSteps));
	delete pNewTimeSteps;

	/// Adds TimeSteps from what users asked to. 
	std::vector<double> AdditionalTimeSteps( itsTimeItersNb );
	ARM_GP_Vector* currentTimeSteps = const_cast< ARM_PDEMethod* >(this)->GetTimeSteps();
	double OtherDatesStep = (*currentTimeSteps)[currentTimeSteps ->size() -1]/(itsTimeItersNb+2);

	if (AdditionalTimeSteps.size())
	{
		AdditionalTimeSteps[0] = OtherDatesStep;
		for(size_t i=1 ; i< itsTimeItersNb ; ++i)
			AdditionalTimeSteps[i] = AdditionalTimeSteps[i-1]+OtherDatesStep;
	}

	std::vector<double> * newTimeSteps = MergeSortedVectorNoDuplicates( &AdditionalTimeSteps, &currentTimeSteps->GetValues() );
	currentTimeSteps = NULL;

	const_cast< ARM_PDEMethod* >(this)->SetTimeSteps( ARM_GP_Vector(*newTimeSteps ));
	delete newTimeSteps;
}

////////////////////////////////////////////////////
///	Class  : ARM_PDEMethod
///	Routine: Induct
///	Returns: 
///	Action : induct from one time to another!
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_PDEMethod::Induct( const ARM_PricingModel& model, ARM_PricingStatesPtr& states,  double toTime)
{
#ifdef __GP_STRICT_VALIDATION
    if(!GetTimeSteps())
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "No schedule initialised" );
#endif

	/// This is to handle case where we call Induct
	/// several times on the same date (Hunt Kennedy Bootstrap calibration)
	if (fabs(toTime)>K_NEW_DOUBLE_TOL && fabs(itsPrevToTime - toTime)<K_NEW_DOUBLE_TOL)
	{
		SetLastTimeIdx(itsPrevLastTimeIdx);
	}

	/// Update previous toTime and previous lastTimeIdx
	itsPrevToTime      = toTime;
	itsPrevLastTimeIdx = GetLastTimeIdx();

    double lastTimeStep = GetLastTimeStep();
	if (lastTimeStep < 0)
		lastTimeStep = 0;

#ifdef __GP_STRICT_VALIDATION
    if(lastTimeStep < toTime - K_NEW_DOUBLE_TOL)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "Inconsistency in backpropagation schedule index");
#endif

	int lastTimeIdx = GetLastTimeIdx();


	if( lastTimeStep <= toTime + K_NEW_DOUBLE_TOL )
		return states; /// nothing to do!

    /// Find ending time index
	int endTimeIdx;
    for( endTimeIdx=lastTimeIdx; endTimeIdx>0; --endTimeIdx )
        if( GetTimeStep(endTimeIdx) <= toTime + K_NEW_DOUBLE_TOL )
            break;

	int timeIdx=lastTimeIdx-1;
    for(;timeIdx>=endTimeIdx;--timeIdx)
		itsPDENumericalScheme->Induct( model, states, timeIdx);

	if (model.NeedStatesEval(endTimeIdx))
	{
		model.TreeStatesToModelStates(states, endTimeIdx);
	}

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

	SetLastTimeIdx(endTimeIdx);
	
    return states;
}



////////////////////////////////////////////////////
///	Class  : ARM_PDEMethod
///	Routine: ComputeExercise
///	Returns: an exercise boundary
///	Action : compute an exercise boundary of an
/// exercise node.
////////////////////////////////////////////////////
ARM_ExerciseBoundary * ARM_PDEMethod::ComputeExercise(const ARM_GP_VectorPtr& payoff, const ARM_GP_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector )
{
	return NULL;
}


////////////////////////////////////////////////////
///	Class  : ARM_PDEMethod
///	Routine: CreatePricerInfo
///	Returns: ARM_PricerInfo*
///	Action : creates the corresponding pricer info
////////////////////////////////////////////////////
ARM_PricerInfo* ARM_PDEMethod::CreatePricerInfo( const ARM_PricingModel& model ) const
{
	return new ARM_PInfo_SingleLoopPDE( itsPDENumericalScheme->getPriceIndex() );
}

////////////////////////////////////////////////////
///	Class  : ARM_PDEMethod
///	Routine: GetSpotProbabilities
///	Returns: ARM_GP_MatrixPtr
///	Action : Return probas to reach a future state from spot date
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_PDEMethod::GetSpotProbabilities(const std::vector<double>& eventTimes) const
{
// FIXMEFRED: mig.vc8 (30/05/2007 18:14:02):cast
	return static_cast<ARM_GP_MatrixPtr>(new ARM_GP_Matrix(0,0));
}


////////////////////////////////////////////////////
///	Class  : ARM_PDEMethod
///	Routine: GetArrowDebreuPrices 
///	Returns: ARM_GP_VectorPtr
///	Action : returns the Arrow Debreu prices at the slice timeIdx
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_PDEMethod::GetArrowDebreuPrices(size_t timeIdx, const ARM_PricingModel& model ) const
{
	return static_cast<ARM_GP_VectorPtr>(NULL);
}


////////////////////////////////////////////////////
///	Class  : ARM_PDEMethod
///	Routine: ControlVariateOnInstrument 
///	Returns: void
///	Action : does a control variate on an instrument
////////////////////////////////////////////////////
void ARM_PDEMethod::ControlVariateOnInstrument( ARM_VectorPtr& numericInstrument, double correctValue,
	const string& curveName, double evalTime, const ARM_PricingStatesPtr& states, const ARM_PricingModel& model ) const
{
//	ControlVariateOnInstrumentAdditive( numericInstrument, correctValue, curveName, evalTime, states, model );
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

