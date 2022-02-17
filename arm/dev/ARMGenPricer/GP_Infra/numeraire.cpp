/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file numeraire.cpp
 *
 *  \brief
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date October 2003
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/comparisonfunctor.h"
#include "gpbase/ostringstream.h"
#include "gpbase/singleton.h"
#include "gpbase/eventviewer.h"

#include "gpinfra/numeraire.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/pricingmodel.h"


#include <algorithm>

#include <iomanip> /// for setprecision()

using std::make_pair;
using std::find;

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_Numeraire
///	Const  : NumeraireTypeTable
/// Meaning: Text corresponding to the type enum
////////////////////////////////////////////////////
const string ARM_Numeraire::NumeraireTypeTable[] =
{
	"Cash",
	"TerminalZc",
	"Terminal Event Zc",
	"Unknown"
};


////////////////////////////////////////////////////
///	Class  : ARM_Numeraire
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_Numeraire::CopyNoCleanUp( const ARM_Numeraire& rhs )
{
    itsTimes	= rhs.itsTimes  ? (ARM_GP_Vector*) rhs.itsTimes->Clone()  : NULL;
	itsValues0	= rhs.itsValues0;
    itsMatIdx	= rhs.itsMatIdx;
}


////////////////////////////////////////////////////
///	Class  : ARM_Numeraire
///	Routine: CleanUp
///	Returns: 
///	Action : Arguments destruction
////////////////////////////////////////////////////
void ARM_Numeraire::CleanUp()
{
    delete itsTimes;
	itsTimes=NULL;
}


////////////////////////////////////////////////////
///	Class  : ARM_Numeraire
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_Numeraire::ARM_Numeraire()
:	itsTimes( new ARM_GP_Vector(1,1.0) ),
	itsValues0(0.0),
	itsMatIdx(0)
{
	CC_ARM_SETNAME(ARM_NUMERAIRE);
}


////////////////////////////////////////////////////
///	Class  : ARM_Numeraire
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_Numeraire::ARM_Numeraire(const ARM_Numeraire& rhs)
:	ARM_RootObject(rhs)
{
    CopyNoCleanUp(rhs);
}


////////////////////////////////////////////////////
///	Class  : ARM_Numeraire
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_Numeraire::~ARM_Numeraire()
{
    CleanUp();
}


////////////////////////////////////////////////////
///	Class  : ARM_Numeraire
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_Numeraire& ARM_Numeraire::operator=(const ARM_Numeraire& rhs)
{
	if(this != &rhs)
	{
		ARM_RootObject::operator=(rhs);
		CleanUp();
		CopyNoCleanUp(rhs);
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_Numeraire
///	Routine: Reset
///	Returns: void
///	Action : Reset the numeraire for each bucket
////////////////////////////////////////////////////
void ARM_Numeraire::Reset(int pricingDir)
{
	ResetNumDiscountMap();
}


////////////////////////////////////////////////////
///	Class  : ARM_Numeraire
///	Routine: ResetLoop
///	Returns: void
///	Action : Reset the numeraire for each loop
////////////////////////////////////////////////////
void ARM_Numeraire::ResetLoop()
{
}

////////////////////////////////////////////////////
///	Class  : ARM_Numeraire
///	Routine: ResetNumDiscountMap
///	Returns: void
///	Action : Reset the numeraire discount map
////////////////////////////////////////////////////
void ARM_Numeraire::ResetNumDiscountMap()
{
	itsNumDiscountMap.clear();
}


////////////////////////////////////////////////////
///	Class   : ARM_Numeraire
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_Numeraire::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;

    os << "\n\n";
    os << indent << "ARM_Numeraire\n";
    os << indent << "-------------\n";
	os << indent << "Type : ";
	os << NumeraireTypeName() << "\n";
	const ARM_NumeraireCash* numCash;
	if((numCash=dynamic_cast<const ARM_NumeraireCash*>(this)))
	{
		os << indent << "Renormalisation : " << (numCash->GetRenormalisation() ? "On\n" : "Off\n");
	}
	os << indent << "Time Lag";
    if(itsTimes)
        os << "\n";
    else
	    os << " : not initialised !\n";

    os << CC_NS(std,fixed);

	if(itsTimes)
	{
		for(int i=0;i<itsTimes->size();++i)
        {
			os << indent << CC_NS(std,setw)(9) << CC_NS(std,setprecision)(3) << (*itsTimes)[i] << "\n";
        }
	}
	os << "Value0 = " << CC_NS(std,setw)(10) << CC_NS(std,setprecision)(5) << itsValues0 << "\n\n";

    return os.str();
}

////////////////////////////////////////////////////
///	Class  : ARM_Numeraire
///	Routine: MoveNumeraire
///	Returns: 
///	Action : Move the numeraire time, does nothing for std numeraire
////////////////////////////////////////////////////

void ARM_Numeraire::MoveNumeraireFwd(){}
void ARM_Numeraire::MoveNumeraireBckwd(){}


////////////////////////////////////////////////////
///	Class   : ARM_Numeraire
///	Routines: DefaultProcessPaidPayoffs
///	Returns : void
///	Action  : Default implementation of paid payoffs
///           "numerisation"
////////////////////////////////////////////////////

void ARM_Numeraire::DefaultProcessPaidPayoffs( const string& payModelName, 
	ARM_VectorPtr& payoffs, 
	double evalTime, 
	const ARM_PricingStatesPtr& states, 
	const ARM_PricingModel& model ) const
{
    /// payModelName is no more used because GetRefModel()
    /// of every model give the payment model

	if( evalTime <  -K_NEW_DOUBLE_TOL )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		ARM_USERNAME + ": the evalTime should never be negative!" );
	
	double finalTime = GetMaturity();
	double zcNum = itsValues0;

	// Use the cache to store numeraire discount factor vector
	ARM_VectorPtr measureChangeInverse;
	NumDiscountMap::iterator it = itsNumDiscountMap.find(evalTime);
	if (it != itsNumDiscountMap.end())
	{
		measureChangeInverse = it->second;
	}
	else
	{
		measureChangeInverse = model.GetRefModel()->GetDiscountFunctor()->DiscountFactor("",evalTime, finalTime, states);
		itsNumDiscountMap.insert(make_pair(evalTime,measureChangeInverse));
	}
	
	
#if defined( __GP_STRICT_VALIDATION )
	if( measureChangeInverse->size() != payoffs->size() )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"numeraire value and payoff have not the same size!");
#endif
	
	for(size_t i=0; i<payoffs->size();++i)
		(*payoffs)[i] *= zcNum/(*measureChangeInverse)[i];
}


////////////////////////////////////////////////////
///	Class   : ARM_Numeraire
///	Routines: ProcessUnPaidPayoffs
///	Returns : void
///	Action  : Default implementation of paid payoffs
///           "numerisation"
////////////////////////////////////////////////////

void ARM_Numeraire::DefaultProcessUnPaidPayoffs( const string& payModelName, 
	ARM_VectorPtr& payoffs, 
	double evalTime, 
	const ARM_PricingStatesPtr& states, 
	const ARM_PricingModel& model ) const
{
    /// payModelName is no more used because GetRefModel()
    /// of every model give the payment model

#if defined( __GP_STRICT_VALIDATION )
	if( evalTime <  -K_NEW_DOUBLE_TOL )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		ARM_USERNAME + ": the evalTime should never be negative!" );
#endif
	
	double finalTime	= GetMaturity();
	double zcNum = GetValue0();

	// Use the cache to store numeraire discount factor vector
	ARM_VectorPtr measureChangeInverse;
	NumDiscountMap::iterator it = itsNumDiscountMap.find(evalTime);
	if (it != itsNumDiscountMap.end())
	{
		measureChangeInverse = it->second;
	}
	else
	{
		measureChangeInverse = model.GetRefModel()->GetDiscountFunctor()->DiscountFactor("",evalTime, finalTime, states);
		itsNumDiscountMap.insert(make_pair(evalTime,measureChangeInverse));
	}
	
#if defined( __GP_STRICT_VALIDATION )
	if( measureChangeInverse->size() != payoffs->size() )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"numeraire value and payoff have not the same size!");
#endif
	
	for(size_t i=0; i<states->size();++i)
		(*payoffs)[i] *= (*measureChangeInverse)[i]/zcNum;
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	ARM_NumeraireFwd
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	ARM_NumeraireTerminalZc
////////////////////////////////////////////////////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_NumeraireTerminalZc
///	Routine: Init
///	Returns: 
///	Action : Schedule initialisation is done to
///          a default schedule if necessary and
///          associated spot values are computed
////////////////////////////////////////////////////
void ARM_NumeraireTerminalZc::Init(const ARM_ZeroCurveFunctor& discFunctor,const string& modelName, const ARM_TimeInfoPtrVector& timeInfos, const ARM_VectorPtr& numeraireTimes )
{
	CleanUp();
	if (numeraireTimes != ARM_VectorPtr(NULL))
	{
#if defined( __GP_STRICT_VALIDATION )
		if( 1 != numeraireTimes->size() )
			ARM_THROW( ERR_INVALID_ARGUMENT, "terminalZc numeraire times should contain only 1 element" );
#endif
		itsTimes = new ARM_GP_Vector( *numeraireTimes );
		itsMatIdx=0;
	}
	else
	{
		size_t nbEvents = timeInfos.size();
		double lastPayTime = 0;

		for(int i=0; i<nbEvents;i++)
		{
			const ARM_GP_Vector& lastPayTimes = timeInfos[i]->GetPayTimes();
			if(lastPayTime<lastPayTimes[lastPayTimes.size()-1])
				lastPayTime =lastPayTimes[lastPayTimes.size()-1] ;
		}
		itsTimes = new ARM_GP_Vector(1,lastPayTime);
	}

	ARM_PricingStatesPtr dumStates(NULL);
	double df = (*(discFunctor.DiscountFactor(modelName,0.0,(*itsTimes)[0],dumStates)))[0];
	itsValues0 = df;
	itsMatIdx=0;

}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	ARM_NumeraireTerminalEventZc
////////////////////////////////////////////////////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_NumeraireTerminalEventZc
///	Routine: Init
///	Returns: 
///	Action : Schedule initialisation is done to
///          a default schedule if necessary and
///          associated spot values are computed
////////////////////////////////////////////////////
void ARM_NumeraireTerminalEventZc::Init(const ARM_ZeroCurveFunctor& discFunctor,const string& modelName,const ARM_TimeInfoPtrVector& timeInfos, const ARM_VectorPtr& numeraireTimes )
{
	CleanUp();
	if (numeraireTimes != ARM_VectorPtr(NULL))
	{
#if defined( __GP_STRICT_VALIDATION )
		if( 1 != numeraireTimes->size() )
			ARM_THROW( ERR_INVALID_ARGUMENT, "TerminalEventZc numeraire times should contain only 1 element" );
#endif
		itsTimes = new ARM_GP_Vector( *numeraireTimes );
	}
	else
	{
		size_t nbEvents = timeInfos.size();
		double lastEventTime = timeInfos[nbEvents-1]->GetEventTime();
		itsTimes = new ARM_GP_Vector(1,lastEventTime);
	}

	ARM_PricingStatesPtr dumStates(NULL);
	double df = (*(discFunctor.DiscountFactor(modelName,0.0,(*itsTimes)[0],dumStates)))[0];
	itsValues0 = df;
	itsMatIdx=0;
}




////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	ARM_NumeraireCash
////////////////////////////////////////////////////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_NumeraireCash
///	Routine: Init
///	Returns: 
///	Action : multiplies by the appropriate final numeraire!
////////////////////////////////////////////////////

void ARM_NumeraireCash::Init(const ARM_ZeroCurveFunctor& discFunctor,const string& modelName,const ARM_TimeInfoPtrVector& timeInfos, const ARM_VectorPtr& numeraireTimes )
{
	CleanUp();
	itsTimes = new ARM_GP_Vector(1,0.0);
	itsMatIdx=0;
	itsValues0 = 1.0;
	itsDiscountMap.erase(itsDiscountMap.begin(),itsDiscountMap.end());
}


////////////////////////////////////////////////////
///	Class  : ARM_NumeraireCash
///	Routine: Reset
///	Returns: void
///	Action : Reset the numeraire for each bucket
////////////////////////////////////////////////////

void ARM_NumeraireCash::Reset(int pricingDir)
{ 
	ARM_Numeraire::Reset(pricingDir);
	if(		ARM_NumMethod::GP_FWDLOOKING == pricingDir 
		||  ARM_NumMethod::GP_FWDBCKWDLOOKING == pricingDir )
	{
		itsDiscountMap.erase(itsDiscountMap.begin(),itsDiscountMap.end());
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_NumeraireCash
///	Routines: FinalizeInduction
///	Returns : void
///	Action  : Calibrate the discountingTerm
////////////////////////////////////////////////////

void ARM_NumeraireCash::FinalizeInduction(
	double evalTime,
	const ARM_PricingStatesPtr& states,
	const ARM_PricingModel& model )
{
	/// renormalisation
	if( itsRenormalisationIsOn )
	{
		ARM_ZeroCurvePtr zc = model.GetZeroCurve();
		double df0 = zc->DiscountPrice(evalTime/K_YEAR_LEN);
		ARM_GP_VectorPtr discount = GetDiscount(evalTime);

		///  renormalisation to be rigorously exact on df
		double sum = 0;
		for(size_t i=0;i<discount->size(); ++i)
			sum += (*discount)[i];
		sum /= discount->size();
		(*discount) *= df0/sum;
	}
}



////////////////////////////////////////////////////
///	Class   : ARM_NumeraireCash
///	Routines: ProcessPaidPayoffs
///	Returns : void
///	Action  : Default implementation of paid payoffs
///           "numerisation"
////////////////////////////////////////////////////

void ARM_NumeraireCash::ProcessPaidPayoffs( const string& payModelName, 
	ARM_VectorPtr& payoffs, 
	double evalTime,
	const ARM_PricingStatesPtr& states, 
	const ARM_PricingModel& model ) const
{
	if (!IsDiscountEmpty())
	{
		ARM_GP_VectorPtr discount = GetDiscount(evalTime);

		for(size_t i=0; i<payoffs->size();++i)
			(*payoffs)[i] *= (*discount)[i];
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_NumeraireCash
///	Routines: ProcessUnPaidPayoffs
///	Returns : void
///	Action  : Default implementation of paid payoffs
///           "numerisation"
////////////////////////////////////////////////////

void ARM_NumeraireCash::ProcessUnPaidPayoffs( const string& payModelName, 
	ARM_VectorPtr& payoffs, 
	double evalTime, 
	const ARM_PricingStatesPtr& states, 
	const ARM_PricingModel& model ) const
{
	if (!IsDiscountEmpty())
	{
		ARM_GP_VectorPtr discount = GetDiscount(evalTime);

		for(size_t i=0; i<payoffs->size();++i)
			(*payoffs)[i] /= (*discount)[i];
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_NumeraireCash
///	Routines: Update
///	Returns : void
///	Action  : Update the numeraire rolligng discount
////////////////////////////////////////////////////
void ARM_NumeraireCash::Update(
		const ARM_PricingModel& model, 
		const ARM_PricingStatesPtr& states,
		size_t timeIdx )
{
	if( model.GetNumMethod() != ARM_NumMethodPtr(NULL) )
	{
		ARM_NumMethod::GP_PricingDirection pricingDirection = model.GetNumMethod()->GetPricingDirection();
		ARM_NumMethod::GP_PricingDirection loopPricingDirection = model.GetNumMethod()->GetPricingDirCurrLoop();
		if( pricingDirection== ARM_NumMethod::GP_BCKWDLOOKING)
		{
			// Do Nothing
		}
		else if(pricingDirection== ARM_NumMethod::GP_FWDLOOKING || pricingDirection==ARM_NumMethod::GP_FWDBCKWDLOOKING)
		{
			if (loopPricingDirection == ARM_NumMethod::GP_FWDLOOKING)
			{
				/// For the first previous
				double time		= model.GetNumMethod()->GetTimeStep(timeIdx);
				double nextTime	= model.GetNumMethod()->GetTimeStep(timeIdx+1);
				if( itsDiscountMap.size() == 0 )
				{
					itsDiscountMap.insert(make_pair(nextTime,model.LocalDiscounts(timeIdx,nextTime-time,states)));
				}
				else
				{
					ARM_GP_VectorPtr discountingTerm = CreateClonedPtr<ARM_GP_Vector>(&*GetDiscount(time));
					(*discountingTerm) *= (*model.LocalDiscounts(timeIdx,nextTime-time,states));
					itsDiscountMap.insert(make_pair(nextTime,discountingTerm));
				}
			}
		}
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, "unknonw pricing direction" );
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_NumeraireCash
///	Routines: GetDiscount
///	Returns : void
///	Action  : Return the discount for a time
////////////////////////////////////////////////////

ARM_GP_VectorPtr ARM_NumeraireCash::GetDiscount(double time) const
{
	DiscountMap::const_iterator it = itsDiscountMap.find(time);

	if (it != itsDiscountMap.end())
		return it->second;
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT,"ARM_NumeraireCash : Discount is missing for the time ");
		return ARM_GP_VectorPtr(NULL);
	}
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	ARM_NumeraireRolling
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_NumeraireRolling
///	Routine: MoveNumeraire
///	Returns: 
///	Action : Move the numeraire time
////////////////////////////////////////////////////

void ARM_NumeraireRolling::MoveNumeraireFwd()
{
	++itsMatIdx;
#if defined( __GP_STRICT_VALIDATION )
	if( itsMatIdx >= itsTimes->size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "itsMatIdx >= itsTimes->size()" );
#endif
}

void ARM_NumeraireRolling::MoveNumeraireBckwd()
{
	--itsMatIdx;
#if defined( __GP_STRICT_VALIDATION )
	if( itsMatIdx <-1 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "itsMatIdx <0" );
#endif
}


////////////////////////////////////////////////////
///	Class  : ARM_NumeraireRolling
///	Routine: Reset
///	Returns: void
///	Action : Reset the numeraire for each bucket
////////////////////////////////////////////////////

void ARM_NumeraireRolling::Reset(int pricingDir)
{ 
	ARM_Numeraire::Reset(pricingDir);
	if(		ARM_NumMethod::GP_FWDLOOKING == pricingDir 
		||  ARM_NumMethod::GP_FWDBCKWDLOOKING == pricingDir )
	{
		itsMatIdx = -1;
		itsRollingDiscount.resize(itsTimes->size()); 
	}
	else if( ARM_NumMethod::GP_BCKWDLOOKING == pricingDir )
		itsMatIdx = itsTimes->size();
}

////////////////////////////////////////////////////
///	Class  : ARM_NumeraireRolling
///	Routine: ResetLoop
///	Returns: void
///	Action : Reset the numeraire for each loop
////////////////////////////////////////////////////
void ARM_NumeraireRolling::ResetLoop()
{
	itsMatIdx = itsTimes->size();
}

////////////////////////////////////////////////////
///	Class  : ARM_NumeraireRolling
///	Routine: Init
///	Returns: 
///	Action : Schedule initialisation is done to
///          a default schedule if necessary and
///          associated spot values are computed
////////////////////////////////////////////////////
void ARM_NumeraireRolling::ProcessPaidPayoffs( const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
		const ARM_PricingStatesPtr& states, const ARM_PricingModel& model ) const
{
	const ARM_NumMethodPtr numMethod = model.GetNumMethod();

	if (numMethod != ARM_NumMethodPtr(NULL) )
	{
		for(size_t statesIdx=0; statesIdx<payoffs->size();++statesIdx )
		{
			double df = (*itsRollingDiscount[itsMatIdx])[statesIdx];
			(*payoffs)[statesIdx]*=df;
		}
	}

	DefaultProcessPaidPayoffs(
		payModelName, 
		payoffs,
		evalTime, 
		states,
		model);
}


////////////////////////////////////////////////////
///	Class  : ARM_NumeraireRolling
///	Routine: Init
///	Returns: 
///	Action : Schedule initialisation is done to
///          a default schedule if necessary and
///          associated spot values are computed
////////////////////////////////////////////////////

void ARM_NumeraireRolling::ProcessUnPaidPayoffs( const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
	const ARM_PricingStatesPtr& states, const ARM_PricingModel& model ) const
{
	const ARM_NumMethodPtr numMethod = model.GetNumMethod();

	if (numMethod != ARM_NumMethodPtr(NULL) )
	{
		for(size_t statesIdx=0; statesIdx<payoffs->size();++statesIdx )
		{
			double df = (*itsRollingDiscount[itsMatIdx])[statesIdx];
			(*payoffs)[statesIdx]/=df;
		}
	}

	DefaultProcessUnPaidPayoffs(
		payModelName, 
		payoffs,
		evalTime, 
		states,
		model);

}

////////////////////////////////////////////////////
///	Class   : ARM_NumeraireRolling
///	Routines: Update
///	Returns : void
///	Action  : Update the numeraire rolligng discount
////////////////////////////////////////////////////
void ARM_NumeraireRolling::Update(
		const ARM_PricingModel& model, 
		const ARM_PricingStatesPtr& states,
		size_t timeIdx )
{
	double evalTime = model.GetNumMethod()->GetTimeStep(timeIdx);
	// For the first previous
	double startTime = 0.0;
	if (itsMatIdx > 0)
		startTime = (*itsTimes)[itsMatIdx-1];
	double endTime = (*itsTimes)[itsMatIdx];

	ARM_NumMethod::GP_PricingDirection pricingDirection = model.GetNumMethod()->GetPricingDirection();
	ARM_NumMethod::GP_PricingDirection loopPricingDirection = model.GetNumMethod()->GetPricingDirCurrLoop();
	if( pricingDirection== ARM_NumMethod::GP_BCKWDLOOKING)
	{
		// Do Nothing
	}
	else if(pricingDirection== ARM_NumMethod::GP_FWDLOOKING || pricingDirection==ARM_NumMethod::GP_FWDBCKWDLOOKING)
	{
		if (loopPricingDirection == ARM_NumMethod::GP_FWDLOOKING)
		{
			if( !itsMatIdx )
				itsRollingDiscount[itsMatIdx] = model.GetRefModel()->GetDiscountFunctor()->DiscountFactor( "",startTime, endTime, states);
			else
			{
				ARM_GP_VectorPtr rollingDiscount = model.GetRefModel()->GetDiscountFunctor()->DiscountFactor("",evalTime, endTime, states);
				(*rollingDiscount) /= (*model.GetRefModel()->GetDiscountFunctor()->DiscountFactor("",evalTime, startTime, states));
				itsRollingDiscount[itsMatIdx] = ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(itsRollingDiscount[itsMatIdx-1]->Clone()));
				(*itsRollingDiscount[itsMatIdx]) *= (*rollingDiscount);
			}
		}
	}
	else ARM_THROW( ERR_INVALID_ARGUMENT, "unknown pricing direction" );
}



////////////////////////////////////////////////////
///	Class  : ARM_NumeraireRollingPayment
///	Routine: Init
///	Returns: 
///	Action : Schedule initialisation is done to
///          a default schedule if necessary and
///          associated spot values are computed
////////////////////////////////////////////////////
void ARM_NumeraireRollingPayment::Init(const ARM_ZeroCurveFunctor& discFunctor,const string& modelName,const ARM_TimeInfoPtrVector& timeInfos, const ARM_VectorPtr& numeraireTimes )
{
	CleanUp();
	if (numeraireTimes != ARM_VectorPtr(NULL))
	{
		itsTimes = new ARM_GP_Vector( *numeraireTimes );
		itsMatIdx=0;
	}
	else
	{
		size_t nbEvents = timeInfos.size();
		itsTimes = new ARM_GP_Vector(nbEvents);
		for (size_t i = 0; i < nbEvents; ++i)
		{
			const ARM_GP_Vector& lastPayTimes = timeInfos[i]->GetPayTimes();
			double lastPayTime = lastPayTimes[lastPayTimes.size()-1];
			(*itsTimes)[i] = lastPayTime;
		}
	}

	itsValues0 = 1.0;
	itsValues0	= 1.0;
	itsMatIdx	= 0;

	itsRollingDiscount.resize(0);
}

////////////////////////////////////////////////////
///	Class  : ARM_NumeraireRollingEvent
///	Routine: Init
///	Returns: 
///	Action : Schedule initialisation is done to
///          a default schedule if necessary and
///          associated spot values are computed
////////////////////////////////////////////////////
void ARM_NumeraireRollingEvent::Init(const ARM_ZeroCurveFunctor& discFunctor,const string& modelName,const ARM_TimeInfoPtrVector& timeInfos, const ARM_VectorPtr& numeraireTimes )
{
	CleanUp();
	if (numeraireTimes != ARM_VectorPtr(NULL))
	{
		itsTimes = new ARM_GP_Vector( *numeraireTimes );
		itsMatIdx=0;
	}
	else
	{
		size_t nbEvents = timeInfos.size();
		itsTimes = new ARM_GP_Vector(nbEvents);
		for (size_t i = 0; i < nbEvents; ++i)
		{
			double lastEventTime = timeInfos[i]->GetEventTime();
			(*itsTimes)[i] = lastEventTime;
		}
	}

	itsValues0 = 1.0;
	itsMatIdx	= 0;
	itsValues0	= 1.0;
	
	itsRollingDiscount.resize(0);
}


////////////////////////////////////////////////////
///	Class   : ARM_NumeraireCash
///	Routines: Update
///	Returns : void
///	Action  : Update the numeraire rolligng discount
////////////////////////////////////////////////////
void ARM_NumeraireRollingCash::Update(
		const ARM_PricingModel& model, 
		const ARM_PricingStatesPtr& states,
		size_t timeIdx )
{
	if( model.GetNumMethod() != ARM_NumMethodPtr(NULL) )
	{
		ARM_NumMethod::GP_PricingDirection pricingDirection = model.GetNumMethod()->GetPricingDirection();
		ARM_NumMethod::GP_PricingDirection loopPricingDirection = model.GetNumMethod()->GetPricingDirCurrLoop();
		if( pricingDirection== ARM_NumMethod::GP_BCKWDLOOKING)
		{
			// Do Nothing
		}
		else if(pricingDirection== ARM_NumMethod::GP_FWDLOOKING || pricingDirection==ARM_NumMethod::GP_FWDBCKWDLOOKING)
		{
			if (loopPricingDirection == ARM_NumMethod::GP_FWDLOOKING)
			{
				/// For the first previous
				double time		= model.GetNumMethod()->GetTimeStep(timeIdx);
				double nextTime	= model.GetNumMethod()->GetTimeStep(timeIdx+1);
				if( itsDiscountMap.size() == 0 )
				{
					itsDiscountMap.insert(make_pair(nextTime,model.GetRefModel()->GetDiscountFunctor()->DiscountFactor("",time, nextTime, states)));
				}
				else
				{
					ARM_GP_VectorPtr discountingTerm = CreateClonedPtr<ARM_GP_Vector>(&*GetDiscount(time));
					(*discountingTerm) *= (*model.GetRefModel()->GetDiscountFunctor()->DiscountFactor("",time, nextTime, states));
					itsDiscountMap.insert(make_pair(nextTime,discountingTerm));
				}
			}
		}
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, "unknonw pricing direction" );
	}
}

CC_END_NAMESPACE()



/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
