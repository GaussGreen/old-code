/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file irfwdmod.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpmodels/IRFwdMod.h"

#include "gpbase/ostringstream.h"
#include "gpbase/utilityport.h"

#include "gpinfra/typedef.h"
#include "gpinfra/lexerdec.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/functorop.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/timeinfo.h"
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/pricingmodeltype.h"


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_IrFwdMod
///	Routine: Constructor
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////
ARM_IrFwdMod::ARM_IrFwdMod(const ARM_ZeroCurvePtr& zc)
:	ARM_PricingModelIR(zc)
{
	CC_ARM_SETNAME(ARM_IRFWDMOD);
}


///////////////////////////////////////////////////
///	Routine: toString
///	Returns: size of the deal description
///	Action : 
////////////////////////////////////////////////////
string ARM_IrFwdMod::toString(const string& indent, const string& nextIndent) const
{ 
	CC_Ostringstream os;
	os << "\n\n===========> ARM_IrFwdMod <===========\n";
	if( GetNumMethod() != ARM_NumMethodPtr(NULL) )
	    os << "with corresponding numerical method:\n" << GetNumMethod()->toString();

	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_IrFwdMod
///	Routine: Clone method
///	Returns: a new copy of this
///	Action : Copy this
////////////////////////////////////////////////////
ARM_Object* ARM_IrFwdMod::Clone() const
{
	return new ARM_IrFwdMod( *this );
}


////////////////////////////////////////////////////
///	Class  : ARM_IrFwdMod
///	Routine: Init
///	Returns: void
///	Action : Initialise before any pricing... for this
///			model, does nothing
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_IrFwdMod::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
    if( GetNumMethod() == ARM_NumMethodPtr(NULL) )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
        "not numerical method set with the model IRFwd.. Please advise!");
    }

    if(GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_FWDBCKWDLOOKING)
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
        "Forward backward method not supported yet");
    }

	size_t nbEvents = timeInfos.size();
    std::vector<double> timeSteps(nbEvents);
	for( int i=0; i<nbEvents; ++i)
		timeSteps[i] = timeInfos[i]->GetEventTime();

    GetNumMethod()->SetTimeSteps( std::vector<double>(timeSteps) );
	if(GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_BCKWDLOOKING)
		GetNumMethod()->SetLastTimeIdx( nbEvents-1 );
	else if(GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_FWDLOOKING)
		GetNumMethod()->SetLastTimeIdx(0);

	/// We put a payoff with zero at first date
	ARM_PricingStatesPtr pricingStates( new ARM_PricingStates(1,1,1) );
    pricingStates->SetPayoff(0,0,0.0);
    return pricingStates;

}



////////////////////////////////////////////////////
///	Class  : ARM_IrFwdMod
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_IrFwdMod::~ARM_IrFwdMod()
{}


////////////////////////////////////////////////////
///	Class  : ARM_IrFwdMod
///	Routine: Copy Constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_IrFwdMod::ARM_IrFwdMod( const ARM_IrFwdMod& rhs)
:	ARM_PricingModelIR( rhs )
{}


////////////////////////////////////////////////////
///	Class  : ARM_IrFwdMod
///	Routine: Assignment operator
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_IrFwdMod& ARM_IrFwdMod::operator=( const ARM_IrFwdMod& rhs )
{
	if( this !=	 &rhs )
	    ARM_PricingModelIR::operator = ( rhs );
	return *this;
}




////////////////////////////////////////////////////
///	Class  : ARM_IrFwdMod
///	Routine: DiscountFactor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_IrFwdMod::DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
        const ARM_PricingStatesPtr& states) const
{

#ifdef __GP_STRICT_VALIDATION

	if( evalTime > maturityTime )
	{
		CC_Ostringstream os;
		os << "Trying to price a discountFactor in the past.\n" 
			<< ARM_USERNAME  << ": please advise\n";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}
#endif

	/// does not use modelStates as this is useless for this model
	/// does not use curveName!
	double fwdDf = GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN)
		/ GetZeroCurve()->DiscountPrice(evalTime/K_YEAR_LEN);
// FIXMEFRED: mig.vc8 (30/05/2007 16:11:16):cast
	return static_cast<ARM_VectorPtr>(new std::vector<double>(1,fwdDf));
}





////////////////////////////////////////////////////
///	Class  : ARM_IrFwdMod
///	Routine: VanillaCaplet
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_IrFwdMod::VanillaCaplet(
		const string& curveName, 
		double evalTime,
		double payTime,			/// not used for convexity correction
		double period,
        double payNotional,
		double fwdResetTime,	/// used for volatility computation
		double fwdStartTime,
        double fwdEndTime,
		double fwdPeriod,
        const std::vector<double>& strikesPerState,
        int capFloor,
        const ARM_PricingStatesPtr& states) const
{
#ifdef __GP_STRICT_VALIDATION

	if( evalTime > fwdStartTime)
	{
		CC_Ostringstream os;
		os << "Trying to price a Libor with a startTime in the past.\n" 
			<< ARM_USERNAME  << ": please advise\n";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}

	if( fwdStartTime > fwdEndTime )
	{
		CC_Ostringstream os;
		os << "Trying to price a Libor with an endTime shorter than its StartTime.\n" 
			<< ARM_USERNAME  << ": please advise\n";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}

	if( fwdStartTime > payTime )
	{
		CC_Ostringstream os;
		os << "Trying to price a Libor with an endTime shorter than its StartTime.\n" 
			<< ARM_USERNAME  << ": please advise\n";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}

	if( period <= 0 )
	{
		CC_Ostringstream os;
		os << "Trying to price a Libor with a non strictly positive period (=" << period << ").\n"
			<< ARM_USERNAME  << ": please advise\n";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}
#endif

	/// does not use modelStates as this is useless for this model
	/// does not use curveName!
    ARM_PricingStatesPtr singleState( new ARM_PricingStates(1,1,1) );

	ARM_VectorPtr dfStart	= GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdStartTime,singleState);

	ARM_VectorPtr  dfEnd	= GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdEndTime,singleState);

	ARM_VectorPtr dfPay	    = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,payTime,singleState);

	double libor	= ((*dfStart)[0]/(*dfEnd)[0]-1.0)/fwdPeriod;
	double result	= capFloor == K_CAP? 
		CC_Max( libor - strikesPerState[0], 0.0 ) * period * (*dfPay)[0] * payNotional
	:	CC_Max( strikesPerState[0] - libor, 0.0 ) * period * (*dfPay)[0] * payNotional;

// FIXMEFRED: mig.vc8 (30/05/2007 16:11:29):cast
	return static_cast<ARM_VectorPtr>(new std::vector<double>(1,result));
}


////////////////////////////////////////////////////
///	Class  : ARM_IrFwdMod
///	Routine: VanillaSwaption
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_IrFwdMod::VanillaSwaption(
		const string& curveName,
		double evalTime,
		double swapResetTime,
		const std::vector<double>& fixNotional,
		const std::vector<double>& floatNotional,
		double floatStartTime,
		double floatEndTime,		
		const std::vector<double>& floatResetTimes,
		const std::vector<double>& floatStartTimes,
		const std::vector<double>& floatEndTimes,
		const std::vector<double>& floatIntTerms,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,  
        const ARM_GP_Matrix& strikesPerState,
        int callPut,
        const ARM_PricingStatesPtr& states,
		bool isConstantNotional,
		bool isConstantSpread,
		bool isConstantStrike) const
{
	/// TO BE UPDATED
	/// Check that the notional is constant
	double swapNotional = fixNotional[0];
	if (!(isConstantNotional&&isConstantSpread&&isConstantStrike))
				ARM_THROW( ERR_INVALID_ARGUMENT, "The Model can not price a swaption with variable notional, Spread or Strike!" );


    if( !GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) )
    {
        /// We need to compute the floating leg by forward method and no more by double notional
        /// but we have not at the moment all floating leg datas => throw an error
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "Swaption pricing not implemented for differents discount & fixing curves" );
    }

#ifdef __GP_STRICT_VALIDATION

	if( evalTime > floatStartTime )
	{
		CC_Ostringstream os;
		os << "Trying to price a SwapRate with a startTime in the past:"
			<< floatStartTime << " days.\n" 
			<< ARM_USERNAME  << ": please advise\n";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}

	if( floatEndTime < floatStartTime )
	{
		CC_Ostringstream os;
		os << "Trying to price a  SwapRate with an endTime shorter than its StartTime: " 
			<< "floatStartTime : " << floatStartTime << "days > "
			<< "floatEndTime : " << floatEndTime << "days.\n"
			<< ARM_USERNAME  << ": please advise\n";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}

	/// test fix paytimes
	if( fixPayTimes[0] < floatStartTime )
	{
		CC_Ostringstream os;
		os << "Trying to price a  SwapRate with a first pay times shorter than its StartTime: " 
			<< "floatStartTime : " << floatStartTime << "days > "
			<< "fixPayTimes[0] : " << fixPayTimes[0] << "days.\n"
			<< ARM_USERNAME  << ": please advise\n";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}

	for(int iPayTimes=1; iPayTimes<fixPayTimes.size(); ++iPayTimes )
		if(fixPayTimes[iPayTimes]<= fixPayTimes[iPayTimes-1])
		{
			CC_Ostringstream os;
			os << "Trying to price a  SwapRate with fix pay times that are non strictly increasing: " 
				<< "fixPayTimes[" << iPayTimes-1 << "]: " << fixPayTimes[iPayTimes-1] << "days > "
				<< "fixPayTimes[" <<  iPayTimes	 << "]: " << fixPayTimes[iPayTimes]   << "days.\n"
				<< ARM_USERNAME  << ": please advise\n";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}


	for(int iPeriod=0; iPeriod<fixPayTimes.size(); ++iPeriod )
	{
		if( fixPayPeriods[iPeriod] <= 0 )
		{
			CC_Ostringstream os;
			os << "Trying to price a SwapRate with a non strictly positive period ["
				<< iPeriod << "] = (" << fixPayPeriods[iPeriod] << ").\n"
				<< ARM_USERNAME  << ": please advise\n";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}
	}
#endif

    double dfEvalTime = GetZeroCurve()->DiscountPrice(evalTime/K_YEAR_LEN);

	double dfStart	= GetZeroCurve()->DiscountPrice(floatStartTime/K_YEAR_LEN) / dfEvalTime;
	double dfEnd	= GetZeroCurve()->DiscountPrice(floatEndTime/K_YEAR_LEN) / dfEvalTime;

	int	i, size = fixPayTimes.size();
	double annuity = 0.0;
	
	for(i=0 ; i<size; ++i)
		annuity += GetZeroCurve()->DiscountPrice(fixPayTimes[i]/K_YEAR_LEN) / dfEvalTime 
		    * fixPayPeriods[i];

	double swap		= (dfStart-dfEnd) / annuity; 
	double result	= callPut == K_CALL ?
		CC_Max( swap - strikesPerState(0,0), 0.0 ) * annuity * swapNotional
	:	CC_Max( strikesPerState(0,0) - swap, 0.0 ) * annuity * swapNotional;

// FIXMEFRED: mig.vc8 (30/05/2007 16:11:38):cast
	return static_cast<ARM_VectorPtr>(new std::vector<double>(1,result));
}



////////////////////////////////////////////////////
///	Class  : ARM_IrFwdMod
///	Routine: Induct
///	Returns: ARM_PricingStatesPtr
///	Action : backward forward induct!
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_IrFwdMod::Induct(ARM_PricingStatesPtr& states,double toTime)
{
    double lastTimeStep = GetNumMethod()->GetLastTimeStep();

	if( lastTimeStep < 0.0 )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"Trying to induct in negative time, Please advise!");

	return GetNumMethod()->Induct( *this, states, toTime );
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: ProcessPaidPayoffs
///	Returns : void
///	Action  : change a paid payoff if necessary (useful for
///				change of measure)
////////////////////////////////////////////////////
void ARM_IrFwdMod::ProcessPaidPayoffs(const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
	const ARM_PricingStatesPtr& states ) const
{
#if defined( __GP_STRICT_VALIDATION )
	if( evalTime <  -K_NEW_DOUBLE_TOL )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			ARM_USERNAME + ": the evalTime should never be negative!" );
#endif

    double dfEvalTime = GetZeroCurve()->DiscountPrice(evalTime/K_YEAR_LEN);
	for(size_t i=0; i<payoffs->size();++i)
        (*payoffs)[i] *= dfEvalTime ;
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: ProcessPaidPayoffs
///	Returns : void
///	Action  : change a paid payoff if necessary (useful for
///				change of measure)
////////////////////////////////////////////////////
void ARM_IrFwdMod::ProcessUnPaidPayoffs(const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
	const ARM_PricingStatesPtr& states ) const
{
	/// is it a closed form?
	if( fabs(evalTime) >  K_NEW_DOUBLE_TOL )
	{
        double dfEvalTime = GetZeroCurve()->DiscountPrice(evalTime/K_YEAR_LEN);
	    for(size_t i=0; i<payoffs->size();++i)
            (*payoffs)[i] /= dfEvalTime ;
    }
}


////////////////////////////////////////////////////
///	Class   : ARM_IrFwdMod
///	Routines: LocalDrifts, 
///		Variances, NumMethodStateLocalGlobalVariances, VarianceToTime,
///		FirstPricingStates, ComputeModelTimes,PostInit
///		MCModelStatesFromToNextTime
///     TreeStatesToModelStates
///	Returns : throw an exception
///	Action  : not implemented because of no use
////////////////////////////////////////////////////
void ARM_IrFwdMod::IntegratedLocalDrifts(
	const std::vector<double>& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts) const
{
    throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                   "Unimplemented <IntegratedLocalDrifts> method");
}

void ARM_IrFwdMod::NumMethodStateLocalVariances( const std::vector<double>& timeSteps,
	ARM_MatrixVector& localVariances ) const
{
    throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                   "Unimplemented <NumMethodStateLocalVariances> method");
}

void ARM_IrFwdMod::ModelStateLocalVariances( const std::vector<double>& timeSteps,
	ARM_MatrixVector& localVariances ) const
{
    throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                   "Unimplemented <ModelStateLocalVariances> method");
}

double ARM_IrFwdMod::VarianceToTime(double var,double minTime,double maxTime) const
{
    throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                   "Unimplemented <VarianceToTime> method");
}

ARM_VectorPtr  ARM_IrFwdMod::VanillaSpreadOptionLet(const string& curveName,
													double evalTime,
													int callPut,
													double startTime,
													double endTime,
													double resetTime,
													double payTime,
													double payPeriod,
													double notional,
													double coeffLong,
													double coeffShort,
													const std::vector<double>& strikes,
													double swapLongFloatStartTime,
													double swapLongFloatEndTime,
													const std::vector<double>& swapLongFixPayTimes,
													const std::vector<double>& swapLongFixPayPeriods,
													double swapShortFloatStartTime,
													double swapShortFloatEndTime,
													const std::vector<double>& swapShortFixPayTimes,
													const std::vector<double>& swapShortFixPayPeriods,
													const ARM_PricingStatesPtr& states) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"VanillaSpreadOption : unimplemented function for ARM_IrFwdMod Model!");
}


ARM_PricingStatesPtr ARM_IrFwdMod::FirstPricingStates( size_t bucketSize ) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
        "Unimplemented <FirstPricingStates> method");
}

std::vector<double>& ARM_IrFwdMod::ComputeModelTimes(  const ARM_TimeInfoPtrVector& timeInfos )
{
	return new std::vector<double>(0);
}

////////////////////////////////////////////////////
///	Class  : ARM_IrFwdMod
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_IrFwdMod::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
        "Unimplemented <MCModelStatesFromToNextTime> method");
}

////////////////////////////////////////////////////
///	Class  : ARM_IrFwdMod
///	Routine: TreeStatesToModelStates
///	Returns: void 
///	Action : 
////////////////////////////////////////////////////
void ARM_IrFwdMod::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
        "Unimplemented <TreeStatesToModelStates> method");
}

////////////////////////////////////////////////////
///	Class   : ARM_IrFwdMod
///	Routines: AdviseBreakPointTimes
///	Returns : void
///	Action  : sets the corresponding suggested break point times to the model param
///////////////////////////////////////
void ARM_IrFwdMod::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, 
							  ARM_ModelParam* inputModelParam, 
							  size_t factorNb )
{
    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
        "Unknown type... an ARM_IrFwdMod does not supports any param" );
}

////////////////////////////////////////////////////
///	Class  : ARM_IrFwdMod
///	Routine: GetType
///	Returns: int
///	Action : tells the type of the model
////////////////////////////////////////////////////

int ARM_IrFwdMod::GetType() const
{
	return MT_NON_STOCHASTIC_MODEL;
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

