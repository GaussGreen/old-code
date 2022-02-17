/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file InfFwdMod.cpp
 *
 *  \brief
 *
 *	\author  A. Schauly
 *	\version 1.0
 *	\date March 2005
 */

/// this header comes firts as it includes some preprocessor constants! (because of modelnamemap)
#include "gpbase/removeidentifiedwarning.h"

/// this header comes firts as it includes some preprocessor constants!
#include "gpmodels/InfFwdMod.h"
#include "gpmodels/multiassets.h"

#include "gpbase/ostringstream.h"
#include "gpbase/utilityport.h"
#include "gpbase/datestrip.h"

#include "gpinfra/typedef.h"
#include "gpinfra/lexerdec.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/functorop.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/timeinfo.h"
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/pricingmodeltype.h"
#include "gpinfra/modelnamemap.h"


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_InfFwdMod
///	Routine: Constructor
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////
ARM_InfFwdMod::ARM_InfFwdMod(const ARM_InfCurvPtr& infc)
:	ARM_PricingModelInflation(infc)
{
}


///////////////////////////////////////////////////
///	Routine: toString
///	Returns: size of the deal description
///	Action : 
////////////////////////////////////////////////////
string ARM_InfFwdMod::toString(const string& indent, const string& nextIndent) const
{ 
	CC_Ostringstream os;
	os << "\n\n===========> ARM_InfFwdMod <===========\n";
	if( GetNumMethod() != ARM_NumMethodPtr(NULL) )
	    os << "with corresponding numerical method:\n" << GetNumMethod()->toString() << "\n";

	const ARM_PricingModelIR* PricingModelIR = getIRModel();

	if( !PricingModelIR )
		os << "with corresponding IRModel:" << PricingModelIR->toString() << "\n";

	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_InfFwdMod
///	Routine: Clone method
///	Returns: a new copy of this
///	Action : Copy this
////////////////////////////////////////////////////
ARM_Object* ARM_InfFwdMod::Clone() const
{
	return new ARM_InfFwdMod( *this );
}

////////////////////////////////////////////////////
///	Class  : ARM_InfFwdMod
///	Routine: UpdateLinks
///	Returns: void
///	Action : maintains links with itsIRModel
////////////////////////////////////////////////////
void ARM_InfFwdMod::UpdateLinks( const ARM_MultiAssetsModel& multiAssetsModel )
{
	const ARM_PricingModelIR* itsRefModel = getIRModel();
	const ARM_ModelNameMap* modelMap = multiAssetsModel.GetModelMap();

	string IRModelName, 
		itsModelName = GetModelName();

	/// Find the name of its IRModel
	if( !itsRefModel )
	{
		/// If there is no IRModel, its name is found through its Other Models in the ModelNameMap
		ARM_IntVector itsOtherModels = (*modelMap)[itsModelName]->OtherModelRefNb();

		if( itsOtherModels.size() != 1 )
			ARM_THROW( ERR_INVALID_ARGUMENT, "Wrong number of otherModels" );

		IRModelName = (*modelMap)[itsOtherModels[0]]->ModelName();
	}
	else
		IRModelName = itsRefModel->GetModelName();

	/// The updated IRModel is set. 
	setIRModel( (*modelMap)[IRModelName]->Model() );
}


////////////////////////////////////////////////////
///	Class  : ARM_InfFwdMod
///	Routine: Init
///	Returns: void
///	Action : Initialise before any pricing... for this
///			model, does nothing
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_InfFwdMod::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
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

	/// Checks if there is an IRModel (cannot work without). 
	if( !getIRModel() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": InfFwdModel works in a multiasset envirronment only. Please advise.");

	ARM_NumMethodPtr numMethod=GetNumMethod();

	size_t nbEvents = timeInfos.size();
    std::vector<double> timeSteps(nbEvents);
	for( int i=0; i<nbEvents; ++i)
		timeSteps[i] = timeInfos[i]->GetEventTime();

    numMethod->SetTimeSteps( std::vector<double>(timeSteps) );
	if(numMethod->GetPricingDirection() == ARM_NumMethod::GP_BCKWDLOOKING)
		numMethod->SetLastTimeIdx( nbEvents-1 );
	else if(numMethod->GetPricingDirection() == ARM_NumMethod::GP_FWDLOOKING)
		numMethod->SetLastTimeIdx(0);

	double firstInductTime = timeInfos[0]->GetEventTime();

	return numMethod->Init(*this,firstInductTime);
}



////////////////////////////////////////////////////
///	Class  : ARM_InfFwdMod
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_InfFwdMod::~ARM_InfFwdMod()
{}


////////////////////////////////////////////////////
///	Class  : ARM_InfFwdMod
///	Routine: Copy Constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_InfFwdMod::ARM_InfFwdMod( const ARM_InfFwdMod& rhs)
:	ARM_PricingModelInflation( rhs )
{}


////////////////////////////////////////////////////
///	Class  : ARM_InfFwdMod
///	Routine: Assignment operator
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_InfFwdMod& ARM_InfFwdMod::operator=( const ARM_InfFwdMod& rhs )
{
	if( this !=	 &rhs )
	    ARM_PricingModelInflation::operator = ( rhs );
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfFwdMod
///	Routine: Induct
///	Returns: ARM_PricingStatesPtr
///	Action : backward forward induct!
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_InfFwdMod::Induct(ARM_PricingStatesPtr& states,double toTime)
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
void ARM_InfFwdMod::ProcessPaidPayoffs(const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
	const ARM_PricingStatesPtr& states ) const
{
#if defined( __GP_STRICT_VALIDATION )
	if( evalTime <  -K_NEW_DOUBLE_TOL )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			ARM_USERNAME + ": the evalTime should never be negative!" );
#endif

    double dfEvalTime = getIRModel()->GetZeroCurve()->DiscountPrice(evalTime/K_YEAR_LEN);
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
void ARM_InfFwdMod::ProcessUnPaidPayoffs(const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, 
	const ARM_PricingStatesPtr& states ) const
{
	/// is it a closed form?
	if( fabs(evalTime) >  K_NEW_DOUBLE_TOL )
	{
        double dfEvalTime = getIRModel()->GetZeroCurve()->DiscountPrice(evalTime/K_YEAR_LEN);
	    for(size_t i=0; i<payoffs->size();++i)
            (*payoffs)[i] /= dfEvalTime ;
    }
}


////////////////////////////////////////////////////
///	Class   : ARM_InfFwdMod
///	Routines: LocalDrifts, 
///		Variances, NumMethodStateLocalGlobalVariances, VarianceToTime,
///		FirstPricingStates, ComputeModelTimes,PostInit
///		MCModelStatesFromToNextTime
///     TreeStatesToModelStates
///	Returns : throw an exception
///	Action  : not implemented because of no use
////////////////////////////////////////////////////
void ARM_InfFwdMod::IntegratedLocalDrifts(
	const std::vector<double>& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts) const
{

}

void ARM_InfFwdMod::NumMethodStateLocalVariances( const std::vector<double>& timeSteps,
	ARM_MatrixVector& localVariances ) const
{

}

void ARM_InfFwdMod::ModelStateLocalVariances( const std::vector<double>& timeSteps,
	ARM_MatrixVector& localVariances ) const
{
}

double ARM_InfFwdMod::VarianceToTime(double var,double minTime,double maxTime) const
{
	return 1.;
}

ARM_PricingStatesPtr ARM_InfFwdMod::FirstPricingStates( size_t bucketSize ) const
{
	return ARM_PricingStatesPtr( new ARM_PricingStates(bucketSize,1,0,1) );
}

std::vector<double>& ARM_InfFwdMod::ComputeModelTimes(  const ARM_TimeInfoPtrVector& timeInfos )
{
	return new std::vector<double>(0);
}

////////////////////////////////////////////////////
///	Class  : ARM_InfFwdMod
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_InfFwdMod::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{

}

////////////////////////////////////////////////////
///	Class  : ARM_InfFwdMod
///	Routine: TreeStatesToModelStates
///	Returns: void 
///	Action : 
////////////////////////////////////////////////////
void ARM_InfFwdMod::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
        "Unimplemented <TreeStatesToModelStates> method");
}

////////////////////////////////////////////////////
///	Class   : ARM_InfFwdMod
///	Routines: AdviseBreakPointTimes
///	Returns : void
///	Action  : sets the corresponding suggested break point times to the model param
///////////////////////////////////////
void ARM_InfFwdMod::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, 
							  ARM_ModelParam* inputModelParam, 
							  size_t factorNb )
{
    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
        "Unknown type... an ARM_InfFwdMod does not supports any param" );
}

////////////////////////////////////////////////////
///	Class  : ARM_InfFwdMod
///	Routine: GetType
///	Returns: int
///	Action : tells the type of the model
////////////////////////////////////////////////////

int ARM_InfFwdMod::GetType() const
{
	return MT_NON_STOCHASTIC_MODEL;
}


////////////////////////////////////////////////////
///	Class  : ARM_InfFwdMod
///	Routine: CPISpot
///	Returns: ARM_GP_VectorPtr
///	Action : returns CPISpot at a certain date and in a certain state of the world
////////////////////////////////////////
ARM_GP_VectorPtr ARM_InfFwdMod::CPISpot( const string& InfcurveName, double evalTime, double CPITime, string DCFLag, long DailyInterp,string ResetLag,const ARM_PricingStatesPtr& states) const
{
	ARM_GP_VectorPtr returnedVector = ARM_GP_VectorPtr( new std::vector<double>(1) );
	(*returnedVector)[0] = getInfCurve()->CPIInterpolate( CPITime,DCFLag, DailyInterp, ResetLag );
	return returnedVector; 
}

////////////////////////////////////////////////////
///	Class  : ARM_InfFwdMod
///	Routine: CPIForward
///	Returns: ARM_GP_VectorPtr
///	Action : returns CPIForward at a certain date and in a certain state of the world
////////////////////////////////////////
ARM_GP_VectorPtr ARM_InfFwdMod::CPIForward(	const string& InfcurveName, double evalTime, double CPITime, double FixingTime, const ARM_PricingStatesPtr& states) const 
{ 
	return ARM_GP_VectorPtr(NULL); 
}

////////////////////////////////////////////////////
///	Class  : ARM_InfFwdMod
///	Routine: ConvexityAdjustment
///	Returns: ARM_GP_VectorPtr
///	Action : returns CPIForward at a certain date and in a certain state of the world
////////////////////////////////////////
ARM_GP_VectorPtr ARM_InfFwdMod::ConvexityAdjustment(const string& InfcurveName, double evalTime,double tenor,double maturityTime, const ARM_PricingStatesPtr& states) const 
{ 
	ARM_GP_VectorPtr result( new std::vector<double>( 1 ) );
	(*result)[0] = 1;

	return result; 
}

////////////////////////////////////////////////////
///	Class  : ARM_InfFwdMod
///	Routine: ForwardCPIRatio
///	Returns: ARM_GP_VectorPtr
///	Action : returns ForwardCPIRatio at a certain date and in a certain state of the world
////////////////////////////////////////
ARM_GP_VectorPtr ARM_InfFwdMod::ForwardCPIRatio(const string& InfcurveName, double evalTime,double tenor,double CPITime, double maturityTime, const ARM_PricingStatesPtr& states) const 
{ 
	ARM_GP_VectorPtr returnedVector = ARM_GP_VectorPtr( new std::vector<double>(1) );
	(*returnedVector)[0] = getInfCurve()->CPIInterpolate( CPITime, CPITime );
	(*returnedVector)[0] /= getInfCurve()->CPIInterpolate( CPITime-tenor, CPITime-tenor );
	return returnedVector; 
}


////////////////////////////////////////////////////
///	Class  : ARM_InfFwdMod
///	Routine: YoYCap
///	Returns: ARM_GP_VectorPtr
///	Action : returns YoYCap Value at a certain date and in a certain state of the world
////////////////////////////////////////
ARM_GP_VectorPtr ARM_InfFwdMod::YoYCapFloor(const string& irCurveName,const string& infCurveName, double evalTime,double Strike,double FloatMargin,int CapFloor, const ARM_DateStripPtr& numDateStrip,const ARM_DateStripPtr& denomDateStrip,double itsSpread,const ARM_PricingStatesPtr& states) const
{
	if( denomDateStrip->GetResetDates()->Elt(0) < evalTime )
		Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
        "Try to evaluate a Cap in the past");

	size_t nbStates = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;

	/// This vector contains the num CPIs fixingDates (and not really the FlowStartDates!)
	std::vector<double> * numCPIDCFDates = numDateStrip->GetFlowStartDates();
	std::vector<double> * numCPIDates = numDateStrip->GetResetDates();
	std::vector<double> * denomCPIDates = denomDateStrip->GetResetDates();
	std::vector<double> * paymentDates = numDateStrip->GetFlowEndDates();
	std::vector<double> * NumInterestTerms = numDateStrip->GetInterestTerms();
	double dayCountCorrection;
	double spread = itsSpread;
	double AsOfDate = GetAsOfDate().GetJulian();

	size_t i, nbFlows = numCPIDates->size();
	ARM_GP_VectorPtr result( new std::vector<double>( nbStates ) );
	ARM_GP_VectorPtr numDFs, CPIRatios;
	std::vector<double>::iterator iterResult, iterDFs, iterCPIRatios; 

	// Compute the floating leg
	for( i = 0 ; i< nbFlows; i++ )
	{
		numDFs = DiscountFactor( irCurveName, evalTime, (*paymentDates)[i] - AsOfDate, states );
		CPIRatios = ForwardCPIRatio( infCurveName, evalTime, (*numCPIDates)[i]-(*denomCPIDates)[i], (*numCPIDates)[i]-AsOfDate, (*numCPIDCFDates)[i]-AsOfDate,states );
		iterDFs = numDFs->begin();
		iterCPIRatios = CPIRatios->begin();
		dayCountCorrection = NumInterestTerms->Elt(i);

		for( iterResult = result->begin() ; iterResult != result->end() ; iterResult++, iterDFs++, iterCPIRatios++ )
			(*iterResult) += MAX(0, ((*iterCPIRatios)-1.+spread-Strike)*(*iterDFs) * dayCountCorrection );
	}

	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfFwdMod
///	Routine: OATCapFloor
///	Returns: ARM_GP_VectorPtr
///	Action : returns YoYCap Value at a certain date and in a certain state of the world
////////////////////////////////////////
ARM_GP_VectorPtr ARM_InfFwdMod::OATCapFloor(const string& irCurveName,const string& infCurveName, double evalTime,double Strike,double FloatMargin,int CapFloor, const ARM_DateStripPtr& numDateStrip,const ARM_DateStripPtr& denomDateStrip,double itsSpread,const ARM_PricingStatesPtr& states) const
{
	if( denomDateStrip->GetResetDates()->Elt(0) < evalTime )
		Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
        "Try to evaluate a Cap in the past");

	size_t nbStates = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;

	/// This vector contains the num CPIs fixingDates (and not really the FlowStartDates!)
	std::vector<double> * numCPIDCFDates = numDateStrip->GetFlowStartDates();
	std::vector<double> * numCPIDates = numDateStrip->GetResetDates();
	std::vector<double> * denomCPIDates = denomDateStrip->GetResetDates();
	std::vector<double> * paymentDates = numDateStrip->GetFlowEndDates();
	std::vector<double> * NumInterestTerms = numDateStrip->GetInterestTerms();
	double dayCountCorrection;
	double spread = itsSpread;
	double AsOfDate = GetAsOfDate().GetJulian();
	double denomCPIDate = (*denomCPIDates)[0];

	size_t i, nbFlows = numCPIDates->size();
	ARM_GP_VectorPtr result( new std::vector<double>( nbStates ) );
	ARM_GP_VectorPtr numDFs, CPIRatios;
	std::vector<double>::iterator iterResult, iterDFs, iterCPIRatios; 

	// Compute the floating leg
	for( i = 0 ; i< nbFlows; i++ )
	{
		numDFs = DiscountFactor( irCurveName, evalTime, (*paymentDates)[i] - AsOfDate, states );
		CPIRatios = ForwardCPIRatio( infCurveName, evalTime, (*numCPIDates)[i]-denomCPIDate, (*numCPIDates)[i]-AsOfDate, (*numCPIDCFDates)[i]-AsOfDate,states );
		iterDFs = numDFs->begin();
		iterCPIRatios = CPIRatios->begin();
		dayCountCorrection = NumInterestTerms->Elt(i);

		for( iterResult = result->begin() ; iterResult != result->end() ; iterResult++, iterDFs++, iterCPIRatios++ )
			(*iterResult) += MAX(0, ((*iterCPIRatios)+spread-Strike)*(*iterDFs) * dayCountCorrection );
	}

	return result;
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

