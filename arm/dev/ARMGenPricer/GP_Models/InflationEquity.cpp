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
 *	\date April 2005
 */

/// this header comes firts as it includes some preprocessor constants! (because of modelnamemap)
#include "gpbase/removeidentifiedwarning.h"

#include "gpbase/ostringstream.h"
#include "gpbase/utilityport.h"
#include "gpbase/datestrip.h"
#include "gpbase/vectormanip.h"

#include "gpmodels/InflationEquity.h"
#include "gpmodels/InflationEquityModelParams.h"
#include "gpmodels/multiassets.h"
#include "gpmodels/ModelParamsHW1F.h"

#include "gpinfra/typedef.h"
#include "gpinfra/lexerdec.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/functorop.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/timeinfo.h"
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/pricingmodeltype.h"
#include "gpinfra/discretisationscheme.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/modelnamemap.h"
#include "gpinfra/curvemodelparam.h"

#include "gpclosedforms/vanilla_bs.h"

/// kernel
#include <inst/portfolio.h>


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: Constructor
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////
ARM_InflationEquityMod::ARM_InflationEquityMod(const ARM_InfCurvPtr& infc, const ARM_PricingModelIRPtr& irMod, double publicationLag, const ARM_ModelParams* params)
:	ARM_PricingModelInflation(infc,params), itsInitValue(0),itsCPIManager(new CPIManager()),itsPublicationLag(publicationLag), itsOtherModelRank(0)
{
}


///////////////////////////////////////////////////
///	Routine: toString
///	Returns: size of the deal description
///	Action : 
////////////////////////////////////////////////////
string ARM_InflationEquityMod::toString(const string& indent, const string& nextIndent) const
{ 
	CC_Ostringstream os;
	os << indent << "ARM_InflationEquityModel\n";
	os << indent << "------------------------\n";
	os << ARM_PricingModel::toString(indent);
	os << "\n\n";

	const ARM_PricingModelIR* PricingModelIR = getIRModel();

	if( !PricingModelIR )
		os << indent <<"with corresponding IRModel: \n" << PricingModelIR->toString() << "\n";

	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: Clone method
///	Returns: a new copy of this
///	Action : Copy this
////////////////////////////////////////////////////
ARM_Object* ARM_InflationEquityMod::Clone() const
{
	return new ARM_InflationEquityMod( *this );
}

////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: PostProcessing
///	Returns: PostProcessing
///	Action : PostProcessing
////////////////////////////////////////////////////
void ARM_InflationEquityMod::PostProcessing(const ARM_ModelFitter& modelFitter)
{
	((ARM_ModelParamsInflationEquity*) GetModelParams())->PostProcessing(modelFitter, this);
}


////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: PricingTimeSteps
///	Returns: ARM_GP_Vector*
///	Action : Computes Time Steps for diffusion
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_InflationEquityMod::PricingTimeSteps(const ARM_TimeInfoPtrVector& timeInfos)
{
	ARM_DiscretisationScheme& discretisationScheme = ARM_EventTime();
	ARM_GP_Vector* ptimeSteps = discretisationScheme.ModelTimesFromTimesInfo(timeInfos, *this );
	itsCPIManager = CPIManagerPtr( new CPIManager );

	ARM_TimeInfoPtrVector::const_iterator timeInfoIterator, timeInfoBegin = timeInfos.begin(), timeInfoEnd = timeInfos.end();

	for( timeInfoIterator = timeInfoBegin ; timeInfoIterator != timeInfoEnd ; timeInfoIterator++ )
	{
		ARM_AdditionalTimeInfoPtrVector AdditionalTimeInfoVector = (*timeInfoIterator)->getAdditionalTimeInfos();
		ARM_AdditionalTimeInfoPtrVector::iterator iterAdditionalInfo, AdditionalInfoBegin = AdditionalTimeInfoVector.begin(), 
			AdditionalInfoEnd = AdditionalTimeInfoVector.end();

		for( iterAdditionalInfo = AdditionalInfoBegin ; iterAdditionalInfo != AdditionalInfoEnd ; iterAdditionalInfo++ )
		{
			CPITimeInfo* cpiTimeInfo = dynamic_cast<CPITimeInfo*> (&**iterAdditionalInfo);

			if( cpiTimeInfo )
			{
				double CPITime = cpiTimeInfo->getCPITime();
				double CPIPublicationTime = CPITime+itsPublicationLag;

				if( CPIPublicationTime >= K_NEW_DOUBLE_TOL )
				{
					itsCPIManager->AddDate( CPITime );
					itsCPIManager->EnqueueCPIInfo( CPITime, CPIInfoPtr( new CPIInfo( cpiTimeInfo->getInfcurveName(), 
						cpiTimeInfo->getDCFLag(), cpiTimeInfo->getResetLag(), cpiTimeInfo->getDailyInterp() ) ) );

					cpiTimeInfo->setPublishTime(CPIPublicationTime);
					ptimeSteps->push_back(CPIPublicationTime);
				}
			}
		}
	}
	ARM_GP_Vector* newptimeSteps = SortSTLBased( *ptimeSteps );
	delete ptimeSteps;
	ptimeSteps = VectorUnique( *newptimeSteps );
	delete newptimeSteps;
	return ptimeSteps;
}


////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: Init
///	Returns: void
///	Action : Initialise before any pricing... for this
///			model, does nothing
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_InflationEquityMod::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
    int nbEvents=timeInfos.size();
    bool isSpotUse = (nbEvents == 0);

	/// Checks if there is an IRModel
	if( !getIRModel() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": InflationEquityModel works in a multiasset envirronment only. Please advise.");
	
    if(!isSpotUse)
    {
        // Numerical method and numeraire are needed to go on
        ARM_NumMethodPtr numMethod=GetNumMethod();
		if( numMethod == ARM_NumMethodPtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numerical method not set in Inflation Equity model!");

        /// test the numeraire and its type!
		ARM_NumerairePtr numeraire=GetNumeraire();
        if( numeraire == ARM_NumerairePtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numeraire not set in the Inflation Equity model!");

		/// creates the model schedule (smart pointor for exception safety!)
		ARM_DiscretisationScheme& discretisationScheme = ARM_EventTime();
		CC_NS(std,auto_ptr)<ARM_GP_Vector> ptimeSteps( discretisationScheme.ModelTimesFromTimesInfo(timeInfos, *this) );

        //// Initialise the numeraire
		numeraire->Init(*(GetDiscountFunctor()),payModelName,timeInfos);

		/// Set the basic schedule in the numerical method and...
		numMethod->SetTimeSteps(*ptimeSteps);

		double firstInductTime = timeInfos[0]->GetEventTime();

		/// ...initialise it
		return numMethod->Init(*this,firstInductTime);
    }
    else
    {
        // Compute a single model states set to (0.0,...,0.0)
        int nbDir = FactorCount();
        ARM_PricingStatesPtr initStates(new ARM_PricingStates(1,nbDir,0));
        for(int i=0;i<nbDir;++i)
            initStates->SetModelState(0,i,0.0);
		
        return initStates;
    }
}



////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_InflationEquityMod::~ARM_InflationEquityMod()
{
}


////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: Copy Constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_InflationEquityMod::ARM_InflationEquityMod( const ARM_InflationEquityMod& rhs)
:	ARM_PricingModelInflation( rhs ), itsInitValue( rhs.itsInitValue ), itsPublicationLag(rhs.itsPublicationLag), itsCPIManager(rhs.itsCPIManager),itsOtherModelRank(rhs.itsOtherModelRank)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: Assignment operator
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_InflationEquityMod& ARM_InflationEquityMod::operator=( const ARM_InflationEquityMod& rhs )
{
	if( this !=	 &rhs )
	    ARM_PricingModelInflation::operator = ( rhs );
	return *this;
}

////////////////////////////////////////////////////
///	Class   : ARM_InflationEquityMod
///	Routines: AdviseCurrentCalib
///	Returns : void
///	Action  : Advises the model that of the calibration  ... 
///           The model advises just the model params
////////////////////////////////////////////////////
void ARM_InflationEquityMod::AdviseCurrentCalib(ARM_ModelFitter& modelFitter)
{
	GetModelParams()->PostProcessing(modelFitter,this); 
}

////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: UpdateLinks
///	Returns: void
///	Action : maintains links with itsIRModel
////////////////////////////////////////////////////
void ARM_InflationEquityMod::UpdateLinks( const ARM_MultiAssetsModel& multiAssetsModel )
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
	itsOtherModelRank = getIRModel()->GetModelRank();

	/// Correlation between EquityInflation and IRModel is updated. 
	ARM_GP_MatrixPtr itsCorrelMatrix = multiAssetsModel.GetCorrelSubMatrix( IRModelName, itsModelName );

	/// Check that the IR/Inflation CorrelMatrix has the right size
	if( itsCorrelMatrix->cols() != 2 || itsCorrelMatrix->rows() != 2 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "Inflation/IR Correl Matrix has the wrong size!" );

	setCorrelWithBonds( itsCorrelMatrix->Elt(0,1) );
}

////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: IntegratedLocalDrifts
///	Returns: A vector saving the local drift of
///          the state variable
///	Action : Compute local relative drifts of the
///          state variable between each time step
////////////////////////////////////////////////////
void ARM_InflationEquityMod::IntegratedLocalDrifts(
	const ARM_GP_Vector& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts) const
{
	/// Dates we want to compute local drifts for
	size_t nbSteps	= timeSteps.size();
    double step		= timeSteps[0],nextStep;

	relativeDrifts = ARM_GP_MatrixPtr(NULL);
	absoluteDrifts	= ARM_GP_MatrixPtr( new ARM_GP_Matrix( nbSteps-1, 1, 0.0 ) );

	for(size_t i=0;i<nbSteps-1;++i)
	{
		/// [i] => local variance from ti->ti+1
		nextStep=timeSteps[i+1];
		(*absoluteDrifts)(i,0) = 0;
		step=nextStep;
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: NumMethodStateLocalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local and global variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_InflationEquityMod::NumMethodStateLocalVariances( const ARM_GP_Vector& timeSteps,
	ARM_MatrixVector& localVariances ) const
{
}

////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: NumMethodStateLocalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local and global variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_InflationEquityMod::ModelStateLocalVariances( const ARM_GP_Vector& timeSteps,
	ARM_MatrixVector& localVariances ) const
{
	size_t nbSteps	= timeSteps.size();
	size_t modelRank= GetModelRank();
    double step		= timeSteps[0], nextStep;
	size_t offsetIndex	= (nbSteps-1)*modelRank;
	double localVariance;
	double rho = getCorrelWithBonds();
	double alpha = ((const ARM_ModelParamsInflationEquity* const) GetModelParams())->getMultiplier();

	const ARM_ModelParamsHW1F* otherModelParams = dynamic_cast<const ARM_ModelParamsHW1F*> (getIRModel()->GetModelParams() ); 
	const ARM_CurveModelParam * cmp = dynamic_cast<const ARM_CurveModelParam *> (&(GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility )));

#if defined(__GP_STRICT_VALIDATION)
	if( localVariances.size()!= offsetIndex ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "localDrifts.size() != offsetIndex" );
#endif
	localVariances.resize((nbSteps-1)*(modelRank+1));
	size_t i;

	for(i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];
		/// [i] => local variance from ti->ti+1

		localVariance = -alpha*2*rho*otherModelParams->HW1FEqFxZcCovariance( *cmp, otherModelParams, step, nextStep, nextStep );
		localVariance += alpha*alpha*otherModelParams->HW1FZcCovariance( otherModelParams, otherModelParams, step, nextStep, nextStep );
		localVariance += ((const ARM_ModelParamsInflationEquity* const) GetModelParams())->StateLocalVariance(step,nextStep);

		localVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(1,localVariance );
		step=nextStep;
	}

}

void ARM_InflationEquityMod::NumMethodStateLocalCovariances( const ARM_GP_Vector& timeSteps, ARM_MatrixVector& localCorrels, const ARM_MultiAssetsModel& multiassets )
{
	size_t modelRank = GetModelRank();
	const ARM_ModelNameMap* const modelMap = multiassets.GetModelMap();
	size_t i,nbSteps = timeSteps.size();
    double step		 = timeSteps[0],nextStep;
	double rho		 = getCorrelWithBonds();
	double cov;
	double denom=0.;
	const ARM_PricingModelIR* itsIRModel = getIRModel();
	size_t otherModelRank = (*modelMap)[itsIRModel->GetModelName()]->Model()->GetModelRank();
	string itsModelName = GetModelName();
	double alpha = ((const ARM_ModelParamsInflationEquity* const) GetModelParams())->getMultiplier();

	const ARM_ModelParamsHW1F* otherModelParams = dynamic_cast<const ARM_ModelParamsHW1F*> (itsIRModel->GetModelParams() );

	if( !otherModelParams )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "inflation equity works with HW1F only !");

	const ARM_CurveModelParam * cmp = dynamic_cast<const ARM_CurveModelParam *> (&(GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility )));
	
	for(i=0;i<nbSteps-1;++i)
	{
		nextStep = timeSteps[i+1];

		cov =  -alpha*otherModelParams->HW1FStateZcCovariance(otherModelParams, otherModelParams, step, nextStep, nextStep, nextStep);
		cov += rho*otherModelParams->HW1FEqFxStateCovariance( *cmp, otherModelParams, step, nextStep, nextStep );

		denom = -alpha*2*rho*otherModelParams->HW1FEqFxZcCovariance( *cmp, otherModelParams, step, nextStep, nextStep );
		denom += alpha*alpha*otherModelParams->HW1FZcCovariance( otherModelParams, otherModelParams, step, nextStep, nextStep );
		denom += ((const ARM_ModelParamsInflationEquity* const) GetModelParams())->StateLocalVariance(step,nextStep);
/*		cov /= sqrt(denom);
		denom = otherModelParams->StateLocalVariance(step,nextStep,nextStep);
		cov /= sqrt(denom);*/
		(*(localCorrels[i]))(modelRank,otherModelRank) = cov;
		(*(localCorrels[i]))(otherModelRank,modelRank) = cov;
		(*(localCorrels[i]))(modelRank,modelRank) = denom;
		step = nextStep;
	}
}

///////////////////////////////////////////////////
///	Class   : ARM_HullWhite1F
///	Routine : FirstPricingStates,
///	Returns : ARM_PricingStatesPtr
///	Action  : create the first pricing states
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_InflationEquityMod::FirstPricingStates( size_t bucketSize ) const
{
	ARM_NumerairePtr numeraire = GetNumeraire();
	const ARM_PricingModelIR* itsIRModel = getIRModel();

	#if defined(__GP_STRICT_VALIDATION)
	/// Checks that the numeraire is not NULL
		if( numeraire == ARM_NumerairePtr(NULL) )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				": Numeraire should not be NULL");
	#endif

	size_t modelNb	= GetModelNb();

	/// Check that the numeraire is TerminalZc or TerminalEventZc
	ARM_Numeraire::NumeraireType type = numeraire->GetType();

	/// Coeffiecient alpha of the model (multiplier) 
	double alpha = ((const ARM_ModelParamsInflationEquity* const) GetModelParams())->getMultiplier();

	/// First model states are set to log( 1/B(0,T)^alpha ) since we diffuse CPISpot(t)/(B(t,T)^alpha * exp( - int_0^t q(s)ds ) * B(0,t)^alpha )
	itsInitValue = log( getInfCurve()->CPIInterpolate( 0,"0M", K_CPILINEAR, "0M" ) );

	/// We finally build the first Pricing states (filled with itsInitValue): 
	/// ARM_PricingStates(nbStates = bucketSize, nbModelStates = 1F , nbPayoffs = 0)
	ARM_PricingStatesPtr pricingStates( new ARM_PricingStates(bucketSize,1,0,1) );

	/// Resetthe CPI Manager
	itsCPIManager->RebuildCPIMap();

	return pricingStates;
}

ARM_GP_Vector* ARM_InflationEquityMod::ComputeModelTimes(  const ARM_TimeInfoPtrVector& timeInfos )
{
	/// No specific ModelTimes to be added
	return new ARM_GP_Vector(0);
}

////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
void ARM_InflationEquityMod::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
#if defined(__GP_STRICT_VALIDATION)
	if(timeIndex<0)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index is negative!");
	if( timeIndex >= GetNumMethod()->GetTimeSteps()->size()-1 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index bigger than max size!");
#endif

	/// ThisTime, ToTime
	double nextTime = GetNumMethod()->GetTimeStep(timeIndex+1);
	double thisTime = GetNumMethod()->GetTimeStep(timeIndex);

	if( thisTime > -K_NEW_DOUBLE_TOL )
	{
		const ARM_MatrixVector& localVar	= GetModelStateLocalVars();
		const ARM_MatrixVector& localStdDev = GetModelStateLocalStdDevs();

		/// useful for computation
		size_t statesNb = states->size();
		size_t modelNb	= GetModelNb();
		double currentState;
		double currentStateCorrection = log( getInfCurve()->CPIInterpolate( 0,"0M", K_CPILINEAR, "0M" ) );

		/// Updates modelStates
		for( size_t i=0;i<statesNb; ++i )
		{
			currentState = states->GetNumMethodState(i,modelNb)+currentStateCorrection;
			states->SetModelState(i,modelNb,currentState);
		}

		double CPITime = nextTime - itsPublicationLag;

		if( itsCPIManager->doesExist(CPITime) )
		{
			CPIInfoPtr itsCPIInfo = itsCPIManager->getCPIInfo( CPITime );
			ComputeCPISpotAndStore( itsCPIInfo->itsInfcurveName, nextTime, CPITime, itsCPIInfo->itsDCFLag,
				itsCPIInfo->itsDailyInterp,itsCPIInfo->itsResetLag, states );
		}

#if defined( __GP_STRICT_VALIDATION )
		if( timeIndex == GetNumMethod()->GetTimeSteps()->size()-1 )
		{
			if( (itsCPIManager->getCPISize() > 1) || ((itsCPIManager->getCPISize() == 1)&& (itsCPIManager->getOldestCPICounter() > 1)) )
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "CPI Manager is not empty!");
		}
#endif
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "trying to diffuse InflationEquity from the past!");
}

////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: TreeStatesToModelStates
///	Returns: void 
///	Action : 
////////////////////////////////////////////////////
void ARM_InflationEquityMod::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
        "Unimplemented <TreeStatesToModelStates> method");
}

////////////////////////////////////////////////////
///	Class   : ARM_InflationEquityMod
///	Routines: AdviseBreakPointTimes
///	Returns : void
///	Action  : sets the corresponding suggested break point times to the model param
///////////////////////////////////////
void ARM_InflationEquityMod::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, 
							  ARM_ModelParam* inputModelParam, 
							  size_t factorNb )
{
	ARM_CurveModelParam* modelParam = dynamic_cast<ARM_CurveModelParam*>(inputModelParam);
	if( !modelParam )
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "expected an ARM_CurveModelParam!");

    double asOfDate = GetAsOfDate().GetJulian();
    int size1       = portfolio->GetSize();  
    ARM_GP_Vector  tmpdates;
    int i;
    
    switch( modelParam->GetType() )
    {
    case ARM_ModelParamType::Volatility:
        {
            double date = portfolio->GetAsset(0)->GetPaymentDates()->Elt(portfolio->GetAsset(0)->GetPaymentDates()->size()-1) - asOfDate;
            tmpdates.push_back(date);
            for(i=1; i<size1; i++) 
            {
                double resetlag = portfolio->GetAsset(i)->GetPaymentDates()->Elt(portfolio->GetAsset(i)->GetPaymentDates()->size()-1) - asOfDate;
                if(fabs (date - resetlag) > FRMVOL_LAG_THRESHOLD)
                {
                    tmpdates.push_back(resetlag);
                    date = resetlag;
                }
				else
				{
					/// ignore this instrument
					portfolio->SetWeight(0.0,i);
				}
			}
			modelParam->UpdateValues(&tmpdates);
        }
        break;
    case ARM_ModelParamType::Multiplier:  break;
        
    default:
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
            "Unknown type... an Inflation Equity model only supports mean Multiplier and Volatility" );
    }
}

////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: VarianceToTime, 
///	Returns: 
///	Action : NOT IMPLEMENTED YET
////////////////////////////////////////////////////
double ARM_InflationEquityMod::VarianceToTime(double var,double minTime,double maxTime) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
          "VarianceToTime : Not Implemented");
	return 1.;
}

////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: GetType
///	Returns: int
///	Action : tells the type of the model
////////////////////////////////////////////////////

int ARM_InflationEquityMod::GetType() const
{
	return MT_EQUITY_MODEL;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// Here starts the Inflation Part of the model //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: ComputeCPISpotAndStore
///	Returns: void
///	Action : computes CPISpot at a certain date and in a certain state of the world
///   This is the implementation of the reconstruction formula. 
////////////////////////////////////////
void ARM_InflationEquityMod::ComputeCPISpotAndStore( const string& InfcurveName, double evalTime, double CPITime, string DCFLag, long DailyInterp,string ResetLag,const ARM_PricingStatesPtr& states) const
{

	ARM_GP_VectorPtr returnedVector = itsCPIManager->GetCPI( CPITime );
#if defined( __GP_STRICT_VALIDATION )
	if( returnedVector == ARM_GP_VectorPtr(NULL) )
	{
#endif	
		/// Multiplier alpha of the numeraire
		double alpha = ((const ARM_ModelParamsInflationEquity* const) GetModelParams())->getMultiplier();
		const ARM_PricingModelIR* itsIRModel = getIRModel();
		const ARM_ModelParamsHW1F* otherModelParams = dynamic_cast<const ARM_ModelParamsHW1F*> (itsIRModel->GetModelParams() ); 
		double rho = getCorrelWithBonds();
		const ARM_CurveModelParam * cmp = dynamic_cast<const ARM_CurveModelParam *> (&(GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility )));
		double resultImprovement = 0;

#if defined( __GP_STRICT_VALIDATION )
		if( !otherModelParams )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "inflation equity works with HW1F only !");
#endif

		ARM_NumerairePtr numeraire=GetNumeraire();

		if( numeraire->GetType() == ARM_Numeraire::TerminalZc || numeraire->GetType() == ARM_Numeraire::TerminalEventZc )
		{
			double probaDate= numeraire->GetMaturity();
			resultImprovement = alpha * otherModelParams->HW1FZcCovariance( otherModelParams, otherModelParams, 0,evalTime, evalTime );
			resultImprovement -= alpha * itsIRModel->IntegratedBondCovariance( 0, evalTime, probaDate, evalTime );

			resultImprovement += alpha*rho*otherModelParams->HW1FEqFxZcCovariance( *cmp, otherModelParams, 0, evalTime, evalTime );
			resultImprovement -= 0.5*alpha*alpha*otherModelParams->HW1FZcCovariance( otherModelParams, otherModelParams, 0, evalTime, evalTime );
			resultImprovement -= 0.5*((const ARM_ModelParamsInflationEquity* const) GetModelParams())->StateLocalVariance(0,evalTime);
			resultImprovement -= rho * getIRModel()->VolatilityScalarProduct( 0, evalTime, evalTime,  GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility) );
			resultImprovement += rho * getIRModel()->VolatilityScalarProduct( 0, evalTime, probaDate,  GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility) );
			resultImprovement = exp(resultImprovement);
		}
		else
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Inflation EquityModel has been implemented for Terminal numeraires only!");
		}

		size_t i,statesNb = states->size();
		size_t modelNb	= GetModelNb();

		double dividend = getInfCurve()->CPIInterpolate( CPITime,DCFLag, DailyInterp, ResetLag )/ getInfCurve()->CPIInterpolate( 0,"0M", DailyInterp, "0M" );

		/// resultVector is filled (at last!)
		returnedVector = ARM_GP_VectorPtr( new ARM_GP_Vector( statesNb ) );

		for( i = 0 ; i < statesNb ; i++ ) 
			(*returnedVector)[i] = exp( states->GetModelState(i, modelNb ) )* dividend*resultImprovement;
		/// CPI = exp( propagated variable ) * integrated dividend

		itsCPIManager->setCPI( CPITime, returnedVector );
#if defined( __GP_STRICT_VALIDATION )
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": CPISpot has already been computed. Please advise.");
#endif

}

////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: IntegratedDividend
///	Returns: double
///	Action : Computes Exp ( - Int_FromTime^ToTime alpha*q(s)ds )
///    Expressed using CPIForwards and volatilities
////////////////////////////////////////
double ARM_InflationEquityMod::IntegratedDividend( double FromTime, double ToTime ) const
{
	if( FromTime < K_DOUBLE_TOL )
	{
		const ARM_PricingModelIR* itsIRModel = getIRModel();
		double alpha = ((const ARM_ModelParamsInflationEquity* const) GetModelParams())->getMultiplier();
		double rho = getCorrelWithBonds();
		double probaDate= GetNumeraire()->GetMaturity();

		double result = 0.5*alpha*getIRModel()->IntegratedBondSquaredVol( FromTime, ToTime, ToTime );
		result -= rho * getIRModel()->VolatilityScalarProduct( FromTime, ToTime, ToTime,  GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility) );
		result = exp( (1-alpha) * result );
		result /= pow( itsIRModel->GetZeroCurve()->DiscountPrice(ToTime/K_YEAR_LEN),alpha);

		return result;
	}
	else
	{
		return IntegratedDividend(0, ToTime)/IntegratedDividend(0, FromTime);
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: CPISpot
///	Returns: ARM_GP_VectorPtr
///	Action : returns CPISpot at a certain date and in a certain state of the world
///   This is the implementation of the reconstruction formula. 
////////////////////////////////////////
ARM_GP_VectorPtr ARM_InflationEquityMod::CPISpot( const string& InfcurveName, double evalTime, double CPITime, string DCFLag, long DailyInterp,string ResetLag,const ARM_PricingStatesPtr& states) const
{
	ARM_GP_VectorPtr returnedVector = itsCPIManager->GetCPI( CPITime );

	if( returnedVector != ARM_GP_VectorPtr(NULL) )
		itsCPIManager->decreaseCounterAndDeleteCPI( CPITime );
	else
		if( CPITime + itsPublicationLag < K_NEW_DOUBLE_TOL )
		{
			returnedVector = ARM_GP_VectorPtr( new ARM_GP_Vector( states->size() ));
			double result = getInfCurve()->CPIInterpolate( CPITime,DCFLag, DailyInterp, ResetLag );
			for( ARM_GP_Vector::iterator iter = returnedVector->begin() ; iter != returnedVector->end() ; iter++ )
				(*iter) = result;
		}
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": CPI did not reset. Please advise.");

	return static_cast<ARM_GP_Vector*> (returnedVector->Clone());
}

////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: CPIForward
///	Returns: ARM_GP_VectorPtr
///	Action : returns CPIForward at a certain date and in a certain state of the world
///  Generally speaking, CPITime - FixingTime = -90
////////////////////////////////////////
ARM_GP_VectorPtr ARM_InflationEquityMod::CPIForward(const string& InfcurveName, double evalTime, double CPITime, double FixingTime, const ARM_PricingStatesPtr& states) const
{ 
	size_t i,statesNb = states->size();

	if( evalTime > K_DOUBLE_TOL )
	{
		/// Get useful params
		double alpha = ((const ARM_ModelParamsInflationEquity* const) GetModelParams())->getMultiplier();
		double rho = getCorrelWithBonds();
		const ARM_ModelParam * vol = &(GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility));

		size_t modelNb	= GetModelNb();

		/// Computes integrated dividend
		double dividend = getInfCurve()->CPIInterpolate( CPITime, CPITime );
		dividend *= IntegratedDividend( 0, FixingTime ) / getInfCurve()->CPIInterpolate( 0,"0M", K_CPILINEAR, "0M" );

		ARM_GP_Vector * result = new ARM_GP_Vector( statesNb );

		/// Changing it into forwards. 
		ARM_GP_VectorPtr DFs = DiscountFactor( InfcurveName, evalTime, FixingTime,states);

		double CorrectionFwd;
		CorrectionFwd =  alpha*(alpha-1)*0.5*getIRModel()->IntegratedBondSquaredVol( evalTime, FixingTime, FixingTime );
		CorrectionFwd += (1-alpha)*rho*getIRModel()->VolatilityScalarProduct( evalTime, FixingTime, FixingTime, *vol );

		for( i = 0 ; i < statesNb ; i++ )
			(*result)[i] = exp( CorrectionFwd + states->GetModelState(i, modelNb ) )* dividend / pow( (*DFs)[i],alpha );
		/// CPI = exp( propagated variable ) * integrated dividend

		return result;
	}
	else
	{
		ARM_GP_Vector * result = new ARM_GP_Vector( statesNb );
		ARM_GP_Vector::iterator iter = result->begin(), iterEnd = result->end();
		double ReturnedCPIForward = getInfCurve()->CPIInterpolate( CPITime, CPITime );
		
		for( ; iter!=iterEnd ; iter++ )
			(*iter) = ReturnedCPIForward;

		return result;
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: ConvexityAdjustment
///	Returns: ARM_GP_VectorPtr
///	Action : returns CPIForward at a certain date and in a certain state of the world
////////////////////////////////////////
ARM_GP_VectorPtr ARM_InflationEquityMod::ConvexityAdjustment(const string& InfcurveName, double evalTime,double tenor,double maturityTime, const ARM_PricingStatesPtr& states) const 
{ 
	size_t statesSize = states->size();
	double alpha = ((const ARM_ModelParamsInflationEquity* const) GetModelParams())->getMultiplier();
	double rho = getCorrelWithBonds();
	const ARM_ModelParam * vol = &(GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility));
	double convexityAdj = 0;

	const ARM_PricingModelIR* itsIRModel = getIRModel();

	/// - int_t^T Gamma(s,T-delta)^2 ds. note that gamma(t,T-delta) = 0 for t>T-delta
	convexityAdj -= itsIRModel->IntegratedBondSquaredVol( evalTime, maturityTime-tenor, maturityTime-tenor );
	/// int_t^T Gamma(s,T)*Gamma(s,T-delta) ds. note that gamma(t,T-delta) = 0 for t>T-delta
	convexityAdj += itsIRModel->IntegratedBondCovariance( evalTime, maturityTime-tenor, maturityTime-tenor, maturityTime);
	convexityAdj *= alpha;

	/// rho * int_t^T Gamma(s,T-delta) * sigma(s) ds 
	convexityAdj += rho * itsIRModel->VolatilityScalarProduct( evalTime, maturityTime-tenor, maturityTime-tenor, *vol );
	/// -rho * int_t^T Gamma(s,T) * sigma(s) ds 
	convexityAdj -= rho * itsIRModel->VolatilityScalarProduct( evalTime, maturityTime-tenor, maturityTime, *vol );
	convexityAdj *= (1-alpha);
	convexityAdj = exp( convexityAdj );


	ARM_GP_VectorPtr result( new ARM_GP_Vector( statesSize ) );

	for( size_t i = 0 ; i < statesSize ; i++ )
		result->Elt(i) = convexityAdj;

	return result; 
}

////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: ForwardCPIRatio
///	Returns: ARM_GP_VectorPtr
///	Action : returns ForwardCPIRatio at a certain date and in a certain state of the world
////////////////////////////////////////
ARM_GP_VectorPtr ARM_InflationEquityMod::ForwardCPIRatio(const string& InfcurveName, double evalTime,double tenor,double CPITime, double maturityTime, const ARM_PricingStatesPtr& states) const 
{ 
	ARM_GP_VectorPtr numVector = CPIForward( InfcurveName, evalTime, CPITime, maturityTime, states );
	ARM_GP_VectorPtr denomVector = CPIForward( InfcurveName, evalTime, CPITime-tenor, maturityTime-tenor, states );
	ARM_GP_VectorPtr Adjustment = ConvexityAdjustment( InfcurveName, evalTime, tenor, maturityTime, states );

	ARM_GP_Vector::iterator iter1,iter2,iter3, iterEnd = numVector->end();

	for( iter1 = numVector->begin(), iter2 = denomVector->begin(), iter3 = Adjustment->begin() ; iter1 != iterEnd ; iter1++,iter2++,iter3++ )
		(*iter1) *= (*iter3)/(*iter2);

	return numVector; 
}

////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: YoYCap
///	Returns: ARM_GP_VectorPtr
///	Action : returns YoYCap at a certain date and in a certain state of the world
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_InflationEquityMod::YoYCapFloor( const string& irCurveName,	const string& infCurveName, double evalTime,double Strike,double FloatMargin, int CapFloor, const ARM_DateStripPtr& numDateStrip,const ARM_DateStripPtr& denomDateStrip,double itsSpread,	const ARM_PricingStatesPtr& states) const
{
	size_t nbStates = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;

	/// This vector contains the num CPIs fixingDates (and not really the FlowStartDates!)
	ARM_GP_Vector * numCPIDCFDates = numDateStrip->GetFlowStartDates();
	ARM_GP_Vector * numCPIDates = numDateStrip->GetResetDates();
	ARM_GP_Vector * denomCPIDates = denomDateStrip->GetResetDates();
	ARM_GP_Vector * paymentDates = numDateStrip->GetFlowEndDates();
	ARM_GP_Vector * NumInterestTerms = numDateStrip->GetInterestTerms();
	double dayCountCorrection;
	double strike = Strike-itsSpread+1.0;
	double AsOfDate = GetAsOfDate().GetJulian();

	size_t i, nbFlows = numCPIDates->size();
	ARM_GP_VectorPtr result( new ARM_GP_Vector( nbStates ) );
	ARM_GP_VectorPtr Caplets;
	ARM_GP_Vector::iterator iterResult, iterCaplets; 

	for( i = 0 ; i< nbFlows; i++ )
	{
		Caplets = YoYCaplet( irCurveName, infCurveName, evalTime, (*numCPIDates)[i]-(*denomCPIDates)[i], (*numCPIDates)[i]-AsOfDate, (*numCPIDCFDates)[i]-AsOfDate, strike, CapFloor, states );
		iterCaplets = Caplets->begin();
		dayCountCorrection = NumInterestTerms->Elt(i);

		for( iterResult = result->begin() ; iterResult != result->end() ; iterResult++, iterCaplets++ )
			(*iterResult) += (*iterCaplets) * dayCountCorrection;
	}

	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: OATCapFloor
///	Returns: ARM_GP_VectorPtr
///	Action : returns OATCap at a certain date and in a certain state of the world
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_InflationEquityMod::OATCapFloor( const string& irCurveName,	const string& infCurveName, double evalTime,double Strike,double FloatMargin, int CapFloor, const ARM_DateStripPtr& numDateStrip,const ARM_DateStripPtr& denomDateStrip,double itsSpread,	const ARM_PricingStatesPtr& states) const
{
	size_t nbStates = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;

	/// This vector contains the num CPIs fixingDates (and not really the FlowStartDates!)
	ARM_GP_Vector * numCPIDCFDates = numDateStrip->GetFlowStartDates();
	ARM_GP_Vector * numCPIDates = numDateStrip->GetResetDates();
	ARM_GP_Vector * denomCPIDates = denomDateStrip->GetResetDates();
	ARM_GP_Vector * paymentDates = numDateStrip->GetFlowEndDates();
	ARM_GP_Vector * NumInterestTerms = numDateStrip->GetInterestTerms();
	double dayCountCorrection;
	double strike = Strike-itsSpread;
	double AsOfDate = GetAsOfDate().GetJulian();

	size_t i, nbFlows = numCPIDates->size();
	ARM_GP_VectorPtr result( new ARM_GP_Vector( nbStates ) );
	ARM_GP_VectorPtr Caplets;
	ARM_GP_Vector::iterator iterResult, iterCaplets;
	double denomCPIDate = (*denomCPIDates)[0];

	for( i = 0 ; i< nbFlows; i++ )
	{
		Caplets = YoYCaplet( irCurveName, infCurveName, evalTime, (*numCPIDates)[i]-denomCPIDate, (*numCPIDates)[i]-AsOfDate, (*numCPIDCFDates)[i]-AsOfDate, strike, CapFloor, states );
		iterCaplets = Caplets->begin();
		dayCountCorrection = NumInterestTerms->Elt(i);

		for( iterResult = result->begin() ; iterResult != result->end() ; iterResult++, iterCaplets++ )
			(*iterResult) += (*iterCaplets) * dayCountCorrection;
	}

	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: integratedVolForBSYoY
///	Returns: double
///	Action : BS Vol to price YoYCaplet
////////////////////////////////////////////////////
double ARM_InflationEquityMod::integratedVolForBSYoY( double FromTime, double ToTime, double tenor ) const
{
	double alpha = ((const ARM_ModelParamsInflationEquity* const) GetModelParams())->getMultiplier();
	double rho = getCorrelWithBonds();
	const ARM_ModelParam * vol = &(GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility));
	const ARM_PricingModelIR* itsIRModel = getIRModel();

	double result = 0;
	result = itsIRModel->IntegratedBondSquaredVol( FromTime, ToTime, ToTime );
	result += itsIRModel->IntegratedBondSquaredVol( FromTime, ToTime-tenor, ToTime-tenor);
	result -= 2*itsIRModel->IntegratedBondCovariance( FromTime, ToTime-tenor, ToTime-tenor, ToTime );
	result *= alpha*alpha;
	result -= 2*alpha*rho*itsIRModel->VolatilityScalarProduct( ToTime-tenor, ToTime, ToTime, *vol );
	result += ((const ARM_ModelParamsInflationEquity* const) GetModelParams())->StateLocalVariance(ToTime-tenor,ToTime);

	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_InflationEquityMod
///	Routine: YoYCap
///	Returns: ARM_GP_VectorPtr
///	Action : returns YoYCap at a certain date and in a certain state of the world
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_InflationEquityMod::YoYCaplet( const string& irCurveName, const string& infCurveName, double EvalTime, double tenor, double CPITime, double maturityTime, double Strike, int CapFloor, const ARM_PricingStatesPtr& states ) const
{
	double payTime = GetAsOfDate().GetJulian();
	payTime = ARM_Date(maturityTime + payTime).GoodBusinessDay(1,const_cast<char*>(irCurveName.c_str())).GetJulian() - payTime;
	ARM_GP_VectorPtr result = ForwardCPIRatio( infCurveName, EvalTime, tenor, CPITime, maturityTime, states );
	ARM_GP_VectorPtr DFs = DiscountFactor( infCurveName, EvalTime, payTime,states); 
	ARM_GP_Vector::iterator iter = result->begin(), iterEnd = result->end(), iterDFs = DFs->begin();
	double vol2 = integratedVolForBSYoY( EvalTime, maturityTime, tenor );
	double vol = sqrt(vol2);

	for( ; iter != iterEnd ; ++iter, ++iterDFs)
		(*iter) = BlackSholes_Formula( (*iter), vol, (*iterDFs), Strike, CapFloor );

	return result;
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

