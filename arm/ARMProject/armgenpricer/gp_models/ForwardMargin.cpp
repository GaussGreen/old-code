/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file ForwardMargin.cpp
 *
 *  \brief
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date April 2004
 */

#include "gpbase/removeidentifiedwarning.h"

/// gpmodels
#include "gpmodels/forwardmargin.h"
#include "gpmodels/multiassets.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/gpvector.h"

/// gpinfra
#include "gpinfra/pricingmodeltype.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelnamemap.h"

/// gpcalib
#include "gpcalib/calibmethod.h"

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_ForwardMargin
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ForwardMargin::ARM_ForwardMargin(
	const ARM_ZeroCurvePtr& shiftZcCurve, 
    ARM_PricingModel* refModel, 
	bool refModelDump)
:   
    ARM_PricingModelIR(shiftZcCurve),
    itsRefModel(refModel),
    itsRefModelDump(refModelDump)
{}


////////////////////////////////////////////////////
///	Class  : 
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ForwardMargin::ARM_ForwardMargin(const ARM_ForwardMargin& rhs)
:
    ARM_PricingModelIR(rhs),
    itsRefModel( rhs.itsRefModel ),
    itsRefModelDump ( rhs.itsRefModelDump)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ForwardMargin
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ForwardMargin::~ARM_ForwardMargin()
{}


////////////////////////////////////////////////////
///	Class  : ARM_ForwardMargin
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_ForwardMargin& ARM_ForwardMargin::operator=(const ARM_ForwardMargin& rhs)
{
	if(this != &rhs)
	{
		ARM_PricingModelIR::operator=(rhs);
        itsRefModel     = rhs.itsRefModel;
        itsRefModelDump = rhs.itsRefModelDump;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_ForwardMargin
///	Routine: DiscountFactor
///	Returns: a vector of Zc(t,T)
///	Action : Compute the forward zero-coupon
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ForwardMargin::DiscountFactor(
	const string& curveName, 
	double evalTime, 
	double maturityTime, 
	const ARM_PricingStatesPtr& states) const
{
    /// Get the deterministic margin between shifted & reference spot yield curves
	double marginRatio = ComputeMarginRatio(evalTime,maturityTime);

    /// Compute df using reference model
    ARM_VectorPtr df= itsRefModel->DiscountFactor( itsRefModel->GetModelName(), evalTime, maturityTime, states );

    for(size_t i=0; i<df->size(); ++i )
        (*df)[i] *= marginRatio;

    return df;
}



////////////////////////////////////////////////////
///	Class   : ARM_ForwardMargin
///	Routines: toString
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
string ARM_ForwardMargin::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    if(itsRefModel)
    {
        if(itsRefModelDump)
            /// The reference model is dumped
            os << indent << "Linked Model : \n" << itsRefModel->toString(indent+nextIndent,nextIndent);
        else
            /// The name of the reference model is only dumped
            os << indent << "Linked Model : \n" << itsRefModel->GetModelName();
    }

    return os.str();
}

////////////////////////////////////////////////////
///	Class   : ARM_ForwardMargin
///	Routine : ValidateCalibMethod
///	Returns :void
///	Action  : call DefaultValidateWithModel
////////////////////////////////////////////////////
void ARM_ForwardMargin::ValidateCalibMethod(ARM_CalibMethod& calibMethod)
{
	calibMethod.DefaultValidateWithModel(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_ForwardMargin
///	Routine: GetType
///	Returns: int
///	Action : tells the type of the model
////////////////////////////////////////////////////

int ARM_ForwardMargin::GetType() const
{
	return MT_NON_STOCHASTIC_MODEL;
}


////////////////////////////////////////////////////
///	Class  : ARM_ForwardMargin
///	Routine: ProcessPaidPayoffs,
///          ProcessUnPaidPayoffs,
///          Init
///	Returns: 
///	Action : Call the reference model but
///          take care of the Discount Functor if
///          it needs to be changed
////////////////////////////////////////////////////
void ARM_ForwardMargin::ProcessPaidPayoffs(const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, const ARM_PricingStatesPtr& states ) const
{
	/// change the discount functor
    ARM_ZeroCurveFunctor* oldDiscountFunctor=itsRefModel->GetDiscountFunctor();
    itsRefModel->SetDiscountFunctor( GetDiscountFunctor() );
    itsRefModel->ProcessPaidPayoffs( payModelName, payoffs, evalTime, states );
    /// restore it
	itsRefModel->SetDiscountFunctor(oldDiscountFunctor);
}

void ARM_ForwardMargin::ProcessUnPaidPayoffs(const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, const ARM_PricingStatesPtr& states ) const
{
	/// change the discount functor
    ARM_ZeroCurveFunctor* oldDiscountFunctor=itsRefModel->GetDiscountFunctor();
    itsRefModel->SetDiscountFunctor( GetDiscountFunctor() );
    itsRefModel->ProcessUnPaidPayoffs( payModelName, payoffs, evalTime, states );
    /// restore it
	itsRefModel->SetDiscountFunctor(oldDiscountFunctor);
}

ARM_PricingStatesPtr ARM_ForwardMargin::Init( const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos )
{
	/// change the discount functor
    ARM_ZeroCurveFunctor* oldDiscountFunctor=itsRefModel->GetDiscountFunctor();
    itsRefModel->SetDiscountFunctor( GetDiscountFunctor() );
    ARM_PricingStatesPtr result = itsRefModel->Init( payModelName, timeInfos );
    /// restore it
	itsRefModel->SetDiscountFunctor(oldDiscountFunctor);
	return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_ForwardMargin
///	Routine: all the calibration method
///	Returns: 
///	Action : delegates to the reference model after checking the existence
///			of a reference model
////////////////////////////////////////////////////
void ARM_ForwardMargin::Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter)
{   
	CheckRefModel();
	itsRefModel->Re_InitialiseCalibParams(modelFitter);
}

void ARM_ForwardMargin::PreProcessing(ARM_ModelFitter& modelFitter)
{
	CheckRefModel();
	itsRefModel->PreProcessing(modelFitter);
}


void ARM_ForwardMargin::PostProcessing(const ARM_ModelFitter& modelFitter)
{
	CheckRefModel();
	itsRefModel->PostProcessing(modelFitter);
}

void ARM_ForwardMargin::AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter)
{
	CheckRefModel();
	itsRefModel->AdviseCurrentCalibSecIndex(index,modelFitter);
}


void ARM_ForwardMargin::AdviseCurrentCalib(ARM_ModelFitter& modelFitter)
{
	CheckRefModel();
	itsRefModel->AdviseCurrentCalib(modelFitter);
}

void ARM_ForwardMargin::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, 
							  ARM_ModelParam* modelParam, 
							  size_t factorNb )
{
	CheckRefModel();
	itsRefModel->AdviseBreakPointTimes(portfolio,modelParam,factorNb );
}


ARM_PricingStatesPtr ARM_ForwardMargin::FirstPricingStates( size_t bucketSize ) const
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "FirstPricingStates : unimplemented function !"); return ARM_PricingStatesPtr(NULL); }

std::vector<double>& ARM_ForwardMargin::ComputeModelTimes( const ARM_TimeInfoPtrVector& timeInfos )
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "ComputeModelTimes : unimplemented function !"); return NULL; }

bool ARM_ForwardMargin::SupportBackwardInduction() const
{   
	CheckRefModel();
	return itsRefModel->SupportBackwardInduction();
}

bool ARM_ForwardMargin::SupportForwardInduction() const
{
	CheckRefModel();
	return itsRefModel->SupportForwardInduction();
}

bool ARM_ForwardMargin::SupportAnalyticMarginal() const
{
	CheckRefModel();
	return itsRefModel->SupportAnalyticMarginal();
}


void ARM_ForwardMargin::IntegratedLocalDrifts(const std::vector<double>& timeSteps,ARM_GP_MatrixPtr& relativeDrifts,ARM_GP_MatrixPtr& absoluteDrifts) const
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "IntegratedLocalDrifts : unimplemented function !"); }

void ARM_ForwardMargin::ModelStateLocalVariances( const std::vector<double>& timeSteps,ARM_MatrixVector& localVariances ) const
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "LocalVariances : unimplemented function !"); }

void ARM_ForwardMargin::NumMethodStateLocalVariances( const std::vector<double>& timeSteps,ARM_MatrixVector& localVariances ) const
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "LocalVariances : unimplemented function !"); }


bool ARM_ForwardMargin::NeedsToCholeskyDecomposeFactors( ) const
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "NeedsToCholeskyDecomposeFactors : unimplemented function !"); return false; }

double ARM_ForwardMargin::VarianceToTime(double var,double minTime,double maxTime) const
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "VarianceToTime : unimplemented function !"); return 0.0; }

ARM_VectorPtr  ARM_ForwardMargin::VanillaSpreadOptionLet(const string& curveName,
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
		"VanillaSpreadOption : unimplemented function for ARM_ForwardMargin Model!");
}

void ARM_ForwardMargin::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{  
		/// change the discount functor
    ARM_ZeroCurveFunctor* oldDiscountFunctor=itsRefModel->GetDiscountFunctor();
    itsRefModel->SetDiscountFunctor( GetDiscountFunctor() );
    itsRefModel->MCModelStatesFromToNextTime( states, timeIndex );
    /// restore it
	itsRefModel->SetDiscountFunctor(oldDiscountFunctor);
}

void ARM_ForwardMargin::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "TreeStatesToModelStates : unimplemented function !"); }



////////////////////////////////////////////////////
///	Class   : ARM_ForwardMargin
///	Routines: SetNumeraire
///	Returns : void
///	Action  : sets the numeraire to the pricing model
////////////////////////////////////////////////////
void ARM_ForwardMargin::SetNumeraire(const ARM_NumerairePtr& numerairePtr)
{
	itsRefModel->SetNumeraire( numerairePtr );
}

////////////////////////////////////////////////////
///	Class   : ARM_ForwardMargin
///	Routines: SetNumMethod
///	Returns : void
///	Action  : sets the numerical method to the pricing model
////////////////////////////////////////////////////
void ARM_ForwardMargin::SetNumMethod(const ARM_NumMethodPtr& numMethodPtr)
{
	itsRefModel->SetNumMethod( numMethodPtr );
}

////////////////////////////////////////////////////
///	Class   : ARM_ForwardMargin
///	Routines: Update the links on the reference model
///	Returns : void
///	Action  : change the ref model
////////////////////////////////////////////////////
void ARM_ForwardMargin::UpdateLinks( const ARM_MultiAssetsModel& multiAssetsModel )
{
	const ARM_ModelNameMap* modelMap = multiAssetsModel.GetModelMap();

	/// find linked models
	ARM_IntVector itsOtherModels = (*modelMap)[GetModelName()]->OtherModelRefNb();
	// check that there is one and only one...
	if( itsOtherModels.size() != 1 )
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ForwardMargin: wrong number of otherModels" );

	itsRefModel = &*(*modelMap)[itsOtherModels[0] ]->Model();
}

////////////////////////////////////////////////////
///	Class   : ARM_ForwardMargin
///	Routines: Update the links on the reference model
///	Returns : void
///	Action  : change the ref model
////////////////////////////////////////////////////
const ARM_ModelParams* const ARM_ForwardMargin::GetModelParams() const
{
	if( !itsRefModel )
		ARM_THROW( ERR_INVALID_ARGUMENT, " there is no ref model, hence cannot get model params" );
	return itsRefModel->GetModelParams();
}



////////////////////////////////////////////////////
///	Class   : ARM_ForwardMargin
///	Routines: Update the links on the reference model
///	Returns : void
///	Action  : change the ref model
////////////////////////////////////////////////////
ARM_ModelParams* ARM_ForwardMargin::GetModelParams()
{
	if( !itsRefModel )
		ARM_THROW( ERR_INVALID_ARGUMENT, " there is no ref model, hence cannot get model params" );
	return itsRefModel->GetModelParams();
}


////////////////////////////////////////////////////
///	Class   : ARM_ForwardMargin
///	Routines: LocalDiscounts
///	Returns : void
///	Action  : Computes the LocalDiscounts
////////////////////////////////////////////////////

ARM_VectorPtr ARM_ForwardMargin::LocalDiscounts( size_t timeIdx, double dt, 
	const ARM_PricingStatesPtr& states) const
{
	CheckRefModel();
	return itsRefModel->LocalDiscounts( timeIdx, dt, states);
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

