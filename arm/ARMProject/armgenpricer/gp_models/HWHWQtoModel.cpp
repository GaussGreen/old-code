/*!
 *
 * Copyright (c) CM CIB January 2005 Paris
 *
 *	\version 1.0
 *	\date February 2007
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/HWHWQtoModel.h"

/// gpmodels
#include "gpmodels/HW1F.h"
#include "gpmodels/ModelParamsHW1F.h"
#include "gpmodels/LN_Fx.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/forwardmargin.h"
#include "gpmodels/typedef.h"
#include "gpmodels/local_normal_model.h"

/// gpbase
#include "gpbase/checkarg.h"
#include "gpbase/stringmanip.h"
#include "gpbase/vectormanip.h"
#include "gpbase/eventviewerfwd.h"

/// gpinfra
#include "gpinfra/modelnamemap.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/discretisationscheme.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/pricingcontext.h"

/// gpnummethods
#include "gpnummethods/pdemethod.h"
#include "gpnummethods/treebase.h"
#include "gpnummethods/sampler.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/tree1d.h"
#include "gpnummethods/mcmethod.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/meanrevertingsampler.h"
#include "gpnummethods/normalcentredsampler.h"

/// gpnumlib
#include "gpnumlib/random.h"

#include <algorithm>
CC_USING_NS (std,find)

CC_BEGIN_NAMESPACE( ARM )

const size_t ARM_HWHWQtoModel::NbLocalModels = 1;

const double MAX_LNVOL_FX = 10.0; // 1000%

////////////////////////////////////////////////////
///	Class   : ARM_HWHWQtoModel
///	Routines: Constructor
///	Returns :
///	Action  : builds the object
///	model types contains the integer on the various model
///	using the correlation matrix
////////////////////////////////////////////////////
ARM_HWHWQtoModel::ARM_HWHWQtoModel(	const ARM_ModelNameMap&	modelNameMap, 
	const ARM_CurveMatrix& correlMatrix,
	bool fxFlag)
:
	ARM_MultiAssetsMeanReverting(&modelNameMap,&correlMatrix),
	itsForLocalModel(false)
{
	/// validate and initialize
	Validate(fxFlag);
}

////////////////////////////////////////////////////
///	Class   : ARM_HWHWQtoModel
///	Routine : InitSubModels
///	Returns : void
///	Action  : intialises the sub models that are linked to a model
////////////////////////////////////////////////////
void ARM_HWHWQtoModel::Validate( bool fxFlag )
{
	/// number of models
	ARM_ModelNameMap* modelNameMap	= GetModelMap();
	if(modelNameMap->size() < NbModels)
    {
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : not enough models to build HW x HW x Quanto");
    }
    else if(modelNameMap->size() > NbModels)
    {
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : too much models to build HW x HW x Quanto");
    }

	/// correlation
	if( !GetCorrelMatrix() )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : please provide a correlation");

	if(	GetCorrelMatrix()->rows() != 3 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : correlation should be 3x3!");	

	/// check model types
	if( !dynamic_cast<ARM_HullWhite1F*>( &*(*modelNameMap)[DomModel]->Model() ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : dom model has to be HW1F IR Model");

	if( !dynamic_cast<ARM_HullWhite1F*>( &*(*modelNameMap)[ForModel]->Model() ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : for model has to be HW1F IR Model");

	ARM_LN_Fx* lnModel = dynamic_cast<ARM_LN_Fx*>(	&*(*modelNameMap)[FxModel]->Model() );
	if( !lnModel)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : fx model has to be LN FX Model");
	
	if( !dynamic_cast<ARM_ForwardMargin*>(&*(*modelNameMap)[DomBasisModel]->Model()) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : dom basis model has to be a Forward Margin Model");

	if( !dynamic_cast<ARM_ForwardMargin*>(&*(*modelNameMap)[ForBasisModel]->Model()) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : for basis model has to be a Forward Margin Model");

	if(modelNameMap->size() + NbLocalModels >= NbModels + 1)
    {
	    if( !dynamic_cast<ARM_Local_Normal_Model*>(&*(*modelNameMap)[ForLocalModel]->Model()) )
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : local model for Foreign has to be a Local Normal Model");
        itsForLocalModel = true;
    }
    
    /// FX currency checks
	ARM_ModelParams* modelParams= (*modelNameMap)[FxModel]->Model()->GetModelParams();
	ARM_ModelParams_Fx* fxModelParams = dynamic_cast<ARM_ModelParams_Fx*>(modelParams);
	string domCcy((*modelNameMap)[DomModel]->Model()->GetZeroCurve()->GetCurrencyUnit()->GetCcyName());
	string fxDomCcy(fxModelParams->GetDomCurve()->GetCurrencyUnit()->GetCcyName());
	if( fxDomCcy != domCcy)
	{
		CC_Ostringstream os;
		os << ARM_USERNAME << " : The cross FXmodel has domestic curve == " << fxDomCcy << " while the IR domestic model ccy is == " << domCcy;
		ARM_THROW( ERR_INVALID_ARGUMENT, os.str() );
	}
	string forCcy((*modelNameMap)[ForModel]->Model()->GetZeroCurve()->GetCurrencyUnit()->GetCcyName());
	string fxForCcy(fxModelParams->GetForCurve()->GetCurrencyUnit()->GetCcyName());
	if( fxForCcy != forCcy)
	{
		CC_Ostringstream os;
		os << ARM_USERNAME << " : The cross FXmodel has foreign curve == " << fxForCcy << " while the IR foreign model ccy is == " << forCcy;
		ARM_THROW( ERR_INVALID_ARGUMENT, os.str() );
	}

	///Initialize the FXModel Correlation.
	ARM_EqFxBase* fxModel = dynamic_cast<ARM_EqFxBase*>(	&*(*modelNameMap)[FxModel]->Model() );
	fxModel->SetCorrelMatrix( *GetCorrelMatrix() );
	// and set factor count for FX
	if (fxFlag)
		lnModel->SetFactorCount(1);
	else
		lnModel->SetFactorCount(0);
	
	
}


////////////////////////////////////////////////////
///	Class   : ARM_HWHWQtoModel
///	Routines: Copy constructor
///	Returns :
///	Action  : constructor
ARM_HWHWQtoModel::ARM_HWHWQtoModel(const ARM_HWHWQtoModel& rhs)
:	ARM_MultiAssetsMeanReverting(rhs),
	itsForLocalModel(rhs.itsForLocalModel)
{
    itsDriftCorrections = rhs.itsDriftCorrections != ARM_GP_MatrixPtr(NULL)
                          ? ARM_GP_MatrixPtr( static_cast< ARM_GP_Matrix* >(rhs.itsDriftCorrections->Clone()) )
                          : ARM_GP_MatrixPtr(NULL);
}


////////////////////////////////////////////////////
///	Class   : ARM_HWHWQtoModel
///	Routine : SupportAnalyticMarginal
///	Returns : void
///	Action  : decide if over sampling between model time
///           step is needed
////////////////////////////////////////////////////
bool ARM_HWHWQtoModel::SupportAnalyticMarginal() const
{
    const ARM_NumerairePtr numeraire = GetNumeraire();
    return  numeraire == ARM_NumerairePtr(NULL) ||
            numeraire->GetType() == ARM_Numeraire::RollingEvent;
}

////////////////////////////////////////////////////
///	Class   : ARM_HWHWQtoModel
///	Routine : NeedMCIntegProcess
///	Returns : ARM_BoolVector
///	Action  : to say how MC method give simulated processes
////////////////////////////////////////////////////
ARM_BoolVector ARM_HWHWQtoModel::NeedMCIntegProcess() const
{
    ARM_BoolVector SimulProcess(FactorCount(),true);
	if (FxModel<FactorCount())
		SimulProcess[FxModel] = false;

	ARM_NumerairePtr numeraire = GetNumeraire();
	if(numeraire->GetType() == ARM_Numeraire::Cash)
	{
		/// Domestic & foreign processes are incremental because
		/// the discretisation scheme is Euler
		SimulProcess[DomModel] = false;
		SimulProcess[ForModel] = false;
	}

    return SimulProcess;
}

///////////////////////////////////////////////////
///	Class   : ARM_HWHWQtoModel
///	Routine : Init
///	Returns : ARM_PricingStatesPtr
///	Action  : intialise the model
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_HWHWQtoModel::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
	/// Validate that the payment currency is the domestic model!
    const ARM_ModelNameMap* const modelMap = GetModelMap();
    const string& domBasisModelName = (*modelMap)[DomBasisModel]->Model()->GetModelName();
    const string& domModelName      = (*modelMap)[DomModel]->Model()->GetModelName();
	if( payModelName != domBasisModelName && payModelName != domModelName)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": HW HW Qto model allow only paimenent in " + domBasisModelName + " or " + domModelName);

	ARM_NumerairePtr numeraire = GetNumeraire();
    ARM_NumMethodPtr numMethod = GetNumMethod();
    if( numMethod != ARM_NumMethodPtr(NULL) && 
        numeraire != ARM_NumerairePtr(NULL) &&
        numMethod->GetPricingDirection() != ARM_NumMethod::GP_AMBIGUOUS &&
        !(numMethod->GetPricingDirection() == ARM_NumMethod::GP_BCKWDLOOKING &&
            numeraire->GetType() == ARM_Numeraire::Cash ||
			numeraire->GetType() == ARM_Numeraire::TerminalZc) &&
        !((numMethod->GetPricingDirection() == ARM_NumMethod::GP_FWDLOOKING ||
			numMethod->GetPricingDirection() == ARM_NumMethod::GP_FWDBCKWDLOOKING) &&
            (numeraire->GetType() == ARM_Numeraire::RollingEvent ||
             numeraire->GetType() == ARM_Numeraire::Cash ||
             numeraire->GetType() == ARM_Numeraire::RollingCash)) )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
        ": only Cash or TerminalZc or RollingEvent numeraire supported by HW HW Qto model");

	if( numMethod != ARM_NumMethodPtr(NULL) && 
        numeraire != ARM_NumerairePtr(NULL) &&
		(numMethod->GetPricingDirection() == ARM_NumMethod::GP_FWDLOOKING 
		|| numMethod->GetPricingDirection() == ARM_NumMethod::GP_FWDBCKWDLOOKING) &&
        numeraire->GetType() == ARM_Numeraire::Cash)
	   dynamic_cast<ARM_NumeraireCash&>((*numeraire)).SetRenormalisation(false);

	/// Delegate to common code
	return ARM_MultiAssetsModel::Init( payModelName, timeInfos );
}


////////////////////////////////////////////////////
///	Class   : ARM_HWHWQtoModel
///	Routine : ComputeIntegratedFwdFxVCV
///	Returns : void
///	Action  : computes the integrated VCV matrixes for
///           domestic & foreign state variables &
///           forward Fx process
////////////////////////////////////////////////////
void ARM_HWHWQtoModel::ComputeIntegratedFwdFxVCV(double step, double nextStep,
        const ARM_ModelParamsHW1FStd* const domModelParams,
        const ARM_ModelParamsHW1FStd* const forModelParams,
        const ARM_ModelParamsHW1FStd* const fxModelParams,
        const ARM_CurveModelParam& fxVol,
        const ARM_GP_Matrix& correlMatrix,
        ARM_GP_Matrix& variances) const
{
    double varFx,varZcDom,varZcFor;
    double covarZcDomZcFor,covarZcDomFx,covarZcForFx;
    double covarXDomZcDom,covarXDomZcFor,covarXDomFx;
    double covarXForZcDom,covarXForZcFor,covarXForFx;

    /// Domestic variances
	variances(DomModel,DomModel)   = domModelParams->StateLocalVariance(step,nextStep,nextStep);

    /// Foreign variances
	variances(ForModel,ForModel)   = forModelParams->StateLocalVariance(step,nextStep,nextStep);

    /// Domestic & Foreign covariances
	variances(DomModel,ForModel) = correlMatrix(DomModel,ForModel) * ARM_ModelParamsHW1F::HW1FStateCovariance(domModelParams, forModelParams,step,nextStep,nextStep);
    variances(ForModel,DomModel) = variances(DomModel,ForModel);

	if (FxModel<FactorCount())
	{
		/// Forex variances
		varFx       = fxModelParams->StateLocalVariance(step,nextStep,nextStep);
		varZcDom    = ARM_ModelParamsHW1F::HW1FZcCovariance( domModelParams, domModelParams, step,nextStep,nextStep );
		varZcFor    = ARM_ModelParamsHW1F::HW1FZcCovariance( forModelParams, forModelParams, step,nextStep,nextStep );

		covarZcDomZcFor = ARM_ModelParamsHW1F::HW1FZcCovariance( domModelParams, forModelParams, step,nextStep,nextStep );
		covarZcDomFx    = ARM_ModelParamsHW1F::HW1FEqFxZcCovariance( fxVol, domModelParams, step,nextStep,nextStep );
		covarZcForFx    = ARM_ModelParamsHW1F::HW1FEqFxZcCovariance( fxVol, forModelParams, step,nextStep,nextStep );

		variances(FxModel,FxModel) =
			varFx + varZcDom + varZcFor
			- 2*correlMatrix(DomModel,ForModel)*covarZcDomZcFor
			- 2*correlMatrix(DomModel,FxModel)*covarZcDomFx
			+ 2*correlMatrix(ForModel,FxModel)*covarZcForFx;

		double varLimit = K_NEW_DOUBLE_TOL*K_NEW_DOUBLE_TOL;
		if(variances(FxModel,FxModel)<varLimit)
			variances(FxModel,FxModel) = varLimit; // numerical noise

		/// Domestic & Forex covariances
		covarXDomZcDom = ARM_ModelParamsHW1F::HW1FStateZcCovariance( domModelParams, domModelParams, step,nextStep,nextStep,nextStep );
		covarXDomZcFor = ARM_ModelParamsHW1F::HW1FStateZcCovariance( domModelParams, forModelParams, step,nextStep,nextStep,nextStep );
		covarXDomFx    = ARM_ModelParamsHW1F::HW1FStateCovariance( domModelParams, fxModelParams, step,nextStep,nextStep );

		variances(DomModel,FxModel) =
			- covarXDomZcDom
			+ correlMatrix(DomModel,ForModel)*covarXDomZcFor
			+ correlMatrix(DomModel,FxModel)*covarXDomFx;

		variances(FxModel,DomModel) = variances(DomModel,FxModel);

		/// Foreign & Forex covariances
		covarXForZcDom = ARM_ModelParamsHW1F::HW1FStateZcCovariance( forModelParams, domModelParams, step,nextStep,nextStep,nextStep );
		covarXForZcFor = ARM_ModelParamsHW1F::HW1FStateZcCovariance( forModelParams, forModelParams, step,nextStep,nextStep,nextStep );
		covarXForFx    = ARM_ModelParamsHW1F::HW1FStateCovariance( forModelParams, fxModelParams, step,nextStep,nextStep );

		variances(ForModel,FxModel) =
			covarXForZcFor
			- correlMatrix(DomModel,ForModel)*covarXForZcDom
			+ correlMatrix(ForModel,FxModel)*covarXForFx;

		variances(FxModel,ForModel) = variances(ForModel,FxModel);
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_HWHWQtoModel
///	Routine : NumMethodStateLocalGlobalVariances
///	Returns : void
///	Action  : computes the integrated covariances
////////////////////////////////////////////////////

void ARM_HWHWQtoModel::NumMethodStateLocalGlobalVariances( const std::vector<double>& timeSteps,
	ARM_MatrixVector& localVariances,
	ARM_MatrixVector& variances ) const
{
	size_t factorNb	= FactorCount();
	const ARM_CurveMatrix* const correlCurveMatrix  = GetCorrelMatrix();
	ARM_GP_Matrix correlMatrix ;
	const ARM_ModelNameMap* const modelMap = GetModelMap();
	size_t nbSteps	= timeSteps.size();
    double step		= timeSteps[0],nextStep;

#if defined(__GP_STRICT_VALIDATION)
	if( localVariances.size()!= 0 ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "localVariances.size() != 0" );
	if( variances.size()!= 0 ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "variances.size() != 0" );
#endif


	localVariances.resize(nbSteps-1);
	variances.resize(nbSteps);

    ARM_MatrixVector fxLocalVariances(nbSteps-1);

	variances[0]=new ARM_GP_TriangularMatrix(factorNb,0.0);

    ARM_PricingModelPtr domModel    = (*modelMap)[DomModel]->Model();
    ARM_PricingModelPtr forModel    = (*modelMap)[ForModel]->Model();
    ARM_PricingModelPtr fxModel     = (*modelMap)[FxModel]->Model();
	const ARM_ModelParamsHW1FStd* const domModelParams  = static_cast<const ARM_ModelParamsHW1FStd* const>( domModel->GetModelParams() );
	const ARM_ModelParamsHW1FStd* const forModelParams  = static_cast<const ARM_ModelParamsHW1FStd* const>( forModel->GetModelParams() );
	const ARM_ModelParamsHW1FStd* const fxModelParams   = static_cast<const ARM_ModelParamsHW1FStd* const>( fxModel->GetModelParams() );
	const ARM_CurveModelParam& fxVol                    = dynamic_cast<const ARM_CurveModelParam&>( fxModelParams->GetModelParam(fxModelParams->GetVolatilityType()) );

    for(size_t i=0;i<nbSteps-1;++i)
	{
		correlMatrix = correlCurveMatrix->Interpolate(step);

		nextStep			= timeSteps[i+1];

		localVariances[i]	= new ARM_GP_TriangularMatrix(factorNb,0.0);

		ComputeIntegratedFwdFxVCV(step,nextStep,
			domModelParams,forModelParams,fxModelParams,fxVol,correlMatrix,
		    *(localVariances[i]));

		variances[i+1]		= new ARM_GP_TriangularMatrix(factorNb,0.0);

		ComputeIntegratedFwdFxVCV(0.0,nextStep,
			domModelParams,forModelParams,fxModelParams,fxVol,correlMatrix,
			*(variances[i+1]));
        
		if (FxModel<FactorCount())
			fxLocalVariances[i] = static_cast< ARM_GP_Matrix* >(localVariances[i]->Clone());

		step=nextStep;
	}

	if (FxModel<FactorCount())
		(*modelMap)[FxModel]->Model()->SetModelStateLocalVars(fxLocalVariances);
}

////////////////////////////////////////////////////
///	Class   : ARM_HWHWQtoModel
///	Routine : NumMethodStateGlobalVariances
///	Returns : void
///	Action  : computes the integrated covariances
////////////////////////////////////////////////////

void ARM_HWHWQtoModel::NumMethodStateGlobalVariances( const std::vector<double>& timeSteps,
	ARM_MatrixVector& variances ) const
{
	ARM_MatrixVector localVariances;

	NumMethodStateLocalGlobalVariances(timeSteps,localVariances,variances );
}

////////////////////////////////////////////////////
///	Class   : ARM_HWHWQtoModel
///	Routine : NumMethodStateLocalVariances
///	Returns : void
///	Action  : computes the integrated covariances
////////////////////////////////////////////////////

void ARM_HWHWQtoModel::NumMethodStateLocalVariances( const std::vector<double>& timeSteps,
	ARM_MatrixVector& localVariances ) const
{
	ARM_MatrixVector variances;

	NumMethodStateLocalGlobalVariances(timeSteps,localVariances,variances );
}

////////////////////////////////////////////////////
///	Class   : ARM_HWHWQtoModel
///	Routine : BackwardLookingInit
///	Returns : ARM_PricingStatesPtr
///	Action  : Initialisation for backward looking
///           numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_HWHWQtoModel::BackwardLookingInit(ARM_NumMethodPtr& numMethod, double firstInductTime)
{
	ARM_PricingStatesPtr initStates(NULL);
    const ARM_ModelNameMap& modelMap = * GetModelMap();

    ARM_TreeBase* treebase = dynamic_cast<ARM_TreeBase*>( &*numMethod );
	ARM_PDEMethod* pdemethod = dynamic_cast<ARM_PDEMethod*>( &*numMethod );
	
	//tree case
	if( treebase )
	{
		initStates = numMethod->Init(*this,firstInductTime);
	}
	//pde case
	else if(pdemethod)
	{
		initStates = numMethod->Init(*this,firstInductTime);
	}
	else
	{
	    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not a ARM_TreeBase or a ARM_PDEmethod numerical method!");
	}
    return initStates;
}

////////////////////////////////////////////////////
///	Class  : ARM_HWHWQtoModel
///	Routine: TreeStatesToModelStates
///	Returns: void
///	Action : Convert the num method states to model states
////////////////////////////////////////////////////
void ARM_HWHWQtoModel::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{   
#if defined(__GP_STRICT_VALIDATION)
	if(timeIndex<0)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index is negative!");
	if( timeIndex > GetNumMethod()->GetTimeSteps()->size()-1 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index bigger than max size!");
#endif

	const ARM_NumMethodPtr& numMethod = GetNumMethod();
	ARM_TreeBase* treebase = dynamic_cast<ARM_TreeBase*>( &*numMethod );
	ARM_PDEMethod* pdemethod = dynamic_cast<ARM_PDEMethod*>( &*numMethod );
	
	//tree case
	if( treebase )
	{
		const ARM_GP_MatrixPtr& modelStates = states->GetNumMethodStates();
		states->SetModelStates(modelStates);
	}
	//other cases imply an error
	else
	{
	    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not a ARM_TreeBase numerical method!");
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_HWHWQtoModel
///	Routine: ComputeNumeraireTimes
///	Returns: 
///	Action : compute the corresonding numeraire times
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWHWQtoModel::ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const
{
	size_t i,nbEvents = timeInfos.size();

    ARM_VectorPtr numTimes(NULL);

	switch( GetNumeraire()->GetType() )
	{
        case ARM_Numeraire::Cash:
		case ARM_Numeraire::TerminalZc:
        case ARM_Numeraire::RollingCash:
            break;

	    case ARM_Numeraire::RollingEvent:
		{
            std::vector<double>& times = new std::vector<double>(nbEvents);
            for(i=0;i<nbEvents;++i)
			    (*times)[i] = timeInfos[i]->GetEventTime();
            numTimes = ARM_GP_VectorPtr(times);
        }
        break;

        default :
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index bigger than max size!");
    }

    return numTimes;
}


////////////////////////////////////////////////////
///	Class   : ARM_HWHWQtoModel
///	Routine : ForwardLookingInit
///	Returns : ARM_PricingStatesPtr
///	Action  : Initialisation for forward looking
///           numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_HWHWQtoModel::ForwardLookingInit(ARM_NumMethodPtr& numMethod, double firstInductTime)
{
    /// Validate numeraire
	ARM_NumerairePtr numeraire = GetNumeraire();
    if( numeraire->GetType() != ARM_Numeraire::RollingCash &&
		numeraire->GetType() != ARM_Numeraire::TerminalZc )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
        ": with forward looking methods only RollingCash or terminal ZC numeraire supported by HWHWQto model");

    const ARM_ModelNameMap* const modelMap = GetModelMap();
    ARM_HullWhite1F* domModel = static_cast< ARM_HullWhite1F* >( &*((*modelMap)[DomModel]->Model()) );
    ARM_HullWhite1F* forModel = static_cast< ARM_HullWhite1F* >( &*((*modelMap)[ForModel]->Model()) );

	ARM_LN_Fx* fxModel = dynamic_cast<ARM_LN_Fx*>(&*((*modelMap)[FxModel]->Model()));
	ARM_MCMethod* mcMethod = dynamic_cast< ARM_MCMethod* >(&*numMethod);
    ARM_MeanRevertingSamplerND* samplerMR = dynamic_cast<ARM_MeanRevertingSamplerND* >(&*(mcMethod->GetSampler()));
    if(!mcMethod || !samplerMR)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": HWHWQto + MC method requires a Mean Reverting Sampler");

    return mcMethod->Init(*this,firstInductTime);
}


////////////////////////////////////////////////////
///	Class  : ARM_HWHWQtoModel
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : Convert the num method states to model states
////////////////////////////////////////////////////
void ARM_HWHWQtoModel::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
    ARM_GP_MatrixPtr numStates = states->GetNumMethodStates();
    ARM_GP_MatrixPtr modelStates = states->GetModelStates();
    
	if (FxModel<FactorCount())
		(*(GetModelMap()))[FxModel]->Model()->MCModelStatesFromToNextTime(states,timeIndex);
	for(size_t j=0;j<numStates->cols();++j)
	{
		(*modelStates)(DomModel,j)  =  (*numStates)(DomModel,j);
		(*modelStates)(ForModel,j)  =  (*numStates)(ForModel,j);
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_HWHWQtoModel
///	Routines: IntegratedLocalDrifts
///	Returns : void
///	Action  : computes the relative and absolute drift
////////////////////////////////////////////////////

void ARM_HWHWQtoModel::IntegratedLocalDrifts( const std::vector<double>& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const
{
	ARM_MultiAssetsModel::IntegratedLocalDrifts( timeSteps, relativeDrifts, absoluteDrifts );

	AddIntegratedLocalCorrections( timeSteps, absoluteDrifts );
    
}


////////////////////////////////////////////////////
///	Class   : ARM_HWHWQtoModel
///	Routines: AddIntegratedLocalCorrections
///	Returns : void
///	Action  : computes integrated local corrections
///           to be added to absolute drifts to diffuse
///           processes w.r.t. numeraire
////////////////////////////////////////////////////
void ARM_HWHWQtoModel::AddIntegratedLocalCorrections( 
	const std::vector<double>& timeSteps,
	ARM_GP_MatrixPtr& absoluteDrifts) const
{
    const ARM_ModelNameMap* const modelMap = GetModelMap();

	ARM_PricingModelPtr domModel = (*modelMap)[DomModel]->Model();
	ARM_PricingModelPtr forModel = (*modelMap)[ForModel]->Model();

    const ARM_ModelParamsHW1FStd* const domModelMRParams    = dynamic_cast<const ARM_ModelParamsHW1FStd* const>( domModel->GetModelParams() );
    const ARM_ModelParamsHW1FStd* const forModelMRParams    = dynamic_cast<const ARM_ModelParamsHW1FStd* const>( forModel->GetModelParams() );

	ARM_PricingModelPtr fxModel = (*modelMap)[FxModel]->Model();
	const ARM_ModelParams_Fx* const fxModelParam = dynamic_cast<const ARM_ModelParams_Fx* const>( fxModel->GetModelParams() );
    const ARM_ModelParamsHW1FStd* const fxModelMRParamsQ1F = dynamic_cast<const ARM_ModelParamsHW1FStd* const>( fxModel->GetModelParams() );
	const ARM_ModelParamsHW1FStd* const fxModelMRParams = dynamic_cast<const ARM_ModelParamsHW1FStd* const>( fxModel->GetModelParams() );

	ARM_CurveModelParam fxVolParam = static_cast< ARM_CurveModelParam >(fxModelMRParams->GetModelParam(ARM_ModelParamType::Volatility).ToCurveModelParam());
	
	size_t i,nbSteps = timeSteps.size();

	double corrDomFor;
    double corrDomFx;
    double corrForFx;
    double lastTime = 0.0,time;

	double convFx;
    double convDom;
    double convFor;
	double betatT;

	ARM_GP_Matrix correlMatrix;

    switch(GetNumeraire()->GetType())
    {
    case ARM_Numeraire::Cash :
        /// diffusion proba is domestic risk free
	    for( i=0; i+1<nbSteps; ++i )
	    {
            time = timeSteps[i+1];
			
			correlMatrix = GetCorrelMatrix()->Interpolate(time);
			corrDomFx    = correlMatrix(DomModel,FxModel);
			corrForFx    = correlMatrix(ForModel,FxModel);

            /// Qfor <-> Qdom (quanto) correction on foreign process
            (*absoluteDrifts)(i,ForModel) -= corrForFx * ARM_ModelParamsHW1F::HW1FStateCovariance(fxModelMRParams,forModelMRParams,lastTime,time,time);

			/// Fwd fx is diffused : Qdom(spot) <-> Qdom correction on forex process
            if (FxModel<FactorCount())
				(*absoluteDrifts)(i,FxModel) -= corrDomFx * ARM_ModelParamsHW1F::HW1FEqFxZcCovariance(fxVolParam,domModelMRParams,lastTime,time,time);
            
            lastTime = time;
	    }
        break;

    case ARM_Numeraire::RollingCash :
        /// Diffusion proba is domestic spot i.e. neutral forward at time = next time step
	    for( i=0; i+1<nbSteps; ++i )
	    {
            time = timeSteps[i+1];

			correlMatrix = GetCorrelMatrix()->Interpolate(time);
			corrDomFor    = correlMatrix(DomModel,ForModel);
			corrDomFx    = correlMatrix(DomModel,FxModel);
			corrForFx    = correlMatrix(ForModel,FxModel);

            /// Qdom <-> Qdom(spot) correction on domestic process
            (*absoluteDrifts)(i,DomModel) += ARM_ModelParamsHW1F::HW1FStateZcCovariance(domModelMRParams,domModelMRParams,lastTime,time,time,time);

            /// Qfor <-> Qdom (quanto) correction and
            /// Qdom <-> Qdom(spot) correction on foreign process
            (*absoluteDrifts)(i,ForModel) += -corrForFx  * ARM_ModelParamsHW1F::HW1FStateCovariance(fxModelMRParams,forModelMRParams,lastTime,time,time)
                                             -corrDomFor * ARM_ModelParamsHW1F::HW1FStateZcCovariance(forModelMRParams,domModelMRParams,lastTime,time,time,time);

            lastTime = time;
	    }
        break;

    case ARM_Numeraire::TerminalZc :
		for( i=0; i+1<nbSteps; ++i )
	    {
			time = timeSteps[i+1];

			correlMatrix = GetCorrelMatrix()->Interpolate(time);
			corrDomFor	= correlMatrix(DomModel,ForModel);
			corrDomFx	= correlMatrix(DomModel,FxModel);
			corrForFx	= correlMatrix(ForModel,FxModel);
			
			double numMatTime = GetNumeraire()->GetMaturity();
			double eps = -1;
			
			(*absoluteDrifts)(i,DomModel) += 0.;

            (*absoluteDrifts)(i,ForModel) +=	
				+ corrDomFor * (-1.) * ARM_ModelParamsHW1F::HW1FStateZcCovariance(forModelMRParams,domModelMRParams,lastTime,time,time,numMatTime) * eps
				- corrForFx  * (+1.) * ARM_ModelParamsHW1F::HW1FEqFxStateCovariance(fxVolParam,forModelMRParams,lastTime,time,time)
				- (1.)		 * (-1.) * ARM_ModelParamsHW1F::HW1FStateZcCovariance(forModelMRParams,forModelMRParams,lastTime,time,time,numMatTime) * eps;

			// pour stocker les réalisation du change sous la forward-neutre
			if (FxModel<FactorCount())
			{
				convFx =	+ corrDomFx  * (+1.) * ARM_ModelParamsHW1F::HW1FEqFxStateCovariance(fxVolParam,domModelMRParams,lastTime,time,time) * eps;
				convFor =	+ corrDomFor * (-1.) * ARM_ModelParamsHW1F::HW1FStateZcCovariance(domModelMRParams,forModelMRParams,lastTime,time,time,time);
				convDom =	- (1.)		 * (-1.) * ARM_ModelParamsHW1F::HW1FStateZcCovariance(domModelMRParams,domModelMRParams,lastTime,time,time,time);
				betatT = domModelMRParams->BetatT(time,numMatTime);
				(*absoluteDrifts)(i,FxModel) += betatT * ( convFx + convFor + convDom);
			}
           
            lastTime = time;
	    }
        break;

	default:
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
        ": only Cash or RollingCash numeraire supported by HWHWQto model");
    }
}

////////////////////////////////////////////////////
///	Class   : ARM_HWHWQtoModel
///	Routine : DiscountFactor
///	Returns : ARM_VectorPtr
///	Action  : computes discount factor for a given model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWHWQtoModel::DiscountFactor( 
	const string& modelName,
    double evalTime, 
	double maturityTime,
    const ARM_PricingStatesPtr& states) const
{
	const ARM_ModelNameMap& modelMap = *GetModelMap();
    ARM_ModelNameMap::const_iterator modelIdx = modelMap[modelName];

	return modelIdx->Model()->DiscountFactor(modelName, evalTime, maturityTime, states);
}

////////////////////////////////////////////////////
///	Class   : ARM_HWHWQtoModel
///	Routine : Forward
///	Returns : ARM_VectorPtr
///	Action  : Compute the forward Fx using the Fx model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWHWQtoModel::Forward(
	const string& modelName, 
    double evalTime,
	double expiryTime,
	double settlementTime,
	double payTime,
    const ARM_PricingStatesPtr& states) const
{
 	const ARM_ModelNameMap& modelMap = *GetModelMap();

    if( modelName != modelMap[FxModel]->ModelName() )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " 2IR + FX model can't compute a forward fx on " + modelName);

    if(expiryTime < evalTime - K_NEW_DOUBLE_TOL)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " expired forward Fx is not supported");

    const ARM_LN_Fx& fxModel = static_cast<ARM_LN_Fx&>( * modelMap[FxModel]->Model() );

    /// Compute the forward Fx
    return fxModel.Forward(fxModel.GetModelName(),evalTime,expiryTime,settlementTime,payTime,states);
}


////////////////////////////////////////////////////
///	Class   : ARM_HWHWQtoModel
///	Routines: toString
///	Returns :
///	Action  : stringify the object
////////////////////////////////////////////////////
string ARM_HWHWQtoModel::toString(const string& indent, const string& nextIndent) const
{
 	const ARM_ModelNameMap& modelMap = *GetModelMap();

    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Hybrid 2 Interest Rates & Forex model (2IR+FX)\n";
    os << indent << "----------------------------------------------\n\n";

	if( GetRefModel() )
	{
		os << "Payment model : " << GetRefModel()->GetModelName() << "\n\n";
	}

	if( GetCorrelMatrix() )
	{
		os << "Correlation matrix\n";
		os << GetCorrelMatrix()->toString(indent,nextIndent);
	}

    os << indent << "\n\n------> Domestic Stochastic IR Model <------\n";
    os << modelMap[DomModel]->Model()->toString(indent,nextIndent);

    os << indent << "\n\n------> Foreign Stochastic IR Model  <------\n";
    os << modelMap[ForModel]->Model()->toString(indent,nextIndent);

    os << indent << "\n\n------>      Stochastic FX Model     <------\n";
    os << modelMap[FxModel]->Model()->toString(indent,nextIndent);

    os << indent << "\n\n------>     Domestic Basis Model     <------\n";
    ARM_ForwardMargin& domBasisModel = static_cast< ARM_ForwardMargin& >(* modelMap[DomBasisModel]->Model());
    bool prevDump = domBasisModel.GetRefModelDump();
    domBasisModel.SetRefModelDump(false);
    os << domBasisModel.toString(indent,nextIndent);
    domBasisModel.SetRefModelDump(prevDump);

    os << indent << "\n\n------>     Foreign Basis Model      <------\n";
    ARM_ForwardMargin& forBasisModel = static_cast< ARM_ForwardMargin& >(* modelMap[ForBasisModel]->Model());
    prevDump = forBasisModel.GetRefModelDump();
    forBasisModel.SetRefModelDump(false);
    os << forBasisModel.toString(indent,nextIndent);
    forBasisModel.SetRefModelDump(prevDump);

	if(itsForLocalModel)
    {
        os << indent << "\n\n------>     Foreign Local Model      <------\n";
        os << modelMap[ForLocalModel]->Model()->toString(indent,nextIndent);
    }

    return os.str();
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

