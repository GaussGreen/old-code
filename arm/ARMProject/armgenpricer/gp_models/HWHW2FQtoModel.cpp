/*!
 *
 * Copyright (c) CM CIB January 2005 Paris
 *
 *	\version 1.0
 *	\date May 2007
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/HWHW2FQtoModel.h"

/// gpmodels
#include "gpmodels/HW1F.h"
#include "gpmodels/ModelParamsHW1F.h"
#include "gpmodels/HW2F.h"
#include "gpmodels/ModelParamsHW2F.h"
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

const size_t ARM_HWHW2FQtoModel::NbLocalModels = 1;

const double MAX_LNVOL_FX = 10.0; // 1000%

////////////////////////////////////////////////////
///	Class   : ARM_HWHW2FQtoModel
///	Routines: Constructor
///	Returns :
///	Action  : builds the object
///	model types contains the integer on the various model
///	using the correlation matrix
////////////////////////////////////////////////////
ARM_HWHW2FQtoModel::ARM_HWHW2FQtoModel(	const ARM_ModelNameMap&	modelNameMap, 
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
///	Class   : ARM_HWHW2FQtoModel
///	Routine : InitSubModels
///	Returns : void
///	Action  : intialises the sub models that are linked to a model
////////////////////////////////////////////////////
void ARM_HWHW2FQtoModel::Validate( bool fxFlag )
{
	/// number of models
	ARM_ModelNameMap* modelNameMap	= GetModelMap();
	if(modelNameMap->size() < NbModels)
    {
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : not enough models to build HW x HW2F x Quanto");
    }
    else if(modelNameMap->size() > NbModels)
    {
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : too much models to build HW x HW2F x Quanto");
    }

	/// correlation
	if( !GetCorrelMatrix() )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : please provide a correlation");

	if(	GetCorrelMatrix()->rows() != 4 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : correlation should be 4x4!");	

	/// check model types
	if( !dynamic_cast<ARM_HullWhite1F*>( &*(*modelNameMap)[DomModel]->Model() ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : dom model has to be HW1F IR Model");

	if( !dynamic_cast<ARM_HullWhite2F*>( &*(*modelNameMap)[ForModel]->Model() ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : for model has to be HW2F IR Model");

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
///	Class   : ARM_HWHW2FQtoModel
///	Routines: Copy constructor
///	Returns :
///	Action  : constructor
ARM_HWHW2FQtoModel::ARM_HWHW2FQtoModel(const ARM_HWHW2FQtoModel& rhs)
:	ARM_MultiAssetsMeanReverting(rhs),
	itsForLocalModel(rhs.itsForLocalModel)
{
    itsDriftCorrections = rhs.itsDriftCorrections != ARM_GP_MatrixPtr(NULL)
                          ? ARM_GP_MatrixPtr( static_cast< ARM_GP_Matrix* >(rhs.itsDriftCorrections->Clone()) )
                          : ARM_GP_MatrixPtr(NULL);
}


////////////////////////////////////////////////////
///	Class   : ARM_HWHW2FQtoModel
///	Routine : SupportAnalyticMarginal
///	Returns : void
///	Action  : decide if over sampling between model time
///           step is needed
////////////////////////////////////////////////////
bool ARM_HWHW2FQtoModel::SupportAnalyticMarginal() const
{
    const ARM_NumerairePtr numeraire = GetNumeraire();
    return  numeraire == ARM_NumerairePtr(NULL) ||
            numeraire->GetType() == ARM_Numeraire::RollingEvent;
}

////////////////////////////////////////////////////
///	Class   : ARM_HWHW2FQtoModel
///	Routine : NeedMCIntegProcess
///	Returns : ARM_BoolVector
///	Action  : to say how MC method give simulated processes
////////////////////////////////////////////////////
ARM_BoolVector ARM_HWHW2FQtoModel::NeedMCIntegProcess() const
{
	// Labels 
	size_t Dom=0, For1= 1, For2=2, Fx = 3; 

    ARM_BoolVector SimulProcess(FactorCount(),true);
	if (Fx<FactorCount())
		SimulProcess[Fx] = false;

	ARM_NumerairePtr numeraire = GetNumeraire();
	if(numeraire->GetType() == ARM_Numeraire::Cash)
	{
		/// Domestic & foreign processes are incremental because
		/// the discretisation scheme is Euler
		SimulProcess[Dom] = false;
		SimulProcess[For1] = false;
		SimulProcess[For2] = false;
	}

    return SimulProcess;
}

///////////////////////////////////////////////////
///	Class   : ARM_HWHW2FQtoModel
///	Routine : Init
///	Returns : ARM_PricingStatesPtr
///	Action  : intialise the model
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_HWHW2FQtoModel::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
	/// Validate that the payment currency is the domestic mod²el!
    const ARM_ModelNameMap* const modelMap = GetModelMap();
    const string& domBasisModelName = (*modelMap)[DomBasisModel]->Model()->GetModelName();
    const string& domModelName      = (*modelMap)[DomModel]->Model()->GetModelName();
	if( payModelName != domBasisModelName && payModelName != domModelName)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + 
		": HW HW2F Qto model allow only paimenent in " + domBasisModelName + " or " + domModelName);

	ARM_NumerairePtr numeraire = GetNumeraire();
    ARM_NumMethodPtr numMethod = GetNumMethod();
    if( numMethod != ARM_NumMethodPtr(NULL) && 
        numeraire != ARM_NumerairePtr(NULL) &&
        numMethod->GetPricingDirection() != ARM_NumMethod::GP_AMBIGUOUS &&
        (numMethod->GetPricingDirection() == ARM_NumMethod::GP_FWDLOOKING ||
		 numMethod->GetPricingDirection() == ARM_NumMethod::GP_FWDBCKWDLOOKING) &&
        (numeraire->GetType() == ARM_Numeraire::RollingEvent||
         numeraire->GetType() == ARM_Numeraire::Cash ) 
	   )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
        ": only TerminalZc or RollingCash numeraire are supported by HW HW2F Qto model");

	/// Delegate to common code
	return ARM_MultiAssetsModel::Init( payModelName, timeInfos );
}


////////////////////////////////////////////////////
///	Class   : ARM_HWHW2FQtoModel
///	Routine : ComputeIntegratedFwdFxVCV
///	Returns : void
///	Action  : computes the integrated VCV matrixes for
///           domestic & foreign state variables &
///           forward Fx process
////////////////////////////////////////////////////
void ARM_HWHW2FQtoModel::ComputeIntegratedFwdFxVCV(double step, double nextStep,
        const ARM_ModelParamsHW1FStd* const domModelParams,
        const ARM_ModelParamsHW2F* const forModelParams,
        const ARM_ModelParamsHW1FStd* const fxModelParams,
        const ARM_CurveModelParam& fxVol,
        const ARM_GP_Matrix& correlMatrix,
        ARM_GP_Matrix& variances) const
{

	// Pointer initialisation
	ARM_GP_TriangularMatrix* varCovarFor =0; 

	ARM_ModelParamsHW1FStd* for1ModelParams=0;
	ARM_ModelParamsHW1FStd* for2ModelParams=0;

	ARM_ModelParamVector modelParams(2);
	modelParams[0] = 0; 
	modelParams[1] = 0;

	// Labels 
	size_t Dom=0, For1= 1, For2=2, Fx = 3; 

	try { // all allocations are encapsulated in a try/catch block
		// Transform ARM_ModelParamsHW2FStd into two ARM_ModelParamsHW1FStd
		modelParams[0] = dynamic_cast<ARM_ModelParam*> (
			forModelParams->GetModelParam(ARM_ModelParamType::Volatility).Clone()
														); 
		modelParams[1] = dynamic_cast<ARM_ModelParam*> (
			forModelParams->GetModelParam(ARM_ModelParamType::MeanReversion).Clone()
														); 
		for1ModelParams  = new ARM_ModelParamsHW1FStd(modelParams); 
 
		// Update parameters
		const ARM_ModelParam& volRatio = forModelParams->GetModelParam(ARM_ModelParamType::VolatilityRatio); 
		ARM_CurveModelParam* volCurve = dynamic_cast<ARM_CurveModelParam*> (modelParams[0]); 
		if(volRatio.size()!=volCurve->size())
			ARM_THROW( ERR_INVALID_ARGUMENT, "Volatility and Volatility ratio should have the same size" );
		size_t i;
		for(i=0; i<volCurve->size(); i++) {
			double t = volCurve->GetTimeAtPoint(i) ; 
			volCurve->SetValueAtPoint(i, volCurve->GetValueAtPoint(i)*volRatio.GetValue(t));
		}
		const ARM_ModelParam& mrSpread = forModelParams->GetModelParam(ARM_ModelParamType::MeanReversionSpread); 
		ARM_CurveModelParam* mrCurve = dynamic_cast<ARM_CurveModelParam*> (modelParams[1]); 
		if(mrSpread.size()!=mrCurve->size())
			ARM_THROW( ERR_INVALID_ARGUMENT, "MeanReversion and MeanReversionSpread should have the same size" );
		for(i=0; i<mrCurve->size(); i++) {
			double t = mrCurve->GetTimeAtPoint(i) ; 
			mrCurve->SetValueAtPoint(i, mrCurve->GetValueAtPoint(i)+mrSpread.GetValue(t));
		}
		for2ModelParams  = new ARM_ModelParamsHW1FStd(modelParams); 
		
		/// Domestic variances
		variances(Dom,Dom)   = domModelParams->StateLocalVariance(step,nextStep,nextStep);

		/// Foreign variances
		varCovarFor = forModelParams->StateLocalVariance(step,nextStep);
		variances(For1,For1)   = (*varCovarFor)(0,0);
		variances(For1,For2)   = (*varCovarFor)(0,1);
		variances(For2,For2)   = (*varCovarFor)(1,1);
		variances(For2,For1)   = variances(For1,For2); //(*varCovarFor)(1,0);
		
		/// Domestic & Foreign covariances
		variances(Dom,For1) = correlMatrix(Dom,For1) * ARM_ModelParamsHW1F::HW1FStateCovariance(domModelParams, for1ModelParams,step,nextStep,nextStep);
		variances(Dom,For2) = correlMatrix(Dom,For2) * ARM_ModelParamsHW1F::HW1FStateCovariance(domModelParams, for2ModelParams,step,nextStep,nextStep);
		variances(For1,Dom) = variances(Dom,For1);
		variances(For2,Dom) = variances(Dom,For2);

		if (Fx<FactorCount())
		{
			/// Forex variances
			double varFx           = fxModelParams->StateLocalVariance(step,nextStep,nextStep);
			double varZcDom        = ARM_ModelParamsHW1F::HW1FZcCovariance( domModelParams, domModelParams, step,nextStep,nextStep );
			double varZcFor1       = ARM_ModelParamsHW1F::HW1FZcCovariance( for1ModelParams, for1ModelParams, step,nextStep,nextStep );
			double varZcFor2       = ARM_ModelParamsHW1F::HW1FZcCovariance( for2ModelParams, for2ModelParams, step,nextStep,nextStep );
			double covarZcFor1For2 = ARM_ModelParamsHW1F::HW1FZcCovariance( for1ModelParams, for2ModelParams, step,nextStep,nextStep );

			double covarZcDomZcFor1 = ARM_ModelParamsHW1F::HW1FZcCovariance( domModelParams, for1ModelParams, step,nextStep,nextStep );
			double covarZcDomZcFor2 = ARM_ModelParamsHW1F::HW1FZcCovariance( domModelParams, for2ModelParams, step,nextStep,nextStep );
			double covarZcDomFx     = ARM_ModelParamsHW1F::HW1FEqFxZcCovariance( fxVol, domModelParams, step,nextStep,nextStep );
			double covarZcFor1Fx    = ARM_ModelParamsHW1F::HW1FEqFxZcCovariance( fxVol, for1ModelParams, step,nextStep,nextStep );
			double covarZcFor2Fx    = ARM_ModelParamsHW1F::HW1FEqFxZcCovariance( fxVol, for2ModelParams, step,nextStep,nextStep );

			variances(Fx,Fx) =  
				varFx + varZcDom + varZcFor1 + varZcFor2 
				+ 2*correlMatrix(For1,For2)*covarZcFor1For2		// ++
				- 2*correlMatrix(Dom,For1)*covarZcDomZcFor1
				- 2*correlMatrix(Dom,For2)*covarZcDomZcFor2
				- 2*correlMatrix(Dom,Fx)*covarZcDomFx
				+ 2*correlMatrix(For1,Fx)*covarZcFor1Fx
				+ 2*correlMatrix(For2,Fx)*covarZcFor2Fx;

			double varLimit = K_NEW_DOUBLE_TOL*K_NEW_DOUBLE_TOL;
			if(variances(Fx,Fx)<varLimit)
				variances(Fx,Fx) = varLimit; // numerical noise

			/// Domestic & Forex covariances
			double covarXDomZcDom  = ARM_ModelParamsHW1F::HW1FStateZcCovariance( domModelParams, domModelParams, step,nextStep,nextStep,nextStep );
			double covarXDomZcFor1 = ARM_ModelParamsHW1F::HW1FStateZcCovariance( domModelParams, for1ModelParams, step,nextStep,nextStep,nextStep );
			double covarXDomZcFor2 = ARM_ModelParamsHW1F::HW1FStateZcCovariance( domModelParams, for2ModelParams, step,nextStep,nextStep,nextStep );
			double covarXDomFx     = ARM_ModelParamsHW1F::HW1FStateCovariance(   domModelParams, fxModelParams, step,nextStep,nextStep );

			variances(Dom,Fx) =	
				- covarXDomZcDom
				+ correlMatrix(Dom,For1)*covarXDomZcFor1
				+ correlMatrix(Dom,For2)*covarXDomZcFor2
				+ correlMatrix(Dom,Fx)*covarXDomFx; 

			variances(Fx,Dom) = variances(Dom,Fx);

			/// Foreign & Forex covariances
			double covarXFor1ZcDom  = ARM_ModelParamsHW1F::HW1FStateZcCovariance( for1ModelParams, domModelParams, step,nextStep,nextStep,nextStep );
			double covarXFor1ZcFor1 = ARM_ModelParamsHW1F::HW1FStateZcCovariance( for1ModelParams, for1ModelParams, step,nextStep,nextStep,nextStep );
			double covarXFor1ZcFor2 = ARM_ModelParamsHW1F::HW1FStateZcCovariance( for1ModelParams, for2ModelParams, step,nextStep,nextStep,nextStep );
			double covarXFor1Fx     = ARM_ModelParamsHW1F::HW1FStateCovariance(   for1ModelParams, fxModelParams, step,nextStep,nextStep );

			variances(For1,Fx) = 
				covarXFor1ZcFor1
				+ correlMatrix(For1,For2)*covarXFor1ZcFor2 
				- correlMatrix(Dom,For1)*covarXFor1ZcDom
				+ correlMatrix(For1,Fx)*covarXFor1Fx; 
				
			variances(Fx,For1) = variances(For1,Fx);

			double covarXFor2ZcDom  = ARM_ModelParamsHW1F::HW1FStateZcCovariance( for2ModelParams, domModelParams, step,nextStep,nextStep,nextStep );
			double covarXFor2ZcFor2 = ARM_ModelParamsHW1F::HW1FStateZcCovariance( for2ModelParams, for2ModelParams, step,nextStep,nextStep,nextStep );
			double covarXFor2ZcFor1 = ARM_ModelParamsHW1F::HW1FStateZcCovariance( for2ModelParams, for1ModelParams, step,nextStep,nextStep,nextStep );
			double covarXFor2Fx     = ARM_ModelParamsHW1F::HW1FStateCovariance(   for2ModelParams, fxModelParams, step,nextStep,nextStep );

			variances(For2,Fx) = 
				covarXFor2ZcFor2
				+ correlMatrix(For1,For2)*covarXFor2ZcFor1 
				- correlMatrix(Dom,For2)*covarXFor2ZcDom
				+ correlMatrix(For2,Fx)*covarXFor2Fx;

			variances(Fx,For2) = variances(For2,Fx); 

		}
	} catch(...) {
			delete varCovarFor;		varCovarFor =0; 

			delete for1ModelParams; for1ModelParams=0;
			delete for2ModelParams; for2ModelParams=0;
			
			delete modelParams[0];	modelParams[0] = 0; 
			delete modelParams[1];	modelParams[1] = 0;
		throw; 
	}

	
	delete varCovarFor;		varCovarFor =0; 

	delete for1ModelParams; for1ModelParams=0;
	delete for2ModelParams; for2ModelParams=0;
	
	delete modelParams[0];	modelParams[0] = 0; 
	delete modelParams[1];	modelParams[1] = 0;
	
}

////////////////////////////////////////////////////
///	Class   : ARM_HWHW2FQtoModel
///	Routine : NumMethodStateLocalGlobalVariances
///	Returns : void
///	Action  : computes the integrated covariances
////////////////////////////////////////////////////

void ARM_HWHW2FQtoModel::NumMethodStateLocalGlobalVariances( const std::vector<double>& timeSteps,
	ARM_MatrixVector& localVariances,
	ARM_MatrixVector& variances ) const
{
	
	// Labels
	size_t Dom=0, For1=1, For2=2, Fx=3; 

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
	const ARM_ModelParamsHW1FStd* const domModelParams  
		= static_cast<const ARM_ModelParamsHW1FStd* const>( domModel->GetModelParams() );
	const ARM_ModelParamsHW2F* const forModelParams  
		= static_cast<const ARM_ModelParamsHW2F* const>( forModel->GetModelParams() );
	const ARM_ModelParamsHW1FStd* const fxModelParams   
		= static_cast<const ARM_ModelParamsHW1FStd* const>( fxModel->GetModelParams() );
	const ARM_CurveModelParam& fxVol                    
		= dynamic_cast<const ARM_CurveModelParam&>( fxModelParams->GetModelParam(fxModelParams->GetVolatilityType()) );
	const ARM_ModelParam& forCorrelParam 
		= forModelParams->GetModelParam(ARM_ModelParamType::Correlation); 

    for(size_t i=0;i<nbSteps-1;++i)
	{
		correlMatrix = correlCurveMatrix->Interpolate(step);
		correlMatrix(For1,For2) = forCorrelParam.GetValue(step);
		correlMatrix(For2,For1) = correlMatrix(For1,For2);

		nextStep			= timeSteps[i+1];

		localVariances[i]	= new ARM_GP_TriangularMatrix(factorNb,0.0);

		ComputeIntegratedFwdFxVCV(step,nextStep,
			domModelParams,forModelParams,fxModelParams,fxVol,correlMatrix,
		    *(localVariances[i]));

		variances[i+1]		= new ARM_GP_TriangularMatrix(factorNb,0.0);

		ComputeIntegratedFwdFxVCV(0.0,nextStep,
			domModelParams,forModelParams,fxModelParams,fxVol,correlMatrix,
			*(variances[i+1]));
        
		if (Fx<FactorCount())
			fxLocalVariances[i] = static_cast< ARM_GP_Matrix* >(localVariances[i]->Clone());

		step=nextStep;
	}

	if (Fx<FactorCount())
		(*modelMap)[FxModel]->Model()->SetModelStateLocalVars(fxLocalVariances);
}

////////////////////////////////////////////////////
///	Class   : ARM_HWHW2FQtoModel
///	Routine : BackwardLookingInit
///	Returns : ARM_PricingStatesPtr
///	Action  : Initialisation for backward looking
///           numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_HWHW2FQtoModel::BackwardLookingInit(ARM_NumMethodPtr& numMethod, double firstInductTime)
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
///	Class  : ARM_HWHW2FQtoModel
///	Routine: TreeStatesToModelStates
///	Returns: void
///	Action : Convert the num method states to model states
////////////////////////////////////////////////////
void ARM_HWHW2FQtoModel::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
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
///	Class  : ARM_HWHW2FQtoModel
///	Routine: ComputeNumeraireTimes
///	Returns: 
///	Action : compute the corresonding numeraire times
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWHW2FQtoModel::ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const
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
///	Class   : ARM_HWHW2FQtoModel
///	Routine : ForwardLookingInit
///	Returns : ARM_PricingStatesPtr
///	Action  : Initialisation for forward looking
///           numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_HWHW2FQtoModel::ForwardLookingInit(ARM_NumMethodPtr& numMethod, double firstInductTime)
{
    /// Validate numeraire
	ARM_NumerairePtr numeraire = GetNumeraire();
    if( numeraire->GetType() != ARM_Numeraire::RollingCash &&
		numeraire->GetType() != ARM_Numeraire::TerminalZc )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
        ": with forward looking methods only RollingCash or terminal ZC numeraire are supported by HWHW2FQto model");

    const ARM_ModelNameMap* const modelMap = GetModelMap();
    ARM_HullWhite1F* domModel = static_cast< ARM_HullWhite1F* >( &*((*modelMap)[DomModel]->Model()) );
    ARM_HullWhite2F* forModel = static_cast< ARM_HullWhite2F* >( &*((*modelMap)[ForModel]->Model()) );

	ARM_LN_Fx* fxModel = dynamic_cast<ARM_LN_Fx*>(&*((*modelMap)[FxModel]->Model()));
	ARM_MCMethod* mcMethod = dynamic_cast< ARM_MCMethod* >(&*numMethod);
    ARM_MeanRevertingSamplerND* samplerMR = dynamic_cast<ARM_MeanRevertingSamplerND* >(&*(mcMethod->GetSampler()));
    if(!mcMethod || !samplerMR)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": HWHW2FQto + MC method requires a Mean Reverting Sampler");

    return mcMethod->Init(*this,firstInductTime);
}


////////////////////////////////////////////////////
///	Class  : ARM_HWHW2FQtoModel
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : Convert the num method states to model states
////////////////////////////////////////////////////
void ARM_HWHW2FQtoModel::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
    ARM_GP_MatrixPtr numStates = states->GetNumMethodStates();
    ARM_GP_MatrixPtr modelStates = states->GetModelStates();
    
	if (3<FactorCount())
		(*(GetModelMap()))[FxModel]->Model()->MCModelStatesFromToNextTime(states,timeIndex);
	for(size_t j=0;j<numStates->cols();++j)
	{
		for(size_t i=0; i<3; ++i) 
			(*modelStates)(i,j)  =  (*numStates)(i,j);
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_HWHW2FQtoModel
///	Routines:  LocalDrifts
///	Returns : void
///	Action  : computes the relative and absolute drift
////////////////////////////////////////////////////

void ARM_HWHW2FQtoModel::IntegratedLocalDrifts( const std::vector<double>& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const
{
	ARM_MultiAssetsModel::IntegratedLocalDrifts( timeSteps, relativeDrifts, absoluteDrifts );

    if( GetNumMethod().IsNull() || (GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_FWDLOOKING) || (GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_FWDBCKWDLOOKING))
    {
	    AddIntegratedLocalCorrections( timeSteps, absoluteDrifts );
    }
}


////////////////////////////////////////////////////
///	Class   : ARM_HWHW2FQtoModel
///	Routines: AddIntegratedLocalCorrections
///	Returns : void
///	Action  : computes integrated local corrections
///           to be added to absolute drifts to diffuse
///           processes w.r.t. numeraire
////////////////////////////////////////////////////
void ARM_HWHW2FQtoModel::AddIntegratedLocalCorrections( 
	const std::vector<double>& timeSteps,
	ARM_GP_MatrixPtr& absoluteDrifts) const
{
	// Labels 
	size_t Dom=0, For1= 1, For2=2, Fx = 3; 

	// Pointer initialisation
	ARM_ModelParamsHW1FStd* for1ModelParams=0;
	ARM_ModelParamsHW1FStd* for2ModelParams=0;
	ARM_ModelParamVector modelParams(2);
	modelParams[0] = 0;
	modelParams[1] = 0;

	try {
		const ARM_ModelNameMap* const modelMap = GetModelMap();

		ARM_PricingModelPtr domModel = (*modelMap)[DomModel]->Model();
		ARM_PricingModelPtr forModel = (*modelMap)[ForModel]->Model();

		const ARM_ModelParamsHW1FStd* const domModelParams    
			= dynamic_cast<const ARM_ModelParamsHW1FStd* const>( domModel->GetModelParams() );
		const ARM_ModelParamsHW2F* const forModelParams    
			= dynamic_cast<const ARM_ModelParamsHW2F* const>( forModel->GetModelParams() );

		ARM_PricingModelPtr fxModel = (*modelMap)[FxModel]->Model();
		const ARM_ModelParamsHW1FStd* const fxModelParams 
			= dynamic_cast<const ARM_ModelParamsHW1FStd* const>( fxModel->GetModelParams() );
		ARM_CurveModelParam fxVolParam = static_cast< ARM_CurveModelParam >(fxModelParams->GetModelParam(ARM_ModelParamType::Volatility).ToCurveModelParam());
		
		size_t i,nbSteps = timeSteps.size();

		// Transform ARM_ModelParamsHW2FStd into two ARM_ModelParamsHW1FStd
		modelParams[0] = dynamic_cast<ARM_ModelParam*> (
			forModelParams->GetModelParam(ARM_ModelParamType::Volatility).Clone()
														); 
		modelParams[1] = dynamic_cast<ARM_ModelParam*> (
			forModelParams->GetModelParam(ARM_ModelParamType::MeanReversion).Clone()
														); 
		for1ModelParams  = new ARM_ModelParamsHW1FStd(modelParams); 

		// Update parameters
		const ARM_ModelParam& volRatio = forModelParams->GetModelParam(ARM_ModelParamType::VolatilityRatio); 
		ARM_CurveModelParam* volCurve = dynamic_cast<ARM_CurveModelParam*> (modelParams[0]); 
		if(volRatio.size()!=volCurve->size())
			ARM_THROW( ERR_INVALID_ARGUMENT, "Volatility and Volatility ratio should have the same size" );
		for(i=0; i<volCurve->size(); i++) {
			double t = volCurve->GetTimeAtPoint(i) ; 
			volCurve->SetValueAtPoint(i, volCurve->GetValueAtPoint(i)*volRatio.GetValue(t));
		}
		const ARM_ModelParam& mrSpread = forModelParams->GetModelParam(ARM_ModelParamType::MeanReversionSpread); 
		ARM_CurveModelParam* mrCurve = dynamic_cast<ARM_CurveModelParam*> (modelParams[1]); 
		if(mrSpread.size()!=mrCurve->size())
			ARM_THROW( ERR_INVALID_ARGUMENT, "MeanReversion and MeanReversionSpread should have the same size" );
		for(i=0; i<mrCurve->size(); i++) {
			double t = mrCurve->GetTimeAtPoint(i) ; 
			mrCurve->SetValueAtPoint(i, mrCurve->GetValueAtPoint(i)+mrSpread.GetValue(t));
		}
		for2ModelParams  = new ARM_ModelParamsHW1FStd(modelParams); 

		const ARM_ModelParam& forCorrelParam 
			= forModelParams->GetModelParam(ARM_ModelParamType::Correlation); 

		double corrDomFor1;
		double corrDomFor2;
		double corrDomFx;
		double corrFor1Fx;
		double corrFor2Fx;
		double corrFor1For2;

		double lastTime = 0.0,time;
		
		double convFx;
		double convDom;
		double convFor;
		double betatT;

		ARM_GP_Matrix correlMatrix;

		switch(GetNumeraire()->GetType())
		{
		case ARM_Numeraire::RollingCash :
			for( i=0; i+1<nbSteps; ++i )
			{
				time = timeSteps[i+1];

				correlMatrix = GetCorrelMatrix()->Interpolate(time);
				corrDomFor1	= correlMatrix(Dom,For1);
				corrDomFor2	= correlMatrix(Dom,For2);
				corrFor1For2= forCorrelParam.GetValue(time);
				corrDomFx	= correlMatrix(Dom,Fx);
				corrFor1Fx	= correlMatrix(For1,Fx);
				corrFor2Fx	= correlMatrix(For2,Fx);
				
				/// Qdom <-> Qdom(spot) correction on domestic process
				(*absoluteDrifts)(i,Dom) += ARM_ModelParamsHW1F::HW1FStateZcCovariance(domModelParams,domModelParams,lastTime,time,time,time);

				/// Qfor <-> Qdom (quanto) correction and
				/// Qdom <-> Qdom(spot) correction on foreign process
				(*absoluteDrifts)(i,For1) += -corrFor1Fx  * ARM_ModelParamsHW1F::HW1FStateCovariance(fxModelParams,for1ModelParams,lastTime,time,time)
											 -corrDomFor1 * ARM_ModelParamsHW1F::HW1FStateZcCovariance(for1ModelParams,domModelParams,lastTime,time,time,time)
										;//	 -corrFor1For2 * ARM_ModelParamsHW1F::HW1FStateZcCovariance(for1ModelParams,for2ModelParams,lastTime,time,time,time);

				(*absoluteDrifts)(i,For2) += -corrFor2Fx  * ARM_ModelParamsHW1F::HW1FStateCovariance(fxModelParams,for2ModelParams,lastTime,time,time)
											 -corrDomFor2 * ARM_ModelParamsHW1F::HW1FStateZcCovariance(for2ModelParams,domModelParams,lastTime,time,time,time)
										;//	 -corrFor1For2 * ARM_ModelParamsHW1F::HW1FStateZcCovariance(for2ModelParams,for1ModelParams,lastTime,time,time,time);

	            lastTime = time;
			}
			break;
		case ARM_Numeraire::TerminalZc :
			for( i=0; i+1<nbSteps; ++i )
			{
				time = timeSteps[i+1];

				correlMatrix = GetCorrelMatrix()->Interpolate(time);
				corrDomFor1	= correlMatrix(Dom,For1);
				corrDomFor2	= correlMatrix(Dom,For2);
				corrDomFx	= correlMatrix(Dom,Fx);
				corrFor1Fx	= correlMatrix(For1,Fx);
				corrFor2Fx	= correlMatrix(For2,Fx);
				corrFor1For2= forCorrelParam.GetValue(time);

				double numMatTime = GetNumeraire()->GetMaturity();
				double eps = -1;
				
				(*absoluteDrifts)(i,Dom) += 0.;

				(*absoluteDrifts)(i,For1) +=	
					+ corrDomFor1 * (-1.) * ARM_ModelParamsHW1F::HW1FStateZcCovariance(for1ModelParams,domModelParams,lastTime,time,time,numMatTime) * eps
					- corrFor1Fx  * (+1.) * ARM_ModelParamsHW1F::HW1FEqFxStateCovariance(fxVolParam,for1ModelParams,lastTime,time,time)
					- (1.)		 * (-1.) * ARM_ModelParamsHW1F::HW1FStateZcCovariance(for1ModelParams,for1ModelParams,lastTime,time,time,numMatTime) * eps
					- (corrFor1For2)* (-1.) * ARM_ModelParamsHW1F::HW1FStateZcCovariance(for1ModelParams,for2ModelParams,lastTime,time,time,numMatTime) * eps
						;

				(*absoluteDrifts)(i,For2) +=	
					+ corrDomFor2 * (-1.) * ARM_ModelParamsHW1F::HW1FStateZcCovariance(for2ModelParams,domModelParams,lastTime,time,time,numMatTime) * eps
					- corrFor2Fx  * (+1.) * ARM_ModelParamsHW1F::HW1FEqFxStateCovariance(fxVolParam,for2ModelParams,lastTime,time,time)
					- (1.)		 * (-1.) * ARM_ModelParamsHW1F::HW1FStateZcCovariance(for2ModelParams,for2ModelParams,lastTime,time,time,numMatTime) * eps
					- (corrFor1For2)* (-1.) * ARM_ModelParamsHW1F::HW1FStateZcCovariance(for2ModelParams,for1ModelParams,lastTime,time,time,numMatTime) * eps
						;

				// pour stocker les réalisations du change sous la forward-neutre
				if (Fx<FactorCount())
				{
					
					convFx =	+ corrDomFx  * (+1.) * ARM_ModelParamsHW1F::HW1FEqFxStateCovariance(fxVolParam,domModelParams,lastTime,time,time) * eps;
					convFor =	+ corrDomFor1 * (-1.) * ARM_ModelParamsHW1F::HW1FStateZcCovariance(domModelParams,for1ModelParams,lastTime,time,time,time);
					convFor +=	 corrDomFor2 * (-1.) * ARM_ModelParamsHW1F::HW1FStateZcCovariance(domModelParams,for2ModelParams,lastTime,time,time,time);
					convDom =	- (1.)		 * (-1.) * ARM_ModelParamsHW1F::HW1FStateZcCovariance(domModelParams,domModelParams,lastTime,time,time,time);
					betatT = domModelParams->BetatT(time,numMatTime);
					(*absoluteDrifts)(i,Fx) += betatT * ( convFx + convFor + convDom);			
				}
           
				lastTime = time;
			}
			break;

		default:
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
			": only TerminalZC numeraire is supported by HWHW2FQto model");
		}
	} catch(...) {
		delete for1ModelParams;		for1ModelParams=0;
		delete for2ModelParams;		for2ModelParams=0;
			
		delete (modelParams[0]);	modelParams[0] = 0; 
		delete (modelParams[1]);	modelParams[1] = 0;
		throw; 
	}

	delete for1ModelParams;		for1ModelParams=0;
	delete for2ModelParams;		for2ModelParams=0;
	
	delete (modelParams[0]);	modelParams[0] = 0; 
	delete (modelParams[1]);	modelParams[1] = 0;
  
	
}

////////////////////////////////////////////////////
///	Class   : ARM_HWHW2FQtoModel
///	Routine : DiscountFactor
///	Returns : ARM_VectorPtr
///	Action  : computes discount factor for a given model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWHW2FQtoModel::DiscountFactor( 
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
///	Class   : ARM_HWHW2FQtoModel
///	Routine : Forward
///	Returns : ARM_VectorPtr
///	Action  : Compute the forward Fx using the Fx model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWHW2FQtoModel::Forward(
	const string& modelName, 
    double evalTime,
	double expiryTime,
	double settlementTime,
	double payTime,
    const ARM_PricingStatesPtr& states) const
{
 	const ARM_ModelNameMap& modelMap = *GetModelMap();

    if( modelName != modelMap[FxModel]->ModelName() )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " HW + HW2F + FX model can't compute a forward fx on " + modelName);

    if(expiryTime < evalTime - K_NEW_DOUBLE_TOL)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " expired forward Fx is not supported");

    const ARM_LN_Fx& fxModel = static_cast<ARM_LN_Fx&>( * modelMap[FxModel]->Model() );

    /// Compute the forward Fx
    return fxModel.Forward(fxModel.GetModelName(),evalTime,expiryTime,settlementTime,payTime,states);
}


////////////////////////////////////////////////////
///	Class   : ARM_HWHW2FQtoModel
///	Routines: toString
///	Returns :
///	Action  : stringify the object
////////////////////////////////////////////////////
string ARM_HWHW2FQtoModel::toString(const string& indent, const string& nextIndent) const
{
 	const ARM_ModelNameMap& modelMap = *GetModelMap();

    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Hybrid 2 Interest Rates & Forex model (HW+HW2F+FX)\n";
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

