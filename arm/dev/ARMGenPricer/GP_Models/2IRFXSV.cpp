/*!
 *
 * Copyright (c) CM CIB January 2005 Paris
 *
 *	\file 1IRFXModel.cpp
 *
 *  \brief N+1 interest rate + N fx multi assets model
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date June 2007
 */

//#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/2IRFXSV.h"

/// gpbase
#include "gpbase/checkarg.h"
#include "gpbase/stringmanip.h"
#include "gpbase/gpmatrixlinalg.h"
#include "gpbase/singleton.h"
#include "gpbase/eventviewerfwd.h"
#include "gpbase/ostringstream.h"

/// gpinfra
#include "gpinfra/modelnamemap.h"
#include "gpinfra/pricingstates.h"

/// gpmodels
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/Heston_Fx.h"
#include "gpmodels/forwardmargin.h"
#include "gpmodels/Local_Functional.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class   : ARM_2IRFXSV
///	Routines: Constructor
///	Returns :
///	Action  : builds the object
///				model types contains the integer on the various model
///	
////////////////////////////////////////////////////
ARM_2IRFXSV::ARM_2IRFXSV(	
	const ARM_ModelNameMap&	modelNameMap, 
	const ARM_CurveMatrix& correlationMatrix)
:
ARM_MultiAssetsMeanReverting(&modelNameMap,&correlationMatrix)
{
	/// validate and initialize
	Validate();
}

	////////////////////////////////////////////////////
///	Class   : ARM_2IRFXSV
///	Routine : Validate
///	Returns : void
///	Action  : intialises the sub models that are linked to a model
////////////////////////////////////////////////////
void ARM_2IRFXSV::Validate( )
{
	ARM_ModelNameMap* modelNameMap	= GetModelMap();

	/// Check linked models consistency

	// Check the size of the model map

	// The model map should contain:
	// _ N IR Model
	// _ N Basis Model
	// _ (N-1) FX Model
	// = 3*N-1 Model


	bool withBasis = false;

	if(modelNameMap->size() == 3)
		withBasis = false;
    else if(modelNameMap->size() == 5)
		withBasis = true;
	else
    {
        ARMTHROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : the model should count 3 or 5 models and not" << modelNameMap->size() << ".");
    }

	
	if( !GetCorrelMatrix() )
        ARMTHROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : please provide a correlation");
	if(	GetCorrelMatrix()->rows() != 3 )
        ARMTHROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : correlation should be 3!");	

	// Check IR models

	if( !dynamic_cast<ARM_QModel1F*>( &*(*modelNameMap)[DomModel]->Model() ) )
		ARMTHROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : the domestic model has to be Q1F IR Model");

	if( !dynamic_cast<ARM_QModel1F*>( &*(*modelNameMap)[ForModel]->Model() ) )
		ARMTHROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : the foreign model has to be Q1F IR Model");

	// Check FX models

	if( !dynamic_cast<ARM_HestonModel_Fx*>( &*(*modelNameMap)[FxModel]->Model() ) )
		ARMTHROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : the FX model has to be Heston FX Model");

	// Check Basis Model

	if (withBasis)
	{
		if( !dynamic_cast<ARM_ForwardMargin*>(&*(*modelNameMap)[DomBasisModel]->Model()) )
			ARMTHROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : the domestic basis model has to be a Forward Margin Model");

		if( !dynamic_cast<ARM_ForwardMargin*>(&*(*modelNameMap)[ForBasisModel]->Model()) )
			ARMTHROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : the foreign basis model has to be a Forward Margin Model");
	}
    
	/// Check the cross currencies
	string domCcy((*modelNameMap)[0]->Model()->GetZeroCurve()->GetCurrencyUnit()->GetCcyName());

	ARM_ModelParams* modelParams= (*modelNameMap)[FxModel]->Model()->GetModelParams();
	ARM_ModelParams_Fx* fxModelParams = dynamic_cast<ARM_ModelParams_Fx*>(modelParams);		
	string fxDomCcy(fxModelParams->GetDomCurve()->GetCurrencyUnit()->GetCcyName());
	if( fxDomCcy != domCcy)
		ARMTHROW( ERR_INVALID_ARGUMENT, ARM_USERNAME << " : The FXModel has domestic curve == " << fxDomCcy << " while the domestic model ccy is == " << domCcy << ".");

	string forCcy((*modelNameMap)[ForModel]->Model()->GetZeroCurve()->GetCurrencyUnit()->GetCcyName());
	string fxForCcy(fxModelParams->GetForCurve()->GetCurrencyUnit()->GetCcyName());
	if( fxForCcy != forCcy)
		ARMTHROW( ERR_INVALID_ARGUMENT, ARM_USERNAME << " : The FXModel has foreign curve == " << fxForCcy << " while the foreign model ccy is == " << forCcy << ".");

	///Initialize the FXModel Correlation.
	ARM_EqFxBase* fxModel = dynamic_cast<ARM_EqFxBase*>(&*(*modelNameMap)[FxModel]->Model());
	fxModel->SetCorrelMatrix( *GetCorrelMatrix() );
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXSV
///	Routines: Copy constructor
///	Returns :
///	Action  : constructor
////////////////////////////////////////////////////

ARM_2IRFXSV::ARM_2IRFXSV(const ARM_2IRFXSV& rhs)
:	ARM_MultiAssetsMeanReverting(rhs)
{
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXSV
///	Routines: destructor
///	Returns :
///	Action  : destructor
////////////////////////////////////////////////////

ARM_2IRFXSV::~ARM_2IRFXSV()
{
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXSV
///	Routines: FactorCount
///	Returns : size_t
///	Action  : Return the number of factors 
/// of the model
////////////////////////////////////////////////////

size_t ARM_2IRFXSV::FactorCount() const
{
	return 5;
}

////////////////////////////////////////////////////
///	Class   : ARM_2IRFXSV
///	Routines: NeedMCIntegProcess
///	Returns : ARM_BoolVector
///	Action  : Which process should be integrated
////////////////////////////////////////////////////

ARM_BoolVector ARM_2IRFXSV::NeedMCIntegProcess() const 
{ 
	ARM_BoolVector ret(FactorCount(), false);
	
	ret[DomModel] = true;
	ret[ForModel] = true;

	return ret;
};


///////////////////////////////////////////////////
///	Class   : ARM_2IRFXSV
///	Routine : Init
///	Returns : ARM_PricingStatesPtr
///	Action  : intialise the model
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_2IRFXSV::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
	/// Validate that the payment currency is the domestic model!
    const string& domBasisModelName = GetDomBasisModel()->GetModelName();
    const string& domModelName      = GetDomModel()->GetModelName();
	if( payModelName != domBasisModelName && payModelName != domModelName)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": NP1IRNFX + FX model allow only paimenent in " + domBasisModelName + " or " + domModelName);

	ARM_NumerairePtr numeraire = GetNumeraire();
    ARM_NumMethodPtr numMethod = GetNumMethod();

    if( numMethod == ARM_NumMethodPtr(NULL) ||
        numeraire->GetType() != ARM_Numeraire::RollingCash)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +": only Rolling Cash numeraire supported by N+1IR + NFX model");

	if( numMethod != ARM_NumMethodPtr(NULL) && 
		(numMethod->GetPricingDirection() != ARM_NumMethod::GP_FWDLOOKING 
		&& numMethod->GetPricingDirection() != ARM_NumMethod::GP_FWDBCKWDLOOKING))

	   ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +": only MC or AMC supported by N+1IR + NFX model");

	/// Delegate to common code
	return ARM_MultiAssetsModel::Init( payModelName, timeInfos );
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXSV
///	Routine : IntegratedLocalDrifts
///	Returns : void
///	Action  : compute integrated relative and absolute drifts
////////////////////////////////////////////////////
void ARM_2IRFXSV::IntegratedLocalDrifts(
	const ARM_GP_Vector& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts) const
{
	size_t i,nbSteps = timeSteps.size();

	ARM_MultiAssetsModel::IntegratedLocalDrifts( timeSteps, relativeDrifts, absoluteDrifts );

	ARM_PricingModelPtr domModel = GetDomModel();
    const ARM_ModelParamsHW1FStd* const domModelMRParams    = dynamic_cast<const ARM_ModelParamsHW1FStd* const>( domModel->GetModelParams() );

	ARM_PricingModelPtr forModel = GetForModel();
	const ARM_ModelParamsHW1FStd* forModelMRParams = dynamic_cast<const ARM_ModelParamsHW1FStd* const>(forModel->GetModelParams());
	ARM_PricingModelPtr fxModel = GetFxModel();
	const ARM_Heston_ModelParams* const fxModelMRParamsHeston = dynamic_cast<const ARM_Heston_ModelParams* const>( fxModel->GetModelParams() );

	ARM_CurveModelParam fxVolParam;
	fxVolParam = fxModelMRParamsHeston->GetVolAlpha();
	
	double corrDomFor;
    double corrForFx;
    double lastTime=0.0,time;

	const ARM_CurveMatrix* const correlCurveMatrix  = GetCorrelMatrix();
	ARM_GP_Matrix correlMatrix;

    switch(GetNumeraire()->GetType())
    {
    case ARM_Numeraire::RollingCash :
        /// Diffusion proba is domestic spot i.e. neutral forward at time = next time step
	    for( i=0; i+1<nbSteps; ++i )
	    {
			correlMatrix = correlCurveMatrix->Interpolate(timeSteps[i+1]);
			corrDomFor = correlMatrix(DomModel,ForModel);
			corrForFx = correlMatrix(ForModel,FxModel);

            time = timeSteps[i+1];


            /// Qdom <-> Qdom(spot) correction on domestic process
            (*absoluteDrifts)(i,0) += ARM_ModelParamsHW1F::HW1FStateZcCovariance(domModelMRParams,domModelMRParams,lastTime,time,time,time);

			/// Qfor <-> Qdom (quanto) correction and
			/// Qdom <-> Qdom(spot) correction on foreign process
			(*absoluteDrifts)(i,1) += -corrForFx  * ARM_ModelParamsHW1F::HW1FEqFxStateCovariance(fxVolParam,forModelMRParams,lastTime,time,time)
									  +corrDomFor * ARM_ModelParamsHW1F::HW1FStateZcCovariance(forModelMRParams,domModelMRParams,lastTime,time,time,time);

            lastTime = time;
	    }
        break;

    default:
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": only RollingCash numeraire supported by N+1 IR + N FX model");
    }
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXSV
///	Routine : ModelStateLocalVariancesAndStdDev
///	Returns : void
///	Action  : computes model covariances
////////////////////////////////////////////////////

void ARM_2IRFXSV::ModelStateLocalVariances( 
const ARM_GP_Vector& timeSteps,
ARM_MatrixVector& localVariances) const
{
	size_t nbSteps	= timeSteps.size();
	size_t factorNb = FactorCount();

	localVariances.resize(nbSteps-1);

	for(size_t i=0;i<nbSteps-1;++i)
		localVariances[i] = new ARM_GP_Matrix(factorNb,factorNb);
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXSV
///	Routines: ModelStateLocalVariancesAndStdDev
///	Returns : void
///	Action  : computes the local variance and std dev 
////////////////////////////////////////////////////
void ARM_2IRFXSV::ModelStateLocalVariancesAndStdDev( const ARM_GP_Vector& timeSteps)
{
	ARM_PricingModel::ModelStateLocalVariancesAndStdDev(timeSteps);
}

////////////////////////////////////////////////////
///	Class   : ARM_2IRFXSV
///	Routines: NumMethodStateLocalVariances
///	Returns : void
///	Action  : NumMethodStateLocalVariances
////////////////////////////////////////////////////
void ARM_2IRFXSV::NumMethodStateLocalVariances( 
		const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const
{
	int i;
	size_t factorNb	= FactorCount();
	const ARM_CurveMatrix* const correlCurveMatrix  = GetCorrelMatrix();
	ARM_GP_Matrix correlMatrix ;
	const ARM_ModelNameMap* const modelMap = GetModelMap();
	size_t nbSteps	= timeSteps.size();
    double step		= timeSteps[0],nextStep;

#if defined(__GP_STRICT_VALIDATION)
	if( localVariances.size()!= 0 ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "localVariances.size() != 0" );
#endif

	localVariances.resize(nbSteps-1);

    ARM_PricingModelPtr domModel    = GetDomModel();
	const ARM_ModelParamsHW1FStd*  domModelParams;
	domModelParams = static_cast<const ARM_ModelParamsHW1FStd* const>( domModel->GetModelParams() );
	
	const ARM_ModelParamsHW1FStd*  fxModelParams;
	const ARM_Heston_ModelParams*  fxParams;
	ARM_CurveModelParam	   fxVol;

	ARM_PricingModelPtr forModel = GetForModel();
	const ARM_ModelParamsHW1FStd*  forModelParams;
	forModelParams = dynamic_cast<const ARM_ModelParamsHW1FStd*>(forModel->GetModelParams());
	ARM_PricingModelPtr fxModel = GetFxModel();
	fxModelParams = dynamic_cast<const ARM_ModelParamsHW1FStd*>( fxModel->GetModelParams() );
	fxParams = dynamic_cast<const ARM_Heston_ModelParams*>( fxModelParams );
	fxVol = fxParams->GetVolAlpha();
	double alpha = fxParams->Alpha();

	double varFx,varZcDom,varZcFor;
    double covarZcDomZcFor,covarZcDomFx,covarZcForFx;
    double covarXDomZcDom,covarXDomZcFor,covarXDomFx;
    double covarXForZcDom,covarXForZcFor,covarXForFx;

	ARM_MatrixVector fxLocalVariances;
	fxLocalVariances.resize(nbSteps-1);

	for(i=0;i<nbSteps-1;++i)
	{
		nextStep			= timeSteps[i+1];

		correlMatrix = correlCurveMatrix->Interpolate(step);

		localVariances[i]	= new ARM_GP_TriangularMatrix(factorNb,0.0);

		/// Domestic variances
		(*localVariances[i])(DomModel,DomModel)   = domModelParams->StateLocalVariance(step,nextStep,nextStep);

		/// Foreign variances
		(*localVariances[i])(ForModel,ForModel)   = forModelParams->StateLocalVariance(step,nextStep,nextStep);

		/// Domestic & Foreign covariances
		(*localVariances[i])(DomModel,ForModel) = correlMatrix(DomModel,ForModel) * ARM_ModelParamsHW1F::HW1FStateCovariance(domModelParams, forModelParams,step,nextStep,nextStep);
		(*localVariances[i])(ForModel,DomModel) = (*localVariances[i])(DomModel,ForModel);

		/// Forex variances
		varFx       = alpha*fxModelParams->StateLocalVariance(step,nextStep,nextStep);
		varZcDom    = ARM_ModelParamsHW1F::HW1FZcCovariance( domModelParams, domModelParams, step,nextStep,nextStep );
		varZcFor    = ARM_ModelParamsHW1F::HW1FZcCovariance( forModelParams, forModelParams, step,nextStep,nextStep );

		covarZcDomZcFor = ARM_ModelParamsHW1F::HW1FZcCovariance( domModelParams, forModelParams, step,nextStep,nextStep );
		covarZcDomFx    = ARM_ModelParamsHW1F::HW1FEqFxZcCovariance( fxVol, domModelParams, step,nextStep,nextStep );
		covarZcForFx    = ARM_ModelParamsHW1F::HW1FEqFxZcCovariance( fxVol, forModelParams, step,nextStep,nextStep );

		(*localVariances[i])(FxModel,FxModel) =
			varFx + varZcDom + varZcFor
			- 2*correlMatrix(DomModel,ForModel)*covarZcDomZcFor
			- 2*correlMatrix(DomModel,FxModel)*covarZcDomFx
			+ 2*correlMatrix(ForModel,FxModel)*covarZcForFx;

		double varLimit = K_NEW_DOUBLE_TOL*K_NEW_DOUBLE_TOL;
		if((*localVariances[i])(FxModel,FxModel)<varLimit)
			(*localVariances[i])(FxModel,FxModel) = varLimit; // numerical noise

		/// Domestic & Forex covariances
		covarXDomZcDom = ARM_ModelParamsHW1F::HW1FStateZcCovariance( domModelParams, domModelParams, step,nextStep,nextStep,nextStep );
		covarXDomZcFor = ARM_ModelParamsHW1F::HW1FStateZcCovariance( domModelParams, forModelParams, step,nextStep,nextStep,nextStep );
		covarXDomFx    = ARM_ModelParamsHW1F::HW1FEqFxStateCovariance( fxVol, domModelParams, step,nextStep,nextStep );

		(*localVariances[i])(DomModel,FxModel) =
			- covarXDomZcDom
			+ correlMatrix(DomModel,ForModel)*covarXDomZcFor
			+ correlMatrix(DomModel,FxModel)*covarXDomFx;

		(*localVariances[i])(FxModel,DomModel) = (*localVariances[i])(DomModel,FxModel);

		/// Foreign & Forex covariances
		covarXForZcDom = ARM_ModelParamsHW1F::HW1FStateZcCovariance( forModelParams, domModelParams, step,nextStep,nextStep,nextStep );
		covarXForZcFor = ARM_ModelParamsHW1F::HW1FStateZcCovariance( forModelParams, forModelParams, step,nextStep,nextStep,nextStep );
		covarXForFx    = ARM_ModelParamsHW1F::HW1FEqFxStateCovariance( fxVol, forModelParams, step,nextStep,nextStep );

		(*localVariances[i])(ForModel,FxModel) =
			covarXForZcFor
			- correlMatrix(DomModel,ForModel)*covarXForZcDom
			+ correlMatrix(ForModel,FxModel)*covarXForFx;

		(*localVariances[i])(FxModel,ForModel) = (*localVariances[i])(ForModel,FxModel);

		(*localVariances[i])(SpotHeston,SpotHeston) = 1.0;
		(*localVariances[i])(VarHeston,VarHeston) = 1.0;

		fxLocalVariances[i] = static_cast< ARM_GP_Matrix* >(localVariances[i]->Clone());

		step=nextStep;
	}

	(*GetModelMap())[FxModel]->Model()->SetModelStateLocalVars(fxLocalVariances);
}

////////////////////////////////////////////////////
///	Class   : ARM_2IRFXSV
///	Routines: NumMethodStateGlobalVariances
///	Returns : void
///	Action  : NumMethodStateGlobalVariances
////////////////////////////////////////////////////
void ARM_2IRFXSV::NumMethodStateGlobalVariances( 
		const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& globalVariances ) const
{
	size_t nbSteps	= timeSteps.size();
	size_t factorNb = FactorCount();

	globalVariances.resize(nbSteps);

	for(size_t i=0;i<nbSteps;++i)
		globalVariances[i] = new ARM_GP_Matrix(factorNb,factorNb);
}


////////////////////////////////////////////////////
///	Class  : ARM_2IRFXSV
///	Routine: TreeStatesToModelStates
///	Returns: void
///	Action : Convert the num method states to model states
////////////////////////////////////////////////////
void ARM_2IRFXSV::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{   
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "ARM_1IRFXModel::TreeStatesToModelStates is not implemented!");
}


////////////////////////////////////////////////////
///	Class  : ARM_2IRFXSV
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : Convert the num method states to model states
////////////////////////////////////////////////////
void ARM_2IRFXSV::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
	int i;

	ARM_GP_MatrixPtr numStates = states->GetNumMethodStates();
    ARM_GP_MatrixPtr modelStates = states->GetModelStates();

	size_t nbStates = numStates->cols();

	(*(GetModelMap()))[FxModel]->Model()->MCModelStatesFromToNextTime(states,timeIndex);

	for(i=0;i<nbStates;++i)
	{
		(*modelStates)(DomModel,i)  =  (*numStates)(DomModel,i);
		(*modelStates)(ForModel,i)  =  (*numStates)(ForModel,i);
	}

}


////////////////////////////////////////////////////
///	Class  : ARM_2IRFXSV
///	Routine: CalibrateFunctional
///	Returns: void
///	Action : Convert the num method states to model states
////////////////////////////////////////////////////

void ARM_2IRFXSV::CalibrateFunctional (
		const ARM_GP_Vector& ResetTimes,
		vector<ARM_DensityFunctor*> densities,
		int sizeGrid,
		double nbStdDev)
{
	int i;
	ARM_LocalFunctionalPtr localFunctional(new ARM_LocalFunctional(sizeGrid,nbStdDev));

	ARM_PricingModelPtr fxModel = GetFxModel();
	ARM_HestonModel_Fx* fxModelHeston = dynamic_cast<ARM_HestonModel_Fx*>(&*fxModel);

	for (i = 0; i < ResetTimes.size(); ++i)
	{	
		double resetTime = ResetTimes[i]/K_YEAR_LEN;
		double fwd = fxModelHeston->ComputeFwdAtTime(ResetTimes[i]);
		localFunctional->CalibrateHestonFx(ResetTimes[i],fwd,fxModelHeston,(*densities[i]));
	}
	fxModelHeston->SetLocalFunctional(localFunctional);
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

