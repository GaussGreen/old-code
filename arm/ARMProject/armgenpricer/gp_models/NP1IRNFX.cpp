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
 *	\date Jully 2006
 */

//#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/NP1IRNFX.h"

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
#include "gpmodels/Q1F_Fx.h"
#include "gpmodels/Q1F.h"
#include "gpmodels/ModelParamsQ1F.h"
#include "gpmodels/forwardmargin.h"
#include "gpmodels/Local_Functional.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class   : ARM_NP1IRNFXModel
///	Routines: Constructor
///	Returns :
///	Action  : builds the object
///				model types contains the integer on the various model
///	
////////////////////////////////////////////////////
ARM_NP1IRNFXModel::ARM_NP1IRNFXModel(	
	const ARM_ModelNameMap&	modelNameMap, 
	const ARM_CurveMatrix& correlationMatrix)
:
ARM_MultiAssetsMeanReverting(&modelNameMap,&correlationMatrix)
{
	/// validate and initialize
	Validate();
}

	////////////////////////////////////////////////////
///	Class   : ARM_NP1IRNFXModel
///	Routine : Validate
///	Returns : void
///	Action  : intialises the sub models that are linked to a model
////////////////////////////////////////////////////
void ARM_NP1IRNFXModel::Validate( )
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

	if((modelNameMap->size()+1)%3 == 0)
	{
		itsNbCcy = (modelNameMap->size()+1)/3;
		withBasis = false;
	}
    else if((modelNameMap->size()+1)%2 == 0)
	{
		itsNbCcy = (modelNameMap->size()+1)/2;
		withBasis = true;
	}
	else
    {
        ARMTHROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : the model should count 3*N-1 model and not" << modelNameMap->size() << ".");
    }

	
	if( !GetCorrelMatrix() )
        ARMTHROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : please provide a correlation");
	if(	GetCorrelMatrix()->rows() != 2*itsNbCcy-1 )
        ARMTHROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : correlation should be "<< 2*itsNbCcy-1 << "x"<< 2*itsNbCcy-1 <<"!");	

	// Check IR models
	
	int i;

	for (i = 0; i < itsNbCcy; ++i)
	{
		if( !dynamic_cast<ARM_QModel1F*>( &*(*modelNameMap)[i]->Model() ) )
			ARMTHROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : the "<< i <<"th model has to be Q1F IR Model");
	}

	// Check FX models

	for (i = 0; i < itsNbCcy-1; ++i)
	{
		if( !dynamic_cast<ARM_QModel1F_Fx*>( &*(*modelNameMap)[itsNbCcy+i]->Model() ) )
			ARMTHROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : the "<< itsNbCcy+i <<"th model has to be Q1F IR Model");
	}

	// Check Basis Model

	if (withBasis)
	{
		for (i = 0; i < itsNbCcy; ++i)	
		{
			if( !dynamic_cast<ARM_ForwardMargin*>(&*(*modelNameMap)[2*itsNbCcy-1+i]->Model()) )
				ARMTHROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : the " << 2*itsNbCcy-1+i <<"th model has to be a Forward Margin Model");
		}
	}
    
	/// Check the cross currencies
	string domCcy((*modelNameMap)[0]->Model()->GetZeroCurve()->GetCurrencyUnit()->GetCcyName());

	for (i = 0; i < itsNbCcy-1; ++i) {
		ARM_ModelParams* modelParams= (*modelNameMap)[itsNbCcy+i]->Model()->GetModelParams();
		ARM_ModelParams_Fx* fxModelParams = dynamic_cast<ARM_ModelParams_Fx*>(modelParams);		
		string fxDomCcy(fxModelParams->GetDomCurve()->GetCurrencyUnit()->GetCcyName());
		if( fxDomCcy != domCcy)
			ARMTHROW( ERR_INVALID_ARGUMENT, ARM_USERNAME << " : The " << i << "th FXModel has domestic curve == " << fxDomCcy << " while the 0th IR model ccy is == " << domCcy << ".");

		string forCcy((*modelNameMap)[1+i]->Model()->GetZeroCurve()->GetCurrencyUnit()->GetCcyName());
		string fxForCcy(fxModelParams->GetForCurve()->GetCurrencyUnit()->GetCcyName());
		if( fxForCcy != forCcy)
		{
			ARMTHROW( ERR_INVALID_ARGUMENT, ARM_USERNAME << " : The " << i << "th FXModel has foreign curve == " << fxForCcy << " while the " << i+1 << "th IR model ccy is == " << forCcy << ".");
		}
	}

	///Initialize the FXModel Correlation.
	for (i = 0; i < itsNbCcy-1; ++i)
	{
		ARM_EqFxBase* fxModel = dynamic_cast<ARM_EqFxBase*>(&*(*modelNameMap)[itsNbCcy+i]->Model());
		fxModel->SetCorrelMatrix( *GetCorrelMatrix() );
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_NP1IRNFXModel
///	Routines: Copy constructor
///	Returns :
///	Action  : constructor
////////////////////////////////////////////////////

ARM_NP1IRNFXModel::ARM_NP1IRNFXModel(const ARM_NP1IRNFXModel& rhs)
:	ARM_MultiAssetsMeanReverting(rhs),
itsNbCcy(rhs.itsNbCcy)
{
}


////////////////////////////////////////////////////
///	Class   : ARM_NP1IRNFXModel
///	Routines: destructor
///	Returns :
///	Action  : destructor
////////////////////////////////////////////////////

ARM_NP1IRNFXModel::~ARM_NP1IRNFXModel()
{
}


////////////////////////////////////////////////////
///	Class   : ARM_NP1IRNFXModel
///	Routines: FactorCount
///	Returns : size_t
///	Action  : Return the number of factors 
/// of the model
////////////////////////////////////////////////////

size_t ARM_NP1IRNFXModel::FactorCount() const
{
	return 2*itsNbCcy-1;
}

////////////////////////////////////////////////////
///	Class   : ARM_NP1IRNFXModel
///	Routines: NeedMCIntegProcess
///	Returns : ARM_BoolVector
///	Action  : Which process should be integrated
////////////////////////////////////////////////////

ARM_BoolVector ARM_NP1IRNFXModel::NeedMCIntegProcess() const 
{ 
	ARM_BoolVector ret(FactorCount(), false);
	
	int i;
	for (i = 0; i < itsNbCcy; ++i)
		ret[i] = true;

	return ret;
};

////////////////////////////////////////////////////
///	Class   : ARM_NP1IRNFXModel
///	Routines: ForModel
///	Returns : ARM_BoolVector
///	Action  : Return the ith foreign model
////////////////////////////////////////////////////

ARM_PricingModelPtr ARM_NP1IRNFXModel::ForModel(int i) const
{
	if ((i < 0) || (i > itsNbCcy-1))
	{
		ARMTHROW( ERR_INVALID_ARGUMENT, "the " << i << "th foreign model is not available in N+1 IR + N Fx.")
	}

	return (*(GetModelMap()))[1+i]->Model();
}

////////////////////////////////////////////////////
///	Class   : ARM_NP1IRNFXModel
///	Routines: FxModel
///	Returns : ARM_BoolVector
///	Action  : Return the ith FX model
////////////////////////////////////////////////////

ARM_PricingModelPtr ARM_NP1IRNFXModel::FxModel(int i) const
{
	if ((i < 0) || (i > itsNbCcy-1))
	{
		ARMTHROW( ERR_INVALID_ARGUMENT, "the " << i << "th Fx model is not available in N+1 IR + N Fx.")
	}

	return (*(GetModelMap()))[itsNbCcy+i]->Model();
}

///////////////////////////////////////////////////
///	Class   : ARM_NP1IRNFXModel
///	Routine : Init
///	Returns : ARM_PricingStatesPtr
///	Action  : intialise the model
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_NP1IRNFXModel::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
	/// Validate that the payment currency is the domestic model!
    const string& domBasisModelName = DomBasisModel()->GetModelName();
    const string& domModelName      = DomModel()->GetModelName();
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

	ARM_QModel1F_Fx* fxModel;

	int i;
	for (i = 0; i < itsNbCcy-1; ++i)
	{
		fxModel = dynamic_cast<ARM_QModel1F_Fx*>(&*((*GetModelMap())[itsNbCcy+i]->Model()));
		fxModel->SetIsForexDiffusion(true);
	}

	/// Delegate to common code
	return ARM_MultiAssetsModel::Init( payModelName, timeInfos );
}


////////////////////////////////////////////////////
///	Class   : ARM_NP1IRNFXModel
///	Routine : IntegratedLocalDrifts
///	Returns : void
///	Action  : compute integrated relative and absolute drifts
////////////////////////////////////////////////////
void ARM_NP1IRNFXModel::IntegratedLocalDrifts(
	const std::vector<double>& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts) const
{
	size_t i,j,nbSteps = timeSteps.size();

	ARM_MultiAssetsModel::IntegratedLocalDrifts( timeSteps, relativeDrifts, absoluteDrifts );

	ARM_PricingModelPtr domModel = DomModel();
    const ARM_ModelParamsHW1FStd* const domModelMRParams    = dynamic_cast<const ARM_ModelParamsHW1FStd* const>( domModel->GetModelParams() );

    vector<const ARM_ModelParamsHW1FStd*>  forModelMRParams;
	vector<const ARM_ModelParamsHW1FStd*>  fxModelMRParams;

	for (i = 0; i < itsNbCcy-1; ++i)
	{
		ARM_PricingModelPtr forModel = ForModel(i);
		forModelMRParams.push_back(dynamic_cast<const ARM_ModelParamsHW1FStd* const>(forModel->GetModelParams()));
		ARM_PricingModelPtr fxModel = FxModel(i);
		fxModelMRParams.push_back(dynamic_cast<const ARM_ModelParamsHW1FStd* const>( fxModel->GetModelParams() ));
	}
	
	vector<double> corrDomFor;
    vector<double> corrForFx;
    double lastTime=0.0,time;

	const ARM_CurveMatrix* const correlCurveMatrix  = GetCorrelMatrix();
	ARM_GP_Matrix correlMatrix;

	// Correlation interpolation
	for( i=0; i+1<nbSteps; ++i )
	{
		correlMatrix = correlCurveMatrix->Interpolate(timeSteps[i+1]);
		for (j=0; j < itsNbCcy-1; ++j)
		{
			corrDomFor.push_back(correlMatrix(0,1+j));
			corrForFx.push_back(correlMatrix(1+j,itsNbCcy+j));
		}
	}

    switch(GetNumeraire()->GetType())
    {
    case ARM_Numeraire::RollingCash :
        /// Diffusion proba is domestic spot i.e. neutral forward at time = next time step
	    for( i=0; i+1<nbSteps; ++i )
	    {
            time = timeSteps[i+1];


            /// Qdom <-> Qdom(spot) correction on domestic process
            (*absoluteDrifts)(i,0) += ARM_ModelParamsHW1F::HW1FStateZcCovariance(domModelMRParams,domModelMRParams,lastTime,time,time,time);

			/// Qfor <-> Qdom (quanto) correction and
			/// Qdom <-> Qdom(spot) correction on foreign process
			for (j = 0; j < itsNbCcy-1; ++j)
			{
				(*absoluteDrifts)(i,1+j) += -corrForFx[j]  * ARM_ModelParamsHW1F::HW1FStateCovariance(fxModelMRParams[j],forModelMRParams[j],lastTime,time,time)
											 +corrDomFor[j] * ARM_ModelParamsHW1F::HW1FStateZcCovariance(forModelMRParams[j],domModelMRParams,lastTime,time,time,time);
			}

            lastTime = time;
	    }
        break;

    default:
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": only RollingCash numeraire supported by N+1 IR + N FX model");
    }
}


////////////////////////////////////////////////////
///	Class   : ARM_NP1IRNFXModel
///	Routine : ModelStateLocalVariancesAndStdDev
///	Returns : void
///	Action  : computes model covariances
////////////////////////////////////////////////////

void ARM_NP1IRNFXModel::ModelStateLocalVariances( 
const std::vector<double>& timeSteps,
ARM_MatrixVector& localVariances) const
{
	size_t nbSteps	= timeSteps.size();
	size_t factorNb = FactorCount();

	localVariances.resize(nbSteps-1);

	for(size_t i=0;i<nbSteps-1;++i)
		localVariances[i] = new ARM_GP_Matrix(factorNb,factorNb);
}


////////////////////////////////////////////////////
///	Class   : ARM_NP1IRNFXModel
///	Routines: ModelStateLocalVariancesAndStdDev
///	Returns : void
///	Action  : computes the local variance and std dev 
////////////////////////////////////////////////////
void ARM_NP1IRNFXModel::ModelStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps)
{
	ARM_PricingModel::ModelStateLocalVariancesAndStdDev(timeSteps);
}

////////////////////////////////////////////////////
///	Class   : ARM_NP1IRNFXModel
///	Routines: NumMethodStateLocalVariances
///	Returns : void
///	Action  : NumMethodStateLocalVariances
////////////////////////////////////////////////////
void ARM_NP1IRNFXModel::NumMethodStateLocalVariances( 
		const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const
{
	int i,j,k;
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

    ARM_PricingModelPtr domModel    = DomModel();
	vector<const ARM_ModelParamsHW1FStd*>  irModelParams;
	irModelParams.push_back(static_cast<const ARM_ModelParamsHW1FStd* const>( domModel->GetModelParams() ));
	
	vector<const ARM_ModelParamsHW1FStd*>  fxModelParams;
	vector<const ARM_ModelParamsQ1F*>	   fxParams;
	vector<const ARM_CurveModelParam*>	   fxVols;

	for (i = 0; i < itsNbCcy-1; ++i)
	{
		ARM_PricingModelPtr forModel = ForModel(i);
		irModelParams.push_back(dynamic_cast<const ARM_ModelParamsHW1FStd* const>(forModel->GetModelParams()));
		ARM_PricingModelPtr fxModel = FxModel(i);
		fxModelParams.push_back(dynamic_cast<const ARM_ModelParamsHW1FStd* const>( fxModel->GetModelParams() ));
		fxParams.push_back(dynamic_cast<const ARM_ModelParamsQ1F* const>( fxModelParams[i] ));
		fxVols.push_back(dynamic_cast<const ARM_CurveModelParam*>( &fxParams[i]->GetModelParam(fxParams[i]->GetVolatilityType())));
	}

	double covarXIRjZcDom0;
	double covarXIRjZcFork;
	double covarXIRjFXk;
	double covarZcDom0ZcDom0;	
	double covarZcForjZcFork;
	double covarFxjFxk;
	double covarZcForjFXk;
	double covarZcForkFXj;
	double covarZcDom0FXj;
	double covarZcDom0FXk;
	double covarZcDom0ZcForj;
	double covarZcDom0ZcFork;

	vector<ARM_MatrixVector> fxLocalVariances(itsNbCcy-1);

	for (i=0;i<itsNbCcy-1;++i)
		fxLocalVariances[i].resize(nbSteps-1);


	for(i=0;i<nbSteps-1;++i)
	{
		CC_Ostringstream os;
		os << "TimeStep = " << i << std::endl;
		ARM_TheEventViewer.Instance()->AddToMessage(os.str());
		
		nextStep			= timeSteps[i+1];
		localVariances[i]	= new ARM_GP_TriangularMatrix(factorNb,0.0);
		correlMatrix = correlCurveMatrix->Interpolate(nextStep);

		// IR/IR
		for (j = 0; j < itsNbCcy; ++j)
			for (k = j; k < itsNbCcy; ++k)
			{
				(*localVariances[i])(j,k) =  (*localVariances[i])(k,j) = ARM_ModelParamsHW1F::HW1FStateCovariance( irModelParams[j], irModelParams[k], step,nextStep,nextStep )*correlMatrix(j,k);

				CC_Ostringstream os;
				os << "IR" << j << "/IR" << k << "=" << (*localVariances[i])(j,k) << std::endl;
				ARM_TheEventViewer.Instance()->AddToMessage(os.str());
			}

		// IR/FX
		for (j = 0; j < itsNbCcy; ++j)
			for(k = 0; k < itsNbCcy-1; ++k)
			{
				covarXIRjZcDom0 = ARM_ModelParamsHW1F::HW1FStateZcCovariance( irModelParams[j], irModelParams[0], step,nextStep,nextStep,nextStep )*correlMatrix(j,0);
				covarXIRjZcFork = ARM_ModelParamsHW1F::HW1FStateZcCovariance( irModelParams[j], irModelParams[1+k], step,nextStep,nextStep,nextStep )*correlMatrix(j,1+k);
				covarXIRjFXk    = ARM_ModelParamsHW1F::HW1FStateCovariance( irModelParams[j], fxModelParams[k], step,nextStep,nextStep )*correlMatrix(j,itsNbCcy+k);

				(*localVariances[i])(j,itsNbCcy+k) = (*localVariances[i])(itsNbCcy+k,j) =
					- covarXIRjZcDom0
					+ covarXIRjZcFork
					+ covarXIRjFXk;

				CC_Ostringstream os;
				os << "IR" << j << "/FX" << k << "=" << (*localVariances[i])(itsNbCcy+k,j) << std::endl;
				ARM_TheEventViewer.Instance()->AddToMessage(os.str());
			}

		// FX/FX
		for(j = 0; j < itsNbCcy-1; ++j)
			for(k = j; k < itsNbCcy-1; ++k)
			{
				covarZcDom0ZcDom0 = ARM_ModelParamsHW1FStd::HW1FZcCovariance(irModelParams[0],irModelParams[0],step,nextStep,nextStep,nextStep);
				covarZcForjZcFork = ARM_ModelParamsHW1FStd::HW1FZcCovariance(irModelParams[1+j],irModelParams[1+k],step,nextStep,nextStep,nextStep)*correlMatrix(1+j,1+k);
				covarFxjFxk = ARM_ModelParamsHW1FStd::HW1FStateCovariance(fxModelParams[j],fxModelParams[k],step,nextStep,nextStep)*correlMatrix(itsNbCcy+j,itsNbCcy+k);

				covarZcForjFXk = ARM_ModelParamsHW1FStd::HW1FEqFxZcCovariance(*(fxVols[k]),irModelParams[1+j],step,nextStep,nextStep)*correlMatrix(1+j,itsNbCcy+k);
				covarZcForkFXj = ARM_ModelParamsHW1FStd::HW1FEqFxZcCovariance(*(fxVols[j]),irModelParams[1+k],step,nextStep,nextStep)*correlMatrix(1+k,itsNbCcy+j);

				covarZcDom0FXj = ARM_ModelParamsHW1FStd::HW1FEqFxZcCovariance(*(fxVols[j]),irModelParams[0],step,nextStep,nextStep)*correlMatrix(0,itsNbCcy+j);
				covarZcDom0FXk = ARM_ModelParamsHW1FStd::HW1FEqFxZcCovariance(*(fxVols[k]),irModelParams[0],step,nextStep,nextStep)*correlMatrix(0,itsNbCcy+k);

				covarZcDom0ZcForj = ARM_ModelParamsHW1FStd::HW1FZcCovariance(irModelParams[0],irModelParams[1+j],step,nextStep,nextStep,nextStep)*correlMatrix(0,1+j);
				covarZcDom0ZcFork = ARM_ModelParamsHW1FStd::HW1FZcCovariance(irModelParams[0],irModelParams[1+k],step,nextStep,nextStep,nextStep)*correlMatrix(0,1+k);

				(*localVariances[i])(itsNbCcy+j,itsNbCcy+k) = (*localVariances[i])(itsNbCcy+k,itsNbCcy+j) =
					covarZcDom0ZcDom0 
					+ covarZcForjZcFork 
					+ covarFxjFxk 

					- covarZcDom0FXj 
					- covarZcDom0FXk 
					
					+ covarZcForjFXk 
					+ covarZcForkFXj 
					
					- covarZcDom0ZcForj 
					- covarZcDom0ZcFork;

				CC_Ostringstream os;
				os << "FX" << j << "/FX" << k << "=" << (*localVariances[i])(itsNbCcy+j,itsNbCcy+k) << std::endl;
				ARM_TheEventViewer.Instance()->AddToMessage(os.str());

			}

		for (j=0;j<itsNbCcy-1;++j)
			fxLocalVariances[j][i] = static_cast< ARM_GP_Matrix* >(localVariances[i]->Clone());

		step=nextStep;
	}

	for (i=0;i<itsNbCcy-1;++i)
		(*GetModelMap())[itsNbCcy+i]->Model()->SetModelStateLocalVars(fxLocalVariances[i]);
}

////////////////////////////////////////////////////
///	Class   : ARM_NP1IRNFXModel
///	Routines: NumMethodStateGlobalVariances
///	Returns : void
///	Action  : NumMethodStateGlobalVariances
////////////////////////////////////////////////////
void ARM_NP1IRNFXModel::NumMethodStateGlobalVariances( 
		const std::vector<double>& timeSteps,
		ARM_MatrixVector& globalVariances ) const
{
	size_t nbSteps	= timeSteps.size();
	size_t factorNb = FactorCount();

	globalVariances.resize(nbSteps);

	for(size_t i=0;i<nbSteps;++i)
	{
	
		globalVariances[i] = new ARM_GP_Matrix(factorNb,factorNb);
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_NP1IRNFXModel
///	Routine: TreeStatesToModelStates
///	Returns: void
///	Action : Convert the num method states to model states
////////////////////////////////////////////////////
void ARM_NP1IRNFXModel::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{   
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "ARM_1IRFXModel::TreeStatesToModelStates is not implemented!");
}


////////////////////////////////////////////////////
///	Class  : ARM_NP1IRNFXModel
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : Convert the num method states to model states
////////////////////////////////////////////////////
void ARM_NP1IRNFXModel::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
	int i, j;

	ARM_GP_MatrixPtr numStates = states->GetNumMethodStates();
    ARM_GP_MatrixPtr modelStates = states->GetModelStates();

	size_t nbStates = numStates->cols();

	for (i = 0; i < itsNbCcy-1; ++i)
		(*(GetModelMap()))[itsNbCcy+i]->Model()->MCModelStatesFromToNextTime(states,timeIndex);

	for (i = 0; i < itsNbCcy; ++i)
		for(j=0;j<nbStates;++j)
			(*modelStates)(i,j)  =  (*numStates)(i,j);

}


////////////////////////////////////////////////////
///	Class  : ARM_NP1IRNFXModel
///	Routine: CalibrateFunctional
///	Returns: void
///	Action : Convert the num method states to model states
////////////////////////////////////////////////////

void ARM_NP1IRNFXModel::CalibrateFunctional (
		const std::vector<double>& ResetTimes,
		vector<ARM_DensityFunctor*> densities,
		int nbRows,
		int nbCols,
		int sizeGrid,
		double nbStdDev,
		bool rescaling)
{
	int i,j;

	if (nbCols != (itsNbCcy-1))
		ARMTHROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + "ARM_NP1IRNFXModel::CalibrateFunctional the number of columns of the density functor matrix should be " << itsNbCcy << ".")

	ARM_PricingModelPtr domModel = DomModel();
	ARM_QModel1F* domModelQ1F = dynamic_cast<ARM_QModel1F*>(&*domModel);

	for (i = 0; i < nbCols; ++i)
	{
		ARM_LocalFunctionalPtr localFunctional(new ARM_LocalFunctional(sizeGrid,nbStdDev));
		ARM_PricingModelPtr forModel = ForModel(i);
		ARM_QModel1F* forModelQ1F = dynamic_cast<ARM_QModel1F*>(&*forModel);
		ARM_PricingModelPtr fxModel = FxModel(i);
		ARM_QModel1F_Fx* fxModelQ1F = dynamic_cast<ARM_QModel1F_Fx*>(&*fxModel);
		ARM_ModelParamsQ1F* fxModelParams = dynamic_cast<ARM_ModelParamsQ1F*>(fxModel->GetModelParams());

		for (j = 0; j < nbRows; ++j)
		{	
			double resetTime = ResetTimes[j]/K_YEAR_LEN;
			double fwd = fxModelQ1F->ComputeFwdAtTime(ResetTimes[j]);
			double fwdZcFwdFxCov;
			double vol= fxModelQ1F->VarianceFwdFx(0,ResetTimes[j],ResetTimes[j],domModelQ1F,forModelQ1F,fxModelParams,false,fwdZcFwdFxCov);
			vol = sqrt(vol/resetTime); 
			localFunctional->Calibrate(ResetTimes[j],fwd,vol,(*densities[j*nbCols+i]),rescaling);
		}

		fxModelQ1F->SetLocalFunctional(localFunctional);
	}
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

