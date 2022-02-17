/*!
 *
 * Copyright (c) CM CIB January 2005 Paris
 *
 *	\file 1IRFXModel.cpp
 *
 *  \brief 1 Interest Rates model + FX hybrid model
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date Jully 2006
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/1IRFXModel.h"

/// gpmodels
#include "gpmodels/ModelParamsQ1F.h"
#include "gpmodels/Q1F.h"
#include "gpmodels/Smiled_FX.h"
#include "gpmodels/forwardmargin.h"
#include "gpmodels/typedef.h"
#include "gpmodels/2IRFXModel.h"

/// gpbase
#include "gpbase/checkarg.h"
#include "gpbase/stringmanip.h"
#include "gpbase/gpmatrixlinalg.h"
#include "gpbase/eventviewerfwd.h"

/// gpinfra
#include "gpinfra/modelnamemap.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/discretisationscheme.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/pricingcontext.h"

/// gpnumlib
#include "gpnumlib/random.h"

/// gpnummethods
#include "gpnummethods/treebase.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class   : ARM_1IRFXModel
///	Routines: Constructor
///	Returns :
///	Action  : builds the object
///				model types contains the integer on the various model
///	
////////////////////////////////////////////////////
ARM_1IRFXModel::ARM_1IRFXModel(	const ARM_ModelNameMap&	modelNameMap, 
	const ARM_CurveMatrix& correlationMatrix,
	ARM_2IRFXModelPtr Model2IRFX )
:
ARM_MultiAssetsMeanReverting(&modelNameMap,&correlationMatrix),
itsTimeSteps(0),
itsEigenValues(0),
itsEigenVectors(0),
itsFwdFxDrifts(0),
itsHWRelDrifts(0),
itsIntAbsDrifts(0),
its2IRFXModel(CreateClonedPtr<ARM_2IRFXModel>(&*Model2IRFX))

{
	/// validate and initialize
	Validate();

	CalibrateCorrelation(its2IRFXModel);
}

	////////////////////////////////////////////////////
///	Class   : ARM_1IRFXModel
///	Routine : InitSubModels
///	Returns : void
///	Action  : intialises the sub models that are linked to a model
////////////////////////////////////////////////////
void ARM_1IRFXModel::Validate( )
{
	ARM_ModelNameMap* modelNameMap	= GetModelMap();
	/// Check linked models consistency
    if(modelNameMap->size() < NbModels)
    {
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : not enough models to built an hybrid 1IR+FX with basis model");
    }
    else if(modelNameMap->size() > NbModels)
    {
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : too much models to built an hybrid 1IR+FX with basis model");
    }

	if( !GetCorrelMatrix() )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : please provide a correlation");
	if(	GetCorrelMatrix()->rows() < 2 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : correlation should be at least 2x2!");

	/// check model types
	if( !dynamic_cast<ARM_QModel1F*>( &*(*modelNameMap)[DomModel]->Model() ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : dom model has to be Q1F IR Model");

	if(!dynamic_cast<ARM_SmiledModel_Fx*>( &*(*modelNameMap)[FxModel]->Model() ))
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : fx model has to be a Smiled Model");
	
	if( !dynamic_cast<ARM_ForwardMargin*>(&*(*modelNameMap)[DomBasisModel]->Model()) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : dom basis model has to be a Forward Margin Model");


	/// check the cross
	ARM_ModelParams* modelParams= (*modelNameMap)[FxModel]->Model()->GetModelParams();
	ARM_ModelParamsSmiled_Fx* fxModelParams = dynamic_cast<ARM_ModelParamsSmiled_Fx*>(modelParams);

	/*string domCcy((*modelNameMap)[DomModel]->Model()->GetZeroCurve()->GetCurrencyUnit()->GetCcyName());
	string fxDomCcy(fxModelParams->GetDomCurve()->GetCurrencyUnit()->GetCcyName());
	if( fxDomCcy != domCcy)
	{
		CC_Ostringstream os;
		os << ARM_USERNAME << " : The cross FXmodel has domestic curve == " << fxDomCcy << " while the IR domestic model ccy is == " << domCcy;
		ARM_THROW( ERR_INVALID_ARGUMENT, os.str() );
	}*/

	///Initialize the FXModel Correlation.
	ARM_EqFxBase* fxModel = dynamic_cast<ARM_EqFxBase*>(	&*(*modelNameMap)[FxModel]->Model() );
	fxModel->SetCorrelMatrix( *GetCorrelMatrix() );
}


////////////////////////////////////////////////////
///	Class   : ARM_1IRFXModel
///	Routines: Copy constructor
///	Returns :
///	Action  : constructor
////////////////////////////////////////////////////

ARM_1IRFXModel::ARM_1IRFXModel(const ARM_1IRFXModel& rhs)
:	ARM_MultiAssetsMeanReverting(rhs),
itsTimeSteps(rhs.itsTimeSteps),
itsFwdFxDrifts(rhs.itsFwdFxDrifts),
itsHWRelDrifts(rhs.itsHWRelDrifts),
itsIntAbsDrifts(rhs.itsIntAbsDrifts),
its2IRFXModel(CreateClonedPtr(&*rhs.its2IRFXModel))
{
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>( rhs.itsEigenValues, itsEigenValues );
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( rhs.itsEigenVectors, itsEigenVectors );
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( rhs.its2IRFXCorrelMatrix, its2IRFXCorrelMatrix );
}


////////////////////////////////////////////////////
///	Class   : ARM_1IRFXModel
///	Routines: destructor
///	Returns :
///	Action  : destructor
////////////////////////////////////////////////////

ARM_1IRFXModel::~ARM_1IRFXModel()
{
	DeletePointorVector<ARM_GP_Vector>( itsEigenValues );
	DeletePointorVector<ARM_GP_Matrix>( itsEigenVectors );
	DeletePointorVector<ARM_GP_Matrix>( its2IRFXCorrelMatrix );
}


size_t ARM_1IRFXModel::FactorCount() const
{
	ARM_SmiledModel_Fx* fxModel = dynamic_cast<ARM_SmiledModel_Fx*>(&*GetModel(FxModel));
	ARM_ModelParams* modelParams= fxModel->GetModelParams();
	ARM_ModelParamsSmiled_Fx* fxModelParams = dynamic_cast<ARM_ModelParamsSmiled_Fx*>(modelParams);

	if (fxModelParams->GetCorrelType() == ARM_ModelParamsSmiled::Fwd)
		return its2IRFXModel->FactorCount();
	else
		return ARM_MultiAssetsMeanReverting::FactorCount();
}

ARM_BoolVector ARM_1IRFXModel::NeedMCIntegProcess() const 
{ 
	ARM_SmiledModel_Fx* fxModel = dynamic_cast<ARM_SmiledModel_Fx*>(&*GetModel(FxModel));
	ARM_ModelParams* modelParams= fxModel->GetModelParams();
	ARM_ModelParamsSmiled_Fx* fxModelParams = dynamic_cast<ARM_ModelParamsSmiled_Fx*>(modelParams);

	if (fxModelParams->GetCorrelType() == ARM_ModelParamsSmiled::Fwd)
		return its2IRFXModel->NeedMCIntegProcess();
	else
		return ARM_BoolVector(FactorCount(), false);
};


///////////////////////////////////////////////////
///	Class   : ARM_1IRFXModel
///	Routine : CalibrateCorrelation
///	Returns : void
///	Action  : calibrate the correlation 
////////////////////////////////////////////////////

void ARM_1IRFXModel::CalibrateCorrelation( ARM_2IRFXModelPtr Model2IRFX )
{
	// Compute the correlMatrix based on the 2IRFX model
	ARM_SmiledModel_Fx* fxModel = dynamic_cast<ARM_SmiledModel_Fx*>(&*GetModel(FxModel));
	ARM_ModelParams* modelParams= fxModel->GetModelParams();
	ARM_ModelParamsSmiled_Fx* fxModelParams = dynamic_cast<ARM_ModelParamsSmiled_Fx*>(modelParams);
	if ((fxModelParams->GetCorrelType() == ARM_ModelParamsSmiled::CorrelMatrix) && !its2IRFXModel.IsNull())
	{
		ARM_SmiledModel_Fx* fxModel = dynamic_cast<ARM_SmiledModel_Fx*>(&*GetModel(FxModel));
		ARM_GP_Vector resetTimes = fxModel->getResetTimes();
		ARM_GP_Vector settlementTimes = fxModel->getSettlementTimes();
		ARM_GP_Matrix hump;
		ARM_GP_Vector totalVol;
		Model2IRFX->ComputeVolatilitiesAndCorrelMatrix(resetTimes,settlementTimes,hump,its2IRFXCorrelMatrix,totalVol, true);
	}
}

///////////////////////////////////////////////////
///	Class   : ARM_1IRFXModel
///	Routine : Init
///	Returns : ARM_PricingStatesPtr
///	Action  : intialise the model
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_1IRFXModel::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
	/// Validate that the payment currency is the domestic model!
    const ARM_ModelNameMap* const modelMap = GetModelMap();
    const string& domBasisModelName = (*modelMap)[DomBasisModel]->Model()->GetModelName();
    const string& domModelName      = (*modelMap)[DomModel]->Model()->GetModelName();
	if( payModelName != domBasisModelName && payModelName != domModelName)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": 1IFR + FX model allow only paimenent in " + domBasisModelName + " or " + domModelName);

	if (!its2IRFXModel.IsNull())
	{
		const ARM_ModelNameMap* const modelMap2IRFX = its2IRFXModel->GetModelMap();
		ARM_PricingModelPtr domModel    = (*modelMap2IRFX)[ARM_2IRFXModel::DomModel]->Model();
		ARM_PricingModelPtr forModel    = (*modelMap2IRFX)[ARM_2IRFXModel::ForModel]->Model();
		ARM_PricingModelPtr fxModel    = (*modelMap2IRFX)[ARM_2IRFXModel::FxModel]->Model();

		domModel->SetModelNb(0);
		forModel->SetModelNb(1);
		fxModel->SetModelNb(2);

		its2IRFXModel->SetRefModel(  &*( *GetModelMap() )[payModelName]->Model() );
		its2IRFXModel->SetZeroCurve( its2IRFXModel->GetRefModel()->GetZeroCurve() );
	}

	ARM_NumerairePtr numeraire = GetNumeraire();
    ARM_NumMethodPtr numMethod = GetNumMethod();

    if( numMethod.IsNull() || 
        numeraire.IsNull() ||
		((numMethod->GetPricingDirection() != ARM_NumMethod::GP_FWDLOOKING) &&
		(numMethod->GetPricingDirection() != ARM_NumMethod::GP_BCKWDLOOKING) &&
		(numMethod->GetPricingDirection() != ARM_NumMethod::GP_FWDBCKWDLOOKING) ) ||
        (numeraire->GetType() != ARM_Numeraire::TerminalZc &&
		numeraire->GetType() != ARM_Numeraire::RollingCash &&
		numeraire->GetType() != ARM_Numeraire::Cash) )
        //ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": only Monte Carlo & Terminal Zc numeraire supported by 1IR + FX model");
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": only Monte Carlo, Tree Terminal & Zc numeraire supported by 1IR + FX model");

	/// Delegate to common code
	//ARM_PricingStatesPtr facticeStates = its2IRFXModel->Init( payModelName, timeInfos );

	ARM_PricingStatesPtr states = ARM_MultiAssetsModel::Init( payModelName, timeInfos );
	//ARM_PricingStatesPtr states = its2IRFXModel->Init( payModelName, timeInfos );
	

	ARM_PricingModelPtr fxModel = (*modelMap)[FxModel]->Model();
	ARM_ModelParams* modelParams= fxModel->GetModelParams();
	ARM_ModelParamsSmiled_Fx* fxModelParams = dynamic_cast<ARM_ModelParamsSmiled_Fx*>(modelParams);

	return states;
}

////////////////////////////////////////////////////
///	Class   : ARM_1IRFXModel
///	Routines: SetNumMethod
///	Returns :
///	Action  : Set numerical method of the model
////////////////////////////////////////////////////
void ARM_1IRFXModel::SetNumMethod(const ARM_NumMethodPtr& numMethodPtr)
{
    ARM_MultiAssetsModel::SetNumMethod(numMethodPtr);
	its2IRFXModel->SetNumMethod(numMethodPtr);
}


////////////////////////////////////////////////////
///	Class   : ARM_1IRFXModel
///	Routines: SetNumeraire
///	Returns :
///	Action  : Set numeraire of the model
////////////////////////////////////////////////////
void ARM_1IRFXModel::SetNumeraire(const ARM_NumerairePtr& numMethodPtr)
{
    ARM_MultiAssetsModel::SetNumeraire(numMethodPtr);
	its2IRFXModel->SetNumeraire(numMethodPtr);
}

////////////////////////////////////////////////////
///	Class   : ARM_1IRFXModel
///	Routine : ComputeIntegratedVCV
///	Returns : void
///	Action  : computes the integrated VCV matrixes for
///           domestic state variables &
///           forward Fx process
////////////////////////////////////////////////////
void ARM_1IRFXModel::ComputeIntegratedVCV(
		double step, double nextStep,
		const ARM_GP_Vector& resetTimes,
		size_t fxModelNb,
        const ARM_ModelParamsQ1F* const domModelParams,
        const ARM_ModelParamsSmiled* const fxModelParams,
        const ARM_GP_Matrix& correlMatrix,
        ARM_GP_Matrix& variances) const
{
	double stdDevDomModel, stdDevFXModel;
	int i, j;

    /// Domestic variances
	variances(DomModel,DomModel)   = domModelParams->StateLocalVariance(step,nextStep,nextStep);
	stdDevDomModel = sqrt(variances(DomModel,DomModel));

	int nbFwds = resetTimes.size();
	int startTimePos = 0;

	size_t correlIndex = -1;
	if (fxModelParams->GetCorrelType() == ARM_ModelParamsSmiled::CorrelMatrix)
		correlIndex = IndexOfFirstHigherEqInVector_DefaultLast(nextStep,resetTimes);

	/// get the first bigger time
	while(startTimePos< resetTimes.size()
		&& resetTimes[startTimePos] < nextStep )
		++startTimePos;

	for(i=0; i<startTimePos; ++i)
	{
		variances(DomModel,fxModelNb+i)=variances(fxModelNb+i,DomModel)=0.0;
		for(j=i;j<nbFwds;++j)
			variances(fxModelNb+i,fxModelNb+j)=variances(fxModelNb+j,fxModelNb+i)=0.;
	}

	for(i=startTimePos; i<nbFwds; ++i)
	{
		for(j=i;j<nbFwds;++j)
			variances(fxModelNb+i,fxModelNb+j)=variances(fxModelNb+j,fxModelNb+i)=fxModelParams->IntegratedCovariance(step,nextStep,i,j);
	}

	for (i = startTimePos; i < nbFwds; ++i)
	{
		if (fxModelParams->GetCorrelType() == ARM_ModelParamsSmiled::CorrelMatrix)
		{
			stdDevFXModel = sqrt(variances(fxModelNb+i,fxModelNb+i));
			variances(DomModel,FxModel+i) = (*its2IRFXCorrelMatrix[correlIndex])(DomModel,FxModel+i)*stdDevDomModel*stdDevFXModel;
		}
		else
		{
			ARM_GP_Vector fxVolTimes = fxModelParams->GetVolCurve(i)->GetAbscisses();
			ARM_GP_Vector fxVolValues(fxVolTimes.size());
			ARM_CurveModelParam fxVolParam(ARM_ModelParamType::Volatility, &fxVolValues, &fxVolTimes);
			for (j = 0; j < fxVolTimes.size(); ++j)
			fxVolValues[j] = fxModelParams->GetVolCurve(i)->Interpolate(fxVolTimes[j]);
			variances(DomModel,FxModel+i) = correlMatrix(DomModel,FxModel)*ARM_ModelParamsHW1F::HW1FEqFxStateCovariance( fxVolParam, domModelParams, step,nextStep,nextStep );
		}
		variances(FxModel+i,DomModel) = variances(DomModel,FxModel+i);
	}

	//CC_Ostringstream os;
	//os << "Reset Time = " << step << "\n";
	//os << "Var Dom = "<< variances(DomModel,DomModel) << "\n";
	//os << "Var Fwd FX = "<< variances(FxModel+startTimePos,FxModel+startTimePos)<< "\n";
	//os << "Cov Dom/Fwd FX = "<< variances(DomModel,FxModel+startTimePos)<< "\n";
	//ARM_TheEventViewer.Instance()->AddToMessage(os.str());
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : IntegratedLocalDrifts
///	Returns : void
///	Action  : compute integrated relative and absolute drifts
////////////////////////////////////////////////////
void ARM_1IRFXModel::IntegratedLocalDrifts(
	const ARM_GP_Vector& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts) const
{
	ARM_SmiledModel_Fx* fxModel = dynamic_cast<ARM_SmiledModel_Fx*>(&*GetModel(FxModel));
	ARM_ModelParams* modelParams= fxModel->GetModelParams();
	ARM_ModelParamsSmiled_Fx* fxModelParams = dynamic_cast<ARM_ModelParamsSmiled_Fx*>(modelParams);

	if (fxModelParams->GetCorrelType() == ARM_ModelParamsSmiled::Fwd)
	{
		its2IRFXModel->IntegratedLocalDrifts(
			timeSteps,
			relativeDrifts,
			absoluteDrifts);

		int nbRows = absoluteDrifts->rows(), nbCols = absoluteDrifts->cols();

		itsIntAbsDrifts.resize(nbRows,nbCols);

		int i, j;

		for (i = 0; i < nbRows; ++i)
			for (j = 0; j < nbCols; ++j)
			{			
				if (i > 0)
					itsIntAbsDrifts(i,j) += itsIntAbsDrifts(i-1,j) + (*absoluteDrifts)(i,j);
				else
					itsIntAbsDrifts(i,j) = (*absoluteDrifts)(i,j);
			}
	}
	else
	{
		/// nothing as we need the full path to compute localDrift
		/// and compute drift on the fly
		size_t factorsNb = FactorCount();
		size_t nbSteps	 = timeSteps.size();
		size_t modelNb	 = GetModelNb();

		/// Initialise a Relative Drift matrix to 1 (no drift)
		relativeDrifts	= ARM_GP_MatrixPtr( new ARM_GP_Matrix( (nbSteps-1)*(modelNb+1), factorsNb, 1.0 ) );
		absoluteDrifts	= ARM_GP_MatrixPtr( NULL);
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_1IRFXModel
///	Routine : ComputeDrifts
///	Returns : void
///	Action  : computes the integrated VCV matrixes for
///           domestic state variables &
///           forward Fx process
////////////////////////////////////////////////////
void ARM_1IRFXModel::ComputeDrifts(
const ARM_GP_Vector& timeSteps,
const ARM_GP_Vector& resetTimes,
const ARM_GP_Vector& settlementTimes,
const ARM_ModelParamsQ1F* const domModelParams,
const ARM_ModelParamsSmiled* const fxModelParams,
const ARM_MatrixVector& localVariances )
{
	/// nothing as we need the full path to compute localDrift
	/// and compute drift on the fly
	size_t factorsNb = FactorCount();
	size_t nbSteps	 = timeSteps.size();
	size_t nbFwds	 = resetTimes.size();
	double numTime	 = settlementTimes[nbFwds-1];

	if (fxModelParams->GetCorrelType() != ARM_ModelParamsSmiled::Fwd)
	{
		itsFwdFxDrifts.resize(nbSteps-1,nbFwds+1);
		itsHWRelDrifts.resize(nbSteps-1);
		const ARM_CurveMatrix* const correlCurveMatrix  = GetCorrelMatrix();
		ARM_GP_Matrix correlMatrix;

		size_t i, j, k, firstFwd;

		double step		= timeSteps[0],nextStep;
		double betaTjTNum;

		ARM_GP_Vector* eigen;

		itsHWRelDrifts[0] = 1.0;

		for (i = 0; i < nbSteps-1; ++i)
		{
			correlMatrix = correlCurveMatrix->Interpolate(step);

			eigen = itsEigenValues[i];
			nextStep = timeSteps[i+1];

			firstFwd = IndexOfFirstHigherEqInVector_DefaultLast(nextStep,resetTimes);

			for (j = firstFwd+1; j < nbFwds+1; ++j)
			{
				itsFwdFxDrifts(i,j) = 0;
				for (k = 0; k < factorsNb; ++k)
					itsFwdFxDrifts(i,j) += (*localVariances[i])(j,k)*(*localVariances[i])(0,k)*(*eigen)[k];

				betaTjTNum = domModelParams->BetatT(resetTimes[j-1],numTime);
				itsFwdFxDrifts(i,j) *= domModelParams->StateLocalDrift(nextStep,resetTimes[j-1]);
				itsFwdFxDrifts(i,j) = -itsFwdFxDrifts(i,j)*betaTjTNum;
			}
			if (i < nbSteps-2)
				itsHWRelDrifts[i+1] = domModelParams->StateLocalDrift(step,nextStep);
			
			step=nextStep;
		}
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_1IRFXModel
///	Routine : ModelStateLocalVariancesAndStdDev
///	Returns : void
///	Action  : computes model covariances
////////////////////////////////////////////////////

void ARM_1IRFXModel::ModelStateLocalVariances( 
const ARM_GP_Vector& timeSteps,
ARM_MatrixVector& localVariances) const
{ 
	size_t factorNb	= FactorCount();
	const ARM_CurveMatrix* const correlCurveMatrix  = GetCorrelMatrix();
	ARM_GP_Matrix correlMatrix;
	const ARM_ModelNameMap* const modelMap = GetModelMap();
	size_t nbSteps	= timeSteps.size();
    double step		= timeSteps[0],nextStep;

	localVariances.resize(nbSteps-1);
	itsEigenValues.resize(nbSteps-1);

#if defined(__GP_STRICT_VALIDATION)
	if( factorNb < 2 ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "1IR+FX factor number should be greater or equal than 2 !" );
#endif

    ARM_PricingModelPtr domModel    = (*modelMap)[DomModel]->Model();
    ARM_SmiledModel_Fx* fxModel = static_cast<ARM_SmiledModel_Fx*>(&*((*modelMap)[FxModel]->Model()));
	size_t fxModelNb = fxModel->GetModelNb();
	const ARM_ModelParamsQ1F* const domModelParams  = static_cast<const ARM_ModelParamsQ1F* const>( domModel->GetModelParams() );
	const ARM_ModelParamsSmiled* const fxModelParams    = static_cast<const ARM_ModelParamsSmiled* const>( fxModel->GetModelParams() );

	const ARM_GP_Vector& resetTimes = fxModel->getResetTimes();
	const ARM_GP_Vector& settlementTimes = fxModel->getResetTimes();
	size_t nbFwds = resetTimes.size();

	ARM_GP_TriangularMatrix auxLocalVariances(nbFwds+1,0.0);

	if (fxModelParams->GetCorrelType() != ARM_ModelParamsSmiled::Fwd)
	{
		itsACPErrors.resize(nbSteps-1);
		itsEigenVectors.resize(nbSteps-1);
	}
	itsTimeSteps = timeSteps;

	for(size_t i=0;i<nbSteps-1;++i)
	{
		if (fxModelParams->GetCorrelType() == ARM_ModelParamsSmiled::Fwd)
		{
			localVariances[i] = new ARM_GP_Matrix(factorNb,factorNb);
		}
		else
		{
			correlMatrix = correlCurveMatrix->Interpolate(step);

			nextStep			= timeSteps[i+1];

			localVariances[i]	= new ARM_GP_Matrix(nbFwds+1,factorNb);
			itsEigenValues[i]	= new ARM_GP_Vector(factorNb);

			ComputeIntegratedVCV(
				step,
				nextStep,
				resetTimes,
				fxModelNb,
				domModelParams,
				fxModelParams,
				correlMatrix,
				auxLocalVariances);

			itsACPErrors[i] = ACPTransformationWithRescalling(&auxLocalVariances,(*itsEigenValues[i]),(*localVariances[i]));
			itsEigenVectors[i] = CreateClone(localVariances[i]);
		}

		step=nextStep;
	}

    const_cast<ARM_1IRFXModel*>(this)->ComputeDrifts(timeSteps,resetTimes,settlementTimes,domModelParams,fxModelParams,localVariances);
}


////////////////////////////////////////////////////
///	Class   : ARM_1IRFXModel
///	Routines: ModelStateLocalVariancesAndStdDev
///	Returns : void
///	Action  : computes the local variance and std dev 
////////////////////////////////////////////////////
void ARM_1IRFXModel::ModelStateLocalVariancesAndStdDev( const ARM_GP_Vector& timeSteps)
{
	ARM_PricingModel::ModelStateLocalVariancesAndStdDev(timeSteps);
}

////////////////////////////////////////////////////
///	Class   : ARM_1IRFXModel
///	Routines: NumMethodStateLocalVariances
///	Returns : void
///	Action  : NumMethodStateLocalVariances
////////////////////////////////////////////////////
void ARM_1IRFXModel::NumMethodStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const
{
	ARM_SmiledModel_Fx* fxModel = dynamic_cast<ARM_SmiledModel_Fx*>(&*GetModel(FxModel));
	ARM_ModelParams* modelParams= fxModel->GetModelParams();
	ARM_ModelParamsSmiled_Fx* fxModelParams = dynamic_cast<ARM_ModelParamsSmiled_Fx*>(modelParams);

	ARM_MatrixVector fxLocalVariances;

	if (fxModelParams->GetCorrelType() == ARM_ModelParamsSmiled::Fwd)
	{
		ARM_MatrixVector variances;
		its2IRFXModel->NumMethodStateLocalGlobalVariances(
			timeSteps,
			localVariances,
			variances);
	}
	else
	{
		size_t factorsNb     = FactorCount(),
		   timeStepsSize = timeSteps.size(),
		   modelNb		 = GetModelNb(),
		   startTimePos  = 0;

		localVariances.resize((timeStepsSize-1)*(modelNb+1));
		for(size_t i=0;i<timeStepsSize-1 ;++i)
		{
			/// initialize everything
			localVariances[i] = new ARM_GP_Matrix(factorsNb,factorsNb,0.);
			for (size_t k=0; k<factorsNb;k++)
				(*localVariances[i])(k,k) = (*itsEigenValues[i])[k];
		}
   }

	fxLocalVariances.resize(localVariances.size());

	for (int i = 0; i < localVariances.size(); ++i)
		fxLocalVariances[i] = CreateClone(localVariances[i]);

	fxModel->SetNumMethodStateLocalVars(fxLocalVariances);
}

////////////////////////////////////////////////////
///	Class   : ARM_1IRFXModel
///	Routines: NumMethodStateGlobalVariances
///	Returns : void
///	Action  : NumMethodStateGlobalVariances
////////////////////////////////////////////////////
void ARM_1IRFXModel::NumMethodStateGlobalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& globalVariances ) const
{
	ARM_SmiledModel_Fx* fxModel = dynamic_cast<ARM_SmiledModel_Fx*>(&*GetModel(FxModel));
	ARM_ModelParams* modelParams= fxModel->GetModelParams();
	ARM_ModelParamsSmiled_Fx* fxModelParams = dynamic_cast<ARM_ModelParamsSmiled_Fx*>(modelParams);

	if (fxModelParams->GetCorrelType() == ARM_ModelParamsSmiled::Fwd)
	{
		ARM_MatrixVector localVariances;
		its2IRFXModel->NumMethodStateLocalGlobalVariances(
			timeSteps,
			localVariances,
			globalVariances);
	}
	else
	{
		// USELESS !!!!! But we have to compute it for the sampler
		size_t factorsNb     = FactorCount(),
			   timeStepsSize = timeSteps.size(),
			   modelNb		 = GetModelNb(),
			   offsetIndex	 = (timeStepsSize-1)*modelNb;

#if defined(__GP_STRICT_VALIDATION)
		if( globalVariances.size()!= offsetIndex ) 
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,NumMethodStateGlobalVariances: globalVariances.size() != offsetIndex" );
#endif
		globalVariances.resize(timeStepsSize*(modelNb+1));
		for (size_t i=0;i<timeStepsSize;i++){
			globalVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(factorsNb,1.0);
		}
	}
	
}


////////////////////////////////////////////////////
///	Class  : ARM_1IRFXModel
///	Routine: TreeStatesToModelStates
///	Returns: void
/// Nabyl
///	Action : Convert the num method states to model states
////////////////////////////////////////////////////
void ARM_1IRFXModel::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{   
#if defined(__GP_STRICT_VALIDATION)
	if(timeIndex<0)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index is negative!");

	if( timeIndex > GetNumMethod()->GetTimeSteps()->size()-1 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index bigger than max size!");
#endif

	const ARM_NumMethodPtr& numMethod = GetNumMethod();
	ARM_TreeBase* treebase = dynamic_cast<ARM_TreeBase*>( &*numMethod );

	if (treebase==NULL)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "Only ARM_TreeBase method is available as backward method!");

	const ARM_GP_MatrixPtr& modelStates = states->GetNumMethodStates();
	states->SetModelStates(modelStates);

	//ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "ARM_1IRFXModel::TreeStatesToModelStates is not implemented!");
}


void ARM_1IRFXModel::From2to1IRFXModelStates(ARM_PricingStatesPtr& states, int timeIndex)
{}

////////////////////////////////////////////////////
///	Class  : ARM_1IRFXModel
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : Convert the num method states to model states
////////////////////////////////////////////////////
void ARM_1IRFXModel::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
    double currTime  = GetNumMethod()->GetTimeStep(timeIndex);
	double nextTime  = GetNumMethod()->GetTimeStep(timeIndex+1);
	size_t factorsNb = FactorCount();
	size_t statesNb  = states->size();
	size_t modelNb	 = GetModelNb();

	const ARM_ModelNameMap* modelMap	= GetModelMap();
	ARM_SmiledModel_Fx* fxModel = static_cast<ARM_SmiledModel_Fx*>(&*((*modelMap)[FxModel]->Model()));
	const ARM_ModelParamsSmiled* params	=	dynamic_cast< const ARM_ModelParamsSmiled*>(fxModel->GetModelParams());
	if( !params )
	   ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM::Induct : ARM_ModelParamsSmiled* needed" );

	if (params->GetCorrelType() == ARM_ModelParamsSmiled::Fwd)
	{
		size_t k;

		const ARM_ModelNameMap* modelMap2IRFX	= its2IRFXModel->GetModelMap();

		ARM_GP_MatrixPtr numStates = states->GetNumMethodStates();
		ARM_GP_MatrixPtr modelStates = states->GetModelStates();

		ARM_PricingModelPtr domModel	= (*modelMap2IRFX)[ARM_2IRFXModel::DomBasisModel]->Model();
		ARM_PricingModelPtr forModel	= (*modelMap2IRFX)[ARM_2IRFXModel::ForBasisModel]->Model();


		const ARM_GP_Matrix& localVar = *(fxModel->GetNumMethodStateLocalVars()[timeIndex]);
        double lnDrift = -0.5*localVar(ARM_2IRFXModel::FxModel,ARM_2IRFXModel::FxModel);

		ARM_GP_VectorPtr domDf(NULL),forDf(NULL);

		domDf = domModel->DiscountFactor( domModel->GetModelName(), currTime, nextTime, states );
		forDf = forModel->DiscountFactor( forModel->GetModelName(), currTime, nextTime, states );

		for(k = 0 ; k < statesNb ; ++k )
		{
			(*modelStates)(ARM_2IRFXModel::FxModel,k)  *= (*forDf)[k]/(*domDf)[k]*exp(lnDrift+(*numStates)(ARM_2IRFXModel::FxModel,k));

			(*modelStates)(ARM_2IRFXModel::DomModel,k)  =  (*numStates)(ARM_2IRFXModel::DomModel,k);
			(*modelStates)(ARM_2IRFXModel::ForModel,k)  =  (*numStates)(ARM_2IRFXModel::ForModel,k);
		}
	}
	else
	{
		const ARM_GP_Vector& resetTimes = fxModel->getResetTimes();

		size_t i,j,k;
		size_t fwdsNb  = resetTimes.size();
		size_t iLast   = fwdsNb-1;

		int iFirst  = 0;
		while( iFirst < fwdsNb 	&& resetTimes[iFirst]<nextTime)
		{
			++iFirst;
		}

		const ARM_MatrixVector& modelLocalVar	= GetModelStateLocalVars();
		double next,modStd;

		ARM_GP_Vector* eigen = itsEigenValues[timeIndex];
		
		for(k = 0 ; k < statesNb ; ++k )
		{
			next		= 0.0;
			for (j=0; j<factorsNb; ++j )
			{
				modStd	=	(*modelLocalVar[timeIndex])(0,j);
				next		+= modStd * (states->GetNumMethodState(k,j+modelNb));
			}
			states->SetModelState(k,modelNb,states->GetModelState(k,0)*itsHWRelDrifts[timeIndex] + next);

			for (i=iFirst+1;i<=iLast+1;++i)
			{
				next		= 0.0;
				for (j=0; j<factorsNb; ++j )
				{
					modStd	=	(*modelLocalVar[timeIndex])(i,j);
					next	+=	modStd * (states->GetNumMethodState(k,j+modelNb));
				}
				next  += states->GetModelState(k,i+modelNb)+itsFwdFxDrifts(timeIndex,i);
				states->SetModelState(k,i+modelNb,next);
			}
		}
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_1IRFXModel
///	Routines: BackwardLookingInit 
///	Returns :
///	Action  : 
////////////////////////////////////////////////////

ARM_PricingStatesPtr ARM_1IRFXModel::BackwardLookingInit(ARM_NumMethodPtr& numMethod, double firstInductTime)
{
	return its2IRFXModel->BackwardLookingInit(numMethod, firstInductTime);
}

////////////////////////////////////////////////////
///	Class   : ARM_1IRFXModel
///	Routines: toString
///	Returns :
///	Action  : stringify the object
////////////////////////////////////////////////////
string ARM_1IRFXModel::toString(const string& indent, const string& nextIndent) const
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

	ARM_SmiledModel_Fx* fxModel = dynamic_cast<ARM_SmiledModel_Fx*>(&*modelMap[FxModel]->Model());
	ARM_GP_Vector resetTimes = fxModel->getResetTimes();
	ARM_ModelParams* modelParams= fxModel->GetModelParams();
	ARM_ModelParamsSmiled_Fx* fxModelParams = dynamic_cast<ARM_ModelParamsSmiled_Fx*>(modelParams);

	size_t i, j, k;

	os << "Correlation matrix\n";

	if( GetCorrelMatrix() )
	{
		if (fxModelParams->GetCorrelType() == ARM_ModelParamsSmiled::CorrelMatrix)
		{
			for (i = 0; i < its2IRFXCorrelMatrix.size(); ++i)
			{
				os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(4)<< CC_NS(std,setprecision)(0) << "ResetTime " << i << " = " <<  CC_Round(resetTimes[i]) << endl;
				os << string("      ");

				for (k = 0; k < its2IRFXCorrelMatrix[i]->cols(); ++k)
				{
					if (k != 0)
						os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(6)<< CC_NS(std,setprecision)(0) << CC_Round(resetTimes[k-1]);
					else
						os << "dom";

					if (k < its2IRFXCorrelMatrix[i]->cols()-1)
						os << string("|");
				}

				os << endl;

				for (j = 0; j < its2IRFXCorrelMatrix[i]->rows(); ++j)
				{
					if (j != 0)
					{
						os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(9)<< CC_NS(std,setprecision)(0) << CC_Round(resetTimes[j-1]);
						os << string("|");
					}
					else
						os << string("      dom|");

					for (k = 0; k < its2IRFXCorrelMatrix[i]->cols(); ++k)
					{
						os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(6)<< CC_NS(std,setprecision)(2) << (*its2IRFXCorrelMatrix[i])(j,k)*100;

						if (k < its2IRFXCorrelMatrix[i]->cols()-1)
							os << string("|");
					}
					os << endl;
				}
				os << endl;
			}
		}
		else
		{
			os << GetCorrelMatrix()->toString(indent,nextIndent);
		}
	}

	if (itsEigenVectors.size())
	{
		os << "ACP Eigen Vectors \n" << endl;
		for (i = 0; i < itsEigenVectors.size(); ++i)
		{
			os << "Reset Date = " << itsTimeSteps[i+1] << endl;
			for (j = 0; j < itsEigenVectors[i]->cols(); ++j)
			{
				CC_Ostringstream osFactor;
				osFactor << "Factor" << j << " = ";
				os << CC_NS(std,setfill(' ')) << CC_NS(std,scientific) <<  CC_NS(std,setw)(9)<< CC_NS(std,setprecision)(0) << osFactor.str();
				os << CC_NS(std,setw)(7)<< CC_NS(std,setprecision)(2) << (*itsEigenValues[i])(j);
				os << string("|");

				for (k = 0; k < itsEigenVectors[i]->rows(); ++k)
				{
					os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(7)<< CC_NS(std,setprecision)(4) << (*itsEigenVectors[i])(k,j);

					if (k <  itsEigenVectors[i]->rows()-1)
						os << string("|");
				}
				os << endl;
			}
			os << endl;
		}
		os << endl;
	}

	if (itsACPErrors.size())
	{
		os << "ACP Errors \n" << endl;

		for (i = 0; i < itsACPErrors.size(); ++i)
			os << CC_NS(std,fixed) << CC_NS(std,setw)(12)<< CC_NS(std,setprecision)(10) << itsACPErrors[i] << endl;
	}

	os << endl;

    os << indent << "\n\n------> Domestic Stochastic IR Model <------\n";
    os << modelMap[DomModel]->Model()->toString(indent,nextIndent);

    os << indent << "\n\n------>      Stochastic FX Model     <------\n";
    os << modelMap[FxModel]->Model()->toString(indent,nextIndent);

    os << indent << "\n\n------>     Domestic Basis Model     <------\n";
    ARM_ForwardMargin& domBasisModel = static_cast< ARM_ForwardMargin& >(* modelMap[DomBasisModel]->Model());
    bool prevDump = domBasisModel.GetRefModelDump();
    domBasisModel.SetRefModelDump(false);
    os << domBasisModel.toString(indent,nextIndent);
    domBasisModel.SetRefModelDump(prevDump);

    return os.str();
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

