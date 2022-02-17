/*!
 *
 * Copyright (c) CM CIB January 2005 Paris
 *
 *	\file 2IRFXModel.cpp
 *
 *  \brief 2 Interest Rates model + FX hybrid model
 *
 *	\author  E. Benhamou, JM Prié
 *	\version 1.0
 *	\date January 2005
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/2IRFXModel.h"

/// gpmodels
#include "gpmodels/ModelParamsHW1F.h"
#include "gpmodels/Q1F_Fx.h"
#include "gpmodels/CEV_Fx.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/Q1F.h"
#include "gpmodels/ModelParamsQ1F.h"
#include "gpmodels/forwardmargin.h"
#include "gpmodels/typedef.h"
#include "gpmodels/local_sln_model.h"

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
#include "gpinfra/surfacelistmodelparam.h"
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

//gpclosedforms
#include "gpclosedforms/normal.h"

#include <algorithm>
CC_USING_NS (std,find)

CC_BEGIN_NAMESPACE( ARM )

const size_t ARM_2IRFXModel::NbLocalModels = 3;

const double MAX_LNVOL_FX = 10.0; // 1000%

////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routines: Constructor
///	Returns :
///	Action  : builds the object
///	model types contains the integer on the various model
///	using the correlation matrix
////////////////////////////////////////////////////
ARM_2IRFXModel::ARM_2IRFXModel(	const ARM_ModelNameMap&	modelNameMap, 
	const ARM_CurveMatrix& correlMatrix )
:
	ARM_MultiAssetsMeanReverting(&modelNameMap,&correlMatrix),
	itsIsFlooredFxLocalModel(false),
	itsIsCappedFxLocalModel(false),
	itsIsRedemptionFxLocalModel(false)
{
	/// validate and initialize
	Validate();
}

////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routine : InitSubModels
///	Returns : void
///	Action  : intialises the sub models that are linked to a model
////////////////////////////////////////////////////
void ARM_2IRFXModel::Validate( )
{
	ARM_ModelNameMap* modelNameMap	= GetModelMap();
	/// Check linked models consistency
    if(modelNameMap->size() + NbLocalModels < NbModels)
    {
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : not enough models to built an hybrid 2IR+FX with basis model");
    }
    else if(modelNameMap->size() > NbModels)
    {
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : too much models to built an hybrid 2IR+FX with basis model");
    }

	if( !GetCorrelMatrix() )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : please provide a correlation");
	if(	GetCorrelMatrix()->rows() != 3 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : correlation should be 3x3!");	

	/// check model types
	if( !dynamic_cast<ARM_QModel1F*>( &*(*modelNameMap)[DomModel]->Model() ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : dom model has to be Q1F IR Model");

	if( !dynamic_cast<ARM_QModel1F*>( &*(*modelNameMap)[ForModel]->Model() ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : for model has to be Q1F IR Model");

	if(!dynamic_cast<ARM_QModel1F_Fx*>( &*(*modelNameMap)[FxModel]->Model() ) && !dynamic_cast<ARM_CEVModel_Fx*>( &*(*modelNameMap)[FxModel]->Model() ))
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : fx model has to be Q1F FX Model or CEV FX Model");
	
	if( !dynamic_cast<ARM_ForwardMargin*>(&*(*modelNameMap)[DomBasisModel]->Model()) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : dom basis model has to be a Forward Margin Model");

	if( !dynamic_cast<ARM_ForwardMargin*>(&*(*modelNameMap)[ForBasisModel]->Model()) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : for basis model has to be a Forward Margin Model");

    
    if(modelNameMap->size() + NbLocalModels >= NbModels + 1)
    {
	    if( !dynamic_cast<ARM_Local_SLN_Model*>(&*(*modelNameMap)[FlooredFxLocalModel]->Model()) )
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : local model for floored Fx coupon has to be a Local SLN Model");
        itsIsFlooredFxLocalModel = true;
    }
    
    if(modelNameMap->size() + NbLocalModels >= NbModels + 2)
    {
	    if( !dynamic_cast<ARM_Local_SLN_Model*>(&*(*modelNameMap)[CappedFxLocalModel]->Model()) )
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : local model for capped Fx coupon has to be a Local SLN Model");
        itsIsCappedFxLocalModel = true;
    }

    if(modelNameMap->size() + NbLocalModels >= NbModels + 3)
    {
	    if( !dynamic_cast<ARM_Local_SLN_Model*>(&*(*modelNameMap)[RedemptionFxLocalModel]->Model()) )
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : local model for redemption Fx coupon has to be a Local SLN Model");
        itsIsRedemptionFxLocalModel = true;
    }

	/// check the cross
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
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routines: Copy constructor
///	Returns :
///	Action  : constructor
ARM_2IRFXModel::ARM_2IRFXModel(const ARM_2IRFXModel& rhs)
:	ARM_MultiAssetsMeanReverting(rhs) 
{
    itsIsFlooredFxLocalModel = rhs.itsIsFlooredFxLocalModel;
    itsIsCappedFxLocalModel=rhs.itsIsCappedFxLocalModel;
    itsIsRedemptionFxLocalModel=rhs.itsIsRedemptionFxLocalModel;

    itsDriftCorrections = rhs.itsDriftCorrections != ARM_GP_MatrixPtr(NULL)
                          ? ARM_GP_MatrixPtr( static_cast< ARM_GP_Matrix* >(rhs.itsDriftCorrections->Clone()) )
                          : ARM_GP_MatrixPtr(NULL);

	ARM_2IRFXModel&	rhsbis =const_cast<ARM_2IRFXModel&>(rhs); // 
	ARM_PricingModel* refModel = rhsbis.GetRefModel();
	SetRefModel(refModel);

}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routine : SupportAnalyticMarginal
///	Returns : void
///	Action  : decide if over sampling between model time
///           step is needed
////////////////////////////////////////////////////
bool ARM_2IRFXModel::SupportAnalyticMarginal() const
{
    const ARM_NumerairePtr numeraire = GetNumeraire();
    return  numeraire == ARM_NumerairePtr(NULL) ||
            numeraire->GetType() == ARM_Numeraire::RollingEvent ||
			numeraire->GetType() == ARM_Numeraire::RollingCash;
}

////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routine : NeedMCIntegProcess
///	Returns : ARM_BoolVector
///	Action  : to say how MC method give simulated
///           processes
////////////////////////////////////////////////////
ARM_BoolVector ARM_2IRFXModel::NeedMCIntegProcess() const
{
    ARM_BoolVector SimulProcess(FactorCount(),true);
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
///	Class   : ARM_2IRFXModel
///	Routine : Init
///	Returns : ARM_PricingStatesPtr
///	Action  : intialise the model
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_2IRFXModel::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
	/// Validate that the payment currency is the domestic model!
    const ARM_ModelNameMap* const modelMap = GetModelMap();
    const string& domBasisModelName = (*modelMap)[DomBasisModel]->Model()->GetModelName();
    const string& domModelName      = (*modelMap)[DomModel]->Model()->GetModelName();
	if( payModelName != domBasisModelName && payModelName != domModelName)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": 2IFR + FX model allow only paimenent in " + domBasisModelName + " or " + domModelName);

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
        ": only Cash or TerminalZc or RollingEvent numeraire supported by 2IR + FX model");

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
///	Class   : ARM_2IRFXModel
///	Routine : ComputeIntegratedFwdFxVCV
///	Returns : void
///	Action  : computes the integrated VCV matrixes for
///           domestic & foreign state variables &
///           forward Fx process
////////////////////////////////////////////////////
void ARM_2IRFXModel::ComputeIntegratedFwdFxVCV(double step, double nextStep,
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
	//CC_Ostringstream os;
	//os << "Reset Time = " << step << "\n";
	//os << "Var Dom = "<< variances(DomModel,DomModel) << "\n";
	//os << "Var Fwd FX = "<< variances(FxModel,FxModel) << "\n";
	//os << "Cov Dom/Fwd FX = "<< variances(DomModel,FxModel) << "\n";
	//ARM_TheEventViewer.Instance()->AddToMessage(os.str());
}

////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routine : ComputeIntegratedSpotFxVCVInQModel
///	Returns : void
///	Action  : computes the integrated VCV matrixes for
///           domestic & foreign state variables &
///           spot Fx process
////////////////////////////////////////////////////
void ARM_2IRFXModel::ComputeIntegratedSpotFxVCVInQModel(
		double step, double nextStep,
        const ARM_ModelParamsHW1FStd* const domModelParams,
        const ARM_ModelParamsHW1FStd* const forModelParams,
        const ARM_ModelParamsQ1F* const fxModelParams,
        const ARM_GP_Matrix& correlMatrix,
        ARM_GP_Matrix& variances) const
{
    double covarXDomFx,covarXForFx;

    /// Domestic variances
	variances(DomModel,DomModel)   = domModelParams->StateLocalVariance(step,nextStep,nextStep);

    /// Foreign variances
	variances(ForModel,ForModel)   = forModelParams->StateLocalVariance(step,nextStep,nextStep);

    /// Domestic & Foreign covariances
	variances(DomModel,ForModel) = correlMatrix(DomModel,ForModel) * ARM_ModelParamsHW1F::HW1FStateCovariance(domModelParams, forModelParams,step,nextStep,nextStep);
    variances(ForModel,DomModel) = variances(DomModel,ForModel);

    /// Forex variances
	variances(FxModel,FxModel) = ARM_ModelParamsQ1F::Q1FStateCovariance(fxModelParams,fxModelParams,step,nextStep,nextStep);


    /// Domestic & Forex covariances
    const ARM_ModelParamsQ1F* const domParams = static_cast< const ARM_ModelParamsQ1F* >(domModelParams);
    covarXDomFx    = ARM_ModelParamsQ1F::Q1FStateCovariance(fxModelParams,domParams,step,nextStep,nextStep,false );
	variances(DomModel,FxModel) = correlMatrix(DomModel,FxModel)*covarXDomFx;
    variances(FxModel,DomModel) = variances(DomModel,FxModel);

    /// Foreign & Forex covariances
    const ARM_ModelParamsQ1F* const forParams = static_cast< const ARM_ModelParamsQ1F* >(forModelParams);
    covarXForFx    = ARM_ModelParamsQ1F::Q1FStateCovariance(fxModelParams,forParams,step,nextStep,nextStep,false);
	variances(ForModel,FxModel) = correlMatrix(ForModel,FxModel)*covarXForFx;
    variances(FxModel,ForModel) = variances(ForModel,FxModel);
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routine : ComputeIntegratedSpotFxVCV
///	Returns : void
///	Action  : computes the integrated VCV matrixes for
///           domestic & foreign state variables &
///           spot Fx process
////////////////////////////////////////////////////
void ARM_2IRFXModel::ComputeIntegratedSpotFxVCV(
		double step, double nextStep,
        const ARM_ModelParamsHW1FStd* const domModelParams,
        const ARM_ModelParamsHW1FStd* const forModelParams,
        const ARM_ModelParamsHW1FStd* const fxModelParams,
        const ARM_GP_Matrix& correlMatrix,
        ARM_GP_Matrix& variances) const
{
    double covarXDomFx,covarXForFx;

    /// Domestic variances
	variances(DomModel,DomModel)   = domModelParams->StateLocalVariance(step,nextStep,nextStep);

    /// Foreign variances
	variances(ForModel,ForModel)   = forModelParams->StateLocalVariance(step,nextStep,nextStep);

    /// Domestic & Foreign covariances
	variances(DomModel,ForModel) = correlMatrix(DomModel,ForModel) * ARM_ModelParamsHW1F::HW1FStateCovariance(domModelParams, forModelParams,step,nextStep,nextStep);
    variances(ForModel,DomModel) = variances(DomModel,ForModel);

    /// Forex variances
	variances(FxModel,FxModel) = fxModelParams->StateLocalVariance(step,nextStep,nextStep);


    /// Domestic & Forex covariances
    const ARM_ModelParamsQ1F* const domParams = static_cast< const ARM_ModelParamsQ1F* >(domModelParams);
    covarXDomFx    = ARM_ModelParamsHW1FStd::HW1FStateCovariance(fxModelParams,domParams,step,nextStep,nextStep);
	variances(DomModel,FxModel) = correlMatrix(DomModel,FxModel)*covarXDomFx;
    variances(FxModel,DomModel) = variances(DomModel,FxModel);

    /// Foreign & Forex covariances
    const ARM_ModelParamsQ1F* const forParams = static_cast< const ARM_ModelParamsQ1F* >(forModelParams);
    covarXForFx    = ARM_ModelParamsHW1FStd::HW1FStateCovariance(fxModelParams,forParams,step,nextStep,nextStep);
	variances(ForModel,FxModel) = correlMatrix(ForModel,FxModel)*covarXForFx;
    variances(FxModel,ForModel) = variances(ForModel,FxModel);
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routine : ComputeCorrelMatrix
///	Returns : void
///	Action  : compute the correlation based 
/// on a 2IRFXModel
////////////////////////////////////////////////////

void ARM_2IRFXModel::ComputeVolatilitiesAndCorrelMatrix(
			const std::vector<double>& resetTimes,
			const std::vector<double>& settlementTimes,
			ARM_GP_Matrix& volatilities,
			ARM_MatrixVector& correlMatrixVector,
			std::vector<double>& totalVolatilities,
			bool withDomesticIR) const
{
	size_t FxStart;

	if (withDomesticIR)
		FxStart = 1;
	else
		FxStart = 0;

	const ARM_CurveMatrix* const correlCurveMatrix  = GetCorrelMatrix();
	ARM_GP_Matrix correlMatrix;
	const ARM_ModelNameMap* const modelMap = GetModelMap();
	ARM_PricingModelPtr domModel    = (*modelMap)[DomModel]->Model();
    ARM_PricingModelPtr forModel    = (*modelMap)[ForModel]->Model();
    ARM_PricingModelPtr q1fFxModel     = (*modelMap)[FxModel]->Model();
	const ARM_ModelParamsHW1FStd* const domModelParams  = static_cast<const ARM_ModelParamsHW1FStd* const>( domModel->GetModelParams() );
	const ARM_ModelParamsHW1FStd* const forModelParams  = static_cast<const ARM_ModelParamsHW1FStd* const>( forModel->GetModelParams() );
	const ARM_ModelParamsHW1FStd* const fxModelParams   = static_cast<const ARM_ModelParamsHW1FStd* const>( q1fFxModel->GetModelParams() );
	const ARM_ModelParamsQ1F* const fxParams            = dynamic_cast<const ARM_ModelParamsQ1F* const>( fxModelParams );
	const ARM_CurveModelParam& fxVol                    = dynamic_cast<const ARM_CurveModelParam&>( fxModelParams->GetModelParam(fxModelParams->GetVolatilityType()) );

	size_t nbFwds	 = resetTimes.size();
	correlMatrixVector.resize(nbFwds);

	size_t factorsNb = nbFwds + FxStart;

	size_t i, j, k;

	// 15 var/covar
	double varDom, varZcDomi, varZcFori, varFX;
	double covarZcDomiZcFori;
	double covarZcForiDom, covarZcDomiDom, covarFXDom;
	double covarZcDomiZcDomj, covarZcForiZcForj;
	double covarZcForiFX, covarZcForjFX;
	double covarZcDomiFX, covarZcDomjFX;
	double covarZcDomiZcForj, covarZcForiZcDomj;

	correlMatrixVector.resize(nbFwds);
	volatilities.resize(nbFwds,nbFwds);

	totalVolatilities.resize(nbFwds);

	double step, nextStep, sqrdt;

	step = 0.0;

	for (k = 0; k <nbFwds; ++k)
	{
		nextStep = resetTimes[k];

		correlMatrix = correlCurveMatrix->Interpolate(step);

		sqrdt = sqrt(((double)nextStep-step)/resetTimes[nbFwds-1]);
		correlMatrixVector[k] = new ARM_GP_Matrix(factorsNb,factorsNb);

		varDom = ARM_ModelParamsHW1FStd::HW1FStateCovariance(domModelParams,domModelParams,step,nextStep,nextStep);
		varFX = fxParams->StateLocalVariance(step,nextStep,nextStep);
		covarFXDom = ARM_ModelParamsHW1FStd::HW1FEqFxStateCovariance(fxVol,domModelParams,step,nextStep,nextStep)*correlMatrix(DomModel,FxModel);

		if (withDomesticIR)
			(*correlMatrixVector[k])(0,0) = varDom;
		
		for (i = 0; i < nbFwds; ++i)
		{
			varZcDomi = ARM_ModelParamsHW1FStd::HW1FZcCovariance(domModelParams,domModelParams,step,nextStep,settlementTimes[i],settlementTimes[i]);
			varZcFori = ARM_ModelParamsHW1FStd::HW1FZcCovariance(forModelParams,forModelParams,step,nextStep,settlementTimes[i],settlementTimes[i]);

			covarZcDomiFX = ARM_ModelParamsHW1FStd::HW1FEqFxZcCovariance(fxVol,domModelParams,step,nextStep,settlementTimes[i])*correlMatrix(DomModel,FxModel);
			covarZcForiFX = ARM_ModelParamsHW1FStd::HW1FEqFxZcCovariance(fxVol,forModelParams,step,nextStep,settlementTimes[i])*correlMatrix(ForModel,FxModel);
			covarZcDomiZcFori = ARM_ModelParamsHW1FStd::HW1FZcCovariance(domModelParams,forModelParams,step,nextStep,settlementTimes[i])*correlMatrix(DomModel,ForModel);

			covarZcForiDom = ARM_ModelParamsHW1FStd::HW1FStateZcCovariance(domModelParams,forModelParams,step,nextStep,nextStep,settlementTimes[i])*correlMatrix(DomModel,ForModel);
			covarZcDomiDom = ARM_ModelParamsHW1FStd::HW1FStateZcCovariance(domModelParams,domModelParams,step,nextStep,nextStep,settlementTimes[i]);
			
			(*correlMatrixVector[k])(i+FxStart,i+FxStart) = varZcDomi + varZcFori + varFX + 2*(covarZcForiFX - covarZcDomiFX - covarZcDomiZcFori);
			if ((*correlMatrixVector[k])(i+FxStart,i+FxStart) < -K_NEW_DOUBLE_TOL)
				(*correlMatrixVector[k])(i+FxStart,i+FxStart) = 0.0;
			if (withDomesticIR)
			{
				(*correlMatrixVector[k])(0,i+FxStart) = covarFXDom + covarZcForiDom - covarZcDomiDom;
				(*correlMatrixVector[k])(i+FxStart,0) = (*correlMatrixVector[k])(0,i+FxStart);
			}

			for (j = 0; j < i; ++j)
			{
				covarZcDomiZcDomj = ARM_ModelParamsHW1FStd::HW1FZcCovariance(domModelParams,domModelParams,step,nextStep,settlementTimes[i],settlementTimes[j]);
				covarZcForiZcForj = ARM_ModelParamsHW1FStd::HW1FZcCovariance(forModelParams,forModelParams,step,nextStep,settlementTimes[i],settlementTimes[j]);
				covarZcDomiZcForj = ARM_ModelParamsHW1FStd::HW1FZcCovariance(domModelParams,forModelParams,step,nextStep,settlementTimes[i],settlementTimes[j])*correlMatrix(DomModel,ForModel);
				covarZcForiZcDomj = ARM_ModelParamsHW1FStd::HW1FZcCovariance(forModelParams,domModelParams,step,nextStep,settlementTimes[i],settlementTimes[j])*correlMatrix(DomModel,ForModel);
				covarZcDomjFX = ARM_ModelParamsHW1FStd::HW1FEqFxZcCovariance(fxVol,domModelParams,step,nextStep,settlementTimes[j])*correlMatrix(DomModel,FxModel);
				covarZcForjFX = ARM_ModelParamsHW1FStd::HW1FEqFxZcCovariance(fxVol,forModelParams,step,nextStep,settlementTimes[j])*correlMatrix(ForModel,FxModel);

				(*correlMatrixVector[k])(i+FxStart,j+FxStart) = (*correlMatrixVector[k])(j+FxStart,i+FxStart) = covarZcDomiZcDomj + covarZcForiZcForj + varFX + covarZcForiFX + covarZcForjFX - covarZcDomiFX - covarZcDomjFX - covarZcDomiZcForj - covarZcForiZcDomj;
			}

		}


		for (i = k; i < nbFwds; ++i)
		{
			if (((*correlMatrixVector[k])(FxStart+i,FxStart+i) > -K_NEW_DOUBLE_TOL) ||((*correlMatrixVector[k])(FxStart+i,FxStart+i) < -K_NEW_DOUBLE_TOL))
				volatilities(k,i) = sqrt((*correlMatrixVector[k])(FxStart+i,FxStart+i)/((double)nextStep-step)*K_YEAR_LEN);
			else
				ARM_THROW( ERR_INVALID_ARGUMENT, "negative variance" );
			totalVolatilities[i] += (*correlMatrixVector[k])(FxStart+i,FxStart+i);
		}

		for (i = 0; i < factorsNb; ++i)
		{
			for (j = 0; j < i; ++j)
			{
				(*correlMatrixVector[k])(i,j) /= sqrt((*correlMatrixVector[k])(i,i));
				(*correlMatrixVector[k])(i,j) /= sqrt((*correlMatrixVector[k])(j,j));

				(*correlMatrixVector[k])(j,i) /= sqrt((*correlMatrixVector[k])(i,i));
				(*correlMatrixVector[k])(j,i) /= sqrt((*correlMatrixVector[k])(j,j));
			}
		}

		for (i = 0; i < factorsNb; ++i)
			(*correlMatrixVector[k])(i,i) = 1.0;
		
		step = nextStep;
	}

	for (i = 0; i < nbFwds; ++i)
		totalVolatilities[i] = sqrt(totalVolatilities[i]/resetTimes[i]*K_YEAR_LEN);


	// Transform the volatilities into hump
	/*for (k = 0; k <nbFwds; ++k)
	{
		for (i = 0; i < nbFwds; ++i)
			volatilities(k,i) /= totalVolatilities[i];
	}*/
}


// Compute the volatilities used for the 1IRFX
void ARM_2IRFXModel::ComputeVolatilities(
		const std::vector<double>& resetTimes,
		const std::vector<double>& settlementTimes,
		ARM_GP_Matrix& volatilities,
		std::vector<double>& totalVar) const
{

	const ARM_CurveMatrix* const correlCurveMatrix  = GetCorrelMatrix();
	ARM_GP_Matrix correlMatrix;
	const ARM_ModelNameMap* const modelMap = GetModelMap();
	ARM_PricingModelPtr domModel    = (*modelMap)[DomModel]->Model();
    ARM_PricingModelPtr forModel    = (*modelMap)[ForModel]->Model();
    ARM_PricingModelPtr q1fFxModel     = (*modelMap)[FxModel]->Model();
	const ARM_ModelParamsHW1FStd* const domModelParams  = static_cast<const ARM_ModelParamsHW1FStd* const>( domModel->GetModelParams() );
	const ARM_ModelParamsHW1FStd* const forModelParams  = static_cast<const ARM_ModelParamsHW1FStd* const>( forModel->GetModelParams() );
	const ARM_ModelParamsHW1FStd* const fxModelParams   = static_cast<const ARM_ModelParamsHW1FStd* const>( q1fFxModel->GetModelParams() );
	const ARM_ModelParamsQ1F* const fxParams            = dynamic_cast<const ARM_ModelParamsQ1F* const>( fxModelParams );
	const ARM_CurveModelParam& fxVol                    = dynamic_cast<const ARM_CurveModelParam&>( fxModelParams->GetModelParam(fxModelParams->GetVolatilityType()) );

	size_t nbFwds	 = resetTimes.size();
	size_t factorsNb = nbFwds;

	size_t i, k;

	double varZcDomi, varZcFori, varFX;
	double covarZcDomiZcFori, covarZcDomiFX, covarZcForiFX;

	volatilities.resize(nbFwds,nbFwds);
	totalVar.resize(nbFwds);
	double step, nextStep;

	step = 0.0;
	// Just for test
	std::vector<double> totalVol(nbFwds, 0.);
	double sum(0.);

	for (k = 0; k <nbFwds; ++k)
	{
		nextStep = resetTimes[k];

		correlMatrix = correlCurveMatrix->Interpolate(step);

		varFX = fxParams->StateLocalVariance(step,nextStep,nextStep);
		
		for (i = k; i < nbFwds; ++i)
		{
			varZcDomi = ARM_ModelParamsHW1FStd::HW1FZcCovariance(domModelParams,domModelParams,step,nextStep,settlementTimes[i],settlementTimes[i]);
			varZcFori = ARM_ModelParamsHW1FStd::HW1FZcCovariance(forModelParams,forModelParams,step,nextStep,settlementTimes[i],settlementTimes[i]);

			covarZcDomiFX = ARM_ModelParamsHW1FStd::HW1FEqFxZcCovariance(fxVol,domModelParams,step,nextStep,settlementTimes[i])*correlMatrix(DomModel,FxModel);
			covarZcForiFX = ARM_ModelParamsHW1FStd::HW1FEqFxZcCovariance(fxVol,forModelParams,step,nextStep,settlementTimes[i])*correlMatrix(ForModel,FxModel);
			covarZcDomiZcFori = ARM_ModelParamsHW1FStd::HW1FZcCovariance(domModelParams,forModelParams,step,nextStep,settlementTimes[i])*correlMatrix(DomModel,ForModel);
			
			volatilities(k,i) = (varZcDomi + varZcFori + varFX + 2*(covarZcForiFX - covarZcDomiFX - covarZcDomiZcFori));
			
			volatilities(k,i) /= ((nextStep-step)/K_YEAR_LEN);

			volatilities(k,i) = sqrt(volatilities(k,i));

			totalVar(i) += varZcDomi + varZcFori + varFX + 2*(covarZcForiFX - covarZcDomiFX - covarZcDomiZcFori);

			
			//totalVol(i) = volatilities(k,i);
			totalVol(i) = sqrt( totalVar(i)/(nextStep/K_YEAR_LEN) );
			
		}	
		//totalVol(k) /= sqrt( nextStep/K_YEAR_LEN );
		step = nextStep;
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routine : ComputeMeanRevertingVariances
///	Returns : void
///	Action  : Compute mean reverting variances
/// Copy Pasted of the multi asset mean reverting !!!
////////////////////////////////////////////////////

void ARM_2IRFXModel::ComputeMeanRevertingVariances( 
	const std::vector<double>& timeSteps,
	ARM_MatrixVector& localVariances,
	ARM_MatrixVector& variances ) const
{
	size_t factorNb	= FactorCount();
	const ARM_CurveMatrix* const correlCurveMatrix  = GetCorrelMatrix();
	ARM_GP_Matrix correlMatrix;
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

	variances[0]=new ARM_GP_TriangularMatrix(factorNb,0.0);

	size_t i,j,k;
	ARM_ModelNameMap::const_iterator iter,iter2;
	double elem,elem2;
	
	for(i=0;i<nbSteps-1;++i)
	{

		nextStep			= timeSteps[i+1];

		correlMatrix = correlCurveMatrix->Interpolate(nextStep);

		localVariances[i]	= new ARM_GP_TriangularMatrix(factorNb,0.0);
		variances[i+1]		= new ARM_GP_TriangularMatrix(factorNb,0.0);

		j = 0;
		for( iter=modelMap->begin(); iter!=modelMap->end(); modelMap->getNextUsedIter(iter), ++j )
		{
			const ARM_ModelParams* const modelParams = (*iter).Model()->GetModelParams();
			const ARM_ModelParamsHW1FStd* const modelParamsHW1FStd = dynamic_cast<const ARM_ModelParamsHW1FStd* const>( modelParams );

			(*(localVariances[i]))(j,j) = modelParamsHW1FStd->StateLocalVariance(step,nextStep,nextStep);
			(*(variances[i+1]))(j,j)    = modelParamsHW1FStd->StateLocalVariance(0,nextStep,nextStep);

			k=0;
			for( iter2=modelMap->begin(); k<j; modelMap->getNextUsedIter(iter2), ++k )
			{
				const ARM_ModelParams* const modelParams_2 = (*iter2).Model()->GetModelParams();
				const ARM_ModelParamsHW1FStd* const modelParamsHW1FStd_2 = dynamic_cast<const ARM_ModelParamsHW1FStd* const>( modelParams_2 );

				elem = 	correlMatrix(j,k) * ARM_ModelParamsHW1F::HW1FStateCovariance(modelParamsHW1FStd, modelParamsHW1FStd_2,step,nextStep,nextStep);
				elem2= 	correlMatrix(j,k) * ARM_ModelParamsHW1F::HW1FStateCovariance(modelParamsHW1FStd, modelParamsHW1FStd_2,0,nextStep,nextStep);
				(*(localVariances[i]))(j,k)	= elem;
				(*(localVariances[i]))(k,j)	= elem;
				(*(variances[i+1]))(j,k)	= elem2;
				(*(variances[i+1]))(k,j)	= elem2;
			}
		}
		step=nextStep;
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routine : NumMethodStateLocalGlobalVariances
///	Returns : void
///	Action  : computes the integrated covariances
////////////////////////////////////////////////////

void ARM_2IRFXModel::NumMethodStateLocalGlobalVariances( const std::vector<double>& timeSteps,
	ARM_MatrixVector& localVariances,
	ARM_MatrixVector& variances ) const
{
	//CC_Ostringstream os1;
	//os1 << "2IRFX : Start Compute VCV\n";
	//ARM_TheEventViewer.Instance()->AddToMessage(os1.str());

    if( GetNumMethod()->GetSampler() && !GetNumMethod()->GetSampler()->IsIntegratedSampling() )
    {
        /// Compute VCV using only OU covariances
        ComputeMeanRevertingVariances(timeSteps,localVariances,variances);
        return;
    }

    /// In case of integrated sampling, FX variances & covariances must use
    /// foward fx ones then domestic and foreign zero-coupon bond volatilities

	size_t factorNb	= FactorCount();
	const ARM_CurveMatrix* const correlCurveMatrix  = GetCorrelMatrix();
	ARM_GP_Matrix correlMatrix ;
	const ARM_ModelNameMap* const modelMap = GetModelMap();
	size_t nbSteps	= timeSteps.size();
    double step		= timeSteps[0],nextStep;

#if defined(__GP_STRICT_VALIDATION)
	if( factorNb != 3 ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "2IR+FX factor number differs from 3 !" );
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
	const ARM_ModelParamsQ1F* const fxParams            = dynamic_cast<const ARM_ModelParamsQ1F* const>( fxModelParams );
    const ARM_CurveModelParam& fxVol                    = dynamic_cast<const ARM_CurveModelParam&>( fxModelParams->GetModelParam(fxModelParams->GetVolatilityType()) );
    bool isLnFwdFx = fxModelParams->IsLn();

    /// In backward looking fwd Fx is diffused, but in forward looking spot Fx is locally diffused
    bool isFwdFxVCV = GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_BCKWDLOOKING || isLnFwdFx;

	ARM_QModel1F_Fx* qFxModel = dynamic_cast<ARM_QModel1F_Fx*>(&*fxModel);
	bool isQIntegratedVersion = !isFwdFxVCV && qFxModel && qFxModel->GetIntegratedVersion();

	double yf,dyf;
	ARM_GP_Matrix* vcv;
	for(size_t i=0;i<nbSteps-1;++i)
	{
		correlMatrix = correlCurveMatrix->Interpolate(step);

		nextStep			= timeSteps[i+1];

		localVariances[i]	= new ARM_GP_TriangularMatrix(factorNb,0.0);

        if(isFwdFxVCV)
        {
            ComputeIntegratedFwdFxVCV(step,nextStep,
                domModelParams,forModelParams,fxModelParams,fxVol,correlMatrix,
                *(localVariances[i]));

		    variances[i+1]		= new ARM_GP_TriangularMatrix(factorNb,0.0);

            ComputeIntegratedFwdFxVCV(0.0,nextStep,
                domModelParams,forModelParams,fxModelParams,fxVol,correlMatrix,
                *(variances[i+1]));
        }
		else if (fxParams)
		{
			variances[i+1] = new ARM_GP_TriangularMatrix(factorNb,0.0);
			if(isQIntegratedVersion)
			{
				ComputeIntegratedSpotFxVCVInQModel(step,nextStep,
					domModelParams,forModelParams,fxParams,correlMatrix,
					*(localVariances[i]));


				ComputeIntegratedSpotFxVCVInQModel(0.0,nextStep,
					domModelParams,forModelParams,fxParams,correlMatrix,
					*(variances[i+1]));
			}
			else
			{
				/// Simple discretized VCV
				dyf = (nextStep-step)/K_YEAR_LEN;
				yf = nextStep/K_YEAR_LEN;
				vcv = localVariances[i];
				(*vcv)(DomModel,DomModel)	= dyf;
				(*vcv)(ForModel,ForModel)	= dyf;
				(*vcv)(FxModel,FxModel)	= dyf;
				(*vcv)(DomModel,ForModel)	=correlMatrix(DomModel,ForModel) * dyf,	(*vcv)(ForModel,DomModel) = (*vcv)(DomModel,ForModel);
				(*vcv)(DomModel,FxModel)	=correlMatrix(DomModel,FxModel) * dyf,	(*vcv)(FxModel,DomModel)  = (*vcv)(DomModel,FxModel);
				(*vcv)(ForModel,FxModel)	=correlMatrix(ForModel,FxModel) * dyf,   (*vcv)(FxModel,ForModel)  = (*vcv)(ForModel,FxModel);

				vcv = variances[i+1];
				(*vcv)(DomModel,DomModel)	= yf;
				(*vcv)(ForModel,ForModel)	= yf;
				(*vcv)(FxModel,FxModel)	= yf;
				(*vcv)(DomModel,ForModel)	=correlMatrix(DomModel,ForModel) * yf,	(*vcv)(ForModel,DomModel) = (*vcv)(DomModel,ForModel);
				(*vcv)(DomModel,FxModel)	=correlMatrix(DomModel,FxModel) * yf,	(*vcv)(FxModel,DomModel)  = (*vcv)(DomModel,FxModel);
				(*vcv)(ForModel,FxModel)	=correlMatrix(ForModel,FxModel) * yf,	(*vcv)(FxModel,ForModel)  = (*vcv)(ForModel,FxModel);
			}
		}
        else
        {
            ComputeIntegratedSpotFxVCV(step,nextStep,
                domModelParams,forModelParams,fxModelParams,correlMatrix,
                *(localVariances[i]));

		    variances[i+1]		= new ARM_GP_TriangularMatrix(factorNb,0.0);

            ComputeIntegratedSpotFxVCV(0.0,nextStep,
                domModelParams,forModelParams,fxModelParams,correlMatrix,
                *(variances[i+1]));
        }

        fxLocalVariances[i] = static_cast< ARM_GP_Matrix* >(localVariances[i]->Clone());

		step=nextStep;
	}


    /// Save local variance for Fx diffusion (take care, Set...() never clone)
	(*modelMap)[FxModel]->Model()->SetModelStateLocalVars(fxLocalVariances);

	//CC_Ostringstream os2;
	//os2 << "2IRFX : End Compute VCV\n";
	//ARM_TheEventViewer.Instance()->AddToMessage(os2.str());
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routine : NumMethodStateGlobalVariances
///	Returns : void
///	Action  : computes the integrated covariances
////////////////////////////////////////////////////

void ARM_2IRFXModel::NumMethodStateGlobalVariances( const std::vector<double>& timeSteps,
	ARM_MatrixVector& variances ) const
{
	ARM_MatrixVector localVariances;

	ComputeMeanRevertingVariances(
		timeSteps,
		localVariances,
		variances );
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routine : AdjustVol
///	Returns : void
///	Action  : interpolate model vols w.r.t. the schedule
///           of the numerical method
////////////////////////////////////////////////////
void ARM_2IRFXModel::AdjustVol()
{
    const ARM_ModelNameMap& modelMap = * GetModelMap();

    /// Only supported by the tree numerical method (at the moment)
    ARM_TreeBase* treebase = dynamic_cast<ARM_TreeBase*>( &* GetNumMethod() );
    if( !treebase )
	    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": the vol adjustment to schedule is only supported by a ARM_TreeBase numerical method");

	ARM_GP_VectorPtr treeTimes = treebase->ComputeTimeSteps( *this );

    ARM_PricingModelPtr domModel = modelMap[DomModel]->Model();
    ARM_PricingModelPtr forModel = modelMap[ForModel]->Model();
	ARM_QModel1F_Fx& fxModel = static_cast<ARM_QModel1F_Fx&>( * modelMap[FxModel]->Model() );

    AdjustVolToSchedule(treeTimes,domModel,forModel,fxModel);
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routine : BackwardLookingInit
///	Returns : ARM_PricingStatesPtr
///	Action  : Initialisation for backward looking
///           numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_2IRFXModel::BackwardLookingInit(ARM_NumMethodPtr& numMethod, double firstInductTime)
{
	ARM_PricingStatesPtr initStates(NULL);
    const ARM_ModelNameMap& modelMap = * GetModelMap();

    ARM_TreeBase* treebase = dynamic_cast<ARM_TreeBase*>( &*numMethod );
	ARM_PDEMethod* pdemethod = dynamic_cast<ARM_PDEMethod*>( &*numMethod );
	
	//tree case
	ARM_NumerairePtr numeraire = GetNumeraire();
	if( treebase)
	{
		if (numeraire->GetType() == ARM_Numeraire::Cash)
		{
			/// Calibrate domestic & foreign curves independently using a 1D tree
			ARM_TruncatorBase* truncator	= treebase->GetTruncator();
			ARM_ReconnectorBase* reconnector= treebase->GetReconnector();
			ARM_SamplerBase* sampler		= treebase->GetSampler();
			ARM_SmootherBase* smoother      = treebase->GetSmoother();

			/// Build the tree schedule
			ARM_GP_VectorPtr timeSteps = sampler->GetScheduler()->ComputeTimeSteps( *this );
			size_t nbSteps = timeSteps->size();

			/// Get a 1D equivalent sampler
			ARM_SamplerBase* sampler1D = sampler->CorrespondingSampler1D();

			/// Get current numeraire
			ARM_NumerairePtr numeraire = GetNumeraire();

			ARM_Truncator1D* truncator1D;

			size_t defaultDriftOffset=0;
			std::vector<double> defaultTargetPayoffs(0);
			ARM_VectorPtrVector defaultDriftCorrection(0);
			size_t timeIdx;
			double time,nextTime;
			ARM_PricingStatesPtr nullStates(NULL);


			/// Domestic yield curve 1D calibration
			/// -----------------------------------

			/// Build a 1D equivalent tree for domestic/coupon curve
			/// Should impose the already built schedule !
			truncator1D = truncator->CorrespondingTruncator1D(DomModel);
			ARM_Tree1D* domTree1D=new ARM_Tree1D(sampler1D,truncator1D,reconnector,smoother,treebase->GetComputeSpotProbas());
			ARM_NumMethodPtr domTree(domTree1D);
			ARM_PricingModelPtr domModel = modelMap[DomModel]->Model();
			domModel->SetNumeraire(numeraire);
			ARM_NumMethodPtr oldNumMethodDomModel = domModel->GetNumMethod();
			domModel->SetNumMethod(domTree);

			/// Set the model number to its one factor & mono-currency value (=0)
			size_t modelNb = domModel->GetModelNb();
			domModel->SetModelNb( 0 );

			/// Compute targetDf on the domestic basis model
			/// Take care that  the short rate calibrated at current
			/// time step fulfils the DF(0,nextSliceTime)
			ARM_PricingModelPtr domBasisModel = modelMap[DomBasisModel]->Model();
			std::vector<double> targetDomDf(nbSteps);

			ARM_GP_VectorPtr df;
			for(timeIdx=0;timeIdx+1<nbSteps;++timeIdx)
			{
				nextTime = (*timeSteps)[timeIdx+1];
				df=GetRefModel()->DiscountFactor("",0.0,nextTime,nullStates);
				targetDomDf[timeIdx] = (*df)[0];
			}
			nextTime = (*timeSteps)[nbSteps-1] + ARM_TreeBase::LastSliceCalibrationTermInDays;
			df=GetRefModel()->DiscountFactor("",0.0,nextTime,nullStates);
			targetDomDf[nbSteps-1] = (*df)[0];

			/// Calibrate domestic drift correction under domestic risk free probability
			initStates = domTree1D->Init(*domModel,firstInductTime,defaultDriftCorrection,defaultDriftOffset,targetDomDf,timeSteps);

			/// Restore hybrid model number
			domModel->SetModelNb( modelNb );

			/// Get drift corrections of domestic curve
			ARM_VectorPtrVector prevDriftCorrections;
			domTree1D->GetDriftCorrections(prevDriftCorrections);

			/// Free memory (1D objects created from ND ones)
			delete truncator1D;

			/// Foreign yield curve 1D calibration
			/// -----------------------------------

			/// Build a 1D equivalent tree for foreign/funding curve
			/// Should impose the already built schedule !
			truncator1D = truncator->CorrespondingTruncator1D(ForModel);
			ARM_Tree1D* forTree1D = new ARM_Tree1D(sampler1D, truncator1D,reconnector,smoother,treebase->GetComputeSpotProbas());
			ARM_NumMethodPtr forTree(forTree1D);
			ARM_PricingModelPtr forModel = modelMap[ForModel]->Model();
			forModel->SetNumeraire(numeraire);
			ARM_NumMethodPtr oldNumMethodForModel = forModel->GetNumMethod();
			forModel->SetNumMethod(forTree);

			/// Set the model number to its one factor & mono-currency value (=0)
			modelNb = forModel->GetModelNb();
			forModel->SetModelNb( 0 );

			/// Compute targetDf on the foreign basis model
			/// Take care that  the short rate calibrated at current
			/// time step fulfils the DF(0,nextSliceTime)
			ARM_PricingModelPtr forBasisModel = modelMap[ForBasisModel]->Model();
			std::vector<double> forBasisDf(nbSteps);
			for(timeIdx=0;timeIdx+1<nbSteps;++timeIdx)
			{
				nextTime = (*timeSteps)[timeIdx+1];
				df=forBasisModel->DiscountFactor("",0.0,nextTime,nullStates);
				forBasisDf[timeIdx] = (*df)[0];
			}
			nextTime = (*timeSteps)[nbSteps-1] + ARM_TreeBase::LastSliceCalibrationTermInDays;
			df=forBasisModel->DiscountFactor("",0.0,nextTime,nullStates);
			forBasisDf[nbSteps-1] = (*df)[0];

			/// Calibrate foreign drift correction under foreign risk free probability
			initStates = forTree1D->Init(*forModel,firstInductTime,defaultDriftCorrection,defaultDriftOffset,forBasisDf,timeSteps);

			/// Restore hybrid model number
			forModel->SetModelNb( modelNb );

			/// Get drift corrections of foreign curve
			ARM_VectorPtrVector forDriftCor;
			forTree1D->GetDriftCorrections(forDriftCor);

			/// Free memory (1D objects created from ND ones)
			delete truncator1D;


			/// Forex 3D calibration using previous domestic & foreign corrections
			/// ------------------------------------------------------------------

			/// Restore domestic and foreign numerical method
			domModel->SetNumMethod(oldNumMethodDomModel);
			forModel->SetNumMethod(oldNumMethodForModel);

			/// Merge yield curve drift corrections
			ARM_QModel1F_Fx& fxModel = static_cast<ARM_QModel1F_Fx&>( * modelMap[FxModel]->Model() );
			const string& fxModelName = fxModel.GetModelName();


			/// The foreign process is centred then the quanto correction is added to
			/// the foreign calibrated drift correction
			const ARM_ModelParamsHW1FStd* const forModelMRParams = dynamic_cast<const ARM_ModelParamsHW1FStd* const>( forModel->GetModelParams() );
			const ARM_ModelParamsHW1FStd* const fxModelMRParams = dynamic_cast<const ARM_ModelParamsHW1FStd* const>( fxModel.GetModelParams() );
			if(!forModelMRParams || !fxModelMRParams)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
			": foreign/funding and Fx models must be mean reverting based models");

			const ARM_CurveMatrix* const correlCurveMatrix  = GetCorrelMatrix();
			ARM_GP_Matrix correlMatrix;
			double domDriftCor,correl,quantoCorr;
			for(timeIdx=0;timeIdx<nbSteps;++timeIdx)
			{
				time = (*timeSteps)[timeIdx];
				correlMatrix = correlCurveMatrix->Interpolate(time);
				domDriftCor = (*(prevDriftCorrections[timeIdx]))[DomModel];
				prevDriftCorrections[timeIdx]->resize(3);
				(*(prevDriftCorrections[timeIdx]))[DomModel]=domDriftCor; /// not sure of resizing effect

				
				correl = correlMatrix(FxModel,ForModel);
				quantoCorr = -correl * ARM_ModelParamsHW1F::HW1FStateCovariance(fxModelMRParams, forModelMRParams,0.0,time,time);

				(*(prevDriftCorrections[timeIdx]))[ForModel] = (*(forDriftCor[timeIdx]))[0] + quantoCorr;

		/****
		FILE* f=fopen("c:\\temp\\dump2IRFX.txt","a");
		double stateFxVar = static_cast< const ARM_ModelParamsQ1F* const>(fxModel.GetModelParams())->StateLocalVariance(0.0,time,time);
		double qFx = fxModel.GetModelParams()->GetModelParam( ARM_ModelParamType::QParameter).ToCurveModelParam().GetCurve()->Interpolate(time);
		double lnFxCorr = -0.5*qFx*stateFxVar; /// valid only for cst qFx otherwise it must be used in the variance integration stateFxVar
		fprintf(f,"#%3d\tt=\t%6.2lf\tQuantoCorr=\t%15.10lf\tLnFxCorr=\t%15.10lf\n",timeIdx,time,quantoCorr,lnFxCorr);
		fclose(f);
		****/
			}

			/// Get the spot(0)=forward(0,0)
			double fxSpot = fxModel.ComputeFwdAtTime(0.0);

			/// Compute target payoffs to be calibrated by the 3D tree
			std::vector<double> targetLocalPayoffs(nbSteps);
			for(timeIdx=0;timeIdx<nbSteps;++timeIdx)
			{
				time = (*timeSteps)[timeIdx];
				df=forBasisModel->DiscountFactor("",0.0,time,nullStates);
				targetLocalPayoffs[timeIdx] = fxSpot * (*df)[0];
			}

			/// Initialise the 3D tree and calibrate forex drift under domestic risk free probability
			initStates = treebase->Init(*this,firstInductTime,prevDriftCorrections,FxModel,targetLocalPayoffs,timeSteps);

			if(!treebase->GetSampler()->IsIntegratedSampling())
			{
				/// Adjust calibrated vol to tree schedule to force interpolation as used to
				/// evaluate residual FX options in Tree3F
				AdjustVolToSchedule(timeSteps,domModel,forModel,fxModel);
			}
			/// Free memory (1D objects created from ND ones)
			delete sampler1D;	
		}
		else
		{
			initStates = numMethod->Init(*this,firstInductTime);
		}
	}
	//pde case
	else if(pdemethod)
	{
		/// Temp code waiting for actual sampler use inside MCMethod
		const ARM_ModelNameMap* const modelMap = GetModelMap();

		/// Only H&W degenerated IR Q models are allowed
		ARM_QModel1F* domModel = static_cast< ARM_QModel1F* >( &*((*modelMap)[DomModel]->Model()) );
		ARM_QModel1F* forModel = static_cast< ARM_QModel1F* >( &*((*modelMap)[ForModel]->Model()) );

		if(!domModel->IsDegenerateInHW() || !forModel->IsDegenerateInHW())
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
			": IR H&W models required in 2IR + FX MC pricing");
		}

		/// Reset structures used for DF calibration (may be used in another life of the model !)
		domModel->ResetTargetPayoffs();
		forModel->ResetTargetPayoffs();

		initStates = numMethod->Init(*this,firstInductTime);
	}
	//other cases imply an error
	else
	{
	    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not a ARM_TreeBase or a ARM_PDEmethod numerical method!");
	}
    return initStates;
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routine : AdjustVolToSchedule
///	Returns : void
///	Action  : Adjust calibrated volatilities curve to a given schedule
////////////////////////////////////////////////////
void ARM_2IRFXModel::AdjustVolToSchedule(const ARM_GP_VectorPtr& timeSteps, ARM_PricingModelPtr& domModel, ARM_PricingModelPtr& forModel, ARM_PricingModel& fxModel)
{
    ARM_ModelParams* domParams = domModel->GetModelParams();
    ARM_CurveModelParam& domInitVol=domParams->GetModelParam(ARM_ModelParamType::QVol).ToCurveModelParam();

    ARM_ModelParams* forParams = forModel->GetModelParams();
    ARM_CurveModelParam& forInitVol=forParams->GetModelParam(ARM_ModelParamType::QVol).ToCurveModelParam();

    ARM_ModelParams* fxParams = fxModel.GetModelParams();
    ARM_CurveModelParam& fxInitVol=fxParams->GetModelParam(ARM_ModelParamType::QVol).ToCurveModelParam();

    size_t nbSteps = timeSteps->size();

    /// Adjust model volatilities to be consistent to special interpolation rule
    /// used in the tree : replace initial vol values by the interpolated ones
    /// at each slice of the tree.

    ARM_GP_MatrixPtr vols,d1Vols,correls;
    VolatilitiesAndCorrelations(*timeSteps,vols,d1Vols,correls);


    /// Shift upward timeSteps (and then discard time=0) because volatility are
    /// cste right interpolated while Tree3F uses them in cst left style
    size_t i,newNbSteps=nbSteps-1;
    std::vector<double> newTimeSteps(newNbSteps);
    std::vector<double> domVols(newNbSteps);
    std::vector<double> forVols(newNbSteps);
    std::vector<double> fxVols(newNbSteps);

    for(i=0;i+1<nbSteps;++i)
    {
        newTimeSteps[i]=(*timeSteps)[i+1];
        domVols[i]  = (*vols)(ARM_2IRFXModel::DomModel,i);
        forVols[i]  = (*vols)(ARM_2IRFXModel::ForModel,i);
        fxVols[i]   = (*vols)(ARM_2IRFXModel::FxModel,i);
    }

    /// Keep initial vol values if beyond last tree step
    double lastStep = newTimeSteps[newNbSteps-1];
    const std::vector<double>& domInitTimes =  domInitVol.GetCurve()->GetAbscisses();
    std::vector<double> domTimeSteps(newTimeSteps);
    for(i=0;i<domInitTimes.size();++i)
    {
        if(domInitTimes[i] > lastStep)
        {
            domTimeSteps.push_back(domInitTimes[i]);
            domVols.push_back(domInitVol.GetCurve()->GetOrdinate(i));
        }
    }

    const std::vector<double>& forInitTimes =  forInitVol.GetCurve()->GetAbscisses();
    std::vector<double> forTimeSteps(newTimeSteps);
    for(i=0;i<forInitTimes.size();++i)
    {
        if(forInitTimes[i] > lastStep)
        {
            forTimeSteps.push_back(forInitTimes[i]);
            forVols.push_back(forInitVol.GetCurve()->GetOrdinate(i));
        }
    }

    const std::vector<double>& fxInitTimes =  fxInitVol.GetCurve()->GetAbscisses();
    std::vector<double> fxTimeSteps(newTimeSteps);
    for(i=0;i<fxInitTimes.size();++i)
    {
        if(fxInitTimes[i] > lastStep)
        {
            fxTimeSteps.push_back(fxInitTimes[i]);
            fxVols.push_back(fxInitVol.GetCurve()->GetOrdinate(i));
        }
    }

    /// Replace new volatility curves in stochastic models
    std::vector<double> domVolLowerBound(domTimeSteps.size(),(*(domInitVol.GetLowerBound()))[0]);
    std::vector<double> domVolUpperBound(domTimeSteps.size(),(*(domInitVol.GetUpperBound()))[0]);
    ARM_CurveModelParam domVol(ARM_ModelParamType::QVol,&domVols,&domTimeSteps,"DOMQVOL","STEPUPRIGHT",
                               &domVolLowerBound,&domVolUpperBound);
    domParams->SetModelParam(&domVol);

    std::vector<double> forVolLowerBound(forTimeSteps.size(),(*(forInitVol.GetLowerBound()))[0]);
    std::vector<double> forVolUpperBound(forTimeSteps.size(),(*(forInitVol.GetUpperBound()))[0]);
    ARM_CurveModelParam forVol(ARM_ModelParamType::QVol,&forVols,&forTimeSteps,"FORQVOL","STEPUPRIGHT",
                               &forVolLowerBound,&forVolUpperBound);
    forParams->SetModelParam(&forVol);

    std::vector<double> fxVolLowerBound(fxTimeSteps.size(),(*(fxInitVol.GetLowerBound()))[0]);
    std::vector<double> fxVolUpperBound(fxTimeSteps.size(),(*(fxInitVol.GetUpperBound()))[0]);
    ARM_CurveModelParam fxVol(ARM_ModelParamType::QVol,&fxVols,&fxTimeSteps,"FXQVOL","STEPUPRIGHT",
                              &fxVolLowerBound,&fxVolUpperBound);
    fxParams->SetModelParam(&fxVol);

/****
FILE* f=fopen("c:\\temp\\dump2IRFX.txt","a");
fprintf(f,"Domestic, Foreign and Fx vols adjusted to diffusion schedule\n");
fprintf(f,"------------------------------------------------------------\n");
for(i=0;i<newNbSteps;++i)
    fprintf(f,"%6.2lf\t%15.10lf\t%15.10lf\t%15.10lf\n",newTimeSteps[i],100*domVols[i],100*forVols[i],100*fxVols[i]);
fprintf(f,"\n");
fclose(f);
****/

}


////////////////////////////////////////////////////
///	Class  : ARM_2IRFXModel
///	Routine: TreeStatesToModelStates
///	Returns: void
///	Action : Convert the num method states to model states
////////////////////////////////////////////////////
void ARM_2IRFXModel::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
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
	else if (pdemethod)
	{
		//we do nothing because the model sates have already been set before with the 
		//function PdeStatesToModelStates
	}
	//other cases imply an error
	else
	{
	    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not a ARM_TreeBase or a ARM_PDEmethod numerical method!");
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_2IRFXModel
///	Routine: PdeStatesToModelStates
///	Returns: void
///	Action : Convert the num method states to model states
////////////////////////////////////////////////////
void ARM_2IRFXModel::PdeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex, double lambda) const
{   
#if defined(__GP_STRICT_VALIDATION)
	if(timeIndex<0)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index is negative!");
	if( timeIndex > GetNumMethod()->GetTimeSteps()->size()-1 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index bigger than max size!");
#endif

	ARM_GP_MatrixPtr modelStates( static_cast<ARM_GP_Matrix*> (states->GetNumMethodStates()->Clone() ));

	size_t N = modelStates->cols();
	const ARM_ModelNameMap* const modelMap	= GetModelMap();

	//models
	ARM_PricingModelPtr domModel     = (*modelMap)[DomModel]->Model();
	ARM_PricingModelPtr forModel     = (*modelMap)[ForModel]->Model();
	ARM_QModel1F_Fx* fxModel		 = dynamic_cast<ARM_QModel1F_Fx*>(&*(*modelMap)[FxModel]->Model());

	ARM_ZeroCurvePtr& DomZcCurve = domModel->GetZeroCurve();
	ARM_ZeroCurvePtr& ForZcCurve = forModel->GetZeroCurve();


	//date 
	double ti = GetNumMethod()->GetTimeStep(timeIndex);

	ARM_GP_Matrix correlMatrix = GetCorrelMatrix()->Interpolate(ti);

	//correlation
	double rhoxy=correlMatrix(0,1);
	double rhoyz=correlMatrix(1,2);
	double rhozx=correlMatrix(2,0);

	if(lambda==-2.5)//Centred LogFXSpot level 2(interest rate totally centred)
	{
		//Cash numeraire case
		if(GetNumeraire()->GetType() == ARM_Numeraire::Cash)
		{
			const ARM_ModelParamsHW1FStd* const domModelMRParams = static_cast<const ARM_ModelParamsHW1FStd* const>( domModel->GetModelParams() );
			const ARM_ModelParamsHW1FStd* const forModelMRParams = static_cast<const ARM_ModelParamsHW1FStd* const>( forModel->GetModelParams() );
			const ARM_ModelParamsHW1FStd* const fxModelMRParams = dynamic_cast<const ARM_ModelParamsHW1FStd* const>( fxModel->GetModelParams() );

			//mean of Xs in lognormal case
			const ARM_ModelParamsQ1F* const fxParams   = static_cast<const ARM_ModelParamsQ1F* const>( fxModel->GetModelParams() );
			const ARM_CurveModelParam& fxVolParam = static_cast< const ARM_CurveModelParam& >(fxParams->GetModelParam(ARM_ModelParamType::QVol));
			double fwdFx0 =fxModel->ComputeFwdAtTime(ti);
	        double domintegPhi = - ARM_ModelParamsHW1F::HW1FStateZcCovariance( domModelMRParams, domModelMRParams, 0.0, ti, ti, ti );
	        double forintegPhi = - ARM_ModelParamsHW1F::HW1FStateZcCovariance( forModelMRParams, forModelMRParams, 0.0, ti, ti, ti );
			double LNmean = log(fwdFx0) + domintegPhi - forintegPhi;/* -0.5*ARM_ModelParamsQ1F::Q1FStateCovariance(fxParams,fxParams,0.0,ti,ti,true);*/
			 
			//loop
			for (int K = 0; K < N; ++K)
			{
				//no correction for the dom factor
				//correction for the for factor

				//correction of the fx factor
				(*modelStates)(2,K)+=LNmean;//for the moment
				//exponnentiel of the third factor equals SpotFX
				(*modelStates)(2,K)=exp((*modelStates)(2,K));
			}

		}
		else if(GetNumeraire()->GetType() == ARM_Numeraire::TerminalZc)
		{

		}
		//other cases imply an error
		else
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not a Cash or TerminalZc numeraire!");
		}
	}

	if(lambda==-2)//Centred LogFXSpot
	{
		//Cash numeraire case
		if(GetNumeraire()->GetType() == ARM_Numeraire::Cash)
		{
			//mean of Xf
			const ARM_ModelParamsHW1FStd* const forModelMRParams = static_cast<const ARM_ModelParamsHW1FStd* const>( forModel->GetModelParams() );
			const ARM_ModelParamsHW1FStd* const fxModelMRParams = dynamic_cast<const ARM_ModelParamsHW1FStd* const>( fxModel->GetModelParams() );

			//mean of Xs in lognormal case
			const ARM_ModelParamsQ1F* const fxParams   = static_cast<const ARM_ModelParamsQ1F* const>( fxModel->GetModelParams() );
			const ARM_CurveModelParam& fxVolParam = static_cast< const ARM_CurveModelParam& >(fxParams->GetModelParam(ARM_ModelParamType::QVol));
			double fwdFx0 =fxModel->ComputeFwdAtTime(ti);
			double LNmean = log(fwdFx0);/* -0.5*ARM_ModelParamsQ1F::Q1FStateCovariance(fxParams,fxParams,0.0,ti,ti,true);*/
			 
			//loop
			for (int K = 0; K < N; ++K)
			{
				//no correction for the dom factor
				//correction for the for factor

				//correction of the fx factor
				(*modelStates)(2,K)+=LNmean;//for the moment
				//exponnentiel of the third factor equals SpotFX
				(*modelStates)(2,K)=exp((*modelStates)(2,K));
			}

		}
		else if(GetNumeraire()->GetType() == ARM_Numeraire::TerminalZc)
		{

		}
		//other cases imply an error
		else
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not a Cash or TerminalZc numeraire!");
		}
	}


	if(lambda==-1)//Centred Qmapping
	{
		//Cash numeraire case
		if(GetNumeraire()->GetType() == ARM_Numeraire::Cash)
		{
			//mean of Xf
			const ARM_ModelParamsHW1FStd* const forModelMRParams = static_cast<const ARM_ModelParamsHW1FStd* const>( forModel->GetModelParams() );
			const ARM_ModelParamsHW1FStd* const fxModelMRParams = dynamic_cast<const ARM_ModelParamsHW1FStd* const>( fxModel->GetModelParams() );
			double quantoCorr = -rhoyz * ARM_ModelParamsHW1F::HW1FStateCovariance(fxModelMRParams, forModelMRParams,0.0,ti,ti);

			//mean of Xs in lognormal case
			const ARM_ModelParamsQ1F* const fxParams   = static_cast<const ARM_ModelParamsQ1F* const>( fxModel->GetModelParams() );
			const ARM_CurveModelParam& fxVolParam = static_cast< const ARM_CurveModelParam& >(fxParams->GetModelParam(ARM_ModelParamType::QVol));
			double LNmean = -0.5*ARM_ModelParamsQ1F::Q1FStateCovariance(fxParams,fxParams,0.0,ti,ti,false);
			 
			//loop
			for (int K = 0; K < N; ++K)
			{
				//no correction for the dom factor
				//correction of the for factor			
				(*modelStates)(1,K)+=quantoCorr;
				//correction of the fx factor
				(*modelStates)(2,K)+=LNmean;//for the moment
			}

		}
		else if(GetNumeraire()->GetType() == ARM_Numeraire::TerminalZc)
		{

		}
		//other cases imply an error
		else
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not a Cash or TerminalZc numeraire!");
		}
	}

	if(lambda==2)//Uncentred LogSpotFX
	{
		for (int K = 0; K < N; ++K)
		{
			//S0*exponnentiel of the third factor equals SpotFX
			double Fx0 =fxModel->ComputeFwdAtTime(0.0);
			(*modelStates)(2,K)=Fx0*exp((*modelStates)(2,K));
		}
	}
	if(lambda==2.5)//Uncentred LogSpotFX with freezing of gamma in the quantodrift
	{
		for (int K = 0; K < N; ++K)
		{
			//S0*exponnentiel of the third factor equals SpotFX
			double Fx0 =fxModel->ComputeFwdAtTime(0.0);
			(*modelStates)(2,K)=Fx0*exp((*modelStates)(2,K));
		}
	}
	if(lambda==3)//Uncentred SpotFX
	{
		for (int K = 0; K < N; ++K)
		{
			//S0*the third factor equals SpotFX
			double Fx0 =fxModel->ComputeFwdAtTime(0.0);
			(*modelStates)(2,K)=Fx0*(*modelStates)(2,K);
		}
	}
	states->SetModelStates(modelStates);
}


////////////////////////////////////////////////////
///	Class  : ARM_2IRFXModel
///	Routine: ComputeNumeraireTimes
///	Returns: 
///	Action : compute the corresonding numeraire times
////////////////////////////////////////////////////
ARM_VectorPtr ARM_2IRFXModel::ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const
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
///	Class   : ARM_2IRFXModel
///	Routine : ForwardLookingInit
///	Returns : ARM_PricingStatesPtr
///	Action  : Initialisation for forward looking
///           numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_2IRFXModel::ForwardLookingInit(ARM_NumMethodPtr& numMethod, double firstInductTime)
{
    /// Validate numeraire
	ARM_NumerairePtr numeraire = GetNumeraire();
    if( numeraire->GetType() != ARM_Numeraire::RollingEvent &&
        numeraire->GetType() != ARM_Numeraire::Cash &&
        numeraire->GetType() != ARM_Numeraire::RollingCash )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
        ": with forward looking methods only Cash or RollingEvent numeraire supported by 2IRFX model");

    /// Temp code waiting for actual sampler use inside MCMethod
    const ARM_ModelNameMap* const modelMap = GetModelMap();

    /// Only H&W degenerated IR Q models are allowed
    ARM_QModel1F* domModel = static_cast< ARM_QModel1F* >( &*((*modelMap)[DomModel]->Model()) );
    ARM_QModel1F* forModel = static_cast< ARM_QModel1F* >( &*((*modelMap)[ForModel]->Model()) );

    if(!domModel->IsDegenerateInHW() || !forModel->IsDegenerateInHW())
    {
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
        ": IR H&W models required in 2IR + FX MC pricing");
    }

    /// Reset structures used for DF calibration (may be used in another life of the model !)
    domModel->ResetTargetPayoffs();
    forModel->ResetTargetPayoffs();

	ARM_QModel1F_Fx* fxModel = dynamic_cast<ARM_QModel1F_Fx*>(&*((*modelMap)[FxModel]->Model()));
	bool isLnFwdFx = (dynamic_cast<const ARM_ModelParamsHW1FStd* const>( fxModel->GetModelParams() ))->IsLn();

    ARM_MCMethod* mcMethod = dynamic_cast< ARM_MCMethod* >(&*numMethod);

    ARM_MeanRevertingSamplerND* samplerMR = dynamic_cast<ARM_MeanRevertingSamplerND* >(&*(mcMethod->GetSampler()));
    ARM_NormalCentredSamplerND* samplerNC = dynamic_cast<ARM_NormalCentredSamplerND* >(&*(mcMethod->GetSampler()));
    if(!mcMethod || (isLnFwdFx && !samplerMR) ||
		(!isLnFwdFx && !samplerMR && !samplerNC))
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
        ": 2IRFX + MC method requires a Mean Reverting Sampler (LN or Skew Fx) or a Normal Centred Sampler (Skew Fx)");

	/// In all cases Forex is diffused in MC (fwd or spot Fx)
	fxModel->SetIsForexDiffusion(true);
	if(!isLnFwdFx)
	{
		if(samplerMR)
		{
			/// Previous MC diffusion method
			fxModel->SetIntegratedVersion(true);
		}
		else
		{
			/// Simple Euler discretisation scheme
			fxModel->SetIntegratedVersion(false);

			/// Merge all model time steps
			std::vector<double> *modelTimeSteps,*timeSteps;

			/// Save the last event time
			double lastEventTime = mcMethod->GetTimeStep(mcMethod->GetTimeSteps()->size()-1);

			/// No computation via dom, for and fx models ComputeModelTimes( timeInfos )
			/// because of side effects if this fct is upgraded in each model...!
			const std::vector<double>& domVolTimes = (static_cast< const ARM_CurveModelParam& >(domModel->GetModelParams()->GetModelParam(ARM_ModelParamType::QVol).ToCurveModelParam())).GetCurve()->GetAbscisses();
			timeSteps = MergeSortedVectorNoDuplicates(*(mcMethod->GetTimeSteps()),domVolTimes);

			const std::vector<double>& forVolTimes = (static_cast< const ARM_CurveModelParam& >(forModel->GetModelParams()->GetModelParam(ARM_ModelParamType::QVol).ToCurveModelParam())).GetCurve()->GetAbscisses();
			modelTimeSteps = MergeSortedVectorNoDuplicates(*timeSteps,forVolTimes);
			delete timeSteps;
			timeSteps = modelTimeSteps;

			const std::vector<double>& fxVolTimes = (static_cast< const ARM_CurveModelParam& >(fxModel->GetModelParams()->GetModelParam(ARM_ModelParamType::QVol).ToCurveModelParam())).GetCurve()->GetAbscisses();
			modelTimeSteps = MergeSortedVectorNoDuplicates(*timeSteps,fxVolTimes);
			delete timeSteps;
			timeSteps = modelTimeSteps;

			const std::vector<double>& fxQTimes = (static_cast< const ARM_CurveModelParam& >(fxModel->GetModelParams()->GetModelParam(ARM_ModelParamType::QParameter).ToCurveModelParam())).GetCurve()->GetAbscisses();
			modelTimeSteps = MergeSortedVectorNoDuplicates(*timeSteps,fxQTimes);
			delete timeSteps;
			timeSteps = modelTimeSteps;

			std::vector<double> correlTimes = GetCorrelMatrix()->GetTimes();
			modelTimeSteps = MergeSortedVectorNoDuplicates(*timeSteps,correlTimes);

			/// Erase model steps beyong last event time
			std::vector<double>::iterator pos=modelTimeSteps->end();
			--pos;
			size_t nbSteps = modelTimeSteps->size();
			while(pos!=modelTimeSteps->begin() && (*pos)!=lastEventTime)
				--pos;
			++pos;
			modelTimeSteps->erase(pos,modelTimeSteps->end());

			numMethod->SetTimeSteps(*modelTimeSteps);

			delete timeSteps;
			delete modelTimeSteps;
		}
	}

    return mcMethod->Init(*this,firstInductTime);
}


////////////////////////////////////////////////////
///	Class  : ARM_2IRFXModel
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : Convert the num method states to model states
////////////////////////////////////////////////////
void ARM_2IRFXModel::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
    ARM_GP_MatrixPtr numStates = states->GetNumMethodStates();
    ARM_GP_MatrixPtr modelStates = states->GetModelStates();
    size_t j,nbStates = numStates->cols();

	double t		= GetNumMethod()->GetTimeStep(timeIndex);
	double nextt	= GetNumMethod()->GetTimeStep(timeIndex+1);
	double yft		= t/K_YEAR_LEN;
	double yfNextt	= nextt/K_YEAR_LEN;
	double dt		= yfNextt-yft;

    const ARM_ModelNameMap* const modelMap = GetModelMap();

	ARM_GP_Matrix correlMatrix = GetCorrelMatrix()->Interpolate(t);
    double corForFx    = correlMatrix(ForModel,FxModel);

	ARM_PricingModelPtr domModel	= (*modelMap)[DomModel]->Model();
	ARM_PricingModelPtr forModel	= (*modelMap)[ForModel]->Model();
	ARM_QModel1F_Fx* fxModel		= dynamic_cast<ARM_QModel1F_Fx*>(&*((*modelMap)[FxModel]->Model()));
	bool isLnFwdFx = (dynamic_cast<const ARM_ModelParamsHW1FStd* const>( fxModel->GetModelParams() ))->IsLn();

	double mrsDom,mrsFor,qVolDom,qVolFor,qVolFx,qFx;
	std::vector<double> quantoDrift;
	if(!isLnFwdFx)
	{
		mrsDom	= domModel->GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion).GetValue(t);
		mrsFor	= forModel->GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion).GetValue(t);
		qVolDom	= domModel->GetModelParams()->GetModelParam(ARM_ModelParamType::QVol).GetValue(t);
		qVolFor	= forModel->GetModelParams()->GetModelParam(ARM_ModelParamType::QVol).GetValue(t);
		qVolFx	= fxModel->GetModelParams()->GetModelParam(ARM_ModelParamType::QVol).GetValue(t);
		qFx	= fxModel->GetModelParams()->GetModelParam(ARM_ModelParamType::QParameter).GetValue(t);

		/// State dependent quanto correction
		quantoDrift.resize(nbStates);
		double fx,norVolFx,lnVolFx,driftFact = -corForFx*qVolFor*dt;
		double volRelShift = qFx * qVolFx;
		double volAbsShift = (qVolFx - volRelShift) * (*(fxModel->GetPrecomputedFwds()))[timeIndex];
		for(j=0;j<nbStates;++j)
		{
			fx = (*modelStates)(FxModel,j);
			norVolFx = volRelShift * fx + volAbsShift;
			if(fabs(norVolFx)>fabs(fx)*MAX_LNVOL_FX)
				lnVolFx = ((norVolFx>0 && fx>=0) || (norVolFx<0 && fx<=0) ? MAX_LNVOL_FX : -MAX_LNVOL_FX);
			else
				lnVolFx = norVolFx/fx;
			quantoDrift[j] = driftFact * lnVolFx;
		}
	}

    /// First step we use IR model states at current time to get short rates
    /// and then S(t) is diffused from ti to ti+1
    (*(GetModelMap()))[FxModel]->Model()->MCModelStatesFromToNextTime(states,timeIndex);


	if(isLnFwdFx || GetNumeraire()->GetType() != ARM_Numeraire::Cash)
	{
		/// Second step for IR processes model states are identical to numerical method states (valid for timeIndex+1)
		for(j=0;j<nbStates;++j)
		{
			(*modelStates)(DomModel,j)  =  (*numStates)(DomModel,j);
			(*modelStates)(ForModel,j)  =  (*numStates)(ForModel,j);
		}
	}
	else
	{
		/// Second step : Euler scheme on domestic & foreign models with stochastic quanto correction
		double xDom,xFor;
		for(j=0;j<nbStates;++j)
		{
			xDom = (*modelStates)(DomModel,j);
			(*modelStates)(DomModel,j) = xDom*(1.0-mrsDom*dt) + qVolDom * (*numStates)(DomModel,j);

			xFor = (*modelStates)(ForModel,j);
			(*modelStates)(ForModel,j) = xFor*(1.0-mrsFor*dt) + quantoDrift[j] + qVolFor * (*numStates)(ForModel,j);

		}
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routines: IntegratedGlobalDrifts
///	Returns : void
///	Action  : computes global drifts of diffused processes
////////////////////////////////////////////////////
void ARM_2IRFXModel::IntegratedGlobalDrifts(
		const std::vector<double>& timeSteps,
        ARM_GP_MatrixPtr& drifts)
{
    /// - domestic process is a centred MR under domestic risk free probability => non global drift
    /// - foreign process is a centred MR under foreign risk free probability => global drift due to quanto
    ///     effect (which is globally added to the calibrated foreign drift correction)
    /// - forex process has a markovian drift that exactly compensate in expectation the forward fx for
    ///     a q=1 model  and a q=0 model on IR sides
    ///     (the Ito correction -1/2.q.vol.vol.t is added global by the Q mapping function).
    ///     The fx process is then assumed centred.
    ///     But this approximation could be revisited for q < 1 on FX model and/or q > 0 on IR sides

	drifts	= ARM_GP_MatrixPtr( NULL ); // to say that all processes are centred
}

////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routines: IntegratedLocalDrifts
///	Returns : void
///	Action  : computes the relative and absolute drift
////////////////////////////////////////////////////

void ARM_2IRFXModel::IntegratedLocalDrifts( const std::vector<double>& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const
{
	ARM_MultiAssetsModel::IntegratedLocalDrifts( timeSteps, relativeDrifts, absoluteDrifts );

    if( GetNumMethod().IsNull() || 
		(GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_FWDLOOKING) || 
		(GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_FWDBCKWDLOOKING))
    {
        /// Add local corrections depending of numeraire
	    AddIntegratedLocalCorrections( timeSteps, absoluteDrifts );
    }
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routines: AddIntegratedLocalCorrections
///	Returns : void
///	Action  : computes integrated local corrections
///           to be added to absolute drifts to diffuse
///           processes w.r.t. numeraire
////////////////////////////////////////////////////
void ARM_2IRFXModel::AddIntegratedLocalCorrections( 
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
    const ARM_ModelParamsQ1F* const fxModelMRParamsQ1F = dynamic_cast<const ARM_ModelParamsQ1F* const>( fxModel->GetModelParams() );
	const ARM_ModelParamsHW1FStd* const fxModelMRParams = dynamic_cast<const ARM_ModelParamsHW1FStd* const>( fxModel->GetModelParams() );

	bool isLnFwdFx = fxModelMRParams->IsLn();
	bool isQ1F = false;
	ARM_CurveModelParam fxVolParam;
	if ((fxModelMRParams->DoesModelParamExist(ARM_ModelParamType::QVol)) && 
		(fxModelMRParams->DoesModelParamExist(ARM_ModelParamType::QParameter)))
	{
		fxVolParam = static_cast< ARM_CurveModelParam >(fxModelMRParams->GetModelParam(ARM_ModelParamType::QVol).ToCurveModelParam());
		isQ1F = true;
	}
	else if ((fxModelMRParams->DoesModelParamExist(ARM_ModelParamType::Volatility)) && 
		(fxModelMRParams->DoesModelParamExist(ARM_ModelParamType::Beta)))
	{
		fxVolParam = static_cast< ARM_CurveModelParam >(fxModelMRParams->GetModelParam(ARM_ModelParamType::Volatility).ToCurveModelParam());
	}
    

	size_t i,nbSteps =timeSteps.size();

	double corrDomFor;
    double corrDomFx;
    double corrForFx;
    double lastTime=0.0,time;

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

            if(isLnFwdFx)
            {
                /// Fwd fx is diffused : Qdom(spot) <-> Qdom correction on forex process
                (*absoluteDrifts)(i,FxModel) -= corrDomFx * ARM_ModelParamsHW1F::HW1FEqFxZcCovariance(fxVolParam,domModelMRParams,lastTime,time,time);
            }
            else
                /// Spot Fx is diffused : non drift under Qdom

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

            if(isLnFwdFx) // Fwd fx is diffused
            {
                /// Qfor <-> Qdom (quanto) correction and
                /// Qdom <-> Qdom(spot) correction on foreign process
                (*absoluteDrifts)(i,ForModel) += -corrForFx  * ARM_ModelParamsHW1F::HW1FStateCovariance(fxModelMRParams,forModelMRParams,lastTime,time,time)
                                                 +corrDomFor * ARM_ModelParamsHW1F::HW1FStateZcCovariance(forModelMRParams,domModelMRParams,lastTime,time,time,time);

                /// No drift on Fwd fx
            }
            else // Spot fx is diffused
            {
                /// For foreign process quanto correction is markovian but here we use
                /// Piterbarg's approx from 0 to ti+1 for spot Fx vol
                (*absoluteDrifts)(i,ForModel) += -corrForFx  * ARM_ModelParamsHW1F::HW1FStateCovariance(fxModelMRParams,forModelMRParams,lastTime,time,time)
                                                 +corrDomFor * ARM_ModelParamsHW1F::HW1FStateZcCovariance(forModelMRParams,domModelMRParams,lastTime,time,time,time);

				if (isQ1F)
				{
					(*absoluteDrifts)(i,FxModel) += corrDomFx  * ARM_ModelParamsQ1F::Q1FStateZcCovariance(fxModelMRParamsQ1F,domModelMRParams,lastTime,time,time,time);
				}
				else
				{
					(*absoluteDrifts)(i,FxModel) += corrDomFx  * ARM_ModelParamsHW1FStd::HW1FStateZcCovariance(fxModelMRParams,domModelMRParams,lastTime,time,time,time);
				}
            }

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
        ": only Cash or RollingCash numeraire supported by 2IR + FX model");
    }
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routines: EulerLocalDrifts
///	Returns : void
///	Action  : computes the relative and absolute drift
////////////////////////////////////////////////////
void ARM_2IRFXModel::EulerLocalDrifts( const std::vector<double>& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const
{
	ARM_MultiAssetsModel::EulerLocalDrifts( timeSteps, relativeDrifts, absoluteDrifts );

    if((GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_FWDLOOKING) && (GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_FWDBCKWDLOOKING))
    {
        /// Domestic cash numeraire is assumed here !

	    size_t nbSteps =timeSteps.size();
	    const ARM_ModelNameMap* const modelMap	= GetModelMap();
        std::vector<double> yfSteps(timeSteps);
        yfSteps /= K_YEAR_LEN;

	    /// Quanto adjustment
	    ARM_ModelParam& modelParamFxVol			= (*modelMap)[FxModel]->Model()->GetModelParams()->GetModelParam(ARM_ModelParamType::QVol);
	    ARM_ModelParam& modelParamForIRVol		= (*modelMap)[ForModel]->Model()->GetModelParams()->GetModelParam(ARM_ModelParamType::QVol);
	    ARM_GP_Matrix correlMatrix;
	    for( size_t j=0; j<nbSteps-1; ++j )
	    {
			correlMatrix = GetCorrelMatrix()->Interpolate(timeSteps[j]);

		    (*absoluteDrifts)(j,ForModel) = -modelParamFxVol.GetValue(timeSteps[j]) 
			    * modelParamForIRVol.GetValue(timeSteps[j]) 
			    * correlMatrix(FxModel,ForModel) 
			    * (yfSteps[j+1]-yfSteps[j]);

	    }
    }
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routines: IntegratedMarkovianDrift
///	Returns : ARM_GP_MatrixPtr
///	Action  : computes the Integrated Markovian drift
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_2IRFXModel::IntegratedMarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates, const ARM_GP_VectorPtr& driftCorrection ) const
{
	const ARM_ModelNameMap* const modelMap	= GetModelMap();
	ARM_QModel1F_Fx& fxModel = static_cast<ARM_QModel1F_Fx&>( * (*modelMap)[FxModel]->Model() );
    const ARM_ModelParamsHW1FStd* const fxModelMRParams = dynamic_cast<const ARM_ModelParamsHW1FStd* const>( fxModel.GetModelParams() );
    const ARM_Curve& qFxCurve = * (fxModelMRParams->GetModelParam(ARM_ModelParamType::QParameter).ToCurveModelParam().GetCurve());

    ARM_GP_MatrixPtr drift;
    if(qFxCurve == 1.0)
        drift = fxModel.IntegratedMarkovianDrift(timeIdx,numMethodStates,driftCorrection);
    else
    {
        size_t i,j,nbStates = numMethodStates->cols(),nbFactors = driftCorrection->size();

#if defined(__GP_STRICT_VALIDATION)
	if(numMethodStates->rows() != nbFactors)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "inconsistency between factor number in state and calibrated drift size");
#endif

        double t = GetNumMethod()->GetTimeStep(timeIdx);
        double T = GetNumMethod()->GetTimeStep(timeIdx+1);

        /// Set drift corrections to ModelStates to prepare calibration of Zc at current slice
        ARM_PricingStatesPtr states(new ARM_PricingStates);
	    for(j=0;j<nbFactors;++j)
            for(i=0;i<nbStates;++i)
                (*numMethodStates)(j,i) += (*driftCorrection)[j];
        states->SetModelStates(numMethodStates);

        /// ZcDom(t,T) & ZcFor(t,T) will be calibrated to compute
        /// if necessary and accurately Fwd(t,T)=Spot(t).ZcFor(t,T)/ZcDom(t,T)
        InitDomCalibration(t,T,states);
        InitForCalibration(t,T,states);

        /// Initialise markovian drifts matrix and compute values for
        /// the underlying forex process
        drift = fxModel.IntegratedMarkovianDrift(timeIdx,numMethodStates,driftCorrection);


/**** No quanto markovian drift correction because already done through the above fct
        ARM_PricingModelPtr forModel = (*modelMap)[ForModel]->Model();
        const ARM_ModelParamsHW1FStd* const forModelMRParams = dynamic_cast<const ARM_ModelParamsHW1FStd* const>( forModel->GetModelParams() );

        /// Compute markovian drifts of the underlying foreign process (quanto effect)
        /// Integral is cut w.r.t. Q parameter schedule (assumed to be stepwise right constant)
        const std::vector<double>& qTimes  = qFxCurve.GetAbscisses();
        const std::vector<double>& qValues = qFxCurve.GetOrdinates();

        /// Find Tn(t)-1 < t <= Tn(t) and Tn(T)-1< T <= Tn(T)
        size_t idx,nbQ = qTimes.size();
        for(idx=0;idx<nbQ;++idx)
            if(t <= qTimes[idx] + K_NEW_DOUBLE_TOL) break;
        size_t nt=idx;
        for(;idx<nbQ;++idx)
            if(T <= qTimes[idx] + K_NEW_DOUBLE_TOL) break;
        size_t nT=idx;

        /// Loop over constant Q parameter interval to compute the markovian quanto drift :
        /// cov(tp-tp+1,XFor,XFx)*(1-Q(tp))*{1-S(0,tp)/S(tp)}
        double lastt=t,nextt,q,driftCoef;
	    const ARM_GP_TriangularMatrix* const correlMatrix  = GetUsedCorrelMatrix();
        double correl = (*correlMatrix)(FxModel,ForModel);
        for(idx=nt;idx<=nT;++idx)
        {
            nextt=(idx<nT ? qTimes[idx] : T);
            q=qValues[idx<nbQ ? idx : nbQ-1];
            driftCoef = (1-q) * correl * ARM_ModelParamsHW1F::HW1FStateCovariance(fxModelMRParams,forModelMRParams,lastt,nextt,T);
            if(fabs(driftCoef) > ARM_NumericConstants::ARM_TOLERENCE)
            {
                for(i=0;i<nbStates;++i)
                    (*drift)(ForModel,i) += driftCoef * (1 - 1.0/(fxModel.MappingFunction((*numMethodStates)(FxModel,i),1.0,q)));
            }
            lastt = nextt;
        }
****/

/****
        /// Restore states without drift corrections
	    for(j=0;j<nbFactors;++j)
            for(i=0;i<nbStates;++i)
                (*numMethodStates)(j,i) -= (*driftCorrection)[j];
****/
    }

    return drift;
}

////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routine : ComputeFwdFXModelParam
///	Returns : ARM_VectorPtr
///	Action  : Compute q(evalTime, settlementTime) and sigma(evalTime, settlementTime) of the Fwd FX
////////////////////////////////////////////////////
ARM_VectorPtr ARM_2IRFXModel::ComputeFwdFXModelParam(
	double evalTime,
	double settlementTime,
	ARM_PricingContext* context) const
{	
 	const ARM_ModelNameMap* const modelMap	= GetModelMap();
	ARM_QModel1F_Fx* fxModel	= dynamic_cast<ARM_QModel1F_Fx*>(&*(*modelMap)[FxModel]->Model());

    if(settlementTime < evalTime - K_NEW_DOUBLE_TOL)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " settlement time is not supported");

    return fxModel->ComputeFwdFXModelParam(evalTime, settlementTime, settlementTime, context);
}

////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routines: MarkovianDrift
///	Returns : ARM_GP_MatrixPtr
///	Action  : computes the Markovian drift
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_2IRFXModel::MarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates ) const
{
	const ARM_ModelNameMap* const modelMap	= GetModelMap();
	ARM_PricingModelPtr model = (*modelMap)[FxModel]->Model();

    return model->MarkovianDrift(timeIdx,numMethodStates );
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routine : LocalDiscounts
///	Returns : ARM_VectorPtr
///	Action  : computes LocalDiscounts
////////////////////////////////////////////////////
ARM_VectorPtr ARM_2IRFXModel::LocalDiscounts( size_t timeIdx, double dt, 
	const ARM_PricingStatesPtr& states) const
{
	return GetRefModel()->LocalDiscounts(timeIdx,dt,states);
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routines: LocalPayoffs
///	Returns : ARM_VectorPtr
///	Action  : Local payoffs for a 3D hybrid calibration
///           refers to spot Fx at t because we want to fit
///           S(0).Bfor(0,t)=Edom[exp(-integ{0->t,rdom(s)ds}).S(t)]
////////////////////////////////////////////////////
ARM_VectorPtr ARM_2IRFXModel::LocalPayoffs(size_t timeIdx, 
	double dt, const ARM_PricingStatesPtr& states) const
{
	const ARM_ModelNameMap* const modelMap	= GetModelMap();

	ARM_QModel1F_Fx* fxModel = dynamic_cast<ARM_QModel1F_Fx*>(	&*(*modelMap)[FxModel]->Model() );
    const string& fxModelName = fxModel->GetModelName();

    /// Compute spot FX at slice time (starting also at slice time)
    double time = GetNumMethod()->GetTimeStep(timeIdx);

    /// Theoricaly settlementTime=payTime, should be equal to time + 2 business days
	return Forward(fxModelName,time,time,time,time,states);
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routines: InitDomCalibration
///	Returns : void
///	Action  : Set calibration datas for domestic Zc
///           computation
////////////////////////////////////////////////////
void ARM_2IRFXModel::InitDomCalibration(double evalTime, double maturityTime, const ARM_PricingStatesPtr& states) const
{
	if (GetNumeraire()->GetType() == ARM_Numeraire::Cash)
	{
		const ARM_ModelNameMap& modelMap = *GetModelMap();
		ARM_QModel1F& domModel = static_cast<ARM_QModel1F&>( * modelMap[DomModel]->Model() );

		/// If reference model is basis calibrate basis Zc
		double basisRatio=1.0,targetDf;
		if(IsBasisRefModel())
		{
			const ARM_ForwardMargin& domBasisModel = static_cast<ARM_ForwardMargin&>( * modelMap[DomBasisModel]->Model() );
			basisRatio = domBasisModel.ComputeMarginRatio(evalTime,maturityTime);
			targetDf = domBasisModel.GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN);
		}
		else
			targetDf = domModel.GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN);

		/// Set datas to domestic stochastic model
		domModel.SetTargetPayoff(maturityTime,targetDf,ARM_GP_VectorPtr( new std::vector<double>(1,basisRatio) ));
	}

}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routines: InitForCalibration
///	Returns : void
///	Action  : Set calibration datas for foreign Zc
///           computation
////////////////////////////////////////////////////
void ARM_2IRFXModel::InitForCalibration(double evalTime, double maturityTime, const ARM_PricingStatesPtr& states) const
{
	if (GetNumeraire()->GetType() == ARM_Numeraire::Cash)
	{
		const ARM_ModelNameMap& modelMap = *GetModelMap();
		ARM_QModel1F& forModel = static_cast<ARM_QModel1F&>( * modelMap[ForModel]->Model() );

		/// Foreign basis yield curve is always used to calibrate Zc
		const ARM_ForwardMargin& forBasisModel = static_cast<ARM_ForwardMargin&>( * modelMap[ForBasisModel]->Model() );
		double basisRatio = forBasisModel.ComputeMarginRatio(evalTime,maturityTime);
		double targetDf = forBasisModel.GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN);

		/// Compute spot Fx at asOfDate & evalTime
		const ARM_QModel1F_Fx& fxModel = static_cast<ARM_QModel1F_Fx&>( * modelMap[FxModel]->Model() );
		ARM_GP_VectorPtr payoffFactor = fxModel.Forward(fxModel.GetModelName(),evalTime,evalTime,evalTime,evalTime,states);
		if(basisRatio != 1.0)
		{
			for(size_t i=0;i<payoffFactor->size();++i)
				(*payoffFactor)[i] *= basisRatio;
		}

		double targetPayoff = fxModel.ComputeFwdAtTime(0.0) * targetDf;

		/// Set datas to foreign stochastic model
		forModel.SetTargetPayoff(maturityTime,targetPayoff,payoffFactor);
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routine : DiscountFactor
///	Returns : ARM_VectorPtr
///	Action  : computes discount factor for a given model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_2IRFXModel::DiscountFactor( 
	const string& modelName,
    double evalTime, 
	double maturityTime,
    const ARM_PricingStatesPtr& states) const
{
	const ARM_ModelNameMap& modelMap = *GetModelMap();
    ARM_ModelNameMap::const_iterator modelIdx = modelMap[modelName];

	bool isTreeMethod = (dynamic_cast<ARM_TreeBase*>(&*GetNumMethod())!=NULL);

    if( evalTime > K_NEW_DOUBLE_TOL && states != ARM_PricingStatesPtr(NULL) &&
        GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_BCKWDLOOKING &&
		isTreeMethod)
    {
        /// Initialise datas for lattice Zc calibration
        modelsAlias modelPos = static_cast<modelsAlias>(modelIdx - modelMap.begin() + DomModel);

		if(modelPos == DomModel || modelPos == DomBasisModel)
			InitDomCalibration(evalTime,maturityTime,states);

		else if(modelPos == ForModel || modelPos == ForBasisModel)
			InitForCalibration(evalTime,maturityTime,states);

		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " 2IR + FX model can't compute a Df on " + modelName);
    }

    /// Compute the Df with internal calibration if necessary
	return modelIdx->Model()->DiscountFactor(modelName, evalTime, maturityTime, states);
}

////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routines: Libor
///	Returns : a vector of libor values
///	Action  : Discount factor calibration datas then
///           default libor computation
////////////////////////////////////////////////////
ARM_VectorPtr ARM_2IRFXModel::Libor( 
		const string& modelName, 
		double evalTime,
		double fwdStartTime,
		double fwdEndTime,
		double period,
		double fwdResetTime,
		double payTime,
		const ARM_PricingStatesPtr& states) const
{
    /// Set datas for discount factor calibration
	const ARM_ModelNameMap& modelMap = *GetModelMap();
    ARM_ModelNameMap::const_iterator modelIdx = modelMap[modelName];

   if( evalTime > K_NEW_DOUBLE_TOL && states != ARM_PricingStatesPtr(NULL) &&
        GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_BCKWDLOOKING &&
		dynamic_cast<ARM_TreeBase*>(&*GetNumMethod()))
    {
        /// Initialise datas for lattice Zc calibration
        modelsAlias modelPos = static_cast<modelsAlias>(modelIdx - modelMap.begin() + DomModel);

        if(modelPos == DomModel || modelPos == DomBasisModel)
        {
            InitDomCalibration(evalTime,fwdStartTime,states);
            InitDomCalibration(evalTime,fwdEndTime,states);
        }

        else if(modelPos == ForModel || modelPos == ForBasisModel)
        {
            InitForCalibration(evalTime,fwdStartTime,states);
            InitForCalibration(evalTime,fwdEndTime,states);
        }
        else
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " 2IR + FX model can't compute a Libor on " + modelName);
    }

    /// Compute a default libor (no convexity adjustment)
    return DefaultLibor(modelName,evalTime,fwdStartTime,fwdEndTime,period,fwdResetTime,payTime,states);
}

////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routine : InitSwapCalibration
///	Returns : void
///	Action  : initialize all Zc calibrations for further lattice
///			  computations
////////////////////////////////////////////////////
void ARM_2IRFXModel::InitSwapCalibration(
	double evalTime,
	double floatStartTime,
	double floatEndTime, 
	const std::vector<double>& fixPayTimes,
	const std::vector<double>& fixPayPeriods,
	const std::vector<double>& fwdStartTimes, 
	const std::vector<double>& fwdEndTimes, 
	const std::vector<double>& fwdPayPeriods, 
	const std::vector<double>& floatPayTimes, 
	const std::vector<double>& floatPayPeriods, 
	const std::vector<double>& margin,
	bool isDbleNotional,
	const ARM_PricingStatesPtr& states,
	const InitCalibrationFunc& initFunc) const
{
	size_t i,nbFloatFlows = floatPayTimes.size();
	size_t nbFixFlows = fixPayTimes.size();

	std::vector<double> alreadyKnownTimes; /// to avoid repetitive calls even if calibration will occur once
	std::vector<double>::iterator found;

	/// Floating leg calibration
	if(isDbleNotional)
	{
		(this->*initFunc)(evalTime,floatStartTime,states);
		(this->*initFunc)(evalTime,floatEndTime,states);
		alreadyKnownTimes.push_back(floatStartTime);
		alreadyKnownTimes.push_back(floatEndTime);
	}

	for(i=0;i<nbFloatFlows;++i)
	{
		if( (margin[i]!=0.0 || !isDbleNotional) &&
			(found = find(alreadyKnownTimes.begin(),alreadyKnownTimes.end(),floatPayTimes[i]))==alreadyKnownTimes.end())
		{
			(this->*initFunc)(evalTime,floatPayTimes[i],states);
			alreadyKnownTimes.push_back(floatPayTimes[i]);
		}
		if(!isDbleNotional)
		{
			if((found = find(alreadyKnownTimes.begin(),alreadyKnownTimes.end(),fwdStartTimes[i]))==alreadyKnownTimes.end())
			{
				(this->*initFunc)(evalTime,fwdStartTimes[i],states);
				alreadyKnownTimes.push_back(fwdStartTimes[i]);
			}
			if((found = find(alreadyKnownTimes.begin(),alreadyKnownTimes.end(),fwdEndTimes[i]))==alreadyKnownTimes.end())
			{
				(this->*initFunc)(evalTime,fwdEndTimes[i],states);
				alreadyKnownTimes.push_back(fwdEndTimes[i]);
			}
		}
	}

	/// Fixed leg calibration
	for(i=0;i<nbFixFlows;++i)
	{

		if((found = find(alreadyKnownTimes.begin(),alreadyKnownTimes.end(),fixPayTimes[i]))==alreadyKnownTimes.end())
		{
			(this->*initFunc)(evalTime,fixPayTimes[i],states);
			alreadyKnownTimes.push_back(fixPayTimes[i]);
		}
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routine : NPVSwap
///	Returns : ARM_VectorPtr
///	Action  : call NPVSwap on the sub model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_2IRFXModel::NPVSwap(
	const string& modelName, 
	double evalTime,
	double floatStartTime,
	double floatEndTime, 
	const std::vector<double>& fixPayTimes,
	const std::vector<double>& fixPayPeriods,
	const std::vector<double>& fwdStartTimes, 
	const std::vector<double>& fwdEndTimes, 
	const std::vector<double>& fwdPayPeriods, 
	const std::vector<double>& floatPayTimes, 
	const std::vector<double>& floatPayPeriods, 
	const std::vector<double>& margin,
	bool isDbleNotional,
	const std::vector<double>& FixNotional,
	const std::vector<double>& FloatNotional,
	const ARM_GP_Matrix& strikesPerState,
	int payRec,
	const ARM_PricingStatesPtr& states) const
{
	const ARM_ModelNameMap& modelMap = *GetModelMap();
    ARM_ModelNameMap::const_iterator modelIdx = modelMap[modelName];

	ARM_PricingFunctionIR* IRFctor = dynamic_cast<ARM_PricingFunctionIR*>(&*modelIdx->Model());
	if( !IRFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + modelName + " does not support IR function NPVSwap!" );

	if( evalTime > K_NEW_DOUBLE_TOL && states != ARM_PricingStatesPtr(NULL) &&
        GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_BCKWDLOOKING &&
		dynamic_cast<ARM_TreeBase*>(&*GetNumMethod()))
		{
        /// Initialise datas for lattice Zc calibration
        modelsAlias modelPos = static_cast<modelsAlias>(modelIdx - modelMap.begin() + DomModel);

        if(modelPos == DomModel || modelPos == DomBasisModel)
        {
			InitSwapCalibration(evalTime,floatStartTime,floatEndTime,fixPayTimes,fixPayPeriods,
				fwdStartTimes,fwdEndTimes,fwdPayPeriods,floatPayTimes,floatPayPeriods,margin,
				isDbleNotional,states,
				&ARM_2IRFXModel::InitDomCalibration);
		}

        else if(modelPos == ForModel || modelPos == ForBasisModel)
        {
			InitSwapCalibration(evalTime,floatStartTime,floatEndTime,fixPayTimes,fixPayPeriods,
				fwdStartTimes,fwdEndTimes,fwdPayPeriods,floatPayTimes,floatPayPeriods,margin,
				isDbleNotional,states,
				&ARM_2IRFXModel::InitForCalibration);
        }
        else
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " 2IR + FX model can't compute a NPSwap on " + modelName);
    }

	return IRFctor->NPVSwap( modelName, evalTime, floatStartTime, floatEndTime, fixPayTimes, fixPayPeriods,
		fwdStartTimes, fwdEndTimes, fwdPayPeriods, floatPayTimes, floatPayPeriods, margin, isDbleNotional,
		FixNotional,FloatNotional, strikesPerState, payRec, states);
}

////////////////////////////////////////////////////
///	Class   : ARM_PricingModelIR
///	Routines: NPVBasisSwapWithNullStrikeAndNotionalExchange
///	Returns : a vector of swap rate values with null strike
///	Action  : Default Swap Rate computation with null strike and exchanging notionl flow by flow
///				WARNING: Only floating leg is computed with a notional exchanging
ARM_VectorPtr ARM_2IRFXModel::NPVBasisSwap( 
		const string& domCurveName,
		const string& forCurveName, 
		const string& fxCurveName, 
		double evalTime,
		double startTime,
		double endTime,
		int		payRec,
		const std::vector<double>& domResetTimes,	    
		const std::vector<double>& domFwdStartTimes,
		const std::vector<double>& domFwdEndTimes,
		const std::vector<double>& domFlowStartTimes,			
		const std::vector<double>& domFlowEndTimes,	
		const std::vector<double>& domFwdPayPeriods,	
		const std::vector<double>& domPayTimes,
		const std::vector<double>& domPayPeriods,
		const std::vector<double>& domMarginVector,
		const std::vector<double>& domNotionalVector,
		bool                 isDomFlottant,
		const std::vector<double>& forResetTimes,       
		const std::vector<double>& forFwdStartTimes,
		const std::vector<double>& forFwdEndTimes,
		const std::vector<double>& forFlowStartTimes,   		
		const std::vector<double>& forFlowEndTimes,	    
		const std::vector<double>& forFwdPayPeriods,	
		const std::vector<double>& forPayTimes,
		const std::vector<double>& forPayPeriods,
		const std::vector<double>& forMarginVector,
		const std::vector<double>& forNotionalVector,
		bool                 isForFlottant,
		const string&        exNotionalType,    
		const std::vector<double>& fxResetTimes,
		const std::vector<double>& fxSettlTimes,  
		const std::vector<double>& fxPayTimes,
		const ARM_GP_Matrix& domStrikePerState,
		const ARM_GP_Matrix& forStrikePerState,
		const ARM_PricingStatesPtr& states ) const
{
	if( evalTime > K_NEW_DOUBLE_TOL && states != ARM_PricingStatesPtr(NULL) &&
        GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_BCKWDLOOKING &&
		dynamic_cast<ARM_TreeBase*>(&*GetNumMethod()))
	{
		InitSwapCalibration(evalTime,startTime,endTime,domPayTimes,domPayPeriods,
			domFwdStartTimes,domFwdEndTimes,domFwdPayPeriods,domPayTimes,domPayPeriods,domMarginVector,
			false,states,&ARM_2IRFXModel::InitDomCalibration);
		
		InitSwapCalibration(evalTime,startTime,endTime,forPayTimes,forPayPeriods, forFwdStartTimes,
			forFwdEndTimes,forFwdPayPeriods,forPayTimes,forPayPeriods,forMarginVector,
			false,states,&ARM_2IRFXModel::InitForCalibration);
	}

	return ARM_MultiAssetsModel::NPVBasisSwap(domCurveName,forCurveName, fxCurveName, evalTime, startTime, endTime,payRec, domResetTimes,domFwdStartTimes,
		 domFwdEndTimes, domFlowStartTimes,	domFlowEndTimes, domFwdPayPeriods, domPayTimes,	 domPayPeriods,domMarginVector, 
		 domNotionalVector,isDomFlottant,forResetTimes,forFwdStartTimes,forFwdEndTimes, forFlowStartTimes,forFlowEndTimes,	    
		 forFwdPayPeriods,forPayTimes,forPayPeriods,forMarginVector,forNotionalVector,isForFlottant,exNotionalType,fxResetTimes,fxSettlTimes,  
		fxPayTimes, domStrikePerState, forStrikePerState, states ); 
}

////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routine : Forward
///	Returns : ARM_VectorPtr
///	Action  : Compute the forward Fx using the Fx model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_2IRFXModel::Forward(
	const string& modelName, 
    double evalTime,
	double expiryTime,
	double settlementTime,
	double payTime,
    const ARM_PricingStatesPtr& states) const
{
 	const ARM_ModelNameMap& modelMap = *GetModelMap();

    if( modelName != modelMap[FxModel]->ModelName() &&
        !(itsIsFlooredFxLocalModel      && modelName == modelMap[FlooredFxLocalModel]->ModelName()) &&
        !(itsIsCappedFxLocalModel       && modelName == modelMap[CappedFxLocalModel]->ModelName()) &&
        !(itsIsRedemptionFxLocalModel   && modelName == modelMap[RedemptionFxLocalModel]->ModelName()) )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " 2IR + FX model can't compute a forward fx on " + modelName);

    if(expiryTime < evalTime - K_NEW_DOUBLE_TOL)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " expired forward Fx is not supported");

    const ARM_QModel1F_Fx& fxModel = static_cast<ARM_QModel1F_Fx&>( * modelMap[FxModel]->Model() );

     if( evalTime > K_NEW_DOUBLE_TOL && states != ARM_PricingStatesPtr(NULL) &&
        GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_BCKWDLOOKING &&
		dynamic_cast<ARM_TreeBase*>(&*GetNumMethod()))
    {
        /// Initialise datas for lattice Zc calibration
        InitDomCalibration(evalTime,settlementTime,states);
        InitForCalibration(evalTime,settlementTime,states);

        /// No calibration is needed for domestic Df at payTime because it is passed to Forward() only
        /// for convexity adjustment
    }

    /// Compute the forward Fx
    return fxModel.Forward(fxModel.GetModelName(),evalTime,expiryTime,settlementTime,payTime,states);
}

////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routine : CallVectorial
///	Returns : ARM_VectorPtr
///	Action  : Compute the vectorial call using the Fx model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_2IRFXModel::CallVectorial(
	const string& modelName,
    double evalTime,
	double expiryTime,
	double settlementTime,
	const std::vector<double>& strikePerState,
	int callPut,
	double payTime,
    const ARM_PricingStatesPtr& states,
    ARM_PricingContext* context) const
{	
 	const ARM_ModelNameMap& modelMap = *GetModelMap();

    if( modelName != modelMap[FxModel]->ModelName() &&
        !(itsIsFlooredFxLocalModel      && modelName == modelMap[FlooredFxLocalModel]->ModelName()) &&
        !(itsIsCappedFxLocalModel       && modelName == modelMap[CappedFxLocalModel]->ModelName()) &&
        !(itsIsRedemptionFxLocalModel   && modelName == modelMap[RedemptionFxLocalModel]->ModelName()) )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " 2IR + FX model can't compute a call/put fx on " + modelName);

    if(expiryTime < evalTime - K_NEW_DOUBLE_TOL)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " expired option is not supported");

     if( evalTime > K_NEW_DOUBLE_TOL && states != ARM_PricingStatesPtr(NULL) &&
        GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_BCKWDLOOKING &&
		dynamic_cast<ARM_TreeBase*>(&*GetNumMethod()))
    {
        /// Initialise datas for lattice Zc calibration
        InitDomCalibration(evalTime,settlementTime,states);
        InitDomCalibration(evalTime,payTime,states);
        InitForCalibration(evalTime,settlementTime,states);
    }

	 return ARM_MultiAssetsModel::CallVectorial(modelName,evalTime,expiryTime,settlementTime,strikePerState,callPut,payTime,states,context);
	const ARM_MultiAssetsModel& multiAssetModel = dynamic_cast<const ARM_MultiAssetsModel&>( *this );
	return multiAssetModel.CallVectorial(modelName,evalTime,expiryTime,settlementTime,strikePerState,callPut,payTime,states,context);
}


////////////////////////////////////////////////////
///	Class  : ARM_MarketHybridModel
///	Routine: HybridCallVectorial
///	Returns: price an hybrid call that is a spread
///			 between a forward equity or FX strip
///			 and an IR swap :
///			 Call = O1(t)*Et[max(EqFxStrip(T)/O1(T) - (SwapRate(T)-Strike),0)] 
///			 Put = O1(t)*Et[max(SwapRate(T)-Strike - EqFxStrip(T)/O1(T),0)] 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_2IRFXModel::HybridCallVectorial(
	const string& modelName,
    double evalTime,
	double expiryTime,
	int callPut,
	const std::vector<double>& strikesPerState,

	/// Strip of forwards FX (or equity)
	const std::vector<double>& fxExpiryTimes,
	const std::vector<double>& fxSettlementTimes,
	const std::vector<double>& fxPayTimes,
	const std::vector<double>& fxNotionals,

	/// IR Swap
	double swapResetTime,
	const std::vector<double>& fixNotionals,
	const std::vector<double>& floatNotionals,
	double floatStartTime,
	double floatEndTime,
	const std::vector<double>& floatResetTimes,
	const std::vector<double>& floatStartTimes,
	const std::vector<double>& floatEndTimes,
	const std::vector<double>& floatIntTerms,
	const std::vector<double>& fixPayTimes,
	const std::vector<double>& fixPayPeriods,

	const ARM_PricingStatesPtr& states,
    ARM_PricingContext* context) const
{
 	const ARM_ModelNameMap& modelMap = *GetModelMap();

    if( modelName != modelMap[FxModel]->ModelName() )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " 2IR + FX model can't compute a hybrid call/put IR/FX on " + modelName);

    if(expiryTime < evalTime - K_NEW_DOUBLE_TOL)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " expired option is not supported");

	const ARM_PricingFunctionEquity& eqFxModel = dynamic_cast< ARM_PricingFunctionEquity& >( * modelMap[modelName]->Model() );

	return eqFxModel.HybridCallVectorial(modelName,evalTime,expiryTime,callPut,strikesPerState,
				fxExpiryTimes,fxSettlementTimes,fxPayTimes,fxNotionals,
				swapResetTime,fixNotionals,floatNotionals,floatStartTime,floatEndTime,
				floatResetTimes,floatStartTimes,floatEndTimes,floatIntTerms,
				fixPayTimes,fixPayPeriods,
				states,context);
}

////////////////////////////////////////////////////
///	Class  : ARM_2IRFXModel
///	Routine: RangeAccrualVectorial
///	Returns: vector ptr
///	Action : compute a range accrual dual index : IR index and FX index
////////////////////////////////////////////////////
ARM_VectorPtr ARM_2IRFXModel::RangeAccrualVectorial(
		const string& curveName,
		double evalTime,
		double startTime,
		double endTime,
		double payTime,
		const  std::vector<double>& fixingTimes,
		int    payIndexType, 
        double payIndexTerm,
		const  string& fxModelName,
		int    irIndexType, 
		const  std::vector<double>& irIndexResetTimes,
		const  std::vector<double>& irIndexStartTimes,
		const  std::vector<double>& irIndexEndTimes,
		const  std::vector<double>& irIndexTerms,
		const  std::vector<double>& fxDownBarriers,
		const  std::vector<double>& fxUpBarriers,
		const  std::vector<double>& irDownBarriers,
		const  std::vector<double>& irUpBarriers,
		const  std::vector<double>& notionals,
		const  ARM_PricingStatesPtr& states,
		ARM_Vector* eachFixingPrices,
        ARM_PricingContext* context) const
{
	const ARM_ModelNameMap* const modelMap	= GetModelMap();

	//models
	ARM_PricingModelPtr domModel     = (*modelMap)[DomModel]->Model();
	ARM_PricingModelPtr forModel     = (*modelMap)[ForModel]->Model();
	ARM_PricingModelPtr fxModel     = (*modelMap)[fxModelName]->Model();

	const ARM_ModelParamsHW1FStd* const domModelParams  = static_cast<const ARM_ModelParamsHW1FStd* const>( domModel->GetModelParams() );
	const ARM_ModelParamsHW1FStd* const forModelParams  = static_cast<const ARM_ModelParamsHW1FStd* const>( forModel->GetModelParams() );
	

	bool useLocalData=false;
	ARM_Local_SLN_Model* localfxModel = dynamic_cast<ARM_Local_SLN_Model*> (&*fxModel);
	ARM_CurveModelParam fxVol = ARM_CurveModelParam();
	ARM_SurfaceListModelParam localVol = ARM_SurfaceListModelParam();
	ARM_SurfaceListModelParam localShift = ARM_SurfaceListModelParam();
	ARM_ModelParamsHW1FStd* fxModelParams = NULL;
	if (localfxModel)
	{
		useLocalData=true;

		const ARM_ModelParamsHW1FStd* const localfxModelParams   = static_cast<const ARM_ModelParamsHW1FStd* const>( fxModel->GetModelParams() );
		const ARM_ModelParam& sfMPVol = localfxModelParams->GetModelParam(ARM_ModelParamType::Volatility);
		ARM_SurfaceListModelParam sfListMPVol = static_cast<const ARM_SurfaceListModelParam&> (sfMPVol);
		ARM_SurfacePtrVector sfListVol = sfListMPVol.GetSurfaceList();
		localVol.SetSurfaceList(sfListVol);

		const ARM_ModelParam& sfMPShift = localfxModelParams->GetModelParam(ARM_ModelParamType::Shift);
		ARM_SurfaceListModelParam sfListMPShift = static_cast<const ARM_SurfaceListModelParam&> (sfMPShift);
		ARM_SurfacePtrVector sfListShift = sfListMPShift.GetSurfaceList();
		localShift.SetSurfaceList(sfListShift);

		ARM_PricingModel* fxnumericalModel = localfxModel->GetNumericalModel();
		fxModelParams   = static_cast<ARM_ModelParamsHW1FStd* >( fxnumericalModel->GetModelParams() );
		fxVol = dynamic_cast<const ARM_CurveModelParam&>( fxModelParams->GetModelParam(fxModelParams->GetVolatilityType()) );
	}
	else 
	{
		fxModelParams   = static_cast<ARM_ModelParamsHW1FStd* >( fxModel->GetModelParams() );
		fxVol = dynamic_cast<const ARM_CurveModelParam&>( fxModelParams->GetModelParam(fxModelParams->GetVolatilityType()) );
	
	}

	size_t T0=0;
	size_t T1=1;
	size_t T2=2;
	size_t T3=3;
	

	if (	(fxDownBarriers.size()	!=	1) || 
			(fxUpBarriers.size()	!=	1) ||
			(irUpBarriers.size()	!=	1) ||
			(notionals.size()		!=	1) ||
			(irDownBarriers.size()	!=	1) )

		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " 2IR + FX : only one barrier allowed by periods");

	double irUpBarrier		= irUpBarriers[0];
	double irDownBarrier	= irDownBarriers[0];
	if(irDownBarrier >= irUpBarrier)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " irDownBarrier should be inferior than irUpBarrier");
	double fxUpBarrier		= fxUpBarriers[0];
	double fxDownBarrier	= fxDownBarriers[0];
	if(fxDownBarrier >= fxUpBarrier)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " fxDownBarrier should be inferior than fxUpBarrier");
	
	double notional			= notionals[0];
	double stddevFX, stddevFX_Divide_2;
	double stddevFXDown, stddevFXDown_Divide_2, stddevFXUp, stddevFXUp_Divide_2;
	double stddevfwdZC, stddevfwdZC_Divide_2;
	double stddevIRDown, stddevIRUp;
	
	double proba1, proba2, proba3, proba4;

	size_t nbStates = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
	ARM_VectorPtr rangeAccrual (new std::vector<double>(nbStates));

	//if (fabs(notional) < K_NEW_DOUBLE_TOL ) 
	//{
		double varLimit = K_NEW_DOUBLE_TOL*K_NEW_DOUBLE_TOL;

		ARM_VectorPtr dfPaytest	= GetRefModel()->DiscountFactor( "", evalTime, payTime, states );
		ARM_VectorPtr dfPay	= DiscountFactor( curveName, evalTime, payTime, states );
		
		size_t nbFixings = fixingTimes.size();
				
		std::vector<double> ExerciseProbability(nbStates,0.0);
		
		double K_ZC_Up_j=100000, K_ZC_Down_j=-100000, K_FX_Up_j=100000, K_FX_Down_j=-100000, rhoj=0.0;
		double K_IR_Up_j=100000, K_IR_Down_j=-100000;
		
		bool computeFXDownK=true, computeFXUpK=true;
		if (fxDownBarrier <= 0.0) computeFXDownK=false;
		if (fxUpBarrier <= 0.0) computeFXUpK=false;

		bool interpol=false;

		for (size_t i=0; i<nbFixings; i++)
		{
			double fixing = fixingTimes[i];
			ARM_GP_Matrix correlMatrix = GetCorrelMatrix()->Interpolate(fixing);

			/// Forex variances
			
				double varspotFx       = fxModelParams->StateLocalVariance(evalTime,fixing,fixing);
				double varZcDom    = ARM_ModelParamsHW1F::HW1FZcCovariance( domModelParams, domModelParams, evalTime,fixing,fixing );
				double varZcFor    = ARM_ModelParamsHW1F::HW1FZcCovariance( forModelParams, forModelParams, evalTime,fixing,fixing );

				double covarZcDomZcFor = ARM_ModelParamsHW1F::HW1FZcCovariance( domModelParams, forModelParams, evalTime,fixing,fixing );
				double covarZcDomFx    = ARM_ModelParamsHW1F::HW1FEqFxZcCovariance( fxVol, domModelParams, evalTime,fixing,fixing );
				double covarZcForFx    = ARM_ModelParamsHW1F::HW1FEqFxZcCovariance( fxVol, forModelParams, evalTime,fixing,fixing );

				double varFX =
					varspotFx + varZcDom + varZcFor
					- 2*correlMatrix(DomModel,ForModel)*covarZcDomZcFor
					- 2*correlMatrix(DomModel,FxModel)*covarZcDomFx
					 + 2*correlMatrix(ForModel,FxModel)*covarZcForFx;
				
				if(varFX<varLimit)
				varFX = varLimit; // numerical noise
				
				stddevFX = sqrt(varFX);
				stddevFX_Divide_2 = 0.5 * stddevFX;

			if (!useLocalData)
			{	
				stddevFXDown = stddevFX;
				stddevFXUp = stddevFX;

				stddevFXDown_Divide_2 = stddevFX_Divide_2;
				stddevFXUp_Divide_2 = stddevFX_Divide_2;

			}
			else 
			{
				double volFXUp = localVol.GetValue(T0, evalTime, fixing);
				double volFXDown = localVol.GetValue(T1, evalTime, fixing);

				if (fabs(volFXUp) < K_NEW_DOUBLE_TOL)	//variance squeeeze
					stddevFXUp = stddevFX;
				else 
					stddevFXUp = volFXUp * sqrt((fixing-evalTime)/K_YEAR_LEN);
				

				if (fabs(volFXDown) < K_NEW_DOUBLE_TOL)	//variance squeeeze
					stddevFXDown = stddevFX;
				else 
					stddevFXDown = volFXDown * sqrt((fixing-evalTime)/K_YEAR_LEN);
				
				stddevFXDown_Divide_2 = stddevFXDown * 0.5;
				stddevFXUp_Divide_2 = stddevFXUp * 0.5;

			}
						
			//IR variance
			double varZCDomEnd		= ARM_ModelParamsHW1FStd::HW1FZcCovariance(domModelParams,domModelParams,evalTime,fixing,irIndexEndTimes[i]);
			double varZcDomStart    = ARM_ModelParamsHW1FStd::HW1FZcCovariance(domModelParams,domModelParams,evalTime,fixing,fixing );
			double covarZCDomEnd = ARM_ModelParamsHW1FStd::HW1FZcCovariance(domModelParams,domModelParams,evalTime,fixing,fixing,irIndexEndTimes[i]);
			double varfwdZC = varZCDomEnd + varZcDomStart - 2*covarZCDomEnd;

			stddevfwdZC=sqrt(varfwdZC);
			stddevfwdZC_Divide_2 = 0.5 * stddevfwdZC;
			
			if (useLocalData) 
			{
				//Si Local model utilisé, on travaille directement avec la vol du libor qui est normale
				//on ne passe donc pas par la vol lognormale du ZC forward.
				double volIRUp = localVol.GetValue(T2, evalTime, fixing);
				double volIRDown = localVol.GetValue(T3, evalTime, fixing);

				if ( (fabs(volIRUp) < K_NEW_DOUBLE_TOL) || (fabs(volIRDown) < K_NEW_DOUBLE_TOL) )
				{
					ARM_QModel1F_Fx* eqFxModel = dynamic_cast< ARM_QModel1F_Fx* >( localfxModel->GetNumericalModel() );
					ARM_VectorPtr fwdZCMP_Eval_Start = eqFxModel->ComputeFwdZCModelParam(evalTime, fixing, fixing, irIndexEndTimes[i]);
					double volirModel = (*fwdZCMP_Eval_Start)[0];
					ARM_VectorPtr dfStart	= GetRefModel()->DiscountFactor( "", 0, fixing, states );
					
					//ARM_VectorPtr check0 = DiscountFactor( curveName, 0, fixing, states );
					//double check00 = (*check0)[0];
	
					double dfStart0 = (*dfStart)[0];
					ARM_VectorPtr dfEnd	= GetRefModel()->DiscountFactor( "", 0, irIndexEndTimes[i], states );
					
					//ARM_VectorPtr check1 = DiscountFactor( curveName, 0, irIndexEndTimes[i], states );
					//double check11 = (*check1)[0];	
					
					double dfEnd0 = (*dfEnd)[0];
					volirModel =  (1./irIndexTerms[i]) * volirModel * pow((dfStart0 / dfEnd0), 2);
					
					if (fabs(volIRUp) < K_NEW_DOUBLE_TOL)	//variance squeeeze
						stddevIRUp = volirModel*sqrt((fixing-evalTime)/K_YEAR_LEN);
					else 
						stddevIRUp = volIRUp*sqrt((fixing-evalTime)/K_YEAR_LEN);

					if (fabs(volIRDown) < K_NEW_DOUBLE_TOL)	//variance squeeeze
						stddevIRDown = volirModel*sqrt((fixing-evalTime)/K_YEAR_LEN);
					else
						stddevIRDown = volIRDown*sqrt((fixing-evalTime)/K_YEAR_LEN);
				}
				else 
				{
					stddevIRUp = volIRUp * sqrt((fixing-evalTime)/K_YEAR_LEN);
					stddevIRDown = volIRDown*sqrt((fixing-evalTime)/K_YEAR_LEN);
				}
			
			}
			

			//Correlation
			double covarZcDomFX1 = ARM_ModelParamsHW1FStd::HW1FEqFxZcCovariance(fxVol,domModelParams,evalTime,fixing,irIndexEndTimes[i])*correlMatrix(DomModel,FxModel);
			double covarZcDomFX2 = ARM_ModelParamsHW1FStd::HW1FEqFxZcCovariance(fxVol,domModelParams,evalTime,fixing,fixing)*correlMatrix(DomModel,FxModel);
			double covar = covarZcDomFX1 - covarZcDomFX2;
			rhoj =  covar / (stddevfwdZC * stddevFX);
									
			//Calcul de la barrière du fwdZC
			bool computeZCDownK=true, computeZCUpK=true;
			double zcDownBarrier=-10000, zcUpBarrier=10000;
			if (!useLocalData)
			{
				double temp1 = 1+irIndexTerms[i] * irUpBarrier;
				double temp2 = 1+irIndexTerms[i] * irDownBarrier;	//temp2 <= temp1 obligatoirement
				
				if ((temp1>0.0) && (temp2>0.0))
				{
					zcDownBarrier = 1. / temp1;
					zcUpBarrier = 1. / temp2;
				}
				else if (temp2<0.0) 
				{
					if (temp1<0.0)
					{
						zcDownBarrier = 1. / temp2;
						zcUpBarrier = 1. / temp1;
					}
					else if (temp1>0.0)
					{
						zcDownBarrier = 1. / temp1;
					}
				
				}

				if (zcDownBarrier <= 0.0) computeZCDownK=false;
				if (zcUpBarrier <= 0.0) computeZCUpK=false;
			}
			else
			{
				if (computeFXUpK)
					fxUpBarrier = localShift.GetValue(T0, evalTime, fixing);
				if (computeFXDownK)
					fxDownBarrier = localShift.GetValue(T1, evalTime, fixing);
				
				irUpBarrier = localShift.GetValue(T2,evalTime, fixing);
				irDownBarrier = localShift.GetValue(T3,evalTime, fixing);

			}		

			//Attention: payTime non utilisé (pas d'ajustement du forward FX) sous la proba QpayTime.
			ARM_VectorPtr forwardFX = Forward(fxModelName, evalTime, fixing, fixing, payTime, states);

			//ARM_VectorPtr dfFix	= GetRefModel()->DiscountFactor( "", evalTime, fixing, states );
			//ARM_VectorPtr dfEnd	= GetRefModel()->DiscountFactor( "", evalTime, irIndexEndTimes[i], states );

			ARM_VectorPtr dfFix	= DiscountFactor( curveName, evalTime, fixing, states );
			ARM_VectorPtr dfEnd	= DiscountFactor( curveName, evalTime, irIndexEndTimes[i], states );

			for (size_t j=0; j<nbStates; j++)
			{
				
				//strikes associés à binaire sur FX
				double fwdFX = (*forwardFX)[j];
				if (computeFXUpK)	
					K_FX_Up_j	= log(fxUpBarrier/fwdFX)	/stddevFXUp + stddevFXUp_Divide_2;
				if (computeFXDownK) 
					K_FX_Down_j = log(fxDownBarrier/fwdFX)	/stddevFXDown + stddevFXDown_Divide_2;

				if (!useLocalData)
				{
					//strikes associés à binaire sur IR
					double fwdZC = (*dfEnd)[j] / (*dfFix)[j] ;

					if (computeZCUpK)
						K_ZC_Up_j	= log(zcUpBarrier/fwdZC)	/stddevfwdZC + stddevfwdZC_Divide_2;	
					if (computeZCDownK)
						K_ZC_Down_j = log(zcDownBarrier/fwdZC)	/stddevfwdZC + stddevfwdZC_Divide_2;

					proba1=NormalCDF(K_FX_Up_j,K_ZC_Up_j,rhoj);
					proba2=NormalCDF(K_FX_Up_j,K_ZC_Down_j,rhoj);
					proba3=NormalCDF(K_FX_Down_j,K_ZC_Up_j,rhoj);
					proba4=NormalCDF(K_FX_Down_j,K_ZC_Down_j,rhoj);
				}
				else 
				{
					double fwdIR = ((*dfFix)[j] / (*dfEnd)[j] - 1.0) / irIndexTerms[i];
					K_IR_Up_j	= (irUpBarrier - fwdIR ) / stddevIRUp;
					K_IR_Down_j = (irDownBarrier - fwdIR ) / stddevIRDown;

					proba1=NormalCDF(K_FX_Up_j,K_IR_Up_j,rhoj);
					proba2=NormalCDF(K_FX_Up_j,K_IR_Down_j,rhoj);
					proba3=NormalCDF(K_FX_Down_j,K_IR_Up_j,rhoj);
					proba4=NormalCDF(K_FX_Down_j,K_IR_Down_j,rhoj);
				}					
				
				//double test = fwdFX * NormalCDF(log(fwdFX/fxUpBarrier)	/stddevFX + stddevFX_Divide_2) 
				//						- fxUpBarrier * NormalCDF(log(fwdFX/fxUpBarrier)	/stddevFX - stddevFX_Divide_2);

				ExerciseProbability[j] += proba1-proba2-proba3+proba4;

				

			}

			/*	//CHECK FX Digitale
			ARM_VectorPtr testforwardFX = Forward(fxModelName, 0, fixing, fixing, payTime, states);
			ARM_GP_Matrix VV(3,3);
			ComputeIntegratedFwdFxVCV(0,fixing, domModelParams, forModelParams, fxModelParams, fxVol, correlMatrix, VV);
			double std=sqrt(VV(2,2));
			double testK=log(fxUpBarrier/(*testforwardFX)[0])	/std + 0.5*std;
			double prix = NormalCDF(testK,10000,0)*notional;
			*/

				//CHECK IR Digitale
			/*ARM_VectorPtr testdfFix	= DiscountFactor( curveName, 0, fixing, states );
			ARM_VectorPtr testdfEnd	= DiscountFactor( curveName, 0, irIndexEndTimes[i], states );
			double testvarZCDomFull = domModelParams->ZcVarianceSpread(fixing,fixing,irIndexEndTimes[i],fixing);
			ARM_GP_Matrix VV(3,3);
			ComputeIntegratedFwdFxVCV(0,fixing, domModelParams, forModelParams, fxModelParams, fxVol, correlMatrix, VV);
			double teststddevXDom = sqrt(VV(0,0));
			double testK	= (log(zcUpBarrier/((*testdfEnd)[0] / (*testdfFix)[0]))	+ 0.5*testvarZCDomFull + stateFact*teststddevXDom) / (stateFact*teststddevXDom);
			double prix = NormalCDF(testK,10000,0)*notional;*/
			
  
		
		}
		
		double notio_divide_nbFix = notional / nbFixings;
		for (size_t j=0; j<nbStates; j++)
		{
			(*rangeAccrual)[j] = notio_divide_nbFix * ExerciseProbability[j] * (*dfPay)[j];
		}
//	}	//end zero notional
//	else
//	{
//		for (size_t j=0; j<nbStates; j++)
//		{
//			(*rangeAccrual)[j] = 0.0;
//		}
//	}



	return rangeAccrual;
}

////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routines: UnderlyingCovariance
///	Returns : double
///	Action  : compute the local covariance between
///           fromTime -> toTime of forward FXs
///           with settlements at startTimes
////////////////////////////////////////////////////
double ARM_2IRFXModel::UnderlyingCovariance(
    string underlyingType,
    double	fromTime,
    double	toTime,
    double	startTime1,
    double  endTime1,
    double	startTime2,
    double  endTime2,
    double	startTime3,
    double  endTime3,
    double	startTime4,
    double  endTime4) const
{
    /// FIX FIX : underlyingType is not used, always FWD_FX
    stringToUpper(underlyingType);
    if(underlyingType != "FX")
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " 2IR + FX only computes local covariance between forward Fx");

    return (*GetModelMap())[FxModel]->Model()->UnderlyingCovariance(underlyingType,fromTime,toTime,startTime1,endTime1,startTime2,endTime2,startTime3,endTime3,startTime4,endTime4);
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: UpdatePDE3DCoeffs
///	Returns : void
///	Action  : Update the PDE 3D coefficients
////////////////////////////////////////////////////
void ARM_2IRFXModel::UpdatePDE3DCoeffs(
		size_t timeIdx,
		const ARM_PricingStatesPtr& states,
		ARM_GP_VectorPtr& qxx,
		ARM_GP_VectorPtr& qyy,
		ARM_GP_VectorPtr& qzz,
		ARM_GP_VectorPtr& qxy,
		ARM_GP_VectorPtr& qyz,
		ARM_GP_VectorPtr& qzx,
		ARM_GP_VectorPtr& px,
		ARM_GP_VectorPtr& py,
		ARM_GP_VectorPtr& pz,
		ARM_GP_VectorPtr& o,
		double lambda,
		bool IsInit
		) const
{
	const ARM_GP_MatrixPtr& numMethodStates = states->GetNumMethodStates();
	size_t N = numMethodStates->cols();
    const ARM_ModelNameMap* const modelMap	= GetModelMap();

	//models
	ARM_PricingModelPtr domModel     = (*modelMap)[DomModel]->Model();
	ARM_PricingModelPtr forModel     = (*modelMap)[ForModel]->Model();
	ARM_QModel1F_Fx* fxModel		 = dynamic_cast<ARM_QModel1F_Fx*>(&*(*modelMap)[FxModel]->Model());

	ARM_ModelParam& modelParamDomMeanReversion	= domModel->GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion);
	ARM_ModelParam& modelParamForMeanReversion	= forModel->GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion);
	ARM_ModelParam& modelParamDomQVol			= domModel->GetModelParams()->GetModelParam(ARM_ModelParamType::QVol);
	ARM_ModelParam& modelParamForQVol			= forModel->GetModelParams()->GetModelParam(ARM_ModelParamType::QVol);
	ARM_ModelParam& modelParamFxQVol			= fxModel->GetModelParams()->GetModelParam(ARM_ModelParamType::QVol);
	ARM_ModelParam& modelParamFxQ				= fxModel->GetModelParams()->GetModelParam(ARM_ModelParamType::QParameter);
	ARM_ZeroCurvePtr& DomZcCurve = domModel->GetZeroCurve();
	ARM_ZeroCurvePtr& ForZcCurve = forModel->GetZeroCurve();


	//dates needed
	double ti;
	double tk;
	double dt;
	if(IsInit)
	{
		ti = GetNumMethod()->GetTimeStep(timeIdx);
		tk = GetNumMethod()->GetTimeStep(timeIdx-1);
		dt = (tk-ti)/K_YEAR_LEN;	
	}
	else
	{
		ti = GetNumMethod()->GetTimeStep(timeIdx);
		tk = GetNumMethod()->GetTimeStep(timeIdx+1);
		dt = (tk-ti)/K_YEAR_LEN;
	}

	//Calculation of forward instaneous rate
    double Bdom0ti = DomZcCurve->DiscountPrice(ti/K_YEAR_LEN);
    double Bfor0ti = ForZcCurve->DiscountPrice(ti/K_YEAR_LEN);

    double Bdom0tk = DomZcCurve->DiscountPrice(tk/K_YEAR_LEN);
    double Bfor0tk = ForZcCurve->DiscountPrice(tk/K_YEAR_LEN);

    double domFwdRate = log(Bdom0ti/Bdom0tk)/dt;
    double forFwdRate = log(Bfor0ti/Bfor0tk)/dt;
	
	//Calculation of the constants present in the coefficient of the PDE
	double volx=modelParamDomQVol.GetValue(ti);
	double voly=modelParamForQVol.GetValue(ti);
	double volz=modelParamFxQVol.GetValue(ti);

	double mrx=modelParamDomMeanReversion.GetValue(ti);
	double mry=modelParamForMeanReversion.GetValue(ti);
	double q=modelParamFxQ.GetValue(ti);

	ARM_GP_Matrix correlMatrix = GetCorrelMatrix()->Interpolate(ti);

	double rhoxy=correlMatrix(0,1);
	double rhoyz=correlMatrix(1,2);
	double rhozx=correlMatrix(2,0);
	
	double gamma=0.0;

	std::vector<double>::iterator qxxit = qxx->begin();
	std::vector<double>::iterator qyyit = qyy->begin();
	std::vector<double>::iterator qzzit = qzz->begin();
	std::vector<double>::iterator qxyit = qxy->begin();
	std::vector<double>::iterator qyzit = qyz->begin();
	std::vector<double>::iterator qzxit = qzx->begin();
	std::vector<double>::iterator pxit  = px->begin();
	std::vector<double>::iterator pyit  = py->begin();

	if(lambda==-3)//Centred SpotFX
	{
		
	}
	else if(lambda==-2.5)//Tollaly Centred for the rates LogSpotFX
	{
		double mS=0.0;
		double gammaS=0.0;
		double S0ti = fxModel->ComputeFwdAtTime(ti);
		const ARM_GP_MatrixPtr& modelStates = states->GetModelStates();
		 //loop
		for (int K = 0; K < N; ++K)
		{
			*qxxit = -pow (volx,2);
			*qyyit = -pow(voly,2);
			mS	   = S0ti*(1/q-1);
			gamma = q*volz*( 1 + mS/(*modelStates)(2,K) );
			*qzzit = -pow(gamma,2);
			*qxyit = -rhoxy*volx*voly;
			*qyzit = -rhoyz*voly*gamma;
			*qzxit = -rhozx*gamma*volx;
			*pxit = mrx*(*numMethodStates)(0,K);
			*pyit = mry*(*numMethodStates)(1,K)+rhoyz*voly*gamma; 

			qxxit++;
			qyyit++;
			qzzit++;
			qxyit++;
			qyzit++;
			qzxit++;
			pxit++;
			pyit++;

		}
		//Coefficients function of the numeraire

		ARM_GP_VectorPtr rd = domModel->RiskNeutralDrift(timeIdx,modelStates,false); 
		ARM_GP_VectorPtr rf = forModel->RiskNeutralDrift(timeIdx,modelStates,false);
		std::vector<double>::iterator pzit  = pz->begin();
		std::vector<double>::iterator rdit  = rd->begin();
		std::vector<double>::iterator rfit  = rf->begin();


		if(GetNumeraire()->GetType() == ARM_Numeraire::TerminalZc)
		{
		
		}
		else if(GetNumeraire()->GetType() == ARM_Numeraire::Cash)
		{
			//markovian drift
			for (int K = 0; K < N; ++K)
			{
				mS	  = S0ti*(1/q-1);
				gamma = q*volz*( 1 + mS/(*modelStates)(2,K) );
				*pzit= -( (*rdit) - (*rfit)-0.5*gamma*gamma );
				pzit++;
				rdit++;
				rfit++;
			}
			//actualisation
			*o = *(domModel->RiskNeutralDrift(timeIdx,numMethodStates,true))  + domFwdRate;

		}
		else
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not Cash or TerminalZC numeraire!");
		}	
	}
	else if(lambda==-2)//Centred (just of the forward rate for the moment) LogSpotFX
	{
		double mS=0.0;
		double gammaS=0.0;
		double S0ti = fxModel->ComputeFwdAtTime(ti);
		const ARM_GP_MatrixPtr& modelStates = states->GetModelStates();

		 //loop
		for (int K = 0; K < N; ++K)
		{
			*qxxit = -pow (volx,2);
			*qyyit = -pow(voly,2);
			mS	   = S0ti*(1/q-1);
			gamma = q*volz*( 1 + mS/(*modelStates)(2,K) );
			*qzzit = -pow(gamma,2);
			*qxyit = -rhoxy*volx*voly;
			*qyzit = -rhoyz*voly*gamma;
			*qzxit = -rhozx*gamma*volx;
			*pxit = mrx*(*numMethodStates)(0,K);
			*pyit = mry*(*numMethodStates)(1,K)+rhoyz*voly*gamma; 

			qxxit++;
			qyyit++;
			qzzit++;
			qxyit++;
			qyzit++;
			qzxit++;
			pxit++;
			pyit++;

		}
		//Coefficients function of the numeraire

		ARM_GP_VectorPtr rd = domModel->RiskNeutralDrift(timeIdx,modelStates,true); 
		ARM_GP_VectorPtr rf = forModel->RiskNeutralDrift(timeIdx,modelStates,true);
		std::vector<double>::iterator pzit  = pz->begin();
		std::vector<double>::iterator rdit  = rd->begin();
		std::vector<double>::iterator rfit  = rf->begin();


		if(GetNumeraire()->GetType() == ARM_Numeraire::TerminalZc)
		{
		
		}
		else if(GetNumeraire()->GetType() == ARM_Numeraire::Cash)
		{
			//markovian drift
			for (int K = 0; K < N; ++K)
			{
				mS	  = S0ti*(1/q-1);
				gamma = q*volz*( 1 + mS/(*modelStates)(2,K) );
				*pzit= -( (*rdit) - (*rfit)-0.5*gamma*gamma );
				pzit++;
				rdit++;
				rfit++;
			}
			//actualisation
			*o = *(domModel->RiskNeutralDrift(timeIdx,numMethodStates,true))  + domFwdRate;

		}
		else
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not Cash or TerminalZC numeraire!");
		}	
	}
	else if(lambda==-1)//Centred Qmapping
	{
		 //loop
		for (int K = 0; K < N; ++K)
		{
			*qxxit = -pow (volx,2);
			*qyyit = -pow(voly,2);
			*qzzit = -pow(volz,2);
			*qxyit = -rhoxy*volx*voly;
			*qyzit = -rhoyz*voly*volz;
			*qzxit = -rhozx*volz*volx;
			*pxit = mrx*(*numMethodStates)(0,K);
			*pyit = mry*(*numMethodStates)(1,K); 

			qxxit++;
			qyyit++;
			qzzit++;
			qxyit++;
			qyzit++;
			qzxit++;
			pxit++;
			pyit++;

		}
		//Coefficients function of the numeraire

		const ARM_GP_MatrixPtr& modelStates = states->GetModelStates();
		fxModel->MarkovianDriftPDE(timeIdx,modelStates,pz);

		if(GetNumeraire()->GetType() == ARM_Numeraire::TerminalZc)
		{
			*pz= 0.0 - *(pz);
			*o = std::vector<double>(N,0.0);
		}
		else if(GetNumeraire()->GetType() == ARM_Numeraire::Cash)
		{
			//markovian drift and actualisation
			*pz= 0.0 - *(pz);
			*o = *(domModel->RiskNeutralDrift(timeIdx,modelStates,true))  + domFwdRate;
		}
		else
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not Cash or TerminalZC numeraire!");
		}
	}
	else if(lambda==1)//Non centred Qmapping
	{
		 //loop
		for (int K = 0; K < N; ++K)
		{
			*qxxit = -pow (volx,2);
			*qyyit = -pow(voly,2);
			*qzzit = -pow(volz,2);
			*qxyit = -rhoxy*volx*voly;
			*qyzit = -rhoyz*voly*volz;
			*qzxit = -rhozx*volz*volx;
			*pxit = mrx*(*numMethodStates)(0,K);
			gamma=q*volz*exp(q*(*numMethodStates)(2,K))/(exp(q*(*numMethodStates)(0,K))-(1-q));
			*pyit = mry*(*numMethodStates)(1,K)+rhoyz*voly*gamma; 

			qxxit++;
			qyyit++;
			qzzit++;
			qxyit++;
			qyzit++;
			qzxit++;
			pxit++;
			pyit++;

		}
		//Coefficients function of the numeraire

		fxModel->MarkovianDriftPDE(timeIdx,numMethodStates,pz);

		if(GetNumeraire()->GetType() == ARM_Numeraire::TerminalZc)
		{
			double T=GetNumeraire()->GetMaturity();
			const ARM_ModelParamsHW1FStd* const modelParamsDom = static_cast<const ARM_ModelParamsHW1FStd* const>( domModel->GetModelParams() );
			double fxTheta = rhozx*volz*volx*modelParamsDom->BetatT(ti,T);
			*pz= fxTheta + 0.5*q*volz*volz  - *(pz);
			*o = std::vector<double>(N,0.0);
		}
		else if(GetNumeraire()->GetType() == ARM_Numeraire::Cash)
		{
			//markovian drift and actualisation
			*pz=+0.5*q*volz*volz - *(pz);
			*o = *(domModel->RiskNeutralDrift(timeIdx,numMethodStates,true))  + domFwdRate;
		}
		else
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not Cash or TerminalZC numeraire!");
		}
	}
	else if(lambda==2)//NonCentredLogSpotFX
	{
		double mS=0.0;
		double gammaS=0.0;
		double S0ti = fxModel->ComputeFwdAtTime(ti);
		double S0	= fxModel->ComputeFwdAtTime(0.0);
		const ARM_GP_MatrixPtr& modelStates = states->GetModelStates();

		 //loop
		for (int K = 0; K < N; ++K)
		{
			*qxxit = -pow (volx,2);
			*qyyit = -pow(voly,2);
			mS	   = S0ti*(1/q-1);
			gamma = q*volz*( 1 + mS/(*modelStates)(2,K) );
			*qzzit = -pow(gamma,2);
			*qxyit = -rhoxy*volx*voly;
			*qyzit = -rhoyz*voly*gamma;
			*qzxit = -rhozx*gamma*volx;
			*pxit = mrx*(*numMethodStates)(0,K);
			*pyit = mry*(*numMethodStates)(1,K)+rhoyz*voly*gamma; 

			qxxit++;
			qyyit++;
			qzzit++;
			qxyit++;
			qyzit++;
			qzxit++;
			pxit++;
			pyit++;

		}
		//Coefficients function of the numeraire

		ARM_GP_VectorPtr rd = domModel->RiskNeutralDrift(timeIdx,modelStates,true); 
		ARM_GP_VectorPtr rf = forModel->RiskNeutralDrift(timeIdx,modelStates,true);
		std::vector<double>::iterator pzit  = pz->begin();
		std::vector<double>::iterator rdit  = rd->begin();
		std::vector<double>::iterator rfit  = rf->begin();


		if(GetNumeraire()->GetType() == ARM_Numeraire::TerminalZc)
		{
		
		}
		else if(GetNumeraire()->GetType() == ARM_Numeraire::Cash)
		{
			//markovian drift
			for (int K = 0; K < N; ++K)
			{
				mS	   = S0ti*(1/q-1);
				gamma = q*volz*( 1 + mS/(*modelStates)(2,K) );
				*pzit= -(domFwdRate + (*rdit) - forFwdRate - (*rfit)-0.5*gamma*gamma);
				pzit++;
				rdit++;
				rfit++;
			}
			//actualisation
			*o = *(domModel->RiskNeutralDrift(timeIdx,numMethodStates,true))  + domFwdRate;

		}
		else
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not Cash or TerminalZC numeraire!");
		}	
	}
	else if(lambda==2.5)//NonCentredLogSpotFX with freezing of gamma in the quanto drift
	{
		double mS=0.0;
		double gammaS=0.0;
		double S0ti = fxModel->ComputeFwdAtTime(ti);
		double S0	= fxModel->ComputeFwdAtTime(0.0);
		const ARM_GP_MatrixPtr& modelStates = states->GetModelStates();

		 //loop
		for (int K = 0; K < N; ++K)
		{
			*qxxit = -pow (volx,2);
			*qyyit = -pow(voly,2);
			mS	   = S0ti*(1/q-1);
			gamma = q*volz*( 1 + mS/(*modelStates)(2,K) );
			*qzzit = -pow(gamma,2);
			*qxyit = -rhoxy*volx*voly;
			*qyzit = -rhoyz*voly*volz;
			*qzxit = -rhozx*volz*volx;
			*pxit = mrx*(*numMethodStates)(0,K);
			*pyit = mry*(*numMethodStates)(1,K)+rhoyz*voly*volz; 

			qxxit++;
			qyyit++;
			qzzit++;
			qxyit++;
			qyzit++;
			qzxit++;
			pxit++;
			pyit++;

		}
		//Coefficients function of the numeraire

		ARM_GP_VectorPtr rd = domModel->RiskNeutralDrift(timeIdx,modelStates,true); 
		ARM_GP_VectorPtr rf = forModel->RiskNeutralDrift(timeIdx,modelStates,true);
		std::vector<double>::iterator pzit  = pz->begin();
		std::vector<double>::iterator rdit  = rd->begin();
		std::vector<double>::iterator rfit  = rf->begin();


		if(GetNumeraire()->GetType() == ARM_Numeraire::TerminalZc)
		{
		
		}
		else if(GetNumeraire()->GetType() == ARM_Numeraire::Cash)
		{
			//markovian drift
			for (int K = 0; K < N; ++K)
			{
				mS	   = S0ti*(1/q-1);
				gamma = q*volz*( 1 + mS/(*modelStates)(2,K) );
				*pzit= -(domFwdRate + (*rdit) - forFwdRate - (*rfit)-0.5*gamma*gamma);
				pzit++;
				rdit++;
				rfit++;
			}
			//actualisation
			*o = *(domModel->RiskNeutralDrift(timeIdx,numMethodStates,true))  + domFwdRate;

		}
		else
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not Cash or TerminalZC numeraire!");
		}	
	}
	else if(lambda==3)//Non Centred SpotFX
	{
		double mS=0.0;
		double gammaS=0.0;
		double S0ti = fxModel->ComputeFwdAtTime(ti);
		double S0 = fxModel->ComputeFwdAtTime(0.0);
		const ARM_GP_MatrixPtr& modelStates = states->GetModelStates();

		 //loop
		for (int K = 0; K < N; ++K)
		{
			*qxxit = -pow (volx,2);
			*qyyit = -pow(voly,2);
			mS	   = S0ti*(1/q-1);
			gammaS = q*volz*((*modelStates)(2,K)+mS);
			*qzzit = -pow(gammaS,2);
			*qxyit = -rhoxy*volx*voly;
			*qyzit = -rhoyz*voly*gammaS;
			*qzxit = -rhozx*gammaS*volx;
			*pxit = mrx*(*numMethodStates)(0,K);
			*pyit = mry*(*numMethodStates)(1,K)+rhoyz*voly*gammaS; 

			qxxit++;
			qyyit++;
			qzzit++;
			qxyit++;
			qyzit++;
			qzxit++;
			pxit++;
			pyit++;

		}
		//Coefficients function of the numeraire
		ARM_GP_VectorPtr rd = domModel->RiskNeutralDrift(timeIdx,modelStates,true); 
		ARM_GP_VectorPtr rf = forModel->RiskNeutralDrift(timeIdx,modelStates,true);
		std::vector<double>::iterator pzit  = pz->begin();
		std::vector<double>::iterator rdit  = rd->begin();
		std::vector<double>::iterator rfit  = rf->begin();


		if(GetNumeraire()->GetType() == ARM_Numeraire::TerminalZc)
		{
		
		}
		else if(GetNumeraire()->GetType() == ARM_Numeraire::Cash)
		{
			//markovian drift
			for (int K = 0; K < N; ++K)
			{
				*pzit= -(domFwdRate + (*rdit) - forFwdRate - (*rfit))*(*modelStates)(2,K);
				pzit++;
				rdit++;
				rfit++;
			}
			//actualisation
			*o = *(domModel->RiskNeutralDrift(timeIdx,numMethodStates,true))  + domFwdRate;

		}
		else
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not Cash or TerminalZC numeraire!");
		}
	}
	else
	{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": value of lambda does not correspond to a PDE");
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routines: toString
///	Returns :
///	Action  : stringify the object
////////////////////////////////////////////////////
string ARM_2IRFXModel::toString(const string& indent, const string& nextIndent) const
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

    if(itsIsFlooredFxLocalModel)
    {
        os << indent << "\n\n------>     Floored Fx Local Model      <------\n";
        os << modelMap[FlooredFxLocalModel]->Model()->toString(indent,nextIndent);
    }

    if(itsIsCappedFxLocalModel)
    {
        os << indent << "\n\n------>      Capped Fx Local Model      <------\n";
        os << modelMap[CappedFxLocalModel]->Model()->toString(indent,nextIndent);
    }

    if(itsIsRedemptionFxLocalModel)
    {
        os << indent << "\n\n------>    Redemption Fx Local Model    <------\n";
        os << modelMap[RedemptionFxLocalModel]->Model()->toString(indent,nextIndent);
    }

    return os.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routines: ComputeTimeLag
///	Returns :
///	Action  : compute the time lag for a 
/// forward FX fixed at expiry time paid and paid
/// at paytime
////////////////////////////////////////////////////

double ARM_2IRFXModel::ComputeTimeLag(double expiryTime, double payTime, bool isDomestic) const
{
	return exp(ComputeAbsTimeLag(expiryTime,payTime, isDomestic));
}

////////////////////////////////////////////////////
///	Class   : ARM_2IRFXModel
///	Routines: ComputeLogTimeLag
///	Returns :
///	Action  : compute the time lag for a 
/// forward FX fixed at expiry time paid and paid
/// at paytime
////////////////////////////////////////////////////

double ARM_2IRFXModel::ComputeAbsTimeLag(double expiryTime, double payTime, bool isDomestic) const
{
	const ARM_ModelNameMap* const modelMap = GetModelMap();
	ARM_PricingModelPtr domModel    = (*modelMap)[DomModel]->Model();
    ARM_PricingModelPtr forModel    = (*modelMap)[ForModel]->Model();
    ARM_PricingModelPtr q1fFxModel     = (*modelMap)[FxModel]->Model();
	const ARM_ModelParamsHW1FStd* const domModelParams  = static_cast<const ARM_ModelParamsHW1FStd* const>( domModel->GetModelParams() );
	const ARM_ModelParamsHW1FStd* const forModelParams  = static_cast<const ARM_ModelParamsHW1FStd* const>( forModel->GetModelParams() );
	const ARM_ModelParamsHW1FStd* const fxModelParams   = static_cast<const ARM_ModelParamsHW1FStd* const>( q1fFxModel->GetModelParams() );
	const ARM_ModelParamsQ1F* const fxParams            = dynamic_cast<const ARM_ModelParamsQ1F* const>( fxModelParams );
	const ARM_CurveModelParam& fxVol                    = dynamic_cast<const ARM_CurveModelParam&>( fxModelParams->GetModelParam(fxModelParams->GetVolatilityType()) );

	ARM_GP_Matrix correlMatrix = GetCorrelMatrix()->Interpolate(expiryTime);
	double corrDomFor    = correlMatrix(DomModel,ForModel);
	double corrDomFx    = correlMatrix(DomModel,FxModel);
	double corrForFx    = correlMatrix(ForModel,FxModel);


	double drift;

	if (isDomestic)
	{
		double betatT;
		double covarFXDom;
		double covarZcForDom;
		double covarZcDomDom;

		betatT = domModelParams->BetatT(expiryTime,payTime);

		covarFXDom = ARM_ModelParamsHW1FStd::HW1FEqFxStateCovariance(fxVol,domModelParams,0.0,expiryTime,expiryTime)*corrDomFx;
		covarZcForDom = ARM_ModelParamsHW1FStd::HW1FStateZcCovariance(domModelParams,forModelParams,0.0,expiryTime,expiryTime,expiryTime)*corrDomFor;
		covarZcDomDom = ARM_ModelParamsHW1FStd::HW1FStateZcCovariance(domModelParams,domModelParams,0.0,expiryTime,expiryTime,expiryTime);

		drift = betatT*(covarZcDomDom-covarZcForDom-covarFXDom);
	}
	else
	{
		double betatT;
		double covarFXFor;
		double covarZcForFor;
		double covarZcDomFor;

		betatT = forModelParams->BetatT(expiryTime,payTime);

		covarFXFor = ARM_ModelParamsHW1FStd::HW1FEqFxStateCovariance(fxVol,forModelParams,0.0,expiryTime,expiryTime)*corrDomFx;
		covarZcForFor = ARM_ModelParamsHW1FStd::HW1FStateZcCovariance(forModelParams,forModelParams,0.0,expiryTime,expiryTime,expiryTime)*corrDomFor;
		covarZcDomFor = ARM_ModelParamsHW1FStd::HW1FStateZcCovariance(domModelParams,forModelParams,0.0,expiryTime,expiryTime,expiryTime);

		drift = betatT*(covarZcDomFor-covarZcForFor-covarFXFor);
	}

	return drift;
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

