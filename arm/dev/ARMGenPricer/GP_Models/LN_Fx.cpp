/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file LN_Fx.cpp
 *
 *  \brief Lognormal 1 factor FX version
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2004
 */

/// gpinfra
#include "gpinfra/nummethod.h"
#include "gpinfra/discretisationscheme.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/modelnamemap.h"
#include "gpinfra/pricingcontext.h"

/// gpmodel
#include "gpmodels/LN_Fx.h"
#include "gpmodels/HW1F.h"
#include "gpmodels/multiassets.h"

/// gpnummethods
#include "gpnummethods/tree1D.h"
#include "gpnummethods/markoviandriftsampler.h"
#include "gpnummethods/meanrevertingsampler.h"

// gpclosedforms
#include "gpclosedforms/vanilla_bs.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_LN_Fx
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_LN_Fx::ARM_LN_Fx(const ARM_ZeroCurvePtr& zc, 
					 ARM_ModelParamsLN_Fx* modelParams,
					 const ARM_CurveMatrix& correlMatrix)
:	ARM_EqFxBase(zc,modelParams,correlMatrix), 
	itsIRDomModel(NULL),
	itsIRForModel(NULL)
{
	if (!zc.IsNull() && modelParams)
		ARM_EqFxBase::Init();
}

void ARM_LN_Fx::SetFactorCount(size_t fc)
{
	(static_cast<ARM_ModelParamsLN_Fx*> (GetModelParams()))->SetFactorCount(fc);
}

////////////////////////////////////////////////////
///	Class  : ARM_LN_Fx
///	Routine: Copy Constructor
///	Returns: 
///	Action : Copy Constructor
////////////////////////////////////////////////////

ARM_LN_Fx::ARM_LN_Fx( const ARM_LN_Fx& rhs )
:	
	ARM_EqFxBase(rhs),
	itsIRDomModel( rhs.itsIRDomModel ),
	itsIRForModel( rhs.itsIRForModel )
{
}

////////////////////////////////////////////////////
///	Class   : ARM_LN_Fx
///	Routine : UpdateLinks
///	Returns : void
///	Action  : 
////////////////////////////////////////////////////
void ARM_LN_Fx::UpdateLinks( const ARM_MultiAssetsModel& multiAssetsModel )
{
	const ARM_ModelNameMap* modelMap = multiAssetsModel.GetModelMap();
	/// find linked models
	ARM_IntVector itsOtherModels = (*modelMap)[GetModelName()]->OtherModelRefNb();

	if( !itsOtherModels.empty() ){
		// check that there is two and only two...
		if( itsOtherModels.size() != 2 )
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_LN_Fx: wrong number of otherModels" );

		itsIRDomModel = (*modelMap)[ itsOtherModels[DomModel] ]->Model();
		itsIRForModel = (*modelMap)[ itsOtherModels[ForModel] ]->Model();
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_LN_Fx
///	Routine: Forward
///	Returns: a vector of forward (t,T)
///	Action : computes the forward of the equity or fx model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_LN_Fx::Forward(
	const string& curveName, 
	double evalTime,
	double expiryTime,
	double settlementTime,
	double payTime,
	const ARM_PricingStatesPtr& states) const
{
/*	double forwardValuefwdTime	= ComputeFwdAtTime( settlementTime );
	size_t modelNb = GetModelNb();
    
	
	if(		evalTime<=K_NEW_DOUBLE_TOL
		 || states == ARM_PricingStatesPtr(NULL) )
	{
		size_t payoffSize = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
		return ARM_VectorPtr( new ARM_GP_Vector(payoffSize,forwardValuefwdTime) );
	}
	else
	{
		size_t stateSize = states->size();
		ARM_GP_VectorPtr fwdVector = ARM_GP_VectorPtr( new ARM_GP_Vector(stateSize, 0.0 ));
		double forward = ComputeFwdAtTime( settlementTime );
		double stateLocalVar = static_cast< const ARM_ModelParamsHW1F* const>(GetModelParams())->StateLocalVariance(0.0,evalTime,evalTime);
		double C = 0.5*stateLocalVar;

		for( size_t i=0; i<stateSize; ++i )
			(*fwdVector )[i] = forward * exp( states->GetModelState(i,modelNb)-C );
		return fwdVector;
	}*/
	//Get itsLambda which enables to choose the PDE and the reconstruction formula(centred, non centred, Qmapping, SpotFX)
	const ARM_NumMethodPtr& numMethod = GetNumMethod();
	if( !IsFwdModels() || GetNumMethod() == ARM_NumMethodPtr(NULL)
        || states == ARM_PricingStatesPtr(NULL) || (evalTime<=K_NEW_DOUBLE_TOL))
	{
		double forwardValuefwdTime	= ComputeFwdAtTime( settlementTime );
		size_t modelNb = GetModelNb();
    
		
		if(		evalTime<=K_NEW_DOUBLE_TOL
			 || states == ARM_PricingStatesPtr(NULL) )
		{
			size_t payoffSize = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
			return ARM_VectorPtr( new ARM_GP_Vector(payoffSize,forwardValuefwdTime) );
		}
		else
		{
			size_t stateSize = states->size();
			ARM_GP_VectorPtr fwdVector = ARM_GP_VectorPtr( new ARM_GP_Vector(stateSize, 0.0 ));
			double forward = ComputeFwdAtTime( settlementTime );
			double stateLocalVar = static_cast< const ARM_ModelParamsHW1F* const>(GetModelParams())->StateLocalVariance(0.0,evalTime,evalTime);
			double C = 0.5*stateLocalVar;

			for( size_t i=0; i<stateSize; ++i )
				(*fwdVector )[i] = forward * exp( states->GetModelState(i,modelNb)-C );
			return fwdVector;
		}
	}

    size_t i,nbStates = states->size();
    size_t modelNb = GetModelNb();
    ARM_GP_VectorPtr domDf,forDf;
    ARM_GP_VectorPtr fwdFxValues;

    if( (GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_FWDLOOKING)||(GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_FWDBCKWDLOOKING))
    {
        /// Spot Fx is path-dependently already diffused at current time
        if(settlementTime <= evalTime + K_NEW_DOUBLE_TOL)
        {
            domDf = ARM_GP_VectorPtr( new ARM_GP_Vector(nbStates,1.0) );
            forDf = domDf;
        }
        else
        {
            domDf = itsIRDomModel->DiscountFactor(itsIRDomModel->GetModelName(),evalTime,settlementTime,states);
            forDf = itsIRForModel->DiscountFactor(itsIRForModel->GetModelName(),evalTime,settlementTime,states);
        }
		fwdFxValues = ARM_GP_VectorPtr(new ARM_GP_Vector(nbStates, 0.0));

		for( i=0; i<nbStates; ++i )
			(*fwdFxValues )[i] = states->GetModelState(i,modelNb) * (*forDf)[i] / (*domDf)[i];
    }
    else
    {
        /// Yield curves are stochastic : spot is computed using the mapping function then
        /// calibrated domestic/foreign zero-coupons give the forward Fx

		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": To DO");
    }
 
    return fwdFxValues;
};


////////////////////////////////////////////////////
///	Class  : ARM_LN_Fx
///	Routine: CallVectorial
///	Returns: a vector of call vectorial
///	Action : computes equity or fx call
////////////////////////////////////////////////////
ARM_VectorPtr ARM_LN_Fx::CallVectorial(
	const string& modelName,
    double evalTime,
	double expiryTime,
	double settlementTime,
	const ARM_GP_Vector& strikePerState,
	int callPut,
	double payTime,
    const ARM_PricingStatesPtr& states,
    ARM_PricingContext* context) const
{
    if(expiryTime < evalTime - K_NEW_DOUBLE_TOL)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": the fx option resets before event time" );

    bool isIVCalc;
    if(isIVCalc = (expiryTime <= evalTime + K_NEW_DOUBLE_TOL))
        expiryTime = evalTime;  // force evaluation at reset time

	size_t i,nbStates = (evalTime > K_NEW_DOUBLE_TOL && states != ARM_PricingStatesPtr(NULL)) ? states->size(): 1;

#if defined( __GP_STRICT_VALIDATION )
	if( strikePerState.size() != nbStates )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " strike vector size != states size" );
#endif

	/*
	*/
    double yfEval   = evalTime/K_YEAR_LEN;
    double yfExpiry = expiryTime/K_YEAR_LEN;

    ARM_VectorPtr fwdFx = Forward(modelName, evalTime, expiryTime, settlementTime, payTime, states );
	ARM_VectorPtr zcPay; 
	
    ARM_GP_Vector* fxOption = new ARM_GP_Vector(nbStates);

    double yfTtoM = yfExpiry-yfEval;

	double fwdFxVol = 0.0;

    if(isIVCalc)
    {
        /// Just compute the option intrinsic value
	    for(i=0; i<nbStates; ++i )
		    (*fxOption)[i] = CC_Max<double>(callPut*((*fwdFx)[i]-strikePerState[i]),0.0)*(*zcPay)[i];
    }
    else
    {
		/// Restore reference IR models
		ARM_HullWhite1F* domRefModel = 0;
		ARM_HullWhite1F* forRefModel = 0;
		if( itsIRForModel != ARM_PricingModelPtr(NULL) && itsIRDomModel != ARM_PricingModelPtr(NULL) )
		{
			domRefModel = dynamic_cast< ARM_HullWhite1F* > (itsIRDomModel->GetRefModel());
			forRefModel = dynamic_cast< ARM_HullWhite1F* > (itsIRForModel->GetRefModel());
		}
		const ARM_ModelParamsHW1FStd* const fxModelParams   = static_cast<const ARM_ModelParamsHW1FStd* const>( GetModelParams() );
		const ARM_CurveModelParam& fxVol                    = dynamic_cast<const ARM_CurveModelParam&>( fxModelParams->GetModelParam(fxModelParams->GetVolatilityType()) );
		
		double fwdFxVarTtoM;
		if (domRefModel && forRefModel)
		{
			zcPay = itsIRDomModel->DiscountFactor(itsIRDomModel->GetModelName(),evalTime,payTime,states);
			// only valid for constant correlations!!!!
			ARM_CurveMatrix curveCorrel = GetCorrelMatrix();
			ARM_GP_Matrix correlMatrix = curveCorrel.Interpolate(0);
		
			const ARM_ModelParamsHW1FStd* const domModelParams  = static_cast<const ARM_ModelParamsHW1FStd* const>( domRefModel->GetModelParams() );
			const ARM_ModelParamsHW1FStd* const forModelParams  = static_cast<const ARM_ModelParamsHW1FStd* const>( forRefModel->GetModelParams() );

			double domForCorr = correlMatrix(DomModel,ForModel);
			double domFxCorr = correlMatrix(DomModel,FxModel);
			double forFxCorr = correlMatrix(ForModel,FxModel);

			double varFx	= fxModelParams->StateLocalVariance(evalTime,expiryTime,expiryTime);
			double varZcDom	= ARM_ModelParamsHW1F::HW1FZcCovariance( domModelParams, domModelParams, evalTime,expiryTime,expiryTime );
			double varZcFor	= ARM_ModelParamsHW1F::HW1FZcCovariance( forModelParams, forModelParams, evalTime,expiryTime,expiryTime );

			double covarZcDomZcFor = ARM_ModelParamsHW1F::HW1FZcCovariance( domModelParams, forModelParams, evalTime,expiryTime,expiryTime);
			double covarZcDomFx    = ARM_ModelParamsHW1F::HW1FEqFxZcCovariance( fxVol, domModelParams, evalTime,expiryTime,expiryTime );
			double covarZcForFx    = ARM_ModelParamsHW1F::HW1FEqFxZcCovariance( fxVol, forModelParams, evalTime,expiryTime,expiryTime );

			fwdFxVarTtoM = varFx + varZcDom + varZcFor
				- 2*domForCorr*covarZcDomZcFor
				- 2*domFxCorr*covarZcDomFx
				+ 2*forFxCorr*covarZcForFx;
		}
		else
		{
			ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
			double zcT	= ZcCurve->DiscountPrice(payTime/K_YEAR_LEN);
			double zct	= ZcCurve->DiscountPrice(evalTime/K_YEAR_LEN);
			zcPay = ARM_VectorPtr(new ARM_GP_Vector(nbStates,zcT/zct));
			fwdFxVarTtoM = static_cast< const ARM_ModelParamsHW1F* const>(GetModelParams())->StateLocalVariance(evalTime,expiryTime,expiryTime);
		}
	
        fwdFxVol = sqrt(fwdFxVarTtoM/yfTtoM);
		for(i=0; i<nbStates; ++i )
			(*fxOption)[i] = BlackSholes_Formula( (*fwdFx)[i], fwdFxVol, (*zcPay)[i], strikePerState[i], yfTtoM, callPut );
			
    }
    if(context)
    {
        /// Save pricing context
        ARM_PricingCallContext* callContext = context->ToCallContext();
		callContext->SetVol(fwdFxVol);
        callContext->SetShift(1.0);
        callContext->SetForward(fwdFx);
    }
	return ARM_VectorPtr(fxOption);
}


////////////////////////////////////////////////////
///	Class  : ARM_LN_Fx
///	Routine: DigitalVectorial
///	Returns: a vector of digital vectorial
///	Action : computes equity or fx digital
////////////////////////////////////////////////////
ARM_VectorPtr ARM_LN_Fx::DigitalVectorial(
	const string& modelName,
    double evalTime,
	double expiryTime,
	double settlementTime,
	const ARM_GP_Vector& strikePerState,
	double notional,
	int callPut,
	double payTime,
	ARM_DigitType digitType,
	double epsilon,
    const ARM_PricingStatesPtr& states,
    ARM_PricingContext* context) const
{
    /// FIX FIX : review formula with payTime !!

	ARM_GP_VectorPtr result = Forward(modelName, evalTime, expiryTime, settlementTime, payTime, states );

#if defined(__GP_STRICT_VALIDATION)
	if( strikePerState.size() != result->size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
                ": strikePerState.size() != forward.size()" );
#endif

	if(	evalTime<=K_NEW_DOUBLE_TOL|| states == ARM_PricingStatesPtr(NULL) )
	{
		/// get the intrinsic value
		size_t payoffSize = result->size();
		for( size_t i=0; i<payoffSize; ++i )
			(*result)[i] = callPut*(*result)[i]>callPut*strikePerState[i];
	}
	else
	{
		ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
		size_t nbStates = states!=ARM_PricingStatesPtr(NULL)? states->size(): 1;

		double totalVolatility = static_cast< const ARM_ModelParamsHW1F* const>(GetModelParams())->StateLocalVariance(evalTime,expiryTime,expiryTime);
		double zcT	= ZcCurve->DiscountPrice(settlementTime/K_YEAR_LEN);
		double zct	= ZcCurve->DiscountPrice(evalTime/K_YEAR_LEN);
		double df	= zcT/zct;

		/// get the intrinsic value
		size_t payoffSize = result->size();
		for( size_t i=0; i<payoffSize; ++i )
			(*result)[i] = DigitalBlackSholes_Formula((*result)[i], totalVolatility, df, strikePerState[i], callPut );
	}

	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_LN_Fx
///	Routine: FirstPricingStates
///	Returns: ARM_PricingStatesPtr
///	Action  : create the first pricing state
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_LN_Fx::FirstPricingStates( size_t bucketSize ) const
{
	/// ARM_PricingStates(nbStates = bucketSize, nbModelStates = 1F , nbPayoffs = 0)
	ARM_PricingStatesPtr initStates( new ARM_PricingStates(bucketSize,1,0,1) );
    double initValue = 0.0;
	if( itsIRForModel != ARM_PricingModelPtr(NULL) && itsIRDomModel != ARM_PricingModelPtr(NULL) )
        initValue = ComputeFwdAtTime(0.0);

    for(size_t i=0;i<bucketSize;++i)
        initStates->SetModelState(i,0,initValue);

    return initStates;
}


////////////////////////////////////////////////////
///	Class   : ARM_LN_Fx
///	Routines: MarkovianDrift
///	Returns : void
///	Action  : Computes the Markovian drift
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_LN_Fx::MarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates ) const
{
	return ARM_GP_MatrixPtr( new ARM_GP_Matrix(numMethodStates->rows(), numMethodStates->cols(), 0.0 ) );
}


////////////////////////////////////////////////////
///	Class   : ARM_LN_Fx
///	Routines: LocalDiscounts
///	Returns : void
///	Action  : Computes the LocalDiscounts
////////////////////////////////////////////////////
ARM_VectorPtr ARM_LN_Fx::LocalDiscounts(
	size_t timeIdx, 
	double dt, 
	const ARM_PricingStatesPtr& states) const
{
	size_t stateSize= states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
	const ARM_GP_Vector* const timeSteps = GetNumMethod()->GetTimeSteps();
	double startTime = (*timeSteps)[timeIdx];

	ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
	double zct	= ZcCurve->DiscountPrice(startTime/K_YEAR_LEN);
	double zcT	= ZcCurve->DiscountPrice((startTime+dt)/K_YEAR_LEN);

	return ARM_VectorPtr( new ARM_GP_Vector( stateSize, zcT/zct ) );
}



////////////////////////////////////////////////////
///	Class  : ARM_LN_Fx
///	Routine: LocalDrifts
///	Returns: void
///	Action : local drifts
////////////////////////////////////////////////////
void ARM_LN_Fx::IntegratedLocalDrifts(
	const ARM_GP_Vector& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts) const
{
	size_t nbSteps	= timeSteps.size();
    double step		= timeSteps[0],nextStep;
	relativeDrifts	= ARM_GP_MatrixPtr( new ARM_GP_Matrix( nbSteps-1, 1, 0.0 ) );
	absoluteDrifts	= ARM_GP_MatrixPtr(NULL);

	for(size_t i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];
		/// [i] => local variance from ti->ti+1
		(*relativeDrifts)(i,0) = static_cast<const ARM_ModelParamsHW1F* const> (GetModelParams())->StateLocalDrift(step,nextStep);
		step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_LN_Fx
///	Routine: ModelStateLocalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local and global variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_LN_Fx::ModelStateLocalVariances(
	const ARM_GP_Vector& timeSteps,
	ARM_MatrixVector& localVariances ) const
{
	size_t nbSteps	= timeSteps.size();
	size_t modelNb	= GetModelNb();
    double step		= timeSteps[0],
		   nextStep;
	size_t offsetIndex	= (nbSteps-1)*modelNb;
	localVariances.resize((nbSteps-1)*(modelNb+1));
	size_t i;

	for(i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];
		
		/// [i] => local variance from ti->ti+1
		localVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(1,((const ARM_ModelParamsHW1F* const) GetModelParams())->StateLocalVariance(step,nextStep,nextStep));
		
		step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_LN_Fx
///	Routine: SetNumMethod
///	Returns: void 
///	Action : Set the num method
////////////////////////////////////////////////////
void ARM_LN_Fx::SetNumMethod(const ARM_NumMethodPtr& numMethodPtr)
{
	if( !GetFromMultiFactor() )
	{
		const ARM_TreeBase* treeBase = dynamic_cast<const ARM_TreeBase*>(&*numMethodPtr);
		if(treeBase)
		{
			if( dynamic_cast<const ARM_Tree1D*>(treeBase) )
			{
				if( !dynamic_cast<const ARM_MarkovianDriftSampler1D*>(treeBase->GetSampler()) &&
					!dynamic_cast<const ARM_MeanRevertingSampler1D*>(treeBase->GetSampler()))
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
					" Q model requires Markovian drift or Mean revering samplers !" );
			}
			else
			{
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " tree is not 1D compatible" );
			}
		}
	}
	ARM_PricingModel::SetNumMethod( numMethodPtr );
}

////////////////////////////////////////////////////
///	Class   : ARM_LN_Fx
///	Routine : ValidateModelParams
///	Returns : bool
///	Action  : validate model params
////////////////////////////////////////////////////
bool ARM_LN_Fx::ValidateModelParams(const ARM_ModelParams& params) const
{
	const ARM_ModelParamsLN_Fx* ModelParamsLN = dynamic_cast<const ARM_ModelParamsLN_Fx*>(&params);
	if( !ModelParamsLN )
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelparams is not of type ARM_ModelParamsLN_Fx" );
	}
	return true;
}

////////////////////////////////////////////////////
///	Class  : ARM_LN_Fx
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Default initialisation of the model and the
///          associated numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_LN_Fx::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
    int nbEvents=timeInfos.size();
    bool isSpotUse = nbEvents == 0 || (nbEvents==1 && timeInfos[0]->GetEventTime() <= K_NEW_DOUBLE_TOL);
	
    if(!isSpotUse)
    {
        // Numerical method and numeraire are needed to go on
        ARM_NumMethodPtr numMethod=GetNumMethod();
		if( numMethod == ARM_NumMethodPtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numerical method not set in Q model!");

        /// test the numeraire and its type!
		ARM_NumerairePtr numeraire=GetNumeraire();
        if( numeraire == ARM_NumerairePtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numeraire not set in the Q model!");

        if( numeraire->GetType() != ARM_Numeraire::Cash )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
            ": only Cash numeraire supported by LN_Eq model at the moment!");
		
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
///	Class  : ARM_LN_Fx
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_LN_Fx::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
#if defined(__GP_STRICT_VALIDATION)
	if(timeIndex<0)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index is negative!");
	if( timeIndex >= GetNumMethod()->GetTimeSteps()->size()-1 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index bigger than max size!");
#endif

	if( itsIRForModel != ARM_PricingModelPtr(NULL) && itsIRDomModel != ARM_PricingModelPtr(NULL) )
	{
		const ARM_NumMethodPtr numMethod= GetNumMethod();
        size_t nbSteps = numMethod->GetTimeSteps()->size();
		double evalTime	= numMethod->GetTimeStep(timeIndex);
		double nextTime	= numMethod->GetTimeStep(timeIndex+1);
		size_t i,nbStates = states->size();
	    size_t modelNb	= GetModelNb();

		ARM_GP_MatrixPtr numMethodStates    = states->GetNumMethodStates();
		ARM_GP_MatrixPtr modelStates        = states->GetModelStates();

        const ARM_ModelParamsHW1F* modelParams = static_cast< const ARM_ModelParamsHW1F*>(GetModelParams());
            
		ARM_GP_VectorPtr domDf(NULL),forDf(NULL);
		ARM_GP_VectorPtr integDriftDom(NULL),integDriftFor(NULL);

        /// Pure lognormal forward FX S(t,ti+1) is diffused from ti to ti+1 and
        /// then S(ti+1)=S(ti+1,ti+1)

		domDf = itsIRDomModel->DiscountFactor( itsIRDomModel->GetModelName(), evalTime, nextTime, states );
		forDf = itsIRForModel->DiscountFactor( itsIRForModel->GetModelName(), evalTime, nextTime, states );
            
        /// Compute lognormal drift of diffused variable (S(t,ti+1) or S(t))
        const ARM_GP_Matrix& localVar = *(GetModelStateLocalVars()[timeIndex]);
        double lnDrift = -0.5*localVar(modelNb,modelNb);

		for(i=0;i<nbStates;++i )
        {
			/// Compute S(ti,ti+1)
			double fwdFx = (*modelStates)(modelNb,i)*(*forDf)[i]/(*domDf)[i];

            /// Diffuse fwd Fx to ti+1
            (*modelStates)(modelNb,i) = fwdFx * exp((*numMethodStates)(modelNb,i) + lnDrift);
		}
	}
	else 
	{
		double nextTime = GetNumMethod()->GetTimeStep(timeIndex+1);
		size_t factorsNb= FactorCount();
		size_t statesNb = states->size();
		double currentState;
		size_t modelNb	= GetModelNb();
		
		for( size_t i=0;i<statesNb; ++i )
		{
			for( size_t j=0;  j<factorsNb; ++j )
			{
				currentState   = 0.0;
				for( size_t k =0; k<=j; ++k )
				{
					double gaussian = states->GetNumMethodState(i,modelNb+k);
					currentState  += gaussian;
				}
				states->SetModelState(i,j+modelNb,currentState);
			}
		}
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_LN_Fx
///	Routine: TreeStatesToModelStates
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_LN_Fx::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{   
	#if defined(__GP_STRICT_VALIDATION)
	if(timeIndex<0)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index is negative!");
	if( timeIndex > GetNumMethod()->GetTimeSteps()->size()-1 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index bigger than max size!");
#endif

	const ARM_GP_MatrixPtr& numMethodStates = states->GetNumMethodStates();
	states->SetModelStates(numMethodStates);
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

