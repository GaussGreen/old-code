/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file CEV_Fx.cpp
 *
 *  \brief CEV model FX version
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date December 2004
 */

/// this header comes firts as it includes some preprocessor constants!

// gpbase
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/vectormanip.h"
#include "gpbase/comparisonfunctor.h"

//gpinfra
#include "gpinfra/modelnamemap.h"

// gpmodel
#include "gpmodels/QModelAnalytics.h"
#include "gpmodels/Q1F.h"
#include "gpmodels/ModelParamsQ1F.h"
#include "gpmodels/CEV_Fx.h"
#include "gpmodels/ForwardMargin.h"
#include "gpmodels/multiassets.h"

/// gpcalib
#include "gpcalib/calibmethod.h"

/// ARM Kernel
#include <inst/portfolio.h>
#include <inst/option.h>

/// gpclosedforms
#include "gpclosedforms/vanille_bs_formula.h"
#include "gpclosedforms/gaussian_integrals.h"

const size_t GL_PY_NBPOINTS         = 2;
const size_t GL_MIN_NBPOINTS        = 2;
const size_t DECAY_NBSTEP_PER_YEAR  = 4;

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_CEVModel_Fx
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_CEVModel_Fx::ARM_CEVModel_Fx(const ARM_ZeroCurvePtr& zc, 
								 ARM_ModelParamsCEV_Fx* modelParams,
								 const ARM_CurveMatrix& correlMatrix)
:	ARM_EqFxBase(zc,modelParams,correlMatrix), 
	itsIRDomModel(NULL),
	itsIRForModel(NULL)
{
	if (!zc.IsNull() && modelParams)
		ARM_EqFxBase::Init();
}


////////////////////////////////////////////////////
///	Class  : ARM_CEVModel_Fx
///	Routine: Copy Constructor
///	Returns: 
///	Action : Copy Constructor
////////////////////////////////////////////////////

ARM_CEVModel_Fx::ARM_CEVModel_Fx( const ARM_CEVModel_Fx& rhs )
:	
	ARM_EqFxBase(rhs),
	itsIRDomModel( rhs.itsIRDomModel ),
	itsIRForModel( rhs.itsIRForModel )
{}

////////////////////////////////////////////////////
///	Class  : ARM_CEVModel_Fx
///	Routine: destructor
///	Returns: 
///	Action : destructor
////////////////////////////////////////////////////
ARM_CEVModel_Fx::~ARM_CEVModel_Fx()
{
}

////////////////////////////////////////////////////
///	Class   : ARM_CEVModel_Fx
///	Routines: toString
///	Returns :
///	Action  : stringify the object
////////////////////////////////////////////////////
string ARM_CEVModel_Fx::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "CEV FX Model\n";
    os << indent << "-------------\n";
    os << ARM_PricingModel::toString(indent);
    return os.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_CEVModel_Fx
///	Routines: ComputeSettlementCalendar
///	Returns : string
///	Action  : get the settlement calendar
////////////////////////////////////////////////////
string ARM_CEVModel_Fx::ComputeSettlementCalendar(const string& modelName) const
{
	const ARM_ModelParams* modelParam= GetModelParams();
	const ARM_ModelParamsCEV_Fx* fxModelParams = (const ARM_ModelParamsCEV_Fx*) modelParam;

#if defined(__GP_STRICT_VALIDATION)
	if( !dynamic_cast<const ARM_ModelParamsCEV_Fx*>(modelParam) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": model params is not ARM_ModelParamsCEV_Fx*" );
#endif
	char FXCal[7];
	strcpy(FXCal, fxModelParams->GetDomCurve()->GetCurrencyUnit()->GetCcyName());
	strcat(FXCal, fxModelParams->GetForCurve()->GetCurrencyUnit()->GetCcyName());
	return string(FXCal);
}


////////////////////////////////////////////////////
///	Class   : ARM_CEVModel_Fx
///	Routines: ComputeSettlementGap
///	Returns : double
///	Action  : get the settlement gap
////////////////////////////////////////////////////
double ARM_CEVModel_Fx::ComputeSettlementGap(const string& modelName) const
{
	const ARM_ModelParams* modelParam= GetModelParams();
	const ARM_ModelParamsCEV_Fx* fxModelParams = (const ARM_ModelParamsCEV_Fx*) modelParam;

#if defined(__GP_STRICT_VALIDATION)
	if( !dynamic_cast<const ARM_ModelParamsCEV_Fx*>(modelParam) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": model params is not ARM_CEVModel_Fx*" );
#endif
	double domGap = fxModelParams->GetDomCurve()->GetCurrencyUnit()->GetSpotDays();
	double forGap = fxModelParams->GetForCurve()->GetCurrencyUnit()->GetSpotDays();
	return CC_Max(domGap,forGap);
}


////////////////////////////////////////////////////
///	Class  : ARM_CEVModel_Fx
///	Routine: GetCEVFXModelParams
///	Returns: ARM_ModelParamsCEV_Fx*
///	Action : get the fxm model param
////////////////////////////////////////////////////

ARM_ModelParamsCEV_Fx* ARM_CEVModel_Fx::GetCEVFXModelParams()
{
	ARM_ModelParams* modelParam= GetModelParams();
	ARM_ModelParamsCEV_Fx* fxModelParams = (ARM_ModelParamsCEV_Fx*) modelParam;

#if defined(__GP_STRICT_VALIDATION)
	if( !dynamic_cast<const ARM_ModelParamsCEV_Fx*>(modelParam) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": model params is not ARM_ModelParamsCEV_Fx*" );
#endif
	return fxModelParams;
}

////////////////////////////////////////////////////
///	Class  : ARM_CEVModel_Fx
///	Routine: SetIRForeignModel
///	Returns: void
///	Action : sets the interest rate foreign model
////////////////////////////////////////////////////

void ARM_CEVModel_Fx::SetIRForeignModel( const ARM_PricingModelPtr& irModel )
{
	ValidateIRModel( irModel );
	itsIRForModel = irModel;
	GetCEVFXModelParams()->SetForCurve( irModel->GetZeroCurve()  );
}


////////////////////////////////////////////////////
///	Class  : ARM_CEVModel_Fx
///	Routine: SetIRForeignModel
///	Returns: void
///	Action : sets the interest rate domestic model
////////////////////////////////////////////////////
void ARM_CEVModel_Fx::SetIRDomesticModel( const ARM_PricingModelPtr& irModel )
{
	ValidateIRModel( irModel );
	itsIRDomModel = irModel;
	GetCEVFXModelParams()->SetDomCurve( irModel->GetZeroCurve() );
}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_Fx
///	Routine: FirstPricingStates
///	Returns: ARM_PricingStatesPtr
///	Action  : create the first pricing state
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_CEVModel_Fx::FirstPricingStates( size_t bucketSize ) const
{
	ARM_PricingStatesPtr initStates( new ARM_PricingStates(bucketSize,1,0,1) );
    double initValue = 0.0;

	if( itsIRForModel != ARM_PricingModelPtr(NULL) && itsIRDomModel != ARM_PricingModelPtr(NULL) )
        initValue = ComputeFwdAtTime(0.0);

	for(size_t i=0;i<bucketSize;++i)
        initStates->SetModelState(i,0,initValue);

    return initStates;
}

////////////////////////////////////////////////////
///	Class  : ARM_CEVModel_Fx
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : move model states for timeIndex to timeIndex+1
////////////////////////////////////////////////////
void ARM_CEVModel_Fx::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
	if( itsIRForModel != ARM_PricingModelPtr(NULL) && itsIRDomModel != ARM_PricingModelPtr(NULL) )
	{
		const ARM_NumMethodPtr numMethod= GetNumMethod();
        size_t nbSteps = numMethod->GetTimeSteps()->size();
		double evalTime	= numMethod->GetTimeStep(timeIndex);
		double nextTime	= numMethod->GetTimeStep(timeIndex+1);
		size_t i,nbStates = states->size();
	    size_t modelNb	= GetModelNb();

        bool isLnFwdFx;

        ARM_GP_MatrixPtr numMethodStates    = states->GetNumMethodStates();
        ARM_GP_MatrixPtr modelStates        = states->GetModelStates();

	    const ARM_ModelParamsHW1FStd* modelParams = static_cast< const ARM_ModelParamsHW1FStd*>(GetModelParams());
        const ARM_Curve& betaFxCurve = *(modelParams->GetModelParam(ARM_ModelParamType::Beta).ToCurveModelParam().GetCurve());
        isLnFwdFx = betaFxCurve == 1.0;

        ARM_GP_VectorPtr domDf(NULL),forDf(NULL);
        ARM_GP_VectorPtr integDriftDom(NULL),integDriftFor(NULL);
        double beta,ti,tk,Bdom0ti,Bfor0ti,Bdom0tk,Bfor0tk;
        double domFwdRate,forFwdRate;

		double fwdFX0=0.0;

        if(isLnFwdFx)
        {
            /// Pure lognormal forward FX S(t,ti+1) is diffused from ti to ti+1 and
            /// then S(ti+1)=S(ti+1,ti+1)

			domDf = itsIRDomModel->DiscountFactor( itsIRDomModel->GetModelName(), evalTime, nextTime, states );
			forDf = itsIRForModel->DiscountFactor( itsIRForModel->GetModelName(), evalTime, nextTime, states );
        }
        else
        {
            ti = GetNumMethod()->GetTimeStep(timeIndex);
			tk = GetNumMethod()->GetTimeStep(timeIndex+1);

			fwdFX0 = ComputeFwdAtTime(ti);

            beta=betaFxCurve.Interpolate(ti);

            domDf = itsIRDomModel->DiscountFactor( itsIRDomModel->GetModelName(), 0.0, ti, states );
            Bdom0ti = (*domDf)[0];
            forDf = itsIRForModel->DiscountFactor( itsIRForModel->GetModelName(), 0.0, ti, states );
            Bfor0ti = (*forDf)[0];

            domDf = itsIRDomModel->DiscountFactor( itsIRDomModel->GetModelName(), 0.0, tk, states );
            Bdom0tk = (*domDf)[0];
            forDf = itsIRForModel->DiscountFactor( itsIRForModel->GetModelName(), 0.0, tk, states );
            Bfor0tk = (*forDf)[0];

            domFwdRate = log(Bdom0ti/Bdom0tk);
            forFwdRate = log(Bfor0ti/Bfor0tk);

            bool isDriftAdded = true;

	        integDriftDom = itsIRDomModel->IntegratedRiskNeutralDrift( timeIndex, modelStates, isDriftAdded );
	        integDriftFor = itsIRForModel->IntegratedRiskNeutralDrift( timeIndex, modelStates, isDriftAdded );
        }

        /// Compute lognormal drift of diffused variable (S(t,ti+1) or S(t))
        double spotFx,spotFx1,markovDrift,spotFxPower;
		for(i=0;i<nbStates;++i )
        {
            if(isLnFwdFx)
            {
				const ARM_GP_Matrix& localVar = *(GetModelStateLocalVars()[timeIndex]);
				double lnDrift = -0.5*localVar(modelNb,modelNb);
                /// Compute S(ti,ti+1)
			    double fwdFx = (*modelStates)(modelNb,i)*(*forDf)[i]/(*domDf)[i];

                /// Diffuse fwd Fx to ti+1
                (*modelStates)(modelNb,i) = fwdFx*exp((*numMethodStates)(modelNb,i) + lnDrift);
            }
            else
            {
                /// Get S(ti)
                spotFx = (*modelStates)(modelNb,i);
				spotFxPower = pow(spotFx/fwdFX0,beta-1);

				markovDrift = domFwdRate - forFwdRate + (*integDriftDom)[i] - (*integDriftFor)[i];
				spotFx1  = spotFx+markovDrift*spotFx + spotFxPower*spotFx*(*numMethodStates)(modelNb,i);
                
                (*modelStates)(modelNb,i) = spotFx1;
            }
        }
	}
	else 
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": MCModelStatesFromToNextTime is not implemented when the interets rates are deterministic." );
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_Fx
///	Routine: SetNumMethod
///	Returns: void 
///	Action : Set the num method
////////////////////////////////////////////////////
void ARM_CEVModel_Fx::SetNumMethod(const ARM_NumMethodPtr& numMethodPtr)
{
	ARM_PricingModel::SetNumMethod( numMethodPtr );
}


////////////////////////////////////////////////////
///	Class  : ARM_CEVModel_Fx
///	Routine: LocalDrifts
///	Returns: void
///	Action : local drifts
////////////////////////////////////////////////////
void ARM_CEVModel_Fx::IntegratedLocalDrifts(
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
		(*relativeDrifts)(i,0) = ((const ARM_ModelParamsQ1F* const) GetModelParams())->StateLocalDrift(step,nextStep);
		step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_CEVModel_Fx
///	Routine : ARM_PricingStatesPtr
///	Returns : Init
///	Action  : initialise the model
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_CEVModel_Fx::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
    int nbEvents=timeInfos.size();
    bool isSpotUse = (nbEvents == 0) || (nbEvents==1 && timeInfos[0]->GetEventTime() <= K_NEW_DOUBLE_TOL);

	// Compute a single model states set to (0.0,...,0.0)
    ARM_PricingStatesPtr initStates(new ARM_PricingStates(1,1,0));
    for(size_t i=0;i<1;++i)
        initStates->SetModelState(0,i,0.0);
    return initStates;
}


////////////////////////////////////////////////////
///	Class   : ARM_CEVModel_Fx
///	Routine	: GetStochasticModel (static function)
///	Returns : ARM_QModel1F*
///	Action  : the stochastic model for the interest rate model
////////////////////////////////////////////////////

ARM_QModel1F* ARM_CEVModel_Fx::GetStochasticModel( const ARM_PricingModelPtr& model )
{
	if( model != ARM_PricingModelPtr(NULL) )
	{
		if( dynamic_cast<ARM_ForwardMargin*>( &*model ) )
			return dynamic_cast<ARM_QModel1F*>( model->GetRefModel() );
		else
			return dynamic_cast<ARM_QModel1F*>( &*model );
	}
	else
		return NULL;
}


////////////////////////////////////////////////////
///	Class  : ARM_CEVModel_Fx
///	Routine: Forward
///	Returns: a vector of forward (t,T)
///	Action : computes the forward of fx model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_CEVModel_Fx::Forward(
		const string& modelName, 
        double evalTime,
	    double expiryTime,
	    double settlementTime,
	    double payTime,
        const ARM_PricingStatesPtr& states) const
{
	if( !IsFwdModels() || GetNumMethod() == ARM_NumMethodPtr(NULL)
        || states == ARM_PricingStatesPtr(NULL) )
	{
		double forwardValueMaturityTime	= ComputeFwdAtTime( settlementTime );
		
		if( evalTime<=K_NEW_DOUBLE_TOL
			|| states == ARM_PricingStatesPtr(NULL) )
		{
			return ARM_VectorPtr( new ARM_GP_Vector(1,forwardValueMaturityTime) );
		}
		else
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Forward can just evalate for the asof date without numerical method in ARM_ModelParamsCEV_Fx." );
			return ARM_VectorPtr(NULL);
		}
	
	}

    size_t i,nbStates = states->size();
    size_t modelNb = GetModelNb();
    ARM_GP_VectorPtr domDf,forDf;
    ARM_GP_VectorPtr fwdFxValues;

    if( GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_FWDLOOKING)
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
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": forward Fx is just implemented for forward looking numerical method." );
    }

    return fwdFxValues;
}


////////////////////////////////////////////////////
///	Class   : ARM_CEVModel_Fx
///	Routine	: VarianceFwdFx
///	Returns : a double
///	Action  : Computes the variance in [a,b] of the forward FX
///           S(.T)=S(.)Bf(.,T)/Bd(.,T) (T = fwd Fx settlement time)
///           S follows a Q model dynamics and the naïve approximation
///           at time 0, S(0,T).Bf(0,t,T)/Bd(0,t,T), is used in the
///           exact volatility formula for forward FX
///	          If needed, it also computes the covariance in [a,b]
///           of the forward Zc Bd(.,T)/Bf(.,T) and the forward FX
////////////////////////////////////////////////////
double ARM_CEVModel_Fx::VarianceFwdFx(double a, double b, double settlementTime,
                                      const ARM_QModel1F* domRefModel,
                                      const ARM_QModel1F* forRefModel,
                                      const ARM_ModelParamsCEV_Fx* fxParams,
                                      bool isFwdZcFwdFxCov,
                                      double& fwdZcFwdFxCov) const
{
    const ARM_ModelParamsQ1F* domZcParams   = static_cast< const ARM_ModelParamsQ1F*>(domRefModel->GetModelParams());
    const ARM_ModelParamsQ1F* forZcParams   = static_cast< const ARM_ModelParamsQ1F* >(forRefModel->GetModelParams());

	double spotFxVar = fxParams->StateLocalVariance(a,b,b);

    ARM_GP_Matrix correlMatrix = GetCorrelMatrix().Interpolate(a);

    double domForCorr	= correlMatrix(DomModel,ForModel);
	double domFxCorr	= correlMatrix(DomModel,FxModel);
	double forFxCorr	= correlMatrix(ForModel,FxModel);

    double fwdFxVar,domZcVar,forZcVar,domForZcCovar,domZcFxCovar,forZcFxCovar;
    double domFor,domFx,forFx;
    fwdZcFwdFxCov = 0.0;
    if(a + K_NEW_DOUBLE_TOL < b)
    {
        /// Restore fxVolCurve (underlying gaussian process volatility in fact)
        const ARM_CurveModelParam& fxVolParam = static_cast< const ARM_CurveModelParam& >(fxParams->GetModelParam(ARM_ModelParamType::Volatility));

        if(domRefModel->IsDegenerateInHW() && forRefModel->IsDegenerateInHW())
        {
            /// Both IR models are degenerated to H&W 1F
            domZcVar        = ARM_ModelParamsHW1F::HW1FZcCovariance(domZcParams,domZcParams,a,b,settlementTime);
            forZcVar        = ARM_ModelParamsHW1F::HW1FZcCovariance(forZcParams,forZcParams,a,b,settlementTime);
            domForZcCovar   = ARM_ModelParamsHW1F::HW1FZcCovariance(domZcParams,forZcParams,a,b,settlementTime);
            domZcFxCovar    = ARM_ModelParamsHW1F::HW1FEqFxZcCovariance(fxVolParam,domZcParams,a,b,settlementTime);
            forZcFxCovar    = ARM_ModelParamsHW1F::HW1FEqFxZcCovariance(fxVolParam,forZcParams,a,b,settlementTime);

            domFor      = domForCorr*domForZcCovar;
            domFx       = domFxCorr*domZcFxCovar;
            forFx       = forFxCorr*forZcFxCovar;

            ///  Var(fwdFx(.,T))
            fwdFxVar    = spotFxVar + forZcVar + domZcVar + 2*(forFx - domFx - domFor);

            /// Covar(Bd(.T)/Bf(.T),fwdFx(.,T))
            if(isFwdZcFwdFxCov)
                fwdZcFwdFxCov = -(forZcVar + domZcVar - 2*domFor + forFx - domFx);
        }
        else
        {
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_CEVModel_Fx: VarianceFwdFx is not implemented when domestic and foreign interest rate model are not gaussian" );
        }
    }
    else
        fwdFxVar = 0.0;

    return fwdFxVar;
}


////////////////////////////////////////////////////
///	Class   : ARM_CEVModel_Fx
///	Routine	: CallVectorial
///	Returns : a vector of call vectorial
///	Action  : computes fx call
////////////////////////////////////////////////////
ARM_VectorPtr ARM_CEVModel_Fx::CallVectorial(
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
    /// Restore reference IR models
	ARM_QModel1F* domRefModel = GetStochasticModel( itsIRDomModel );
	ARM_QModel1F* forRefModel = GetStochasticModel( itsIRForModel );

    if(expiryTime < evalTime - K_NEW_DOUBLE_TOL)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": the fx option resets before event time" );

    bool isIVCalc;
    if(isIVCalc = (expiryTime <= evalTime + K_NEW_DOUBLE_TOL))
        expiryTime = evalTime;  // force evaluation at reset time

	size_t i,nbStates = (evalTime > K_NEW_DOUBLE_TOL && states != ARM_PricingStatesPtr(NULL))
                        ? states->size(): 1;

#if defined( __GP_STRICT_VALIDATION )
	if( strikePerState.size() != nbStates )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " strike vector size != states size" );
#endif

    double yfEval   = evalTime/K_YEAR_LEN;
    double yfExpiry = expiryTime/K_YEAR_LEN;

    const ARM_ModelParamsCEV_Fx* fxParams = static_cast< const ARM_ModelParamsCEV_Fx*>(GetModelParams());

    double fwdFxQ,zcFxCov;
    double fwdFxVarTtoM = VarianceFwdFx(evalTime,expiryTime,settlementTime,domRefModel,forRefModel,fxParams,false,zcFxCov);

	const ARM_Curve& spotFxBetaCurve = *(fxParams->GetModelParam( ARM_ModelParamType::Beta).ToCurveModelParam().GetCurve());
    if(spotFxBetaCurve == 1.0)
    {
        /// Standard case : pure lognormal forward FX vol
        fwdFxQ = 1.0;
    }
    else if(!isIVCalc)
    {
        /// Compute a correction to get an approximated displaced forward FX volatility
        /// and shift (Piterbarg's "FX volatility skew" 2005 paper)

        /// Merge IR & FX vol schedules
        const ARM_ModelParamsQ1F* domParams   = static_cast< const ARM_ModelParamsQ1F*>(domRefModel->GetModelParams());
        const ARM_GP_Vector& domTimes     = ((ARM_CurveModelParam&) domParams->GetModelParam(domParams->GetVolatilityType())).GetCurve()->GetAbscisses();

        const ARM_ModelParamsQ1F* forParams   = static_cast< const ARM_ModelParamsQ1F* >(forRefModel->GetModelParams());
        const ARM_GP_Vector& forTimes     = ((ARM_CurveModelParam&) forParams->GetModelParam(forParams->GetVolatilityType())).GetCurve()->GetAbscisses();

        const ARM_GP_Vector& fxTimes     = ((ARM_CurveModelParam&) fxParams->GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetAbscisses();

		const ARM_GP_Vector correlTimes = GetCorrelMatrix().GetTimes();

        ARM_GP_VectorPtr times1(MergeSortedVectorNoDuplicates(domTimes,forTimes));
		ARM_GP_VectorPtr times2(MergeSortedVectorNoDuplicates(*times1,fxTimes));
        ARM_GP_VectorPtr times(MergeSortedVectorNoDuplicates(*times2,correlTimes));

        /// Find n(a) & n(b) such that Tn(a)-1 < evalTime <= Tn(a) and Tn(b)-1< expiryTime <= Tn(b)
	    int na = lower_boundPosWithPrecision(*times,evalTime);
	    int nb = lower_boundPosWithPrecision(*times,expiryTime);

	    ARM_GP_Matrix correlMatrix;

	    double domFxCorr	= 0.0;
	    double forFxCorr	= 0.0;

        double fwdFxVarTtot,fwdZcFwdFxCovTtot;
        double lastt=evalTime,yfLastt=yfEval,t,yft,tt,lasttt,zcFxCov,w;
        double volCEVt,khiRatiot,spotFxBetat,etat,volDomZct,volForZct;
        fwdFxQ = 0.0;
        fwdFxVarTtot = 0.0;
        fwdZcFwdFxCovTtot = 0.0;

        size_t n,nbTimes=times->size();
        size_t i,nbPoints;
        for(n=na;n<=nb;++n)
        {
            t   = (n<nbTimes ? ((*times)[n]>expiryTime ? expiryTime : (*times)[n]): expiryTime);
            yft = t/K_YEAR_LEN;

			correlMatrix = GetCorrelMatrix().Interpolate(t);

            /// A numerical integration is needed to get the equivalent shift on forward FX
            nbPoints=static_cast<int>(ceil(yft-yfLastt)*GL_PY_NBPOINTS);
            if(nbPoints<GL_MIN_NBPOINTS)
                nbPoints = GL_MIN_NBPOINTS;

            GaussLegendre_Coefficients GLPoints(nbPoints,yfLastt,yft);

            /// Gauss-Legendre summation between each [Ti,Ti+1]
            lasttt = lastt;
            for(i=0;i<nbPoints;++i)
            {
                tt= GLPoints.get_point(i)*K_YEAR_LEN;
                w = GLPoints.get_weight(i);

                fwdFxVarTtot        += VarianceFwdFx(lasttt,tt,settlementTime,domRefModel,forRefModel,fxParams,true,zcFxCov);
                fwdZcFwdFxCovTtot   += zcFxCov;
                khiRatiot           = fwdZcFwdFxCovTtot/fwdFxVarTtot;

	            spotFxBetat        = fxParams->GetModelParam( ARM_ModelParamType::Beta).ToCurveModelParam().GetCurve()->Interpolate(tt);
	            volCEVt          = fxParams->GetModelParam(ARM_ModelParamType::Volatility).ToCurveModelParam().GetCurve()->Interpolate(tt);
                etat            = volCEVt*(1.0+khiRatiot)*(spotFxBetat-1.0);

                volDomZct       = domRefModel->VolZc(tt,settlementTime);
                volForZct       = forRefModel->VolZc(tt,settlementTime);

                fwdFxQ          += w * fwdFxVarTtot*etat*(volCEVt + volForZct*forFxCorr - volDomZct*domFxCorr);

                lasttt = tt;
            }

            if(n<nb)
            {
                fwdFxVarTtot        += VarianceFwdFx(lasttt,t,settlementTime,domRefModel,forRefModel,fxParams,true,zcFxCov);
                fwdZcFwdFxCovTtot   += zcFxCov;
            }
            else
                fwdFxVarTtot        += VarianceFwdFx(lasttt,t,settlementTime,domRefModel,forRefModel,fxParams,false,zcFxCov);

            lastt   = t;
            yfLastt = yft;
        }

        if(fabs(fwdFxVarTtoM-fwdFxVarTtot)>K_NEW_DOUBLE_TOL)
    		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
            " : inconsistency in variances computation to get an equivalent forward FX Q parameter" );

        fwdFxQ = 1.0 + 2.0*fwdFxQ/(fwdFxVarTtoM*fwdFxVarTtoM);
    }

	/// Compute option values without memory allocation by replacing fwd by option values
    ARM_VectorPtr zcPay = itsIRDomModel->DiscountFactor(itsIRDomModel->GetModelName(),evalTime,payTime,states);
	ARM_VectorPtr fwdFx = Forward(modelName, evalTime, expiryTime, settlementTime, payTime, states );

    ARM_GP_Vector* fxOption = new ARM_GP_Vector(nbStates);

    double yfTtoM = yfExpiry-yfEval;
	double fwdFxVol=0.0;
    if(isIVCalc)
    {
        /// Just compute the option intrinsic value
	    for(i=0; i<nbStates; ++i )
		    (*fxOption)[i] = CC_Max<double>(callPut*((*fwdFx)[i]-strikePerState[i]),0.0)*(*zcPay)[i];

    }
    else
    {
        /// Compute a Q model price formula
        fwdFxVol = sqrt(fwdFxVarTtoM/yfTtoM);
	    for(i=0; i<nbStates; ++i )
		    (*fxOption)[i] = QModelAnalytics::BSQFunction( (*fwdFx)[i], strikePerState[i], fwdFxVol, yfTtoM, fwdFxQ, (*zcPay)[i], callPut, (*fwdFx)[i] );
    }

	return ARM_VectorPtr(fxOption);
}


////////////////////////////////////////////////////
///	Class  : ARM_CEVModel_Fx
///	Routine: DigitalVectorial
///	Returns: a vector of digital vectorial
///	Action : computes equity or fx digital
////////////////////////////////////////////////////
ARM_VectorPtr ARM_CEVModel_Fx::DigitalVectorial(
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
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "ARM_CEVModel_Fx::DigitalVectorial is not implemented");
}


////////////////////////////////////////////////////
///	Class   : ARM_CEVModel_Fx
///	Routine : UpdateLinks
///	Returns : void
///	Action  : update links between FX and interest rate model
////////////////////////////////////////////////////
void ARM_CEVModel_Fx::UpdateLinks( const ARM_MultiAssetsModel& multiAssetsModel )
{
	const ARM_ModelNameMap* modelMap = multiAssetsModel.GetModelMap();
	/// find linked models
	ARM_IntVector itsOtherModels = (*modelMap)[GetModelName()]->OtherModelRefNb();

	if( !itsOtherModels.empty() ){
		// check that there is two and only two...
		if( itsOtherModels.size() != 2 )
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_CEVModel_Fx: wrong number of otherModels" );

		itsIRDomModel = (*modelMap)[ itsOtherModels[DomModel] ]->Model();
		itsIRForModel = (*modelMap)[ itsOtherModels[ForModel] ]->Model();
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_CEVModel_Fx
///	Routine: AdviseBreakPointTimes
///	Returns: void 
///	Action : sets the corresponding suggested break point times to the model param
////////////////////////////////////////////////////
void ARM_CEVModel_Fx::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* inputModelParam, size_t factorNb )
{
	ARM_CurveModelParam* modelParam = dynamic_cast<ARM_CurveModelParam*>(inputModelParam);
	if( !modelParam )
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "expected an ARM_CurveModelParam!");

    double asOfDate			= GetAsOfDate().GetJulian();
    size_t portfolioSize	= portfolio->size();  
    ARM_GP_Vector  tmpdates;
    size_t i;
	
	/// FX option
	bool isFxOption;
	double date;
	isFxOption = (portfolio->GetAsset(0)->GetName() == ARM_OPTION);
	if(!isFxOption)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
			": portfolio can only contain FX options or IR/FX hybrid swaptions" );
	}
	else
	{
		date = static_cast<ARM_Option*>(portfolio->GetAsset(0))->GetExpiry().GetJulian() - asOfDate;
		tmpdates.push_back(date);
	}
	for( i=1; i<portfolioSize; ++i )
	{
		isFxOption = (portfolio->GetAsset(i)->GetName() == ARM_OPTION);
		if(!isFxOption)
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
				": portfolio can only contain FX options or IR/FX hybrid swaptions" );
		}
		else
		{
			double resetlag = static_cast<ARM_Option*>(portfolio->GetAsset(i))->GetExpiry().GetJulian() - asOfDate;
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
}


////////////////////////////////////////////////////
///	Class   : ARM_CEVModel_Fx
///	Routine : ValidateModelParams
///	Returns : bool
///	Action  : validate model params
////////////////////////////////////////////////////
bool ARM_CEVModel_Fx::ValidateModelParams(const ARM_ModelParams& params) const
{
	if( !dynamic_cast<const ARM_ModelParamsCEV_Fx*>(&params))
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelparams is not of type ARM_ModelParamsHeston_EqFx" );
	return true;
}


CC_END_NAMESPACE()