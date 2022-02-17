/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Heston_Fx.cpp
 *
 *  \brief Lognormal 1 factor FX version
 *
 *	\author  E. Ezzine
 *	\version 1.0
 *	\date May 2006
 */

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelnamemap.h"

/// gpmodel
#include "gpmodels/Heston_Fx.h"
#include "gpmodels/ModelParamsQ1F.h"
#include "gpmodels/Q1F.h"
#include "gpmodels/multiassets.h"
#include "gpmodels/ForwardMargin.h"

/// gpclosedforms
#include "gpclosedforms/inverse.h"
#include "gpclosedforms/heston_interface.h"
#include "gpclosedforms/heston_pricer.h"
#include "gpclosedforms/normal.h"

#include <inst/portfolio.h>

CC_BEGIN_NAMESPACE( ARM )

///  define for code clarity
#if defined(FIRST_STATE_VARIABLE)
	#undef	FIRST_STATE_VARIABLE
#endif
#define FIRST_STATE_VARIABLE 0

const double CALL_SPREAD_SHIFT_FX  =	0.0001;
const double MINSTRIKE = 1e-6;
const double MAXSTRIKE = 1e5;
const double NBLEVELS = 3;

////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Fx
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_HestonModel_Fx::ARM_HestonModel_Fx(const ARM_ZeroCurvePtr& zc, 
					 ARM_ModelParamsHeston_Fx* modelParam,
					 MCScheme mcScheme)
:	
ARM_EqFxBase(zc,modelParam,ARM_GP_Matrix(NULL),
new ARM_HestonDensityFunctor()),
itsIRDomModel(NULL),
itsIRForModel(NULL),
itsMCScheme(mcScheme),
itsLocalFunctional(NULL)
{
	if (!zc.IsNull() && modelParam)
		ARM_EqFxBase::Init();
	if(	!dynamic_cast<const ARM_ModelParamsHeston_Fx*>( modelParam ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": accepted model params are only heston Fx model param" );
}

////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Fx
///	Routine: Forward
///	Returns: a vector of forward (t,T)
///	Action : computes the forward of the equity or fx model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HestonModel_Fx::Forward(
	const string& curveName, 
	double evalTime,
	double expiryTime,
	double settlementTime,
	double payTime,
	const ARM_PricingStatesPtr& states) const
{

	double forwardValueMaturityTime	= ComputeFwdAtTime( settlementTime );
		
	if( evalTime<=K_NEW_DOUBLE_TOL || states == ARM_PricingStatesPtr(NULL) )
	{
		return ARM_VectorPtr( new std::vector<double>(1,forwardValueMaturityTime) );
	}
	else
    {
		size_t i,nbStates = states->size();
		size_t modelNb = GetModelNb();
		ARM_GP_VectorPtr domDf,forDf;
		ARM_GP_VectorPtr fwdFxValues;

        /// Spot Fx is path-dependently already diffused at current time
        if(settlementTime <= evalTime + K_NEW_DOUBLE_TOL)
        {
            domDf = ARM_GP_VectorPtr( new std::vector<double>(nbStates,1.0) );
            forDf = domDf;
        }
        else
        {
            domDf = itsIRDomModel->DiscountFactor(itsIRDomModel->GetModelName(),evalTime,settlementTime,states);
            forDf = itsIRForModel->DiscountFactor(itsIRForModel->GetModelName(),evalTime,settlementTime,states);
        }
		fwdFxValues = ARM_GP_VectorPtr(new ARM_GP_Vector(nbStates, 0.0));

		double drift = 1.0;

		if (GetNumeraire()->GetType() == ARM_Numeraire::TerminalZc)
		{
			double numTime = GetNumeraire()->GetMaturity();

			ARM_GP_Matrix correlMatrix = GetCorrelMatrix().Interpolate(evalTime);
			double corrDomFor    = correlMatrix(DomModel,ForModel);
			double corrDomFx    = correlMatrix(DomModel,FxModel);
			double corrForFx    = correlMatrix(ForModel,FxModel);

			const ARM_ModelParamsHW1FStd* const domParam   = static_cast<const ARM_ModelParamsHW1FStd*>( itsIRDomModel->GetModelParams() );
			const ARM_ModelParamsHW1FStd* const forParam   = static_cast<const ARM_ModelParamsHW1FStd*>( itsIRForModel->GetModelParams() );
			const ARM_CurveModelParam& fxVolParam = static_cast< const ARM_CurveModelParam& >(GetModelParams()->GetModelParam(ARM_ModelParamType::QVol));

			double betatT = domParam->BetatT(expiryTime,numTime);

			double covarFXDom = ARM_ModelParamsHW1FStd::HW1FEqFxStateCovariance(fxVolParam,domParam,0.0,expiryTime,expiryTime)*corrDomFor;
			double covarZcForDom = ARM_ModelParamsHW1FStd::HW1FStateZcCovariance(domParam,forParam,0.0,expiryTime,expiryTime,expiryTime)*corrDomFx;
			double covarZcDomDom = ARM_ModelParamsHW1FStd::HW1FStateZcCovariance(domParam,domParam,0.0,expiryTime,expiryTime,expiryTime);

			drift = exp(betatT*(covarZcDomDom-covarZcForDom+covarFXDom));
		}


		for( i=0; i<nbStates; ++i )
			(*fwdFxValues )[i] = states->GetModelState(i,modelNb) * states->GetModelState(i,modelNb+SpotHeston) * (*forDf)[i] / (*domDf)[i] * drift;

		if (!itsLocalFunctional.IsNull())
			fwdFxValues = itsLocalFunctional->Func(evalTime,fwdFxValues);

		return fwdFxValues;
    }


};

void ComputeAvgLevel(
int n,
double evalTime,
double expiryTime,
const ARM_Heston_ModelParams* fxParams,
double alpha,
std::vector<double>& avgTimes,
std::vector<double>& avgLevels
)
{
	const_cast<ARM_Heston_ModelParams*>(fxParams)->SetVolatilityType(ARM_ModelParamType::Volatility);
	avgTimes.resize(n);
	avgLevels.resize(n);

	double deltat = (expiryTime-evalTime)/n;
	double prevTime = evalTime;
	double curTime = evalTime+deltat;

	for (int i = 0; i < n; ++i)
	{
		avgTimes[i] = curTime/K_YEAR_LEN;
		avgLevels[i] = sqrt(fxParams->StateLocalVariance(prevTime,curTime,curTime)/(curTime-prevTime)*K_YEAR_LEN)*(1.0-alpha);
		prevTime = curTime;
		curTime = curTime+deltat;
	}
	const_cast<ARM_Heston_ModelParams*>(fxParams)->SetVolatilityType(ARM_ModelParamType::Sigma);


}


////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Fx
///	Routine: CallVectorial
///	Returns: a vector of call vectorial
///	Action : computes equity or fx call
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HestonModel_Fx::CallVectorial(
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
	ARM_QModel1F* domModel = GetStochasticModel( itsIRDomModel );
	ARM_QModel1F* forModel = GetStochasticModel( itsIRForModel );

	if(		evalTime<=K_NEW_DOUBLE_TOL
		 || states == ARM_PricingStatesPtr(NULL) )
	{
		ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
		size_t nbStates = states!=ARM_PricingStatesPtr(NULL)? states->size(): 1;
		ARM_GP_VectorPtr fwd = Forward(modelName, evalTime, expiryTime, settlementTime, payTime, states );
		double strike = strikePerState[0];
		double level = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(expiryTime);
		double kappa = GetModelParams()->GetModelParam(ARM_ModelParamType::VolMeanReversion).GetValue(expiryTime);
		double v0 = GetModelParams()->GetModelParam(ARM_ModelParamType::InitialVol).GetValue(expiryTime);
		double theta = GetModelParams()->GetModelParam(ARM_ModelParamType::LongTermVol).GetValue(expiryTime);
		double rho = GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).GetValue(expiryTime);
		double nu = GetModelParams()->GetModelParam(ARM_ModelParamType::VolOfVol).GetValue(expiryTime);
		double shift = GetModelParams()->GetModelParam(ARM_ModelParamType::Beta).GetValue(expiryTime);

		std::vector<double> times = ((ARM_CurveModelParam&)GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetAbscisses();
		std::vector<double> levels = ((ARM_CurveModelParam&)GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetOrdinates();

		const ARM_Heston_ModelParams* fxParams = static_cast< const ARM_Heston_ModelParams*>(GetModelParams());
		double alpha = fxParams->Alpha();

		times /= K_YEAR_LEN;
		levels *= (1.0-alpha);

		std::vector<double> avgTimes;
		std::vector<double> avgLevels;
		ComputeAvgLevel(NBLEVELS,evalTime,expiryTime,fxParams,alpha,avgTimes,avgLevels);

		double sigma;
		if (domModel && forModel)
		{
			double cov;
			sigma = VarianceFwdFx(evalTime,expiryTime,settlementTime,domModel,forModel,fxParams,false,cov);
		}
		else		
		{
			sigma = alpha*fxParams->StateLocalVariance(evalTime,expiryTime,expiryTime);
		}

		sigma = sqrt(sigma/expiryTime*K_YEAR_LEN);

		double zcT	= ZcCurve->DiscountPrice(payTime/K_YEAR_LEN);

		double price;

		if (domModel && forModel)
		{
			ARM_MixteHestonOptionPricer pricer(
				expiryTime/K_YEAR_LEN,
				(*fwd)[0], 
				strike, 
				callPut,
				sigma,
				v0, 
				kappa, 
				theta, 
				rho, 
				nu, 
				shift,
				/*times,
				levels,*/
				avgTimes,
				avgLevels);

			price = pricer.price();
		}
		else
		{
			ARM_MixteHestonOptionPricer pricer(
				expiryTime/K_YEAR_LEN,
				(*fwd)[0], 
				strike, 
				callPut,
				sigma,
				v0, 
				kappa, 
				theta, 
				rho, 
				nu,
				shift,
				level);
			price = pricer.price();
		}

		return ARM_VectorPtr( new std::vector<double>(nbStates,price*zcT) );
	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Heston FX can be used just for its closed forms." );

		ARM_GP_VectorPtr dumyFwdVector;
		return dumyFwdVector;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Fx
///	Routine: DigitalVectorial
///	Returns: a vector of digital vectorial
///	Action : computes equity or fx digital by CALLSPREAD
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HestonModel_Fx::DigitalVectorial(
	const string& modelName,
    double evalTime,
	double expiryTime,
	double settlementTime,
	const std::vector<double>& strikePerState,
	double notional,
	int callPut,
	double payTime,
	ARM_DigitType digitType,
	double epsilon,
    const ARM_PricingStatesPtr& states,
    ARM_PricingContext* context) const
{
	if(		evalTime<=K_NEW_DOUBLE_TOL
		 || states == ARM_PricingStatesPtr(NULL) )
	{
		ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
		size_t nbStates = states!=ARM_PricingStatesPtr(NULL)? states->size(): 1;
		ARM_GP_VectorPtr fwd = Forward(modelName, evalTime, expiryTime, settlementTime, payTime, states );
		double strike = strikePerState[0];
		double level = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(expiryTime);
		double kappa = GetModelParams()->GetModelParam(ARM_ModelParamType::VolMeanReversion).GetValue(expiryTime);
		double v0 = GetModelParams()->GetModelParam(ARM_ModelParamType::InitialVol).GetValue(expiryTime);
		double theta = GetModelParams()->GetModelParam(ARM_ModelParamType::LongTermVol).GetValue(expiryTime);
		double rho = GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).GetValue(expiryTime);
		double nu = GetModelParams()->GetModelParam(ARM_ModelParamType::VolOfVol).GetValue(expiryTime);
		double shift = GetModelParams()->GetModelParam(ARM_ModelParamType::Beta).GetValue(expiryTime);

		double zcT	= ZcCurve->DiscountPrice(payTime/K_YEAR_LEN);

		double eps = CALL_SPREAD_SHIFT_FX;
		double signed_eps = callPut*eps;
		ARM_HestonOptionPricer pricerUp(
			expiryTime/K_YEAR_LEN,
			(*fwd)[0], 
			strike + signed_eps, 
			callPut, 
			v0, 
			kappa, 
			theta, 
			rho, 
			nu, 
			shift,
			level);

		ARM_HestonOptionPricer pricerDown(
			expiryTime/K_YEAR_LEN,
			(*fwd)[0], 
			strike - signed_eps, 
			callPut, 
			v0, 
			kappa, 
			theta, 
			rho, 
			nu, 
			shift,
			level);

		double priceUp = pricerUp.price();
		double priceDown = pricerDown.price();

		double digit = zcT*(priceDown - priceUp)/eps;

		return ARM_VectorPtr( new std::vector<double>(nbStates,digit) );
	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Heston FX can be used just for its closed forms." );

		ARM_GP_VectorPtr dumyFwdVector;
		return dumyFwdVector;
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Fx
///	Routine: UpdateDensityFunctor
///	Returns: void
///	Action : Update the density functor at expiryTime		
////////////////////////////////////////////////////

void ARM_HestonModel_Fx::UpdateDensityFunctor(double forward, double expiryTime) 
{
	ARM_HestonDensityFunctor* densityfunctor = dynamic_cast<ARM_HestonDensityFunctor*>(&*GetDensityFunctor());
	if(!densityfunctor)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : the density functor must be of the type ARM_HestonDensityFunctor");

	double vol0 = GetModelParams()->GetModelParam(ARM_ModelParamType::InitialVol).GetValue(expiryTime);
	double kappa = GetModelParams()->GetModelParam(ARM_ModelParamType::VolMeanReversion).GetValue(expiryTime);
	double theta = GetModelParams()->GetModelParam(ARM_ModelParamType::LongTermVol).GetValue(expiryTime);
	double rho = GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).GetValue(expiryTime);
	double nu = GetModelParams()->GetModelParam(ARM_ModelParamType::VolOfVol).GetValue(expiryTime);
	double shift = GetModelParams()->GetModelParam(ARM_ModelParamType::Beta).GetValue(expiryTime);
	double level = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(expiryTime);

	densityfunctor->setV0( vol0 );
	densityfunctor->setKappa( kappa );
	densityfunctor->setTheta( theta );
	densityfunctor->setRho( rho );
	densityfunctor->setNu( nu );
	densityfunctor->setShift( shift );
	densityfunctor->setLevel( level );
}

////////////////////////////////////////////////////
///	Class   : ARM_HestonModel_Fx
///	Routines: ComputeSettlementCalendar
///	Returns : string
///	Action  : get the settlement calendar
////////////////////////////////////////////////////
string ARM_HestonModel_Fx::ComputeSettlementCalendar(const string& modelName) const
{
	const ARM_ModelParams* modelParam= GetModelParams();
	const ARM_ModelParamsHeston_Fx* fxModelParams = (const ARM_ModelParamsHeston_Fx*) modelParam;

#if defined(__GP_STRICT_VALIDATION)
	if( !dynamic_cast<const ARM_ModelParamsHeston_Fx*>(modelParam) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": model params is not ARM_ModelParamsHeston_Fx*" );
#endif
	char FXCal[7];
	strcpy(FXCal, fxModelParams->GetDomCurve()->GetCurrencyUnit()->GetCcyName());
	strcat(FXCal, fxModelParams->GetForCurve()->GetCurrencyUnit()->GetCcyName());
	return string(FXCal);
}


////////////////////////////////////////////////////
///	Class   : ARM_HestonModel_Fx
///	Routines: ComputeSettlementGap
///	Returns : double
///	Action  : get the settlement gap
////////////////////////////////////////////////////
double ARM_HestonModel_Fx::ComputeSettlementGap(const string& modelName) const
{
	const ARM_ModelParams* modelParam= GetModelParams();
	const ARM_ModelParamsHeston_Fx* fxModelParams = (const ARM_ModelParamsHeston_Fx*) modelParam;

#if defined(__GP_STRICT_VALIDATION)
	if( !dynamic_cast<const ARM_ModelParamsHeston_Fx*>(modelParam) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": model params is not ARM_HestonModel_Fx*" );
#endif
	double domGap = fxModelParams->GetDomCurve()->GetCurrencyUnit()->GetSpotDays();
	double forGap = fxModelParams->GetForCurve()->GetCurrencyUnit()->GetSpotDays();
	return CC_Max(domGap,forGap);
}


////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Fx
///	Routine: GetHestonFXModelParams
///	Returns: ARM_ModelParamsHeston_Fx*
///	Action : get the fxm model param
////////////////////////////////////////////////////

ARM_ModelParamsHeston_Fx* ARM_HestonModel_Fx::GetHestonFXModelParams()
{
	ARM_ModelParams* modelParam= GetModelParams();
	ARM_ModelParamsHeston_Fx* fxModelParams = (ARM_ModelParamsHeston_Fx*) modelParam;

#if defined(__GP_STRICT_VALIDATION)
	if( !dynamic_cast<const ARM_ModelParamsHeston_Fx*>(modelParam) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": model params is not ARM_ModelParamsHeston_Fx*" );
#endif
	return fxModelParams;
}

////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Fx
///	Routine: SetIRForeignModel
///	Returns: void
///	Action : sets the interest rate foreign model
////////////////////////////////////////////////////

void ARM_HestonModel_Fx::SetIRForeignModel( const ARM_PricingModelPtr& irModel )
{
	ValidateIRModel( irModel );
	itsIRForModel = irModel;
	GetHestonFXModelParams()->SetForCurve( irModel->GetZeroCurve()  );
}


////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Fx
///	Routine: SetIRForeignModel
///	Returns: void
///	Action : sets the interest rate domestic model
////////////////////////////////////////////////////
void ARM_HestonModel_Fx::SetIRDomesticModel( const ARM_PricingModelPtr& irModel )
{
	ValidateIRModel( irModel );
	itsIRDomModel = irModel;
	GetHestonFXModelParams()->SetDomCurve( irModel->GetZeroCurve() );
}


////////////////////////////////////////////////////
///	Class   : ARM_HestonModel_Fx
///	Routine	: GetStochasticModel (static function)
///	Returns : ARM_QModel1F*
///	Action  : the stochastic model for the interest rate model
////////////////////////////////////////////////////

ARM_QModel1F* ARM_HestonModel_Fx::GetStochasticModel( const ARM_PricingModelPtr& model )
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
///	Class  : ARM_HestonModel_Fx
///	Routine: FirstPricingStates
///	Returns: ARM_PricingStatesPtr
///	Action  : create the first pricing state
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_HestonModel_Fx::FirstPricingStates( size_t bucketSize ) const
{
	ARM_PricingStatesPtr initStates( new ARM_PricingStates(bucketSize,3,0,3) );
    double initValue = 0.0;

	if( itsIRForModel != ARM_PricingModelPtr(NULL) && itsIRDomModel != ARM_PricingModelPtr(NULL) )
        initValue = ComputeFwdAtTime(0.0);

	const ARM_Heston_ModelParams*	modelParamHeston = dynamic_cast<const ARM_Heston_ModelParams*>(GetModelParams());

	double initVar = modelParamHeston->InitialVol();

	for(size_t i=0;i<bucketSize;++i)
	{
        initStates->SetModelState(i,0,initValue);
		initStates->SetModelState(i,SpotHeston,1.0);
		initStates->SetModelState(i,VarHeston,initVar);
	}

    return initStates;
}

////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Fx
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Default initialisation of the model and the
///          associated numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_HestonModel_Fx::Init(const string& payModelName, 
											   const ARM_TimeInfoPtrVector& timeInfos)
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
///	Class  : ARM_HestonModel_Fx
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_HestonModel_Fx::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
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

        ARM_GP_VectorPtr domDf(NULL),forDf(NULL);

		double fwdFX0=0.0;

		domDf = itsIRDomModel->DiscountFactor( itsIRDomModel->GetModelName(), evalTime, nextTime, states );
		forDf = itsIRForModel->DiscountFactor( itsIRForModel->GetModelName(), evalTime, nextTime, states );

		const ARM_Heston_ModelParams*	modelParamHeston = dynamic_cast<const ARM_Heston_ModelParams*>(GetModelParams());

		double dt = (nextTime-evalTime)/K_YEAR_LEN;
		double racdt = sqrt(dt);
		double vvol = modelParamHeston->VolOfVol(nextTime);
		double theta = modelParamHeston->LongTermVol(nextTime);
		double kappa = modelParamHeston->VolMeanReversion(nextTime);
		double rho = modelParamHeston->Correlation(nextTime);
		double rho_ = sqrt(1.0-rho*rho);
		double alpha_ = 1.0-modelParamHeston->Alpha();
		double scale = modelParamHeston->Scaling(nextTime)*alpha_;
		double adjt = exp(- kappa * dt);
		double wesp2 = 1. + 0.5 * vvol * vvol / (kappa * theta);
		double w1E2_var = 0.5 * theta * vvol * vvol * SQR(1. - exp(-kappa*dt));
		double w2E2_var = vvol * vvol * exp(-kappa*dt)*(1.-exp(-kappa*dt))/kappa;
		double E_var, E2_var, psi, tipmo, a_var, b_var, invbeta, p_var, u_var;
		double spot, var, racvar, vvar, varf;
		double K0, K1, K2, K3, K4, w1 = 0.5, w2 = 0.5;
		

		for(i=0;i<nbStates;++i )
        {
          	const ARM_GP_Matrix& localVar = *(GetModelStateLocalVars()[timeIndex]);
			double lnDrift = -0.5*localVar(modelNb,modelNb);
            /// Compute S(ti,ti+1)
			double fwdFx = (*modelStates)(modelNb,i)*(*forDf)[i]/(*domDf)[i];

			/// Diffuse fwd Fx to ti+1
			(*modelStates)(modelNb,i) = fwdFx*exp((*numMethodStates)(modelNb,i) + lnDrift);

			spot = (*modelStates)(modelNb+SpotHeston,i);
			var = (*modelStates)(modelNb+VarHeston,i);

			if (itsMCScheme == Euler)
			{
				racvar = sqrt(var);
				E_var		= var * adjt + theta * (1. - adjt);
				E2_var		= w1E2_var + w2E2_var * var + E_var*E_var;
				vvar	= log(E2_var / (E_var * E_var));
				varf		= vvar > 0. ? E_var * exp(-0.5 * vvar + sqrt(vvar) * (*numMethodStates)(modelNb+VarHeston,i)) : E_var;
				spot *= exp(-0.5*var*scale*scale*dt+racvar*scale*racdt*(rho*(*numMethodStates)(modelNb+VarHeston,i)+rho_*(*numMethodStates)(modelNb+SpotHeston,i)));
			}
			else
			{
				E_var		= var * adjt + theta * (1. - adjt);
				E2_var		= w1E2_var + w2E2_var * var;

				psi			= E2_var / (E_var*E_var);
				
				if(psi < 1.5)
				{
					tipmo	= 2./psi - 1.;
					b_var	= sqrt(tipmo + sqrt(tipmo * tipmo + tipmo));
					a_var	= E_var / (1. + b_var*b_var);
					varf	= a_var * (b_var + (*numMethodStates)(modelNb+VarHeston,i)) * (b_var + (*numMethodStates)(modelNb+VarHeston,i));
				}
				else
				{
					invbeta	= 0.5 * E_var * (psi + 1.);
					p_var	= (psi - 1.)/(psi + 1.);
					u_var	= NormalCDF((*numMethodStates)(modelNb+VarHeston,i));
					varf	= u_var < p_var ? 0. : invbeta * log((1.- p_var)/(1.- u_var));
				}

				K0	= - rho * kappa * theta * dt * scale / vvol;
				K1	= w1 * dt * (kappa * rho * scale / vvol - 0.5 * scale * scale) - rho * scale/ vvol;
				K2	= w1 * dt * (kappa * rho * scale / vvol - 0.5 * scale * scale) + rho * scale/ vvol;
				K3	= w1 * dt * (1. - rho*rho) * scale * scale;
				K4	= w2 * dt * (1. - rho*rho) * scale * scale;

				spot	*= exp(K0 + K1 * var + K2 * varf + sqrt(K3 * var + K4 * varf) * (*numMethodStates)(modelNb+SpotHeston,i));

				
			}

			(*modelStates)(modelNb+SpotHeston,i) = spot;
			(*modelStates)(modelNb+VarHeston,i) = varf;
        }
	}
	else 
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": MCModelStatesFromToNextTime is not implemented when the interets rates are deterministic." );
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_HestonModel_Fx
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
double ARM_HestonModel_Fx::VarianceFwdFx(double a, double b, double settlementTime,
                                      const ARM_QModel1F* domRefModel,
                                      const ARM_QModel1F* forRefModel,
                                      const ARM_Heston_ModelParams* fxParams,
                                      bool isFwdZcFwdFxCov,
                                      double& fwdZcFwdFxCov) const
{
    const ARM_ModelParamsQ1F* domZcParams   = static_cast< const ARM_ModelParamsQ1F*>(domRefModel->GetModelParams());
    const ARM_ModelParamsQ1F* forZcParams   = static_cast< const ARM_ModelParamsQ1F* >(forRefModel->GetModelParams());

	double alpha = fxParams->Alpha();
	double spotFxVar = alpha*fxParams->StateLocalVariance(a,b,b);

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
        ARM_CurveModelParam fxVolParam = fxParams->GetVolAlpha();

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
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_HestonModel_Fx: VarianceFwdFx is not implemented when domestic and foreign interest rate model are not gaussian" );
        }
    }
    else
        fwdFxVar = 0.0;

    return fwdFxVar;
}

////////////////////////////////////////////////////
///	Class   : ARM_HestonModel_Fx
///	Routine : DigitalVectorial
///	Returns : a vector of digital vectorial
///	Action  : computes equity or fx digital
////////////////////////////////////////////////////
void ARM_HestonModel_Fx::UpdateLinks( const ARM_MultiAssetsModel& multiAssetsModel )
{
	const ARM_ModelNameMap* modelMap = multiAssetsModel.GetModelMap();
	/// find linked models
	ARM_IntVector itsOtherModels = (*modelMap)[GetModelName()]->OtherModelRefNb();

	if( !itsOtherModels.empty() ){
		// check that there is two and only two...
		if( itsOtherModels.size() != 2 )
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_QModel1F_Fx: wrong number of otherModels" );

		itsIRDomModel = (*modelMap)[ itsOtherModels[DomModel] ]->Model();
		itsIRForModel = (*modelMap)[ itsOtherModels[ForModel] ]->Model();
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Fx
///	Routine: LocalDrifts
///	Returns: void
///	Action : local drifts
////////////////////////////////////////////////////
void ARM_HestonModel_Fx::IntegratedLocalDrifts(
	const std::vector<double>& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts) const
{
	size_t nbSteps	= timeSteps.size();
    double step		= timeSteps[0],nextStep;
	relativeDrifts	= ARM_GP_MatrixPtr( new ARM_GP_Matrix( nbSteps-1, 3, 0.0 ) );
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
///	Class  : ARM_HestonModel_Fx
///	Routine: AdviseBreakPointTimes
///	Returns: void 
///	Action : sets the corresponding suggested break point times to the model param
////////////////////////////////////////////////////
void ARM_HestonModel_Fx::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* inputModelParam, size_t factorNb )
{
	ARM_CurveModelParam* modelParam = dynamic_cast<ARM_CurveModelParam*>(inputModelParam);
    double asOfDate			= GetAsOfDate().GetJulian();
    size_t portfolioSize	= portfolio->size();  
    std::vector<double>  tmpdates;
    size_t i;
	
	for( i=0; i<portfolioSize; ++i )
	{
		double expiry = portfolio->GetAsset(i)->GetExpiryDate().GetJulian() - asOfDate;
		tmpdates.push_back(expiry);
	}

	if (tmpdates.size())
		modelParam->UpdateValues(&tmpdates);
}


////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Fx
///	Routine: Proba
///	Returns: double 
///	Action : compute the cumulative
////////////////////////////////////////////////////

double ARM_HestonModel_Fx::Proba(
		double settlementTime,
		double strike) const
{
	if (strike < MINSTRIKE)
		return 0.0;

	ARM_QModel1F* domModel = GetStochasticModel( itsIRDomModel );
	ARM_QModel1F* forModel = GetStochasticModel( itsIRForModel );

	double level = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(settlementTime);
	double kappa = GetModelParams()->GetModelParam(ARM_ModelParamType::VolMeanReversion).GetValue(settlementTime);
	double v0 = GetModelParams()->GetModelParam(ARM_ModelParamType::InitialVol).GetValue(settlementTime);
	double theta = GetModelParams()->GetModelParam(ARM_ModelParamType::LongTermVol).GetValue(settlementTime);
	double rho = GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).GetValue(settlementTime);
	double nu = GetModelParams()->GetModelParam(ARM_ModelParamType::VolOfVol).GetValue(settlementTime);
	double shift = GetModelParams()->GetModelParam(ARM_ModelParamType::Beta).GetValue(settlementTime);

	std::vector<double> times = ((ARM_CurveModelParam&)GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetAbscisses();
	std::vector<double> levels = ((ARM_CurveModelParam&)GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetOrdinates();

	const ARM_Heston_ModelParams* fxParams = static_cast< const ARM_Heston_ModelParams*>(GetModelParams());
	double alpha = fxParams->Alpha();

	times /= K_YEAR_LEN;
	levels *= (1.0-alpha);

	double sigma = 0.0;

	if (domModel && forModel)
	{
		double cov;
		sigma = VarianceFwdFx(0.0,settlementTime,settlementTime,domModel,forModel,fxParams,false,cov);
	}
	else		
	{
		sigma = alpha*fxParams->StateLocalVariance(0.0,settlementTime,settlementTime);
	}

	sigma = sqrt(sigma/settlementTime*K_YEAR_LEN);

	double fwd = ComputeFwdAtTime( settlementTime );

	std::vector<double> avgTimes;
	std::vector<double> avgLevels;
	ComputeAvgLevel(NBLEVELS,0,settlementTime,fxParams,alpha,avgTimes,avgLevels);

	ARM_MixteHestonOptionPricer digital(
		settlementTime/K_YEAR_LEN,
		fwd, 
		strike, 
		K_PUT,
		sigma,
		v0, 
		kappa, 
		theta, 
		rho, 
		nu, 
		shift,
		avgTimes,
		avgLevels,
		true);

	return digital.price();
}

////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Fx
///	Routine: VolATM
///	Returns: double 
///	Action : compute the ATM volatility
////////////////////////////////////////////////////

double ARM_HestonModel_Fx::VolATM(
		double settlementTime) const
{
	ARM_QModel1F* domModel = GetStochasticModel( itsIRDomModel );
	ARM_QModel1F* forModel = GetStochasticModel( itsIRForModel );

	double level = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(settlementTime);
	double kappa = GetModelParams()->GetModelParam(ARM_ModelParamType::VolMeanReversion).GetValue(settlementTime);
	double v0 = GetModelParams()->GetModelParam(ARM_ModelParamType::InitialVol).GetValue(settlementTime);
	double theta = GetModelParams()->GetModelParam(ARM_ModelParamType::LongTermVol).GetValue(settlementTime);
	double rho = GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).GetValue(settlementTime);
	double nu = GetModelParams()->GetModelParam(ARM_ModelParamType::VolOfVol).GetValue(settlementTime);
	double shift = GetModelParams()->GetModelParam(ARM_ModelParamType::Beta).GetValue(settlementTime);

	std::vector<double> times = ((ARM_CurveModelParam&)GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetAbscisses();
	std::vector<double> levels = ((ARM_CurveModelParam&)GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetOrdinates();

	const ARM_Heston_ModelParams* fxParams = static_cast< const ARM_Heston_ModelParams*>(GetModelParams());
	double alpha = fxParams->Alpha();

	times /= K_YEAR_LEN;
	levels *= (1.0-alpha);

	double sigma = 0.0;

	if (domModel && forModel)
	{
		double cov;
		sigma = VarianceFwdFx(0.0,settlementTime,settlementTime,domModel,forModel,fxParams,false,cov);
	}
	else		
	{
		sigma = alpha*fxParams->StateLocalVariance(0.0,settlementTime,settlementTime);
	}

	sigma = sqrt(sigma/settlementTime*K_YEAR_LEN);

	double fwd = ComputeFwdAtTime( settlementTime );

	ARM_MixteHestonOptionPricer pricer(
		settlementTime/K_YEAR_LEN,
		fwd, 
		fwd, 
		K_CALL,
		sigma,
		v0, 
		kappa, 
		theta, 
		rho, 
		nu, 
		shift,
		times,
		levels);

	double price = pricer.price();

	return BSImpliedVol(fwd,fwd,settlementTime/K_YEAR_LEN,price,K_CALL);
}


////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Fx
///	Routine: HestonParameter
///	Returns: void
///	Action : compute the ATM volatility
////////////////////////////////////////////////////

void ARM_HestonModel_Fx::HestonParameter(
		const std::vector<double>& resetTimes,
		std::vector<double>& fwds,
		std::vector<double>& volATMs,
		std::vector<double>& times,
		std::vector<double>& levels,
		std::vector<double>& kappas,
		std::vector<double>& v0s,
		std::vector<double>& thetas,
		std::vector<double>& rhos,
		std::vector<double>& nus,
		std::vector<double>& shifts,
		std::vector<double>& sigmas
		) const
{

	fwds.resize(resetTimes.size());
	volATMs.resize(resetTimes.size());
	kappas.resize(resetTimes.size());
	v0s.resize(resetTimes.size());
	thetas.resize(resetTimes.size());
	rhos.resize(resetTimes.size());
	nus.resize(resetTimes.size());
	shifts.resize(resetTimes.size());
	sigmas.resize(resetTimes.size());

	ARM_QModel1F* domModel = GetStochasticModel( itsIRDomModel );
	ARM_QModel1F* forModel = GetStochasticModel( itsIRForModel );

	const ARM_Heston_ModelParams* fxParams = static_cast< const ARM_Heston_ModelParams*>(GetModelParams());
	double alpha = fxParams->Alpha();

	times = ((ARM_CurveModelParam&)GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetAbscisses();
	levels = ((ARM_CurveModelParam&)GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetOrdinates();

	times /= K_YEAR_LEN;
	levels *= (1.0-alpha);

	for (int i = 0; i < resetTimes.size(); ++i)
	{
		kappas[i] = GetModelParams()->GetModelParam(ARM_ModelParamType::VolMeanReversion).GetValue(resetTimes[i]);
		v0s[i] = GetModelParams()->GetModelParam(ARM_ModelParamType::InitialVol).GetValue(resetTimes[i]);
		thetas[i] = GetModelParams()->GetModelParam(ARM_ModelParamType::LongTermVol).GetValue(resetTimes[i]);
		rhos[i] = GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).GetValue(resetTimes[i]);
		nus[i] = GetModelParams()->GetModelParam(ARM_ModelParamType::VolOfVol).GetValue(resetTimes[i]);
		shifts[i] = GetModelParams()->GetModelParam(ARM_ModelParamType::Beta).GetValue(resetTimes[i]);
		if (domModel && forModel)
		{
			double cov;
			sigmas[i] = VarianceFwdFx(0.0,resetTimes[i],resetTimes[i],domModel,forModel,fxParams,false,cov);
		}
		else		
		{
			sigmas[i] = alpha*fxParams->StateLocalVariance(0.0,resetTimes[i],resetTimes[i]);
		}

		sigmas[i] = sqrt(sigmas[i]/resetTimes[i]*K_YEAR_LEN);

		fwds[i] = ComputeFwdAtTime( resetTimes[i] );
		volATMs[i] = VolATM(resetTimes[i]);
	}
}

void ARM_HestonModel_Fx::Probas(
		double settlementTime,
		const std::vector<double>& strikes,
		std::vector<double>& probasVec) const
{
	ARM_QModel1F* domModel = GetStochasticModel( itsIRDomModel );
	ARM_QModel1F* forModel = GetStochasticModel( itsIRForModel );

	double level = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(settlementTime);
	double kappa = GetModelParams()->GetModelParam(ARM_ModelParamType::VolMeanReversion).GetValue(settlementTime);
	double v0 = GetModelParams()->GetModelParam(ARM_ModelParamType::InitialVol).GetValue(settlementTime);
	double theta = GetModelParams()->GetModelParam(ARM_ModelParamType::LongTermVol).GetValue(settlementTime);
	double rho = GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).GetValue(settlementTime);
	double nu = GetModelParams()->GetModelParam(ARM_ModelParamType::VolOfVol).GetValue(settlementTime);
	double shift = GetModelParams()->GetModelParam(ARM_ModelParamType::Beta).GetValue(settlementTime);

	std::vector<double> times = ((ARM_CurveModelParam&)GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetAbscisses();
	std::vector<double> levels = ((ARM_CurveModelParam&)GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetOrdinates();

	const ARM_Heston_ModelParams* fxParams = static_cast< const ARM_Heston_ModelParams*>(GetModelParams());
	double alpha = fxParams->Alpha();

	times /= K_YEAR_LEN;
	levels *= (1.0-alpha);

	double sigma = 0.0;

	if (domModel && forModel)
	{
		double cov;
		sigma = VarianceFwdFx(0.0,settlementTime,settlementTime,domModel,forModel,fxParams,false,cov);
	}
	else		
	{
		sigma = alpha*fxParams->StateLocalVariance(0.0,settlementTime,settlementTime);
	}

	sigma = sqrt(sigma/settlementTime*K_YEAR_LEN);

	double fwd = ComputeFwdAtTime( settlementTime );

	ARM_MixteHestonOptionPricer digital(
		settlementTime/K_YEAR_LEN,
		fwd, 
		fwd, 
		K_PUT,
		sigma,
		v0, 
		kappa, 
		theta, 
		rho, 
		nu, 
		shift,
		times,
		levels,
		true);

	digital.setStrikes(strikes);
	digital.prices(probasVec);
}

class QuantileFunc : public DoubleToDoubleFunc
{
public: 
	QuantileFunc(
		double settlementTime,
		const ARM_HestonModel_Fx* model)
		:
	itsSettlementTime(settlementTime),
	itsModel(model)
	{
	};
	
	virtual double operator()(double strike) const 
	{
		return itsModel->Proba(itsSettlementTime,strike);
	}
private:
	double itsSettlementTime;
	const ARM_HestonModel_Fx* itsModel;
};



////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Fx
///	Routine: Quantile
///	Returns: void 
///	Action : compute the quantile
////////////////////////////////////////////////////

double ARM_HestonModel_Fx::Quantile(
	double settlementTime,
	double proba) const
{
	double fwd = ComputeFwdAtTime( settlementTime );
	QuantileFunc func(settlementTime, this);
	
	double root = Inverse(func,Inverse::REAL)(proba,fwd,fwd/5.,1e-12);
	
	return root;
}


////////////////////////////////////////////////////
///	Class   : ARM_QModel1F_Fx
///	Routines: toString
///	Returns :
///	Action  : stringify the object
////////////////////////////////////////////////////
string ARM_HestonModel_Fx::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Heston FX Model\n";
    os << indent << "-------------\n";
    os << ARM_PricingModel::toString(indent);

	if (!itsLocalFunctional.IsNull())
		os << itsLocalFunctional->toString();

    return os.str();
}


#undef FIRST_STATE_VARIABLE

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

