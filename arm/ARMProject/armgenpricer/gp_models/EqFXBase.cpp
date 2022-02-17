/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file EqFxBase.cpp
 *
 *  \brief base class for equity and FX model
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date May 2006
 */

/// to remove identified warning
#include "gpbase/removeidentifiedwarning.h"

/// include of the header file
#include "gpmodels/EqFxBase.h"

/// gpbase
#include "gpbase/singleton.h" /// for numeraire factory
#include "gpbase/interpolatorvector.h"
#include "gpbase/datestrip.h"
#include "gpbase/cloneutilityfunc.h"
#include "gpbase/globalconstant.h"


/// gpinfra
#include "gpinfra/nummethod.h"
#include "gpinfra/discretisationscheme.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/pricingcontext.h"

/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/densityfunctors.h"

/// gpmodels
#include "gpmodels/fxname.h"
#include "gpmodels/ModelParamsHW1F.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/BS_Model.h"
#include "gpmodels/SABR_Model.h"

/// gpnummethods
#include "gpnummethods/cfmethod.h"
#include "gpnummethods/tree1D.h"

/// gpclosedforms
#include "gpclosedforms/vanille_bs_formula.h"

/// gpcalculators
#include "gpcalculators/forexvanilla.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_EqFxBase
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_EqFxBase::ARM_EqFxBase(
	const ARM_ZeroCurvePtr& zc, 
	ARM_ModelParams* modelParam, 
	const ARM_CurveMatrix& correlMatrix,
	ARM_DensityFunctor* densityFct,
	ARM_PricingModelIR* convAdjustModel)
:	ARM_PricingModel(zc,modelParam, densityFct),
	itsSettlementCalendar( "UNITIALIZED"), 
	itsSettlementGap(-1111111),
	itsCorrelMatrix(correlMatrix),
	itsCallType(ARM_EqFxBase::ClosedFormula),
	itsConvAdjustModel(NULL),
	itsRho(ARM_FlatCurve(0.7)) // to be consistent with PRDC
{
}

////////////////////////////////////////////////////
///	Class  : ARM_EqFxBase
///	Routine: Init
///	Returns: 
///	Action : Initialize
////////////////////////////////////////////////////

void ARM_EqFxBase::Init(){}/*

{

	const ARM_ModelParams_Fx* fxModelParams = dynamic_cast<const ARM_ModelParams_Fx*>(GetModelParams());

#if defined(__GP_STRICT_VALIDATION)
	if( !dynamic_cast<const ARM_ModelParams_Fx*>(GetModelParams()) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": model params is not ARM_ModelParams_Fx*" );
#endif

	ARM_ZeroCurvePtr domCurve = fxModelParams->GetDomCurve();
	ARM_ZeroCurvePtr forCurve = fxModelParams->GetForCurve();

	string domCcy = domCurve->GetCurrencyUnit()->GetCcyName();
	string forCcy = forCurve->GetCurrencyUnit()->GetCcyName();

	SetModelName(forCcy+domCcy);

	itsSettlementCalendar = ComputeSettlementCalendar( GetModelName() );
	itsSettlementGap = ComputeSettlementGap( GetModelName() );
}*/

////////////////////////////////////////////////////
///	Class  : ARM_EqFxBase
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_EqFxBase::ARM_EqFxBase(const ARM_EqFxBase& rhs)
:	ARM_PricingModel(rhs),
	itsSettlementCalendar(rhs.itsSettlementCalendar), 
	itsSettlementGap(rhs.itsSettlementGap),
	itsCorrelMatrix(rhs.itsCorrelMatrix),
	itsCallType(rhs.itsCallType),
	itsConvAdjustModel(CreateClone(rhs.itsConvAdjustModel)),
	itsRho(rhs.itsRho)
{}

////////////////////////////////////////////////////
///	Class   : ARM_EqFxBase
///	Routines: ComputeSettlementCalendar
///	Returns : string
///	Action  : get the settlement calendar
////////////////////////////////////////////////////
string ARM_EqFxBase::ComputeSettlementCalendar(const string& modelName) const
{
	const ARM_ModelParams_Fx* fxModelParams = dynamic_cast<const ARM_ModelParams_Fx*>(GetModelParams());

#if defined(__GP_STRICT_VALIDATION)
	if( !dynamic_cast<const ARM_ModelParams_Fx*>(GetModelParams()) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": model params is not ARM_ModelParams_Fx*" );
#endif
	char FXCal[7];
	strcpy(FXCal, fxModelParams->GetDomCurve()->GetCurrencyUnit()->GetCcyName());
	strcat(FXCal, fxModelParams->GetForCurve()->GetCurrencyUnit()->GetCcyName());
	return string(FXCal);
}


////////////////////////////////////////////////////
///	Class   : ARM_EqFxBase
///	Routines: ComputeSettlementGap
///	Returns : double
///	Action  : get the settlement gap
////////////////////////////////////////////////////
double ARM_EqFxBase::ComputeSettlementGap(const string& modelName) const
{
	const ARM_ModelParams_Fx* fxModelParams = dynamic_cast<const ARM_ModelParams_Fx*>(GetModelParams());

#if defined(__GP_STRICT_VALIDATION)
	if( !dynamic_cast<const ARM_ModelParams_Fx*>(GetModelParams()) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": model params is not ARM_ModelParams_Fx*" );
#endif

	ARM_ZeroCurvePtr domCurve = fxModelParams->GetDomCurve();
	ARM_ZeroCurvePtr forCurve = fxModelParams->GetForCurve();

	double domGap = fxModelParams->GetDomCurve()->GetCurrencyUnit()->GetSpotDays();
	double forGap = fxModelParams->GetForCurve()->GetCurrencyUnit()->GetSpotDays();
	return CC_Max(domGap,forGap);
}

////////////////////////////////////////////////////
///	Class  : ARM_EqFxBase
///	Routine: DiscountFactor
///	Returns: a vector of Zc(t,T)
///	Action : Closed form formula for DF (currently only forward value!)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_EqFxBase::DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
		const ARM_PricingStatesPtr& states) const
{
    /// Waiting for the access to the yield curve with curveName
    ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
	size_t nbStates = states!=ARM_PricingStatesPtr(NULL)? states->size(): 1;

	double zcT = ZcCurve->DiscountPrice(maturityTime/K_YEAR_LEN);
	double zct = ZcCurve->DiscountPrice(evalTime/K_YEAR_LEN);
	
	return ARM_VectorPtr( new std::vector<double>(nbStates,zcT/zct) );
}

////////////////////////////////////////////////////
///	Class  : ARM_EqFxBase
///	Routine: CallVectorial
///	Returns: a vector of call vectorial
///	Action : computes equity or fx call
////////////////////////////////////////////////////
ARM_VectorPtr ARM_EqFxBase::CallVectorial(
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
	ARM_FXName fxName(modelName);
	bool IsInvMkt = fxName.GetIsInvMkt();
	string rightCcy = modelName.substr(3,3);
	string payCcy	= GetPayModelName();
	if(payCcy == "NoName") payCcy = rightCcy;

	ARM_EqFxBase::CallType calltype;
	if(!IsInvMkt)
		calltype = payCcy==rightCcy ? ARM_EqFxBase::ClosedFormula : ARM_EqFxBase::Quanto;
	else
		calltype = payCcy==rightCcy ? ARM_EqFxBase::InvClosedFormula : ARM_EqFxBase::InvQuanto;

	if(	evalTime<=K_NEW_DOUBLE_TOL || states == ARM_PricingStatesPtr(NULL) )
	{
		ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
		size_t nbStates = states!=ARM_PricingStatesPtr(NULL)? states->size(): 1;
		ARM_GP_VectorPtr fwdVect = Forward(modelName, evalTime, expiryTime, settlementTime, payTime, states );
		double strike = strikePerState[0];
		double fwd = (*fwdVect)[0];
		double zcT	= ZcCurve->DiscountPrice(settlementTime/K_YEAR_LEN);

		if ( (calltype == Quanto) || (calltype == InvQuanto) )
		{
			ARM_CFMethod* cfnumMethod = dynamic_cast<ARM_CFMethod*>(&*GetNumMethod());
			if( !cfnumMethod)
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : For the moment only the numerical method Closed Form is supported");
			ARM_GP_Matrix IntegParameters = cfnumMethod->GetCFParameters();
			//Normal Mapping Fwd_payoff
			std::vector<double> glparams = IntegParameters.empty()? std::vector<double>() : IntegParameters.GetColumns(0);
			ARM_GP_Matrix glmatrix= Forward_Mapping(fwd, expiryTime,glparams);

			//Construction of the ARM_FX_Vanilla
			ARM_FXVanilla2D::FXVanilla2D  fxType  = calltype == Quanto ? ARM_FXVanilla2D::Fx1_Fx2 : ARM_FXVanilla2D::InvFx1_InvFx2;
			ARM_FXCall call(strike,callPut,fxType);
			//Construction of the ARM_GaussReplic
			ARM_GaussReplic1D::QuantoType  quantoType = calltype == Quanto ? ARM_GaussReplic1D::InvQuanto : ARM_GaussReplic1D::Quanto;
			ARM_GaussReplic1D gReplic(&call,glmatrix,quantoType);

			//Expectation
			double price= gReplic.Price();
		
			//normalisation
			const ARM_ModelParams_Fx* fxmodelparams = dynamic_cast< const ARM_ModelParams_Fx* >(GetModelParams());
			double spot = calltype == Quanto ? fxmodelparams->GetSpot():1.0;
			double callprice = zcT/spot*price;
			
			return  ARM_VectorPtr( new std::vector<double>(1.0,callprice) );
		}
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": EqFxBase cannot price this kind of call." );
	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": EqFxBase can be used just for its closed forms." );

		ARM_GP_VectorPtr dumyFwdVector;
		return dumyFwdVector;
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_EqFxBase
///	Routine: Forward
///	Returns: a vector of forward (t,T)
///	Action : computes the forward of the equity or fx model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_EqFxBase::Forward(
		const string& modelName, 
        double evalTime,
	    double expiryTime,
	    double settlementTime,
	    double payTime,
        const ARM_PricingStatesPtr& states) const
{
	double forward	= ComputeFwdAtTime( settlementTime );

	double payLagAdjst = PaymentLagAdjst(expiryTime, settlementTime, payTime);
	
	if(	evalTime<=K_NEW_DOUBLE_TOL || states == ARM_PricingStatesPtr(NULL) )
	{
		size_t payoffSize = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
		return ARM_VectorPtr( new std::vector<double>(payoffSize,forward*payLagAdjst) );
	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Mixture FX can be used just for its closed forms." );

		ARM_GP_VectorPtr dumyFwdVector;
		return dumyFwdVector;
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_EqFxBase
///	Routine: UpdateConvAdjTerm
///	Returns: Nothing
///	Action : Set in fxmodelparams the multiplicative
///		convexity adjustment to do to the FX forward
///     when the settlement is inferior to the payment
////////////////////////////////////////////////////
double ARM_EqFxBase::PaymentLagAdjst(
		double expiryTime,
		double settlementTime,
		double payTime) const
{
	if(!itsConvAdjustModel || fabs(payTime - settlementTime) < ARM_GlobalConstant::ARM_SEVENDAYS_LAG)
		return 1.0;
	//validation
	ARM_PricingModelIR* model = dynamic_cast<  ARM_SABR_Model* >( itsConvAdjustModel );
	if(!model)
	{
		model = dynamic_cast<  ARM_BS_Model* >( itsConvAdjustModel );
		if(!model)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : PaymentLag :Convexity adjustment needs either SABR IR or BS IR, please advise");
	}

	//models params
	const ARM_ModelParam& fxModelParamVol = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility);
	const ARM_ModelParam& irModelParamVol = model->GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility);

	//volatilties
	double fxVol = fxModelParamVol.GetValue(expiryTime);
	double tenor = (payTime - settlementTime)/K_YEAR_LEN;
	double irVol = irModelParamVol.GetValue(expiryTime,tenor);
	ARM_ZeroCurvePtr zcCurve = model->GetZeroCurve();
	double ZcStart	= zcCurve->DiscountPrice(settlementTime/K_YEAR_LEN);
	double ZcEnd	= zcCurve->DiscountPrice(payTime/K_YEAR_LEN);
	double irFwd = (1.0 - ZcEnd/ZcStart);

	//Set  in modelparams;
	double mat = expiryTime/K_YEAR_LEN;
	double rho = itsRho.Interpolate(expiryTime);
	double value = exp(-rho*irFwd*irVol*fxVol*mat);

	return value;
}

////////////////////////////////////////////////////
///	Class  : ARM_EqFxBase
///	Routine: NormalDistribution
///	Returns: ARM_GP_VectorPtr
///	Action : compute the function which describes the
///		fwd as a function of a normal centred 
///     variable			
////////////////////////////////////////////////////

ARM_GP_VectorPtr ARM_EqFxBase::NormalDistribution(GaussLegendre_Coefficients glc, double forward, double expiryTime) const
{
	const_cast <ARM_EqFxBase*>(this)->UpdateDensityFunctor(forward,expiryTime);
	int NbPoints = glc.get_order();
	ARM_GP_VectorPtr normaldist(new std::vector<double>(NbPoints,0.0));
	double xi, proba_i;
	for( int i=0; i<NbPoints; ++i)
	{
		xi = glc.get_point(i);
		proba_i = ARM_GaussianAnalytics::cdfNormal(xi); 
		(*normaldist)[i] = GetDensityFunctor()->Quantile(proba_i, forward, expiryTime/K_YEAR_LEN);
	}
	return normaldist;
}

////////////////////////////////////////////////////
///	Class  : ARM_EqFxBase
///	Routine: Black_Forward_Mapping
///	Returns: ARM_VectorPtr
///	Action : Compute the closed form of the Normal
/// mapping of the fwd or its inverse in the Black world
///		
//////////////////////////////////////////////////// 

ARM_VectorPtr ARM_EqFxBase::Black_Forward_Mapping(
		double	fwd,
		double	vol,
		double	expiryTime,
		GaussLegendre_Coefficients glc,
		bool IsInv) const
{
	/// init of call vectorial
	int NbPoints = glc.get_order();
	ARM_GP_VectorPtr fwd_mat(new std::vector<double>(NbPoints,0.0));
	double xi;
	double mat = expiryTime/K_YEAR_LEN;
	if(IsInv)
	{
		for(int i=0;i<NbPoints;i++)
		{
			xi = glc.get_point(i);
			(*fwd_mat)[i]=1.0/fwd*exp(+0.5*vol*vol*mat+sqrt(mat)*vol*xi);
		}
	}
	else
	{
		for(int i=0;i<NbPoints;i++)
		{
			xi = glc.get_point(i);
			(*fwd_mat)[i]=fwd*exp(-0.5*vol*vol*mat+sqrt(mat)*vol*xi);
		}
	}
	return fwd_mat;
}


////////////////////////////////////////////////////
///	Class  : ARM_EqFxBase
///	Routine: Forward_Mapping
///	Returns: ARM_VectorPtr
///	Action : Compute the Normal
/// mapping of the fwd or its inverse with complete Smile 
///		or on a flatwings Smile
//////////////////////////////////////////////////// 

ARM_GP_Matrix ARM_EqFxBase::Forward_Mapping(
		double	fwd,
		double	expiryTime,
		const std::vector<double>& glparams) const
{

	//xmin and xmax
	double x_inf=-6.0;
	double x_sup=+6.0;
	double xmin = -6.0;
	double xmax = +6.0;

	double Kmin= 1.0e-4;
	double Kmax= 1.0e4;
	int NbPoints= 121;

	if( !glparams.empty())
	{
		Kmin=glparams[0];
		Kmax=glparams[1];
		NbPoints = glparams[2];
	}
	
	if ( Kmin > 0.001 )
		xmax=ARM_GaussianAnalytics::cdfNormal_Inv(GetDensityFunctor()->Proba(Kmax,fwd,expiryTime) ); 
	if ( Kmax<1000 )
		xmin=ARM_GaussianAnalytics::cdfNormal_Inv( GetDensityFunctor()->Proba(Kmin,fwd, expiryTime) );
	if(xmax < xmin) CC_NS(std,swap) (xmax, xmin);
	xmax = CC_Min (xmax, x_sup);
	xmin = CC_Max (xmin, x_inf);

	//Fwd mappingy
	GaussLegendre_Coefficients glc( NbPoints, xmin, xmax);
	ARM_GP_VectorPtr normal_fwd_at_mat = NormalDistribution(glc, fwd, expiryTime);

	//fill
	std::vector<double> x_w_fwd= glc.get_points();
	x_w_fwd.fill( glc.get_weights() );
	x_w_fwd.fill( *normal_fwd_at_mat );
	//return

	return ARM_GP_Matrix(3,NbPoints,x_w_fwd);;
}


////////////////////////////////////////////////////
///	Class  : ARM_EqFxBase
///	Routine: DigitalVectorial
///	Returns: a vector of digital vectorial
///	Action : computes equity or fx digital
////////////////////////////////////////////////////
ARM_VectorPtr ARM_EqFxBase::DigitalVectorial(
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
	ARM_GP_VectorPtr callDown, callUp;
	size_t payoffSize = states!=ARM_PricingStatesPtr(NULL)? states->size(): 1;
	ARM_GP_VectorPtr digit(new std::vector<double>(payoffSize,0.0));

	double signed_eps = callPut*epsilon;
	double norm = 1/epsilon;
	std::vector<double> strikeDownPerState, strikeUpPerState;

	switch (digitType)
	{
	case ARM_FXDigitType::backward:
		{
			strikeDownPerState = strikePerState - signed_eps;
			strikeUpPerState = strikePerState;
			break;
		}
	case ARM_FXDigitType::forward:
		{
			strikeDownPerState = strikePerState;
			strikeUpPerState = strikePerState + signed_eps;
			break;
		}
	case ARM_FXDigitType::centred:
		{
			norm = 1/(2*epsilon);
			strikeDownPerState = strikePerState - signed_eps;
			strikeUpPerState = strikePerState + signed_eps;
			break;
		}
	case ARM_FXDigitType::analytic:
		{
			ARM_FXName fxName(modelName);
			bool IsInvMkt = fxName.GetIsInvMkt();
			string rightCcy = modelName.substr(3,3);
			string payCcy	= GetPayModelName();
			if(payCcy == "NoName") payCcy = rightCcy;
		
			ARM_EqFxBase::CallType calltype;
			if(!IsInvMkt)
				calltype = payCcy==rightCcy ? ARM_EqFxBase::ClosedFormula : ARM_EqFxBase::Quanto;
			else
				calltype = payCcy==rightCcy ? ARM_EqFxBase::InvClosedFormula : ARM_EqFxBase::InvQuanto;
			if ( ( calltype == InvQuanto ) || ( calltype == ClosedFormula) )
			{
				ARM_THROW( ERR_INVALID_ARGUMENT, "Analytic DigitalVectorial unimplemented function for ARM_EqFxBase model!");
			}
			else
			{
				const ARM_ModelParams_Fx* fxmodelparams = dynamic_cast< const ARM_ModelParams_Fx* >(GetModelParams());
				double spot = fxmodelparams->GetSpot();
				double newstrike = 1.0/strikePerState[0];
				std::vector<double> newstrikePerState = std::vector<double>(strikePerState.size(),newstrike);
				int newcallPut = -callPut;
				string mktName = fxName.GetMktName();
				const_cast <ARM_EqFxBase*>(this)->SetPayModelName(mktName.substr(3,3));//because DigitCall(DOMFOR,K) in FOR = 1/FWD*(Put(FORDOM,1/K) in DOM + 1/K*DigitPut(FORDOM,1/K) in DOM)
				//To avoid convexification of the forward, we call the digital and the call with payTime = settlemenTime
				callDown = CallVectorial(mktName, evalTime, expiryTime, settlementTime, newstrikePerState, newcallPut, settlementTime, states, context);
				callUp = DigitalVectorial(mktName, evalTime, expiryTime, settlementTime, newstrikePerState, 1.0, newcallPut, settlementTime, digitType, epsilon, states, context);
				*digit = *callDown + newcallPut*newstrike*(*callUp);
				*digit = newcallPut*1.0/spot*notional*(*digit);
				double forZcTpay = fxmodelparams->GetForCurve()->DiscountPrice(payTime/K_YEAR_LEN);
				double forZcTset = fxmodelparams->GetForCurve()->DiscountPrice(settlementTime/K_YEAR_LEN);
				*digit = forZcTpay/forZcTset*(*digit);
				const_cast <ARM_EqFxBase*>(this)->SetPayModelName(payCcy);//we have to reset the payment currency to its initial value
				return digit;
			}
			break;
		}
	}
	callDown = CallVectorial(modelName, evalTime, expiryTime, settlementTime, strikeDownPerState, callPut, payTime, states, context);
	callUp   = CallVectorial(modelName, evalTime, expiryTime, settlementTime, strikeUpPerState, callPut, payTime, states, context);
	*digit = norm*notional*(*callDown-*callUp);
	return digit;
}

////////////////////////////////////////////////////
///	Class   : ARM_EqFxBase
///	Routines: VolatilitiesAndCorrelations
///	Returns :
///	Action  : computes the volatilities its derivatives and the correlation
////////////////////////////////////////////////////
void ARM_EqFxBase::VolatilitiesAndCorrelations( const std::vector<double>& timeSteps, 
			ARM_GP_MatrixPtr& vols,
			ARM_GP_MatrixPtr& d1Vols,
			ARM_GP_MatrixPtr& correls,
			bool linearVol) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, " unimplemented function for ARM_EqFxBase model!");
}



////////////////////////////////////////////////////
///	Class   : ARM_EqFxBase
///	Routines: EulerLocalDrifts
///	Returns :
///	Action  : computes the relative and absolute drift
////////////////////////////////////////////////////

void ARM_EqFxBase::EulerLocalDrifts( const std::vector<double>& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "EulerLocalDrifts is not implemented for ARM_EqFxBase model!");
}


////////////////////////////////////////////////////
///	Class   : ARM_EqFxBase
///	Routines: MarkovianDrift
///	Returns : void
///	Action  : Computes the Markovian drift
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_EqFxBase::MarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates ) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "MarkovianDrift is not implemented for ARM_EqFxBase model!");
}


////////////////////////////////////////////////////
///	Class   : ARM_EqFxBase
///	Routines: LocalDiscounts
///	Returns : void
///	Action  : Computes the LocalDiscounts
////////////////////////////////////////////////////
ARM_VectorPtr ARM_EqFxBase::LocalDiscounts(
	size_t timeIdx, 
	double dt, 
	const ARM_PricingStatesPtr& states) const
{
		ARM_THROW( ERR_INVALID_ARGUMENT, " unimplemented function LocalDiscounts for ARM_EqFxBase model!");
}



////////////////////////////////////////////////////
///	Class  : ARM_EqFxBase
///	Routine: LocalDrifts
///	Returns: void
///	Action : local drifts
////////////////////////////////////////////////////
void ARM_EqFxBase::IntegratedLocalDrifts(
	const std::vector<double>& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, " unimplemented function IntegratedLocalDrifts for ARM_EqFxBase model!");
}


////////////////////////////////////////////////////
///	Class  : ARM_EqFxBase
///	Routine: ModelStateLocalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local variances
///          of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_EqFxBase::ModelStateLocalVariances(
    const std::vector<double>& timeSteps,
    ARM_MatrixVector& localVariances ) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, " unimplemented function IntegratedLocalDrifts for ModelStateLocalVariances model!");
}


////////////////////////////////////////////////////
///	Class  : ARM_EqFxBase
///	Routine: NumMethodStateLocalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local and global variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_EqFxBase::NumMethodStateLocalVariances(
    const std::vector<double>& timeSteps,
    ARM_MatrixVector& localVariances ) const
{
	ModelStateLocalVariances( timeSteps, localVariances );
}


////////////////////////////////////////////////////
///	Class  : ARM_EqFxBase
///	Routine: NumMethodStateGlobalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute global variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_EqFxBase::NumMethodStateGlobalVariances(
    const std::vector<double>& timeSteps,
    ARM_MatrixVector& variances ) const
{
	size_t nbSteps	= timeSteps.size();
	size_t modelNb	= GetModelNb();
    double step		= timeSteps[0],
		   nextStep;
	size_t offsetIndex1 = (nbSteps-1)*modelNb;
	size_t offsetIndex2 = nbSteps*modelNb;
#if defined(__GP_STRICT_VALIDATION)
	if( variances.size()!= offsetIndex2 ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "localDrifts.size() != offsetIndex" );
#endif

	variances.resize(nbSteps*(modelNb+1));

    variances[offsetIndex2+0]=new ARM_GP_TriangularMatrix;
    for(size_t i=0;i<nbSteps-1;++i)
    {
        nextStep=timeSteps[i+1];
        
		/// [i+1] => variance from 0 -> ti+1
        /// we can't sum up local variance !
        variances[offsetIndex2+i+1] = new ARM_GP_TriangularMatrix(1,((const ARM_ModelParamsHW1FStd* const) GetModelParams())->StateLocalVariance(0.0,nextStep,nextStep));
        step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_EqFxBase
///	Routine: ComputeFwdAtTime
///	Returns: double
///	Action : 
////////////////////////////////////////////////////
double ARM_EqFxBase::ComputeFwdAtTime( double evalTime  ) const
{
	const ARM_ModelParams* params = GetModelParams();
	const ARM_ModelParams_Fx* eqFxFunc = dynamic_cast<const ARM_ModelParams_Fx*>( params );

#if defined(__GP_STRICT_VALIDATION)
	if(  !eqFxFunc )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": eqFxFunc is a null pointor" );
#endif

	return eqFxFunc->Forward( evalTime );
}


////////////////////////////////////////////////////
///	Class   : ARM_EqFxBase
///	Routine : ValidateCalibMethod
///	Returns :void
///	Action  : call DefaultValidateWithModel
////////////////////////////////////////////////////
void ARM_EqFxBase::ValidateCalibMethod(ARM_CalibMethod& calibMethod)
{
	calibMethod.DefaultValidateWithModel(*this);
}


////////////////////////////////////////////////////
///	Class   : ARM_EqFxBase
///	Routines: VanillaSpreadOptionLet
///	Returns : void
///	Action  : No default implementation
////////////////////////////////////////////////////

ARM_VectorPtr  ARM_EqFxBase::VanillaSpreadOptionLet(const string& curveName,
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
		"VanillaSpreadOption : unimplemented function for ARM_EqFxBase Model!");
}

////////////////////////////////////////////////////
///	Class  : ARM_EqFxBase
///	Routine: RangeAccrualVectorial
///	Returns: a vector of rangeAccrual vectorial
///	Action : computes rangeAccrual ie 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_EqFxBase::RangeAccrualVectorial(
		const string& modelName,
		double evalTime,
		double startTime,
		double endTime,
		const std::vector<double>& fixingTimes,
		double payTime,
		const std::vector<double>& downBarrierVect,
		const std::vector<double>& upBarrierVect,
		const std::vector<double>& notionalVect,
		const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context) const
{	
	/// init of rangeAccrual vectorial
	size_t payoffSize = states!=ARM_PricingStatesPtr(NULL)? states->size(): 1;
	
	//density
	string payCcy = GetPayModelName();
	string domMarketFXName = modelName.substr(3,3); 
	string forMarketFXName = modelName.substr(0,3); 
	bool direct; 
	double norm = GetZeroCurve()->DiscountPrice(payTime/K_YEAR_LEN);
	double downBarrier = downBarrierVect[0];
	double upBarrier   = upBarrierVect[0];
	int callPut = 1;
	if(domMarketFXName==payCcy)
	{
		direct = true;
	}
	else if (forMarketFXName==payCcy)
	{
		direct = false;
		const ARM_ModelParams_Fx* fxmodelparams = dynamic_cast< const ARM_ModelParams_Fx* >(GetModelParams());
		double spot = fxmodelparams->GetSpot();
		ARM_GP_VectorPtr fwdVect= Forward(modelName, evalTime, payTime, payTime, payTime, states );
		double fwd = (*fwdVect)[0];
		norm *= fwd/spot;
		downBarrier = 1.0/downBarrier;
		upBarrier   = 1.0/upBarrier;
		callPut =-1;
	}
	else 
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The pure rate quanto range accrual is not ready yet");
	}
	ARM_DensityFunctor* fxdensity = dynamic_cast<ARM_DensityFunctor*>(GetDensityFunctor());
	fxdensity->SetIsDirect(direct);

	//commom arguments
	double expiryTime;  
	double notional;

	//range accrual coupon calculation
	double cpnDown, cpnUp;
	int nbFixing = fixingTimes.size();
	double rangeAccrual = 0.0;
	
	for( int k = 0; k < nbFixing; ++k)
	{
		expiryTime = fixingTimes[k];
		notional = notionalVect[0];
		ARM_GP_VectorPtr fwdVect= Forward(modelName, evalTime, expiryTime, expiryTime, expiryTime, states );
		double fwd = (*fwdVect)[0];
		const_cast <ARM_EqFxBase*>(this)->UpdateDensityFunctor(fwd,expiryTime);
		cpnDown	= fxdensity->Proba(downBarrier,fwd,expiryTime/K_YEAR_LEN);
		cpnUp   = fxdensity->Proba(upBarrier,fwd,expiryTime/K_YEAR_LEN);
		rangeAccrual += callPut*(cpnUp - cpnDown);
	}
	rangeAccrual *= notional;
	rangeAccrual /= nbFixing;
	//discounting
	rangeAccrual*=norm;
	return ARM_GP_VectorPtr( new std::vector<double>(payoffSize,rangeAccrual) );
}

////////////////////////////////////////////////////
///	Class  : ARM_EqFxBase
///	Routine: TreeStatesToModelStates
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_EqFxBase::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{
#if defined(__GP_STRICT_VALIDATION)
	if(timeIndex<0)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index is negative!");
	if( timeIndex > GetNumMethod()->GetTimeSteps()->size()-1 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index bigger than max size!");
#endif

	ARM_NumerairePtr numeraire= GetNumeraire();
    if( numeraire->GetType() == ARM_Numeraire::Cash )
		numeraire->Update( *this, states, timeIndex );

	const ARM_GP_MatrixPtr& numMethodStates = states->GetNumMethodStates();
	states->SetModelStates(numMethodStates);
}


////////////////////////////////////////////////////
///	Class  : ARM_EqFXBase
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_EqFxBase::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states, int timeIndex) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "TreeStatesToModelStates : unimplemented function for ARM_EqFXBase Model!");
}


////////////////////////////////////////////////////
///	Class  : ARM_EqFxBase
///	Routine: FirstPricingStates
///	Returns: ARM_PricingStatesPtr
///	Action  : create the first pricing state
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_EqFxBase::FirstPricingStates( size_t bucketSize ) const
{
	/// ARM_PricingStates(nbStates = bucketSize, nbModelStates = 1F , nbPayoffs = 0)
// FIXMEFRED: mig.vc8 (25/05/2007 15:19:11):cast
	return static_cast<ARM_PricingStatesPtr>(new ARM_PricingStates(bucketSize,1,0,1));
}

////////////////////////////////////////////////////
///	Class  : ARM_EqFxBase
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Default initialisation of the model and the
///          associated numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_EqFxBase::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "Init : unimplemented function for ARM_EqFXBase Model!");

    return ARM_PricingStatesPtr(NULL);
}

////////////////////////////////////////////////////
///	Class  : ARM_EqFxBase
///	Routine: AdviseBreakPointTimes
///	Returns: void 
///	Action : sets the corresponding suggested break point times to the model param
////////////////////////////////////////////////////
void ARM_EqFxBase::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, 
							  ARM_ModelParam* inputModelParam, 
							  size_t factorNb )
{
	/// in any case, use the reset date of the portfolio (use a reference to avoid creating a copy, hence the change are kept in the input model Param!)
	ARM_CurveModelParam& modelParam = inputModelParam->ToCurveModelParam();

    double asOfDate			= GetAsOfDate().GetJulian();
    size_t portfolioSize	= portfolio->size();  
    std::vector<double>  tmpdates;
    size_t i;
	
	/// can only accept option security
	for( i=0; i<portfolioSize; ++i )
		if( portfolio->GetAsset(i)->GetName() != ARM_OPTION )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": portfolio can only contain option!" );

	double date = portfolio->GetAsset(0)->GetMaturity().GetJulian()- asOfDate;
    tmpdates.push_back(date);
    for(i=1; i<portfolioSize; i++) 
    {
        double resetlag = portfolio->GetAsset(i)->GetMaturity().GetJulian() - asOfDate;
        if(fabs (date - resetlag) > FRMVOL_LAG_THRESHOLD)
        {
            tmpdates.push_back(resetlag);
            date = resetlag;
        }
    }
	modelParam.UpdateValues(&tmpdates);
}

////////////////////////////////////////////////////
///	Class  : ARM_EqFxBase
///	Routine: ImpliedVol
///	Returns: int
///	Action : tells the type of the model
////////////////////////////////////////////////////
double ARM_EqFxBase::ImpliedVol(const ARM_VanillaArg& arg) const
{
    return DefaultImpliedVol(arg);
}

////////////////////////////////////////////////////
///	Class   : ARM_EqFxBase
///	Routine : ValidateModelParams
///	Returns : bool
///	Action  : validate model params
////////////////////////////////////////////////////
bool ARM_EqFxBase::ValidateModelParams(const ARM_ModelParams& params) const
{
	if( !dynamic_cast<const ARM_ModelParams_Fx*>(&params))
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelparams is not of type ARM_ModelParams_Fx" );
	return true;
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

