/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file BS_Model.cpp
 *
 *  \brief base class for Black Scholes model
 *	\author  E. Benhamou, O. Croissant
 *	\version 1.0
 *	\date October 2004
 */

/// remove identified warning
#include "gpbase/removeidentifiedwarning.h"

#include "gpmodels/BS_Model.h"

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelparamtype.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/surfacelistmodelparam.h"

/// gpbase
#include "gpbase/surface.h"
#include "gpbase/datestrip.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_bs.h"

/// gpmodels
#include "gpmodels/BS_ModelParams.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_BS_Model
///	Routine: Constructor
///	Returns: 
///	Action : builds the ARM_BS_Model
////////////////////////////////////////////////////

ARM_BS_Model::ARM_BS_Model(const ARM_ZeroCurvePtr& zc, const ARM_BS_ModelParams& params)
:	ARM_AnalyticIRModel( zc, params )
{}


////////////////////////////////////////////////////
///	Class   : ARM_BS_Model
///	Routines: Copy constructor
///	Returns :
///	Action  : 
////////////////////////////////////////////////////

ARM_BS_Model::ARM_BS_Model(const ARM_BS_Model& rhs)
:	ARM_AnalyticIRModel( rhs )
{}


////////////////////////////////////////////////////
///	Class  : ARM_BS_Model
///	Routine: operator=
///	Returns: 
///	Action : Assignment operator
////////////////////////////////////////////////////
ARM_BS_Model& ARM_BS_Model::operator=(const ARM_BS_Model& rhs)
{
	if( this != &rhs )
	{
		ARM_AnalyticIRModel::operator =(rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_BS_Model
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_BS_Model::~ARM_BS_Model()
{}


////////////////////////////////////////////////////
///	Class  : ARM_BS_Model
///	Routine: VanillaCaplet
///	Returns: ARM_VectorPtr
///	Action : computes the price of a caplet
////////////////////////////////////////////////////
ARM_VectorPtr ARM_BS_Model::VanillaCaplet(
	const string& curveName, 
	double evalTime,
	double payTime, 
	double period,
    double payNotional,
	double fwdResetTime, 
	double fwdStartTime,
    double fwdEndTime,
	double fwdPeriod,
	const std::vector<double>& strikesPerState,
    int capFloor,
	const ARM_PricingStatesPtr& states) const
{
	
	double tenor = ARM_AnalyticIRModel::GetMatchingTenor((fwdEndTime-fwdStartTime)/K_YEAR_LEN);

	if (tenor <= 1.0)
		return LiborCaplet(curveName, evalTime,payTime, period,payNotional,fwdResetTime,fwdStartTime,
							fwdEndTime,fwdPeriod,strikesPerState,capFloor,states);
	else
		return CMSCaplet(curveName, evalTime,payTime, period,payNotional,fwdResetTime,fwdStartTime,
							fwdEndTime,fwdPeriod, 0., 1, strikesPerState,capFloor,states);
}


ARM_VectorPtr ARM_BS_Model::LiborCaplet(
	const string& curveName, 
	double evalTime,
	double payTime, 
	double period,
    double payNotional,
	double fwdResetTime, 
	double fwdStartTime,
    double fwdEndTime,
	double fwdPeriod,
	const std::vector<double>& strikesPerState,
    int capFloor,
	const ARM_PricingStatesPtr& states) const
{
	ARM_VectorPtr liborRate	= Libor(curveName, evalTime,fwdStartTime, 
		fwdEndTime, period, fwdResetTime, payTime, states );
	ARM_VectorPtr DF		= DiscountFactor( curveName, evalTime, payTime, states );
	const ARM_SurfaceListModelParam* SmileModelParam;
	
	double time				= fwdResetTime-evalTime;
	double tenor			= ARM_AnalyticIRModel::GetMatchingTenor((fwdEndTime-fwdStartTime)/K_YEAR_LEN);
	double volatility		= GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(time,tenor);
	
	if (GetModelParams()->DoesModelParamExist(ARM_ModelParamType::Smile))
	{
		SmileModelParam = dynamic_cast < const ARM_SurfaceListModelParam* > (&(GetModelParams()->GetModelParam(ARM_ModelParamType::Smile)));
		
		if (SmileModelParam)
		{
			double K_F = strikesPerState[0] - (*liborRate)[0];
			double smile = SmileModelParam->GetValue(time,tenor, K_F);
			volatility = volatility + smile;
		}
	}

	double optionValue = payNotional * period * BlackSholes_Formula(
		(*liborRate)[0],volatility*sqrt(time/K_YEAR_LEN), (*DF)[0],strikesPerState[0], capFloor);

	return ARM_VectorPtr( new std::vector<double>(1,optionValue));
}



////////////////////////////////////////////////////
///	Class  : ARM_BS_Model
///	Routine: CMSCaplet
///	Returns: ARM_VectorPtr
///	Action : computes the price of a CMS caplet
////////////////////////////////////////////////////
ARM_VectorPtr ARM_BS_Model::CMSCaplet(
	const string& curveName, 
	double evalTime,
	double payTime,
	double period,
    double payNotional,
	double fwdResetTime, 
	double fwdStartTime,
    double fwdEndTime,
	double fwdPeriod,
	double spread,
	int decapFreq,
	const std::vector<double>& strikesPerState,
    int capFloor,
	const ARM_PricingStatesPtr& states) const
{
	/// not necesary to use  fwdStartTimes, fwdEndTimes, fwdPayPeriods, floatPayTimes, floatPayPeriods
	/// as we are in the case of a vanilla swap!
	std::vector<double> fixPayTimes, fixPayPeriods;
	std::vector<double> dummyFwdStartTimes, dummyFwdEndTimes, dummyFwdPayPeriods, dummyFwdPayTimes, dummyFloatPayTimes, dummyFloatPayPeriods;

	ARM_Currency* ccy			= GetCurrency( curveName );
	char* ccyName				= ccy->GetCcyName();
	
	ARM_Date asofDate = GetAsOfDate();
	double time				= fwdResetTime-evalTime;
	double tenor			= ARM_AnalyticIRModel::GetMatchingTenor((fwdEndTime-fwdStartTime)/K_YEAR_LEN);

	char fixCalendar[100];
    ccy->CalcFixPayCal(fixCalendar);

	/// the second argument can either be a date or a maturity
	ARM_Date startDate = GetAsOfDate();
	startDate.AddDays(fwdStartTime);
	ARM_Date endDate = startDate;
	endDate.AddYears(tenor);

    /// Fixed leg data extraction
	/// get the fixed leg convention
	int fixFreq	    = ccy->GetFixedPayFreq();
	int fixDayCount = ccy->GetFixedDayCount();

	/// fix leg datestrip
	ARM_DateStrip fixDateStrip( startDate, endDate, fixFreq, fixDayCount, fixCalendar,
								K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE,
								fixCalendar );

	/// not cloned hence no need to delete this!
	std::vector<double>& pfixPayTimes   = fixDateStrip.GetPaymentDates();
	std::vector<double>& pfixPayPeriods = fixDateStrip.GetInterestTerms();

	

	/// copy constructor
	fixPayTimes	= std::vector<double>(*pfixPayTimes);
	fixPayPeriods  = std::vector<double>(*pfixPayPeriods);
	std::vector<double> margin = std::vector<double>(1,0.0);
	/// get time from date! 
	for(int i=0; i<pfixPayTimes->size(); ++i )
	{
		fixPayTimes[i] = GetTimeFromDate(fixPayTimes[i]);
	}

	// overload fwdEndTime from adjusted endDate
	endDate.AdjustToBusDate(fixCalendar, K_MOD_FOLLOWING);
	fwdEndTime = endDate - GetAsOfDate();



	ARM_VectorPtr swapRate = SwapRate(
		curveName, 
		evalTime,
		fwdStartTime, 
		fwdEndTime, 
		fixPayTimes,
		fixPayPeriods,
		dummyFwdStartTimes,
        dummyFwdEndTimes,
        dummyFwdPayTimes,
        dummyFloatPayTimes,
        dummyFloatPayPeriods,
		margin,	/// margin
		true,	/// isDbleNotional to avoid computing the float cash flows piece by piece
		ARM_PricingStatesPtr(NULL) );

	double rawSwap = (*swapRate)[0];
	double strike = strikesPerState[0];


	ARM_VectorPtr DF		= DiscountFactor( curveName, evalTime, payTime, states );
	const ARM_SurfaceListModelParam* SmileModelParam;
	const ARM_SurfaceModelParam* CorrelModelParam;
	
	
	double volidx			= GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(time,tenor);
	double lag				= (payTime - fwdResetTime)/K_YEAR_LEN;
	double voldisc			= GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(time,lag);
	double volatility		= volidx;

	double cvxadj = CptCvxAdj(rawSwap, tenor, time/K_YEAR_LEN, volidx, fixFreq);
	
	double swapRateAdj = rawSwap * cvxadj;
	

	if (GetModelParams()->DoesModelParamExist(ARM_ModelParamType::Correlation))
	{
		CorrelModelParam = dynamic_cast < const ARM_SurfaceModelParam* > (&(GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation)));
		
		if (CorrelModelParam)
		{
			double theta = (payTime - fwdResetTime)/360;

			ARM_VectorPtr liborRate	= Libor(curveName, evalTime, fwdResetTime, 
									payTime, theta, fwdResetTime, payTime, states );

			double R0 = (*liborRate)[0];

			double Rho = CorrelModelParam->GetValue(lag,tenor);

			double lagadjmul = CptTimeLagAdj(time/K_YEAR_LEN, lag, volidx, R0, voldisc, Rho);

			swapRateAdj = swapRateAdj * lagadjmul;
		}
	}

	double indexAdj = swapRateAdj;
	double strikeAdj = 0.;
	double K_F = 0.;


	if (decapFreq ==1)
	{
		strikeAdj = strike; // instead of (strike - spread) because strike has already been adjusted from swapleg
		K_F = strikeAdj - swapRateAdj;
	}
	else
	{
		strike = strike + spread; // to find the original strike because strike has already been adjusted from swapleg
		double idxAdjDecapSpread = ForwardDecap(indexAdj, spread, volatility, decapFreq, time/K_YEAR_LEN);
		double idxAdjDecapNoSpread = ForwardDecap(indexAdj, 0., volatility, decapFreq, time/K_YEAR_LEN);
		strikeAdj = DecompRate(strike, decapFreq);
		strikeAdj = strikeAdj - idxAdjDecapSpread + idxAdjDecapNoSpread;
		double strikeAdjNodecap = InverseDecompRate(strikeAdj, decapFreq);
		K_F = strikeAdjNodecap - swapRateAdj;
		indexAdj = idxAdjDecapNoSpread;	
	}

		

	if (GetModelParams()->DoesModelParamExist(ARM_ModelParamType::Smile))
	{
		SmileModelParam = dynamic_cast < const ARM_SurfaceListModelParam* > (&(GetModelParams()->GetModelParam(ARM_ModelParamType::Smile)));
		
		if (SmileModelParam)
		{
			double smile = SmileModelParam->GetValue(time,tenor, K_F);
			
			volatility = volatility + smile;
		}
	}

	if (decapFreq >1)
	{
		double voldecap = VolDecap(swapRateAdj, volatility, decapFreq);
		volatility = voldecap;
	}


	double optionValue = payNotional * period * BlackSholes_Formula(
		indexAdj, volatility*sqrt(time/K_YEAR_LEN), (*DF)[0],strikeAdj, capFloor);

	return ARM_VectorPtr( new std::vector<double>(1,optionValue));
}



////////////////////////////////////////////////////
///	Class   : ARM_BS_Model
///	Routines: VanillaSwaption
///	Returns :
///	Action  : computes a vanilla swaption
////////////////////////////////////////////////////

ARM_VectorPtr ARM_BS_Model::VanillaSwaption(
	const string& curveName,
	double evalTime,
	double swapResetTime,
	const std::vector<double>& fixNotional,
	const std::vector<double>& floatNotional,
	double floatStartTime,
	double floatEndTime,
	const std::vector<double>& floatResetTimes,
	const std::vector<double>& floatStartTimes,
	const std::vector<double>& floatEndTimes,
	const std::vector<double>& floatIntTerms,	
	const std::vector<double>& fixPayTimes,
	const std::vector<double>& fixPayPeriods,
	const ARM_GP_Matrix& strikesPerState,
    int callPut,
	const ARM_PricingStatesPtr& states,
	bool isConstantNotional,
	bool isConstantSpread,
	bool isConstantStrike) const
{
	/// TO BE UPDATED
	/// Check that the notional is constant
	double swapNotional = fixNotional[0];
	if (!(isConstantNotional&&isConstantSpread&&isConstantStrike))
				ARM_THROW( ERR_INVALID_ARGUMENT, "The Model can not price a swaption with variable notional, Spread or Strike!" );
			

	/// not necesary to use  fwdStartTimes, fwdEndTimes, fwdPayPeriods, floatPayTimes, floatPayPeriods
	/// as we are in the case of a vanilla swap!
	std::vector<double> dummyFwdStartTimes, dummyFwdEndTimes, dummyFwdPayPeriods, dummyFwdPayTimes, dummyFloatPayTimes, dummyFloatPayPeriods;


	std::vector<double> margin = std::vector<double>(1,0.0);
	ARM_VectorPtr swapRate = SwapRate(
		curveName, 
		evalTime,
		floatStartTime, 
		floatEndTime, 
		fixPayTimes,
		fixPayPeriods,
		dummyFwdStartTimes,
        dummyFwdEndTimes,
        dummyFwdPayTimes,
        dummyFloatPayTimes,
        dummyFloatPayPeriods,
		margin,	/// margin
		true,	/// isDbleNotional to avoid computing the float cash flows piece by piece
		ARM_PricingStatesPtr(NULL) );

	ARM_VectorPtr annuity = Annuity(
		curveName, 
        evalTime,
		fixPayTimes,
        fixPayPeriods,
		ARM_PricingStatesPtr(NULL) );

	double time			= swapResetTime-evalTime;
	double tenor		= ARM_AnalyticIRModel::GetMatchingTenor((floatEndTime-floatStartTime)/K_YEAR_LEN );
	double volatility		= GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(time,tenor);

	double optionValue = swapNotional * BlackSholes_Formula(
		(*swapRate)[0],volatility*sqrt(time/K_YEAR_LEN),(*annuity)[0],strikesPerState(0,0), callPut);

	return ARM_VectorPtr( new std::vector<double>(1,optionValue));
}


    
////////////////////////////////////////////////////
///	Class   : ARM_BS_Model
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_BS_Model::Clone() const
{
	return new ARM_BS_Model(*this);
}

////////////////////////////////////////////////////
///	Class   : ARM_BS_Model
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_BS_Model::toString(const string& indent,const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Black Scholes Model \n";
    os << indent << "---------------------\n";
    os << ARM_PricingModel::toString(indent);
    return os.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_BS_Model
///	Routine : ValidateModelParams
///	Returns : bool
///	Action  : validate model params
////////////////////////////////////////////////////
bool ARM_BS_Model::ValidateModelParams(const ARM_ModelParams& params) const
{
	const ARM_BS_ModelParams* BS_ModelParams = dynamic_cast<const ARM_BS_ModelParams*>(&params);
	if( !BS_ModelParams )
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelparams is not of type ARM_BS_ModelParams" );
	return true;
}


////////////////////////////////////////////////////
///	Class   : ARM_BS_Model
///	Routine : CptCvxAdj
///	Returns : double
///	Action  : compute natural adj like Summit
////////////////////////////////////////////////////
double ARM_BS_Model::CptCvxAdj(double cmRate, int tenor, double expiry, double IndexVol, int payFreq) const
{
	double t = 1./payFreq;
	double duration  = cDuration( cmRate, tenor, payFreq, t );
	double convexity = cConvexity( cmRate, tenor, payFreq, t );

	// Adjustment with Exp
	double adjmul = exp(IndexVol*IndexVol*expiry)-1.;
	adjmul = 1.+adjmul*cmRate*convexity/duration/2.;

	return(adjmul);
}


////////////////////////////////////
///	Class   : ARM_BS_Model
///	Routine : CptTimeLagAdj
///	Returns : double
///	Action  : compute Time Lag adj
////////////////////////////////////
double ARM_BS_Model::CptTimeLagAdj(double expiry, double lag, double IndexVol, double Forward, double DiscountVol, double Rho, int methUsed) const
{	
	double adjmul = 1.;

	switch (methUsed)
	{
		case 1: // Multi linear used in Summit analytic formulae
		{
			adjmul = Rho * sqrt(exp(pow( lag * DiscountVol * Forward * sqrt(expiry) /
					(1 + lag * Forward ),2)) - 1) *
					sqrt(exp(pow(IndexVol,2)*expiry)-1);

			adjmul = 1. - adjmul;
			break;
		}
		case 2 : // Multi exponential used in kernel analytic formulae
		{
			adjmul = Rho * lag * Forward * IndexVol * DiscountVol * expiry;
			adjmul = adjmul / (1. + lag*Forward);
			adjmul = exp(-adjmul);
			break;
		}
	}

	return(adjmul);
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

