/*!
 *
 * Copyright(c) CDC IXIS CM July 2003 Paris
 *
 *	\file Heston_Fx.cpp
 *
 *  \brief Mixture FX version
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date August 2006
 */

// nag
#include "nag.h"
#include "nage04.h"

/// gpbase
#include "gpbase/port.h"
/// gpmodel
#include "gpmodels/Mixture_Fx.h"
#include "gpmodels/fxname.h"

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelparams.h"

// gpnumlib
#include "gpnumlib/newtonraphson.h"
#include "gpnumlib/dichotomy.h"
#include "gpnumlib/numfunction.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/inverse.h"
#include "gpclosedforms/vanille_bs_formula.h"

/// gp calib
#include "gpcalib/vanillaeqfxoption.h"

/// gp nummethods
#include "gpnummethods/cfmethod.h"

#include "gpcalculators\forexvanilla.h"


CC_BEGIN_NAMESPACE(ARM)


///  define for code clarity
#if defined(FIRST_STATE_VARIABLE)
	#undef	FIRST_STATE_VARIABLE
#endif
#define FIRST_STATE_VARIABLE 0

const int NBVARS = 3;
const int NBPRICES = 5;

const int DECVOLPOS = 0;
const int ALPHAPOS = 1;
const int LAMBDAPOS = 2;

const int INITPOS = 0;
const int LBOUNDPOS = 1;
const int UBOUNDPOS = 2;

const int ATMPOS = 2;

const double VEGACOEFF = 0.0001;


////////////////////////////////////////////////////
///	Class  : ARM_ParamsMixture_Fx
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_ParamsMixture_Fx::ARM_ParamsMixture_Fx(
		const std::vector<double>& lags,
		const std::vector<double>& volATM,
		const std::vector<double>& decVol,
		const std::vector<double>& shift,
		const std::vector<double>& lambda,
		const string& interpolName)
		:
itsLags(lags),
itsVolATM(volATM),
itsDecVol(decVol),
itsShift(shift),
itsLambda(lambda),
itsInterpolName(interpolName)
{
	if ((lags.size() != volATM.size()) &&
		(volATM.size() != decVol.size()) &&
		(decVol.size() != shift.size()) &&
		(shift.size() != lambda.size()))
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + ": lags, vol ATM, decVol, shift and lambda should have the same size.");
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_ParamsMixture_Fx
///	Routines: toString
///	Returns :
///	Action  : stringify the object
////////////////////////////////////////////////////
string ARM_ParamsMixture_Fx::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Parameter FX Mixture Model\n";
    os << indent << "-------------\n";
	
	os << std::setiosflags(std::ios::left);
	
    os << indent << CC_NS(std, setfill(' ')) << CC_NS(std, setw)(9) << "Lags";
	os << indent << CC_NS(std, setfill(' ')) << CC_NS(std, setw)(9) << "VolATM";
	os << indent << CC_NS(std, setfill(' ')) << CC_NS(std, setw)(9) << "DecVol";
	os << indent << CC_NS(std, setfill(' ')) << CC_NS(std, setw)(9) << "Shift";
	os << indent << CC_NS(std, setfill(' ')) << CC_NS(std, setw)(9) << "Lambda" << endl;
	
	size_t i;
	
	for (i = 0; i < itsLags.size(); ++i)
	{
		os << indent << CC_NS(std, setfill(' ')) << CC_NS(std, fixed) << CC_NS(std, setw)(9)  << CC_NS(std, setprecision)(0) << itsLags[i];
		os << indent << CC_NS(std, setfill(' ')) << CC_NS(std, fixed) << CC_NS(std, setw)(9)  << CC_NS(std, setprecision)(2) << itsVolATM[i];
		os << indent << CC_NS(std, setfill(' ')) << CC_NS(std, fixed) << CC_NS(std, setw)(9)  << CC_NS(std, setprecision)(2) << itsDecVol[i];
		os << indent << CC_NS(std, setfill(' ')) << CC_NS(std, fixed) << CC_NS(std, setw)(9)  << CC_NS(std, setprecision)(2) << itsShift[i];
		os << indent << CC_NS(std, setfill(' ')) << CC_NS(std, fixed) << CC_NS(std, setw)(9)  << CC_NS(std, setprecision)(2) << itsLambda[i] << endl;
	}
    return os.str();
}

////////////////////////////////////////////////////
///	Class  : ARM_MixtureModel_Fx
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_MixtureModel_Fx::ARM_MixtureModel_Fx(const ARM_ZeroCurvePtr& zc, 
					 ARM_ModelParamsMixture_Fx* modelParam)
:	ARM_EqFxBase(zc, modelParam, ARM_GP_Matrix(NULL), new ARM_MixtureDensityFunctor()),
	itsIsSigmaInputed(false)
{
	if (!zc.IsNull() && modelParam)
		ARM_EqFxBase::Init();
	Validate();
}

////////////////////////////////////////////////////
///	Class  : ARM_MixtureModel_Fx
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_MixtureModel_Fx::ARM_MixtureModel_Fx(
	const ARM_ZeroCurvePtr& domZc,
	const ARM_ZeroCurvePtr& forZc,
	double spot,
	ARM_ParamsMixture_Fx* mixtureParams)
	:
ARM_EqFxBase(domZc, NULL, ARM_GP_Matrix(NULL), new ARM_MixtureDensityFunctor()),
itsIsSigmaInputed(false)
{
	ARM_CurveModelParam volATM(
		ARM_ModelParamType::Volatility,
		const_cast<ARM_GP_Vector*>(&mixtureParams->GetVolATM()),
		const_cast<ARM_GP_Vector*>(&mixtureParams->GetLags()), 
		"",
        mixtureParams->GetInterpolName());
	
	ARM_CurveModelParam decVol(
		ARM_ModelParamType::Smile,
		const_cast<ARM_GP_Vector*>(&mixtureParams->GetDecVol()),
		const_cast<ARM_GP_Vector*>(&mixtureParams->GetLags()), 
		"",
        mixtureParams->GetInterpolName());
	
	ARM_CurveModelParam shift(
		ARM_ModelParamType::Shift,
		const_cast<ARM_GP_Vector*>(&mixtureParams->GetShift()),
		const_cast<ARM_GP_Vector*>(&mixtureParams->GetLags()), 
		"",
        mixtureParams->GetInterpolName());
	
	ARM_CurveModelParam lambda(
		ARM_ModelParamType::QParameter,
		const_cast<ARM_GP_Vector*>(&mixtureParams->GetLambda()),
		const_cast<ARM_GP_Vector*>(&mixtureParams->GetLags()), 
		"",
        mixtureParams->GetInterpolName());
	
	ARM_ModelParamVector empty;
	ARM_ModelParamsMixture_Fx modelParams(empty, domZc, forZc, spot);
	
	modelParams.SetModelParam(&volATM);
	modelParams.SetModelParam(&decVol);
	modelParams.SetModelParam(&shift);
	modelParams.SetModelParam(&lambda);
	
	SetModelParams(modelParams);
	
	ARM_EqFxBase::Init();
	
	Validate();
}


///////////////////////////////////////////////////
///	Class  : ARM_MixtureModel_Fx
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
void ARM_MixtureModel_Fx::Validate() const
{
	const ARM_ModelParamsMixture_Fx* modelparams =  dynamic_cast<const ARM_ModelParamsMixture_Fx*>(GetModelParams());
	if (!modelparams)
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + ": accepted model params are only Mixture Fx model param");
	/// checks that the model has the following model param :( 4 Model Param in total)
	static const string modelParamsName("Fx Mixture  Model Param");
	
	/// 1- Volatility at the money
	if (!modelparams->DoesModelParamExist(ARM_ModelParamType::Volatility))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": " + 
		modelParamsName + " requires " + ARM_ModelParamType::GetTypeString(ARM_ModelParamType::Volatility));
	/// 2- Skew
	if (!modelparams->DoesModelParamExist(ARM_ModelParamType::Smile))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": " + 
		modelParamsName + " requires " + ARM_ModelParamType::GetTypeString(ARM_ModelParamType::Smile));
	/// 3- shift
	if (!modelparams->DoesModelParamExist(ARM_ModelParamType::Shift))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": " + 
		modelParamsName + " requires " + ARM_ModelParamType::GetTypeString(ARM_ModelParamType::Shift));
	/// 4- Q
	if (!modelparams->DoesModelParamExist(ARM_ModelParamType::QParameter))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": " + 
		modelParamsName + " requires " + ARM_ModelParamType::GetTypeString(ARM_ModelParamType::QParameter));
	
	/// 5- Sigma ????
	if (modelparams->DoesModelParamExist(ARM_ModelParamType::Sigma))
		CC_MUTABLE(ARM_MixtureModel_Fx, itsIsSigmaInputed) = true;
}

////////////////////////////////////////////////////
///	Class  : ARM_MixtureModel_Fx
///	Routine: CallVectorial
///	Returns: a vector of call vectorial
///	Action : computes equity or fx call
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MixtureModel_Fx::CallVectorial(
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
	string rightCcy = modelName.substr(3, 3);
	string payCcy	= GetPayModelName();
	if (payCcy == "NoName")
		payCcy = rightCcy;
	
	ARM_EqFxBase::CallType calltype;
	if (!IsInvMkt)
		calltype = payCcy == rightCcy ? ARM_EqFxBase::ClosedFormula : ARM_EqFxBase::Quanto;
	else
		calltype = payCcy == rightCcy ? ARM_EqFxBase::InvClosedFormula : ARM_EqFxBase::InvQuanto;
	
	if (evalTime <= K_NEW_DOUBLE_TOL || states == ARM_PricingStatesPtr(NULL))
	{
		ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
		size_t nbStates = states != ARM_PricingStatesPtr(NULL)? states->size(): 1;
		double strike = strikePerState[0];
		double zcT	= ZcCurve->DiscountPrice(payTime/K_YEAR_LEN);
		ARM_GP_VectorPtr fwdVect = Forward(modelName, evalTime, expiryTime, settlementTime, payTime, states);
		double fwd = (*fwdVect)[0];
		if ((calltype == Quanto) || (calltype == InvQuanto))
		{
			return ARM_EqFxBase::CallVectorial(modelName, evalTime, expiryTime, settlementTime, strikePerState, callPut, payTime, states, context);	
		}
		else if (calltype == ClosedFormula || calltype == InvClosedFormula)// Closed Formula
		{
			double norm = 1.0;
			if (calltype == InvClosedFormula)
			{
				zcT	= ZcCurve->DiscountPrice(expiryTime/K_YEAR_LEN);
				const ARM_ModelParams_Fx* fxmodelparams = dynamic_cast< const ARM_ModelParams_Fx* >(GetModelParams());
				double spot = fxmodelparams->GetSpot();
				double forZcTpay = fxmodelparams->GetForCurve()->DiscountPrice(payTime/K_YEAR_LEN);
				double forZcTset = fxmodelparams->GetForCurve()->DiscountPrice(settlementTime/K_YEAR_LEN);
				norm = forZcTpay/forZcTset*strike*1.0/spot;
				strike = 1.0/strike;
				callPut = -callPut;
				//to avoid the convexification of the fwd
				fwdVect = Forward(modelName, evalTime, expiryTime, settlementTime, payTime, states);
				fwd = (*fwdVect)[0];
			}
			double callPrice = 0.0;
			if (expiryTime >= K_NEW_DOUBLE_TOL)
			{
				double volATM = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(expiryTime);
				double decVol = GetModelParams()->GetModelParam(ARM_ModelParamType::Smile).GetValue(expiryTime);
				double alpha = GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(expiryTime);
				double lambda = GetModelParams()->GetModelParam(ARM_ModelParamType::QParameter).GetValue(expiryTime);
				double sqrtT = sqrt(expiryTime/K_YEAR_LEN);
				double stdDev1 = (volATM - decVol)*sqrtT;
				double stdDevATM = volATM*sqrtT;
				double stdDev2 = itsIsSigmaInputed ? GetModelParams()->GetModelParam(ARM_ModelParamType::Sigma).GetValue(expiryTime) * sqrtT :
				CalibVol2(fwd, fwd, callPut, stdDev1, stdDevATM, alpha, lambda);
				
				/// get the intrinsic value
				
				callPrice = CallMixturePrice(fwd, strike, callPut, stdDev1, stdDev2, alpha, lambda)*zcT;
			}
			else
			{
				if (callPut*fwd > callPut*strike)
				{
					callPrice = callPut* (fwd - strike)*zcT;
				}
			}
			
			return ARM_VectorPtr(new std::vector<double>(1, norm*callPrice));
		}
		else
			ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Mixture FX can price this kind of call.");
	}
	else
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Mixture FX can be used just for its closed forms.");
		
		ARM_GP_VectorPtr dumyFwdVector;
		return dumyFwdVector;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_MixtureModel_Fx
///	Routine: DigitalVectorial
///	Returns: a vector of digital vectorial
///	Action : computes equity or fx digital
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MixtureModel_Fx::DigitalVectorial(
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
	if (digitType == ARM_FXDigitType::analytic)
	{
		ARM_FXName fxName(modelName);
		bool IsInvMkt = fxName.GetIsInvMkt();
		string rightCcy = modelName.substr(3, 3);
		string payCcy	= GetPayModelName();
		if (payCcy == "NoName")
			payCcy = rightCcy;
		
		ARM_EqFxBase::CallType calltype;
		if (!IsInvMkt)
			calltype = payCcy == rightCcy ? ARM_EqFxBase::ClosedFormula : ARM_EqFxBase::Quanto;
		else
			calltype = payCcy == rightCcy ? ARM_EqFxBase::InvClosedFormula : ARM_EqFxBase::InvQuanto;
		
		if (calltype == InvClosedFormula)
		{
			return ARM_EqFxBase::DigitalVectorial(modelName, evalTime, expiryTime, settlementTime, 
				strikePerState, notional, callPut, payTime, digitType, epsilon, states, context);
		}
		else if (calltype == Quanto)
		{
			double strike = 1.0/strikePerState[0];
			std::vector<double> strikes(strikePerState);
			return ARM_EqFxBase::DigitalVectorial(modelName, evalTime, expiryTime, settlementTime, 
				1.0/strikePerState, notional, -callPut, payTime, digitType, epsilon, states, context);
		}
		else
		{
			if (evalTime <= K_NEW_DOUBLE_TOL || states == ARM_PricingStatesPtr(NULL))
			{
				double digitalPrice = 0.0;
				ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
				double zcT	= ZcCurve->DiscountPrice(payTime/K_YEAR_LEN);
				ARM_GP_VectorPtr fwd = Forward(modelName, evalTime, expiryTime, settlementTime, payTime, states);
				double strike = strikePerState[0];
				if (calltype == InvQuanto)
				{
					strike  = 1/strike;
					callPut = -callPut;
				}
				if (expiryTime >= K_NEW_DOUBLE_TOL)
				{
					size_t nbStates = states != ARM_PricingStatesPtr(NULL)? states->size(): 1;
					double volATM = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(expiryTime);
					double decVol = GetModelParams()->GetModelParam(ARM_ModelParamType::Smile).GetValue(expiryTime);
					double alpha = GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(expiryTime);
					double lambda = GetModelParams()->GetModelParam(ARM_ModelParamType::QParameter).GetValue(expiryTime);
					
					double sqrtT = sqrt(expiryTime/K_YEAR_LEN);
					
					double stdDev1 = (volATM - decVol)*sqrtT;
					double stdDevATM = volATM*sqrtT;
					
					double stdDev2 = itsIsSigmaInputed ? GetModelParams()->GetModelParam(ARM_ModelParamType::Sigma).GetValue(expiryTime) * sqrtT :
					CalibVol2((*fwd)[0], (*fwd)[0], callPut, stdDev1, stdDevATM, alpha, lambda);
					/// get the intrinsic value
					digitalPrice = notional*DigitalMixturePrice((*fwd)[0], strike, callPut, stdDev1, stdDev2, alpha, lambda)*zcT;
				}
				else
				{
					if (callPut* (*fwd)[0] > callPut*strike)
						digitalPrice = zcT;
				}
				
				size_t payoffSize = states != ARM_PricingStatesPtr(NULL)? states->size(): 1;
				return ARM_VectorPtr(new std::vector<double>(payoffSize, digitalPrice));
			}
			else
			{
				ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Mixture FX can be used just for its closed forms.");
				
				ARM_GP_VectorPtr dumyFwdVector;
				return dumyFwdVector;
			}
		}
	}
	else
	{
		return ARM_EqFxBase::DigitalVectorial(modelName, evalTime, expiryTime, settlementTime, strikePerState, notional, callPut, payTime, digitType, epsilon, states, context);
	}
}


///////////////////////////////////////////////////
///	Class  : ARM_SABR_Model
///	Routine: ImpliedVol
///	Returns: double
///	Action : To Calculate the Implied Volatility
///  By defaut using BS Formula
////////////////////////////////////////////////////
double ARM_MixtureModel_Fx::ImpliedVol(const ARM_VanillaArg& arg) const
{
	if (!dynamic_cast < const ARM_VanillaEqFxOption *>(& arg))
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + ": ARM_MixtureModel_Fx: VanillaArg is not a good type to compte implid vol");
	
	
	double expiryTime = ((ARM_VanillaEqFxOption&)arg).GetExpiry();
	double settlementTime = ((ARM_VanillaEqFxOption&)arg).GetFwdTime();
	double strike = ((ARM_VanillaEqFxOption&)arg).GetStrike();
	int callPut = ((ARM_VanillaEqFxOption&)arg).GetCallPut();
	double payTime = ((ARM_VanillaEqFxOption&)arg).GetPayTime();
	
	ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
	ARM_GP_VectorPtr fwd = Forward(arg.GetCurveName(), arg.GetEvalTime(), expiryTime, settlementTime, payTime, ARM_PricingStatesPtr(NULL));
	double volATM = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(expiryTime);
	double decVol = GetModelParams()->GetModelParam(ARM_ModelParamType::Smile).GetValue(expiryTime);
	double alpha = GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(expiryTime);
	double lambda = GetModelParams()->GetModelParam(ARM_ModelParamType::QParameter).GetValue(expiryTime);
	
	double sqrtT = sqrt(expiryTime/K_YEAR_LEN);
	double stdDev1 = (volATM - decVol)*sqrtT;
	double stdDevATM = volATM*sqrtT;
	double stdDev2 = itsIsSigmaInputed ? GetModelParams()->GetModelParam(ARM_ModelParamType::Sigma).GetValue(expiryTime) * sqrtT :
	CalibVol2((*fwd)[0], (*fwd)[0], callPut, stdDev1, stdDevATM, alpha, lambda);
	
	/// get the intrinsic value
	double price = CallMixturePrice((*fwd)[0], strike, callPut, stdDev1, stdDev2, alpha, lambda);
	double zcT	= ZcCurve->DiscountPrice(settlementTime/K_YEAR_LEN);
	double vol = VanillaImpliedVol_BS((*fwd)[0], strike, expiryTime, price*zcT, callPut, &stdDev2);
	
	return vol;
}


////////////////////////////////////////////////////
///	Class  : ARM_MixtureModel_Fx
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Default initialisation of the model and the
///          associated numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_MixtureModel_Fx::Init(const string& payModelName, 
											   const ARM_TimeInfoPtrVector& timeInfos)
{
	return ARM_PricingStatesPtr(new ARM_PricingStates(1, 1));
}


////////////////////////////////////////////////////
///	Class  : ARM_MixtureModel_Fx
///	Routine: ModelStateLocalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local and global variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_MixtureModel_Fx::ModelStateLocalVariances(
	const std::vector<double>& timeSteps,
	ARM_MatrixVector& localVariances) const
{
	ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Mixture FX can be used just for its closed forms.");
}

////////////////////////////////////////////////////
///	Class  : ARM_MixtureModel_Fx
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_MixtureModel_Fx::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states, int timeIndex) const
{
	ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Mixture FX can be used just for its closed forms.");
}

////////////////////////////////////////////////////
///	Class  : ARM_MixtureModel_Fx
///	Routine: UpdateDensityFunctor
///	Returns: void
///	Action : Update the density functor at expiryTime		
////////////////////////////////////////////////////

void ARM_MixtureModel_Fx::UpdateDensityFunctor(double forward, double expiryTime, double tenor) 
{
	ARM_MixtureDensityFunctor* densityfunctor = dynamic_cast<ARM_MixtureDensityFunctor*>(&*GetDensityFunctor());
	if (!densityfunctor)
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + " : the density functor must be of the type ARM_MixtureDensityFunctor");
	
	double volATM = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(expiryTime);
	double decVol = GetModelParams()->GetModelParam(ARM_ModelParamType::Smile).GetValue(expiryTime);
	double alpha = GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(expiryTime);
	double lambda = GetModelParams()->GetModelParam(ARM_ModelParamType::QParameter).GetValue(expiryTime);
	
	double vol1 = volATM - decVol;
	double sqrtT = sqrt(expiryTime/K_YEAR_LEN);
	double stdDevATM = volATM*sqrtT;
	double stdDev1   = vol1*sqrtT;
	double callPut = 1;
	double stdDev2 = itsIsSigmaInputed ? GetModelParams()->GetModelParam(ARM_ModelParamType::Sigma).GetValue(expiryTime) * sqrtT :	
	CalibVol2(forward, forward, callPut, stdDev1, stdDevATM, alpha, lambda);
	double vol2 = stdDev2/sqrtT;
	
	densityfunctor->setVol1(vol1);
	densityfunctor->setVol2(vol2);
	densityfunctor->setAlpha(alpha);
	densityfunctor->setLambda(lambda);
}

////////////////////////////////////////////////////
///	Class  : ARM_MixtureModel_Fx
///	Routine: CallMixturePrice
///	Returns: double
///	Action : compute the call price in the Mixture
/// model
////////////////////////////////////////////////////

double ARM_MixtureModel_Fx::CallMixturePrice(double fwd, 
	double strike, 
	int callPut, 
	double stdDev1, 
	double stdDev2, 
	double alpha, 
	double lambda)
{
	double price1 = lambda*BlackSholes_Formula(fwd - alpha, stdDev1, 1.0, strike, callPut);
	double price2 = (1.0 - lambda)*BlackSholes_Formula(fwd + alpha*lambda/ (1.0 - lambda), stdDev2, 1.0, strike, callPut);
	
	return price1 + price2;
}


////////////////////////////////////////////////////
///	Class  : ARM_MixtureModel_Fx
///	Routine: DigitalMixturePrice
///	Returns: double
///	Action : compute the digital price in the Mixture
/// model
////////////////////////////////////////////////////

double ARM_MixtureModel_Fx::DigitalMixturePrice(double fwd, 
	double strike, 
	int callPut, 
	double stdDev1, 
	double stdDev2,
	double alpha,
	double lambda)
{
	double digit1 =  lambda*DigitalBlackSholes_Formula(fwd - alpha, stdDev1, 1.0, strike, callPut);
	double digit2 = (1.0 - lambda)*DigitalBlackSholes_Formula(fwd + alpha*lambda/ (1.0 - lambda), stdDev2, 1.0, strike, callPut);
	
	double sum = digit1 + digit2;
	
	return sum;
}

class CalibVol2Func : public ARM_GP::UnaryFunc<double, double> 
{
public: 
	CalibVol2Func(
		double fwd,
		double strike,
		int callPut,
		double stdDev1,
		double alpha,
		double lambda)
		:
	itsFwd(fwd),
		itsStrike(strike),
		itsCallPut(callPut),
		itsStdDev1(stdDev1),
		itsAlpha(alpha),
		itsLambda(lambda)
	{
	};
	
	virtual double operator()(double stdDev2) const 
	{
		return 	ARM_MixtureModel_Fx::CallMixturePrice(itsFwd, itsStrike, itsCallPut, itsStdDev1, stdDev2, itsAlpha, itsLambda);
	}
private:
	double itsFwd;
	double itsStrike;
	int itsCallPut;
	double itsStdDev1;
	double itsAlpha;
	double itsLambda;
};

////////////////////////////////////////////////////
///	Class  : ARM_MixtureModel_Fx
///	Routine: CalibVol2
///	Returns: double
///	Action : Calibrate the seconde volatility
/// of the mixture model based on the ATM price
////////////////////////////////////////////////////


double ARM_MixtureModel_Fx::CalibVol2(double fwd, double strike, int callPut, double stdDev1, double stdDevATM, double alpha, double lambda)
{
	double ATMPrice = BlackSholes_Formula(fwd, stdDevATM, 1.0, strike, 1.0, callPut);
	
	CalibVol2Func func(fwd, strike, callPut, stdDev1, alpha, lambda);
	
	UnaryFuncWithNumDerivative<double> funcWithDeriv(func);
	
	T_DichotomySolver<UnaryFuncWithNumDerivative<double> > dicho_solver(funcWithDeriv, ATMPrice);
	
	dicho_solver.setInitialGuess(2*stdDevATM - stdDev1);
	
	double root = dicho_solver.Solve();
	
	T_SmoothNewtonRaphsonSolver<UnaryFuncWithNumDerivative<double> > nr_solver(funcWithDeriv, ATMPrice);
	
	nr_solver.setInitialGuess(root);
	
	root = nr_solver.Solve();
	
	return root;
}

struct MixtureData
{
	double fwd;
	int callPut;
	double T;
	double sqrtT;
	double K[NBPRICES];
	double vols[NBPRICES];
	double vegas[NBPRICES];
};



void NAG_CALL Err_BSMixture(Integer n,
	double* x, 
    double* fx,
	double* g,
    Nag_Comm *comm)
{
	double err = 0.;
	double stdDev1, stdDev2;
	double price, vol;
	
	MixtureData* data = (MixtureData*) comm->p;
	
	
	stdDev1 = (data->vols[ATMPOS] - x[DECVOLPOS])*data->sqrtT;
	stdDev2 = ARM_MixtureModel_Fx::CalibVol2(data->fwd, data->K[ATMPOS], data->callPut, stdDev1, data->vols[ATMPOS]*data->sqrtT, x[ALPHAPOS], x[LAMBDAPOS]);
	
	for (int i = 0; i < NBPRICES; i++)
	{		
		price = ARM_MixtureModel_Fx::CallMixturePrice(data->fwd, data->K[i], data->callPut, stdDev1, stdDev2, x[ALPHAPOS], x[LAMBDAPOS]);
		vol = data->vols[i];
		vol = VanillaImpliedVol_BS(data->fwd, data->K[i], data->T, price, data->callPut, &vol);
		err +=(vol - data->vols[i])* (vol - data->vols[i])*data->vegas[i]*VEGACOEFF;
	}
	
	*fx = sqrt(err);
}


////////////////////////////////////////////////////
///	Class  : ARM_MixtureModel_Fx
///	Routine: CalibMixture
///	Returns: void
///	Action : Calibrate the mixture model
/// with 5 prices
////////////////////////////////////////////////////

void ARM_MixtureModel_Fx::CalibMixture(
		double fwd, 
		double T,
		int callPut, 
		const vector<double>& K, 
		const vector<double>& vols,
		const vector<double>& decVol,
		const vector<double>& alpha,
		const vector<double>& lambda,
		vector<double>& outParams)
{
	double err;
    Nag_BoundType bound;
	bound = Nag_Bounds;
    double x[NBVARS];
	double g[NBVARS];
	double bl[NBVARS];
	double bu[NBVARS];
	
	MixtureData data;
	
	x[DECVOLPOS] = decVol[INITPOS];
	x[ALPHAPOS] = alpha[INITPOS];
	x[LAMBDAPOS] = lambda[INITPOS];
	bl[DECVOLPOS] = decVol[LBOUNDPOS];
	bu[DECVOLPOS] = decVol[UBOUNDPOS];
	bl[ALPHAPOS] = alpha[LBOUNDPOS];
	bu[ALPHAPOS] = alpha[UBOUNDPOS];
	bl[LAMBDAPOS] = lambda[LBOUNDPOS];
	bu[LAMBDAPOS] = lambda[UBOUNDPOS];
	
    Nag_E04_Opt options;	
	nag_opt_init(&options);
	
	data.fwd = fwd;
	data.callPut = callPut;
	data.T = T;
	data.sqrtT = sqrt(T);
	
	for (size_t i = 0; i < NBPRICES; ++i)
	{
		data.K[i] = K[i];
		data.vols[i] = vols[i];
		double vega = 0.0;
		BS(
			fwd,
			K[i],
			T,
			vols[i],
			NULL,
			NULL,
			&vega);
		
		data.vegas[i] = vega;
	}
	
    Nag_Comm comm;
    comm.p = &data;
	
	
    static NagError fail, fail2;
    fail.print = FALSE;
	options.list				= FALSE;
    options.print_level			= Nag_NoPrint;
    options.output_level		= Nag_NoOutput;
	options.minor_print_level	= Nag_NoPrint;
	options.print_deriv			= Nag_D_NoPrint;
	
	/// optimisation with boundary using no derivatives (e04jbc)
    nag_opt_bounds_no_deriv(NBVARS, &Err_BSMixture, bound, bl, bu, x, &err, g, &options, &comm, &fail);
	
	/// free the option (e04xzc)
    nag_opt_free(&options, "all", &fail2);
	
	double stdDev1 = (vols[ATMPOS] - x[DECVOLPOS])*data.sqrtT;
	double stdDev2 = ARM_MixtureModel_Fx::CalibVol2(fwd, K[ATMPOS], callPut, stdDev1, vols[ATMPOS]*data.sqrtT, x[ALPHAPOS], x[LAMBDAPOS]);
	
	double priceATM = ARM_MixtureModel_Fx::CallMixturePrice(fwd, fwd, callPut, stdDev1, stdDev2, x[ALPHAPOS], x[LAMBDAPOS]);
	double volATM = vols[ATMPOS];
	volATM = VanillaImpliedVol_BS(fwd, fwd, T, priceATM, callPut, &volATM);
	
	outParams.resize(NBVARS + 2);
	
	outParams[0] = volATM;
	outParams[1 + DECVOLPOS] = x[DECVOLPOS] + volATM - vols[ATMPOS];
	outParams[1 + ALPHAPOS] = x[ALPHAPOS];
	outParams[1 + LAMBDAPOS] = x[LAMBDAPOS];
	outParams[NBVARS + 1] = err;
}

#undef FIRST_STATE_VARIABLE

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

