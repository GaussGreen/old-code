/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file SABR_Model.cpp
 *
 *  \brief base class for SABR Model
 *	\author  E. Benhamou, O. Croissant
 *	\version 1.0
 *	\date October 2004
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/SABR_Model.h"

/// gpclosedforms
#include "gpclosedforms/extendedsabrformula.h"
#include "gpclosedforms/sabrimpliedvol.h"
#include "gpclosedforms/vanilla_bs.h"

///gpbase
#include "gpbase/numericconstant.h"
#include "gpbase/globalconstant.h"
#include "gpbase/surface.h"
#include "gpbase/surfacetypedef.h"

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelparamtype.h"
#include "gpinfra/surfacelistmodelparam.h"

/// gpmodels
#include "gpmodels/SABR_ModelParams.h"

/// gpcalib
#include "gpcalib/vanillaarg.h"
#include "gpcalib/vanillaswaption.h"
#include "gpcalib/vanillacap.h"
#include "gpcalib/modelfitter.h"

/// Kernel
#include <inst/portfolio.h>


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_SABR_Model
///	Routine: Constructor
///	Returns: 
///	Action : builds the ARM_SABR_Model
////////////////////////////////////////////////////

ARM_SABR_Model::ARM_SABR_Model(const ARM_ZeroCurvePtr& zc,
	const ARM_SABR_ModelParams& params, 
	const long& TypeofModel,
	size_t IntegrationStep  )
:
	ARM_AnalyticIRModel( zc, params, new ARM_SABRDensityFunctor() ),
	itsIntegrationStep(IntegrationStep),
	itsTypeOfModel(TypeofModel),
    itsCurrentArg(NULL)
{}


////////////////////////////////////////////////////
///	Class   : ARM_SABR_Model
///	Routines: Copy constructor
///	Returns :
///	Action  : 
////////////////////////////////////////////////////

ARM_SABR_Model::ARM_SABR_Model(const ARM_SABR_Model& rhs)
:
	ARM_AnalyticIRModel( rhs ),
	itsIntegrationStep( rhs.itsIntegrationStep),
	itsTypeOfModel(rhs.itsTypeOfModel),
    itsCurrentArg(rhs.itsCurrentArg)
{}


////////////////////////////////////////////////////
///	Class  : ARM_SABR_Model
///	Routine: operator=
///	Returns: 
///	Action : Assignment operator
////////////////////////////////////////////////////
ARM_SABR_Model& ARM_SABR_Model::operator=(const ARM_SABR_Model& rhs)
{
	if( this != &rhs )
	{
        delete itsCurrentArg;
		ARM_AnalyticIRModel::operator =(rhs);
		itsIntegrationStep = rhs.itsIntegrationStep;
		itsTypeOfModel = rhs.itsTypeOfModel;
        itsCurrentArg=rhs.itsCurrentArg;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_SABR_Model
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_SABR_Model::~ARM_SABR_Model()
{
    delete itsCurrentArg;
}


////////////////////////////////////////////////////
///	Class  : ARM_SABR_Model
///	Routine: UpdateDensityFunctor
///	Returns: void
///	Action : Update the density functor at expiryTime and tenor		
////////////////////////////////////////////////////
void ARM_SABR_Model::UpdateDensityFunctor(double forward, double expiryTime, double tenor) 
{
	ARM_SABRDensityFunctor* densityfunctor = dynamic_cast<ARM_SABRDensityFunctor*>(&*GetDensityFunctor());
	if(!densityfunctor)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : the density functor must be of the type ARM_SABRDensityFunctor");
	double sigma = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(expiryTime,tenor);
	double rho = GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).GetValue(expiryTime,tenor);
	double nu = GetModelParams()->GetModelParam(ARM_ModelParamType::VolOfVol).GetValue(expiryTime,tenor);
	double beta = GetModelParams()->GetModelParam(ARM_ModelParamType::Beta).GetValue(expiryTime,tenor);
	int sabrType = ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACT;

	//alpha to sigma
	double expiry = expiryTime/K_YEAR_LEN;
	double alpha = SABR_ComputeAlphaFromSigmaATM(forward,forward,expiry,sigma,beta,rho,nu,sabrType);

	densityfunctor->SetAlpha( alpha );
	densityfunctor->SetBeta( beta );
	densityfunctor->SetRho( rho );
	densityfunctor->SetNu( nu );
	densityfunctor->SetSabrType( sabrType );
}
////////////////////////////////////////////////////
///	Class  : ARM_SABR_Model
///	Routine: VanillaCaplet
///	Returns: ARM_VectorPtr
///	Action : computes the price of a caplet
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SABR_Model::VanillaCaplet(
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
	
	double time				= fwdResetTime-evalTime;
	double tenor			= ARM_AnalyticIRModel::GetMatchingTenor((fwdEndTime-fwdStartTime)/K_YEAR_LEN);
	double expiry			= time/K_YEAR_LEN;

	double impliedVol = ((ARM_SABR_ModelParams*)GetModelParams())->ImpliedVol((*liborRate)[0],strikesPerState[0],time,tenor,itsTypeOfModel,itsIntegrationStep);
    double optionValue = payNotional*period* BlackSholes_Formula((*liborRate)[0],impliedVol*sqrt(expiry),(*DF)[0],strikesPerState[0],capFloor);
	
	return ARM_VectorPtr( new std::vector<double>(1,optionValue));
}

////////////////////////////////////////////////////
///	Class   : ARM_SABR_Model
///	Routines: VanillaSwaption
///	Returns :
///	Action  : computes a vanilla swaption
////////////////////////////////////////////////////

ARM_VectorPtr ARM_SABR_Model::VanillaSwaption(
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
	double expiry		= time/K_YEAR_LEN;

    double impliedVol = ((ARM_SABR_ModelParams*)GetModelParams())->ImpliedVol((*swapRate)[0],strikesPerState(0,0),time,tenor,itsTypeOfModel,itsIntegrationStep);
    double optionValue = swapNotional  * BlackSholes_Formula((*swapRate)[0],impliedVol*sqrt(expiry),(*annuity)[0],strikesPerState(0,0),callPut);

	return ARM_VectorPtr( new std::vector<double>(1,optionValue));
}

///	Class   : ARM_SFRM
///	Routines: AdviseCurrentCalibSecIndex
///	Returns : void
///	Action  : Advises the model that the current index of the calibration security
///             is the index given ... The model advises just the model params
////////////////////////////////////////////////////
void ARM_SABR_Model::AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter)
{
    //itsCurrentArg = modelFitter.GetVanillaArgVector(index);
}
    
////////////////////////////////////////////////////
///	Class   : ARM_SABR_Model
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_SABR_Model::Clone() const
{
	return new ARM_SABR_Model(*this);
}

////////////////////////////////////////////////////
///	Class   : ARM_SABR_Model
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_SABR_Model::toString(const string& indent,const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Normal Model \n";
    os << indent << "---------------------\n";
    os << ARM_PricingModel::toString(indent);
    return os.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_SABR_Model
///	Routine : ValidateModelParams
///	Returns : bool
///	Action  : validate model params
////////////////////////////////////////////////////
bool ARM_SABR_Model::ValidateModelParams(const ARM_ModelParams& params) const
{
	const ARM_SABR_ModelParams* SABR_ModelParams = dynamic_cast<const ARM_SABR_ModelParams*>(&params);
	if( !SABR_ModelParams )
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelparams is not of type ARM_SABR_ModelParams" );
	return true;
}

///////////////////////////////////////////////////
///	Class  : ARM_SABR_Model
///	Routine: ImpliedVol
///	Returns: double
///	Action : To Calculate the Implied Volatility
///  By defaut using BS Formula
////////////////////////////////////////////////////
double ARM_SABR_Model::ImpliedVol(const ARM_VanillaArg& arg) const
{
    double underlyingValue,time, tenor, strike;
	std::vector<double> margin = std::vector<double>(1,0.0);

    switch( arg.GetType() )
	{
    case ARM_VanillaArg::VANILLA_CAP:
        {
            if(((ARM_VanillaCapDigitalArg&)arg).GetPayPeriods()->size() != 1)
                ARM_THROW( ERR_INVALID_ARGUMENT, "ImpliedVol: Only one caplet is avalaible" );
            
            double startTime = (*((ARM_VanillaCapDigitalArg&)arg).GetStartTimes())[0];
            double endTime   = (*((ARM_VanillaCapDigitalArg&)arg).GetEndTimes())[0];
            double resetTime = (*((ARM_VanillaCapDigitalArg&)arg).GetResetTimes())[0];            
            double period = (*(((ARM_VanillaCapDigitalArg&)arg).GetPayPeriods()))[0];
            time = (*((ARM_VanillaCapDigitalArg&)arg).GetResetTimes())[0]-arg.GetEvalTime();
            strike = (*((ARM_VanillaCapDigitalArg&)arg).GetStrikes())[0];                       
            tenor = ARM_AnalyticIRModel::GetMatchingTenor((endTime - startTime)/K_YEAR_LEN);            
            
            if(fabs(endTime - (*(((ARM_VanillaCapDigitalArg&)arg).GetPayTimes()))[0]) > ARM_GlobalConstant::ARM_SEVENDAYS_LAG)
                ARM_THROW( ERR_INVALID_ARGUMENT, "ImpliedVol: pricing caplet with EndDate <> PayDate isn't avalaible, please advise!" );
           
            ARM_VectorPtr liborRate	= Libor(arg.GetCurveName(), 
                arg.GetEvalTime(),
                startTime, 
                endTime,
                period,
                resetTime, 
                endTime, 
                ARM_PricingStatesPtr(NULL));

            underlyingValue =(*liborRate)[0];
        }
		break;
	case ARM_VanillaArg::VANILLA_SWAPTION:
        {
            double startTime = ((ARM_VanillaSwaptionArg&)arg).GetStartTime();
            double endTime   = ((ARM_VanillaSwaptionArg&)arg).GetEndTime();
            std::vector<double> dummytTimes;       
            ARM_VectorPtr swapRate = SwapRate(arg.GetCurveName(),
            arg.GetEvalTime(),
		    startTime, 
		    endTime, 
		    *((ARM_VanillaSwaptionArg&)arg).GetFixPayTimes(),
		    *((ARM_VanillaSwaptionArg&)arg).GetFixPayPeriods(),
		    dummytTimes,
            dummytTimes,
            dummytTimes,
            dummytTimes,
            dummytTimes,
		    margin,	/// margin
		    true,	/// isDbleNotional to avoid computing the float cash flows piece by piece
		    ARM_PricingStatesPtr(NULL) );

            underlyingValue =(*swapRate)[0];
            tenor = ARM_AnalyticIRModel::GetMatchingTenor((endTime - startTime)/K_YEAR_LEN);
            time = ((ARM_VanillaSwaptionArg&)arg).GetResetTime()-arg.GetEvalTime();
            strike = (*((ARM_VanillaSwaptionArg&)arg).GetStrikes())[0];
        }
        break;
    default:
        ARM_THROW( ERR_INVALID_ARGUMENT, "ImpliedVol: Only either std caplet or std swaption is avalaible" );
        break;
    }
    double value = ((ARM_SABR_ModelParams*)GetModelParams())->ImpliedVol(underlyingValue,
                    strike,time,tenor,itsTypeOfModel,itsIntegrationStep);

    return value;
}

////////////////////////////////////////////////////
///	Class   : ARM_SABR_Model
///	Routine : ComputeDerivatives
///	Returns : double
///	Action  : computes the derivatives of the SABR model
////////////////////////////////////////////////////
double ARM_SABR_Model::PartialDerivative( const ARM_ModelParam& modelParam, 
       size_t number, 
       size_t factorNb,
       const ARM_VanillaArg& arg ,
       ARM_MktTargetType targetFuncType)
{
    double value;
	switch( arg.GetType() )
	{
	case ARM_VanillaArg::VANILLA_CAP:
        {
            if(((ARM_VanillaCapDigitalArg&)arg).GetPayPeriods()->size() != 1)
                ARM_THROW( ERR_INVALID_ARGUMENT, "ImpliedVol: Only one caplet is avalaible" );
            
            double startTime = (*((ARM_VanillaCapDigitalArg&)arg).GetStartTimes())[0];
            double endTime   = (*((ARM_VanillaCapDigitalArg&)arg).GetEndTimes())[0];
            double resetTime = (*((ARM_VanillaCapDigitalArg&)arg).GetResetTimes())[0];            
            double period = (*(((ARM_VanillaCapDigitalArg&)arg).GetPayPeriods()))[0];
            double time = (*((ARM_VanillaCapDigitalArg&)arg).GetResetTimes())[0]-arg.GetEvalTime();
            double strike = (*((ARM_VanillaCapDigitalArg&)arg).GetStrikes())[0];                       
            double tenor = ARM_AnalyticIRModel::GetMatchingTenor((endTime - startTime)/K_YEAR_LEN);   
			
			if(fabs(endTime - (*(((ARM_VanillaCapDigitalArg&)arg).GetPayTimes()))[0]) > ARM_GlobalConstant::ARM_SEVENDAYS_LAG)
                ARM_THROW( ERR_INVALID_ARGUMENT, "ImpliedVol: pricing caplet with EndDate <> PayDate isn't avalaible, please advise!" );
           
            ARM_VectorPtr liborRate	= Libor(arg.GetCurveName(), 
                arg.GetEvalTime(),
                startTime, 
                endTime,
                period,
                resetTime, 
                endTime, 
                ARM_PricingStatesPtr(NULL));

            double underlyingValue =(*liborRate)[0];
            value = ((ARM_SABR_ModelParams*)GetModelParams())->PartialDerivative(underlyingValue,
                        strike,time,tenor,itsTypeOfModel,itsIntegrationStep, modelParam.GetType());

            if(targetFuncType == ARM_CalibrationTarget::PriceTarget)
            {
                ARM_VectorPtr DF		= DiscountFactor( arg.GetCurveName(),
                    arg.GetEvalTime(), 
                    endTime,
                    ARM_PricingStatesPtr(NULL)  );
	
                double impliedvol = ((ARM_SABR_ModelParams*)GetModelParams())->ImpliedVol(underlyingValue,
                                       strike,time,tenor,itsTypeOfModel,itsIntegrationStep);
                double notional   = (*((ARM_VanillaCapDigitalArg&)arg).GetNotionals())[0];
                double expiry     = time/K_YEAR_LEN;
                int callput       = arg.GetCallPut();

                double vega = notional*BlackSholes_Derivative_2(underlyingValue,
                    impliedvol*sqrt(expiry),
                    (*DF)[0],
                    strike,
                    callput);
                value*=vega;
            }
        }
		break;
	case ARM_VanillaArg::VANILLA_SWAPTION:
        {
			std::vector<double> margin = std::vector<double>(1,0.0);
            std::vector<double> dummyTimes;
            double startTime =((ARM_VanillaSwaptionArg&)arg).GetStartTime();
            double endTime   =((ARM_VanillaSwaptionArg&)arg).GetEndTime();
            double strike = (*((ARM_VanillaSwaptionArg&)arg).GetStrikes())[0];
            double tenor = ARM_AnalyticIRModel::GetMatchingTenor((endTime - startTime)/K_YEAR_LEN);
            double time = ((ARM_VanillaSwaptionArg&)arg).GetResetTime()-arg.GetEvalTime();

            ARM_VectorPtr swapRate = SwapRate(
            arg.GetCurveName(), 
            arg.GetEvalTime(),
            startTime, 
            endTime, 
            *((ARM_VanillaSwaptionArg&)arg).GetFixPayTimes(),
            *((ARM_VanillaSwaptionArg&)arg).GetFixPayPeriods(),
            dummyTimes,
            dummyTimes,
            dummyTimes,
            dummyTimes,
            dummyTimes,
            margin,	/// margin
            true,	/// isDbleNotional to avoid computing the float cash flows piece by piece
            ARM_PricingStatesPtr(NULL) );

            double underlyingValue =(*swapRate)[0];

            value = ((ARM_SABR_ModelParams*)GetModelParams())->PartialDerivative(underlyingValue,
                        strike,time,tenor,itsTypeOfModel,itsIntegrationStep, modelParam.GetType());
            
            if(targetFuncType == ARM_CalibrationTarget::PriceTarget)
            {
                ARM_VectorPtr annuity = Annuity(arg.GetCurveName(), 
                arg.GetEvalTime(),
                *((ARM_VanillaSwaptionArg&)arg).GetFixPayTimes(),
                *((ARM_VanillaSwaptionArg&)arg).GetFixPayPeriods(),
                ARM_PricingStatesPtr(NULL) );

                double impliedvol = ((ARM_SABR_ModelParams*)GetModelParams())->ImpliedVol(underlyingValue,
                                       strike,time,tenor,itsTypeOfModel,itsIntegrationStep);
                double notional   = ((ARM_VanillaSwaptionArg&)arg).GetNotional();
                double expiry     = time/K_YEAR_LEN;
                int callput       = arg.GetCallPut();

                double vega = notional*BlackSholes_Derivative_2(underlyingValue,
                    impliedvol*sqrt(expiry),
                    (*annuity)[0],
                    strike,
                    callput);
                value*=vega;
            }
        }
		break;
	default:
        ARM_THROW( ERR_INVALID_ARGUMENT, "PartialDerivative: Only either std caplet or std swaption is avalaible" );
		break;
	}

    return value;
}

////////////////////////////////////////////////////
///	Class   : ARM_SABR_Model
///	Routine : PostProcessing
///	Returns : double
///	Action  : update the final surafec after calibration
////////////////////////////////////////////////////
void ARM_SABR_Model::PostProcessing(const ARM_ModelFitter& modelFitter)
{
    ((ARM_SABR_ModelParams*)GetModelParams())->PostProcessing(modelFitter, this);
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

