/*!
 *
 * Copyright (c) CM CIB January 2005 Paris
 *
 *	\file HybridIRFX.cpp
 *
 *  \brief 1 FXmodel  + 1 IR or FX model + 1 Payment model
 *
 *	\author  K.BELKHEIR
 * 	\version 1.0
 *	\date March 2007
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/HybridIRFX.h"

/// gpmodels
#include "gpmodels/EqFxBase.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/typedef.h"
#include "gpmodels/fxname.h"
#include "gpmodels/BS_Model.h"
#include "gpmodels/SABR_Model.h"

/// gpbase
#include "gpbase/singleton.h"

//gpclosedforms
#include "gpnumlib/gaussiananalytics.h"
#include "gpclosedforms/normal.h"

/// gpinfra
#include "gpinfra/modelnamemap.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelparams.h"

/// gpnummethods
#include "gpnummethods/cfmethod.h"

/// gpcalculators
#include "gpcalculators/forexvanilla.h"
#include "gpcalculators\fxvanillafactory.h"

#include <algorithm>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class   : ARM_HybridIRFX
///	Routines: Constructor
///	Returns :
///	Action  : builds the object
///				model types contains the integer on the various model
///	
////////////////////////////////////////////////////
ARM_HybridIRFX::ARM_HybridIRFX(	const ARM_ModelNameMap&	modelNameMap, 
	const ARM_CurveMatrix& correlationMatrix )
:
ARM_MultiAssetsModel(&modelNameMap,correlationMatrix.empty() ? NULL:&correlationMatrix)
{
}
///	Class   : ARM_HybridIRFX
///	Routines: Copy constructor
///	Returns :
///	Action  : constructor
ARM_HybridIRFX::ARM_HybridIRFX(const ARM_HybridIRFX& rhs)
:	ARM_MultiAssetsModel(rhs)
{
}

///////////////////////////////////////////////////
///	Class   : ARM_HybridIRFX
///	Routine : Init
///	Returns : ARM_PricingStatesPtr
///	Action  : intialise the model
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_HybridIRFX::Init(const string& payModelName,
			const ARM_TimeInfoPtrVector& timeInfos)
{
	///to complete
	/// Delegate to common code
	return ARM_MultiAssetsModel::Init( payModelName, timeInfos );
}



////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFX
///	Routine: CallVectorial
///	Returns: a vector of call vectorial
///	Action : computes fx call which can be paid
///		in one of the three currencies present 
///     the HybridIRFX model with 2FX as input
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HybridIRFX::DigitalVectorial(
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
	if(	evalTime<=K_NEW_DOUBLE_TOL || states == ARM_PricingStatesPtr(NULL) )
	{
		// Marketmodel corresponding to the FX of the payoff
		ARM_FXName  fxName(modelName);
		string marketPayoffName = fxName.GetMktName();
		string payCcy = GetPayModelName();
		const ARM_ModelNameMap& modelMap = *GetModelMap();
		string forMarketPayoffName = marketPayoffName.substr(0,3); 

		// Marketmodel corresponding to the FX of the quanto, with the convention
		// that marketQuantoName=marketPayoffName without quanto
		string quantoName(forMarketPayoffName+payCcy);
		string marketQuantoName = marketPayoffName;
		if ( forMarketPayoffName != payCcy)
		{
			string quantoName(forMarketPayoffName+payCcy);
			fxName = ARM_FXName(quantoName);
			marketQuantoName = fxName.GetMktName();
		}
			
		// check if the 2FX model has the good model
		if(!(modelMap.TestIfModelExisting(marketQuantoName)))
		{
			quantoName = modelName.substr(3,3)+ payCcy;
			fxName = ARM_FXName(quantoName);
			string marketQuantoName2 = fxName.GetMktName();
			if(!modelMap.TestIfModelExisting(marketQuantoName2))
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+": HybridIRFX model needs either" + marketQuantoName + "or" + marketQuantoName2 +" to quanto adjustment, please advise");
			marketQuantoName = marketQuantoName2;
		}

		// check if we have only one FX
		if(marketPayoffName==marketQuantoName)
		{
			ARM_ModelNameMap::const_iterator payOffIter = modelMap[marketPayoffName];
			ARM_EqFxBase* payoffModel = dynamic_cast<ARM_EqFxBase*>(&*payOffIter->Model());
			if(!payoffModel)
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The submodel corresponding to the payoff FX must be an FX model");
			if(modelMap.TestIfModelExisting(payCcy) ){
				ARM_ModelNameMap::const_iterator fxIter = modelMap[payCcy];
				ARM_PricingModelIR* irModel = dynamic_cast<ARM_PricingModelIR*>(&*fxIter->Model());
				if(!irModel)
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : RangeAccrualVectorial: IR model is required");
				payoffModel->SetConvAdjustModel(irModel);
			}
			ARM_CurveMatrix* correlMatrix = const_cast<ARM_HybridIRFX*>(this)->GetCorrelMatrix();
			payoffModel->SetCorrelMatrix(*correlMatrix);
			return payoffModel->DigitalVectorial(modelName, evalTime, expiryTime, settlementTime, strikePerState, notional, callPut, payTime, digitType, epsilon, states, context);
		}

		/// callspread for the pure quanto digital
		size_t payoffSize = states!=ARM_PricingStatesPtr(NULL)? states->size(): 1;
		ARM_GP_VectorPtr digit(new std::vector<double>(payoffSize,0.0));
		double signed_eps = callPut*epsilon;
		double norm = 1/epsilon;
		std::vector<double> strikeDownPerState, strikeUpPerState;
		ARM_GP_VectorPtr callDown, callUp;

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
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": No analytic method for pure quanto digital." );
				break;
			}
		}
		callDown = CallVectorial(modelName, evalTime, expiryTime, settlementTime, strikeDownPerState, callPut, payTime, states, context);
		callUp   = CallVectorial(modelName, evalTime, expiryTime, settlementTime, strikeUpPerState, callPut, payTime, states, context);
		*digit = norm*notional*(*callDown-*callUp);
		return digit;
	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Hybrid IRFX can be used just for its closed forms." );
		ARM_GP_VectorPtr dumyDigitVector;
		return dumyDigitVector;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFX
///	Routine: SnowBall Price
///	Returns: a simple snowball version
///	Action : computes SnowBall ie 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HybridIRFX::SnowBallVectorial(
		const string& modelName,
		double evalTime,
		double payTime,
		const  std::vector<double>& settelmentTimes,
		const  std::vector<double>& fixingTimes,
		int	   callPut,
		const  std::vector<double>& strikes,
		const  std::vector<double>& coeffs,
		const  std::vector<double>& leverages,
		const  std::vector<double>& couponMin,
		const  std::vector<double>& couponMax,
		const  ARM_PricingStatesPtr& states,
		bool isPerformOrNo,
		ARM_PricingContext* context) const
{
	// Eval time consistency
	if(	evalTime<=K_NEW_DOUBLE_TOL || states == ARM_PricingStatesPtr(NULL) )
	{	
		// For Df computation
		const ARM_ModelNameMap& modelMap = *GetModelMap();
		ARM_FXName  fxName(modelName);
		ARM_ModelNameMap::const_iterator fxIter = modelMap[fxName.GetMktName()];
		ARM_EqFxBase* fxModel = dynamic_cast<ARM_EqFxBase*>(&*fxIter->Model());
		if(!fxModel)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The submodel corresponding to the first FX must be an FX model");
		
		const ARM_ModelParams_Fx* fxModelParams = dynamic_cast< const ARM_ModelParams_Fx* >( fxModel->GetModelParams() );
		if (!fxModelParams)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The model param should be FX.");

		ARM_ZeroCurvePtr domCurve = fxModelParams->GetDomCurve(), forCurve = fxModelParams->GetForCurve();
		//double spot = fxModelParams->GetSpot();

		// Global variables
		size_t couponSize = fixingTimes.size(), stateSize = ( !states.IsNull() ) ? states->size() : 1;
		double globCoeff(1.), strike_i(0.);
		std::vector<double> optionPrice_i(stateSize, 0.), value(stateSize, 0.);
		double payDf = (isPerformOrNo) ? forCurve->DiscountPrice(payTime/K_YEAR_LEN) : domCurve->DiscountPrice(payTime/K_YEAR_LEN);

		string invFXName = modelName.substr(3,3)+modelName.substr(0,3);

		int i,j;
		for(i=0; i<couponSize; i++)
		{
			globCoeff=1.;
			for(j=i+1; j<couponSize; j++)
				globCoeff *= coeffs[j];
			// Limit case 
			bool limitPerfCase = ( (isPerformOrNo) && (strikes[i]==0.) && ( (couponMin[i]!=couponMax[i]) || ( (couponMin[i]!=0.) && (couponMax[i]!=ARM_NumericConstants::ARM_INFINITY) ) ) );
			if (limitPerfCase)
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : No null strike allowed with couponMin/Max FXBallPerf mode!");
				
			if (isPerformOrNo)
			{
				if ( fabs(couponMin[i])<1. )
				{
					strike_i		= strikes[i]/(1.-callPut*couponMin[i]);
					optionPrice_i	= (*(this)->CallVectorial (invFXName, evalTime, fixingTimes[i], settelmentTimes[i], std::vector<double> (1, 1.0/strike_i), -callPut, payTime, states))*strikes[i];
				} else
					optionPrice_i	= std::vector<double>(stateSize, 0.);
				
			} else {
				strike_i		= strikes[i]+callPut*couponMin[i];

				optionPrice_i	= (*(this)->CallVectorial (modelName, evalTime, fixingTimes[i], settelmentTimes[i], std::vector<double> (1, strike_i), callPut, payTime, states));
			}
			
			if (isPerformOrNo)
			{
				if ( fabs(couponMax[i])<1. )
				{
					strike_i		= strikes[i]/(1.-callPut*couponMax[i]);
					optionPrice_i	-= (*(this)->CallVectorial (invFXName, evalTime, fixingTimes[i], settelmentTimes[i], std::vector<double> (1, 1.0/strike_i), -callPut, payTime, states))*strikes[i];
				} else
					optionPrice_i	-= std::vector<double>(stateSize, 0.);
			} else {
				strike_i		= strikes[i]+callPut*couponMax[i];
				optionPrice_i	-= (*(this)->CallVectorial (modelName, evalTime, fixingTimes[i], settelmentTimes[i], std::vector<double> (1, strike_i), callPut, payTime, states));
			}

			optionPrice_i	+= payDf*couponMin[i];

			if (i>0)
				value *= coeffs[i];
			value += leverages[i]*optionPrice_i;
			//value += globCoeff*(leverages[i]*optionPrice_i);
		}

		return ARM_GP_VectorPtr( new std::vector<double>(value) );
	} else {
		ARM_THROW( ERR_PAYOFF_CALCULATION_PB, ARM_USERNAME + " : Forward pricing is not allowed with this FX model!");
		ARM_GP_VectorPtr dumyFwdVector;
		return dumyFwdVector;
	}
} 

////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFX
///	Routine: CallVectorial
///	Returns: a vector of call vectorial
///	Action : computes fx call which can be paid
///		in one of the three currencies present 
///     the HybridIRFX model with 2FX as input
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HybridIRFX::CallVectorial(
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
	if(	evalTime<=K_NEW_DOUBLE_TOL || states == ARM_PricingStatesPtr(NULL) )
	{
		// Marketmodel corresponding to the FX of the payoff
		ARM_FXName  fxName(modelName);
		string marketPayoffName = fxName.GetMktName();
		string rightCcy = marketPayoffName.substr(3,3);
		string payCcy = GetPayModelName();
		if(payCcy == "NoName") payCcy = rightCcy;
		const ARM_ModelNameMap& modelMap = *GetModelMap();
		string forMarketPayoffName = marketPayoffName.substr(0,3); 

		// Marketmodel corresponding to the FX of the quanto, with the convention
		// that marketQuantoName=marketPayoffName without quanto
		string quantoName(forMarketPayoffName+payCcy);
		string marketQuantoName = marketPayoffName;
		if ( forMarketPayoffName != payCcy)
		{
			string quantoName(forMarketPayoffName+payCcy);
			fxName = ARM_FXName(quantoName);
			marketQuantoName = fxName.GetMktName();
		}
			
		// check if the 2FX model has the good model
		if(!(modelMap.TestIfModelExisting(marketQuantoName)))
		{
			quantoName = modelName.substr(3,3)+ payCcy;
			fxName = ARM_FXName(quantoName);
			string marketQuantoName2 = fxName.GetMktName();
			if(!modelMap.TestIfModelExisting(marketQuantoName2))
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+": HybridIRFX model needs either" + marketQuantoName + "or" + marketQuantoName2 +" to quanto adjustment, please advise");
			marketQuantoName = marketQuantoName2;
		}
		
		// payoff fxModel
		ARM_ModelNameMap::const_iterator payOffIter = modelMap[marketPayoffName];
		ARM_EqFxBase* payoffModel = dynamic_cast<ARM_EqFxBase*>(&*payOffIter->Model());
		if(!payoffModel)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The submodel corresponding to the payoff FX must be an FX model");

		//quanto fxModel
		ARM_ModelNameMap::const_iterator quantoIter = modelMap[marketQuantoName];
		ARM_EqFxBase* quantoModel = dynamic_cast<ARM_EqFxBase*>(&*quantoIter->Model());
		if(!quantoModel)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The submodel corresponding to the quanto FX must be an FX model");

		// pay Model
		if(modelMap.TestIfModelExisting(payCcy) ){
			ARM_ModelNameMap::const_iterator payIter = modelMap[payCcy];
			ARM_PricingModelIR* irModel = dynamic_cast<ARM_PricingModelIR*>(&*payIter->Model());
			if(!irModel)
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The submodel corresponding to the payment curve must be an IR model");
			// set de ir modelin the EqFXBase
			payoffModel->SetConvAdjustModel(irModel);
			quantoModel->SetConvAdjustModel(irModel);
		}

		//correls
		ARM_GP_Matrix matrix(GetCorrelMatrix()->Interpolate(expiryTime));

		// check if we have only one FX
		if(marketPayoffName==marketQuantoName)
			return payoffModel->CallVectorial(modelName, evalTime, expiryTime, settlementTime, strikePerState, callPut, payTime, states, context);

		/// init of call vectorial
		size_t payoffSize = states!=ARM_PricingStatesPtr(NULL)? states->size(): 1;
		ARM_GP_VectorPtr callVect(new std::vector<double>(payoffSize,0.0));
	
		//NumMethod Validation
		ARM_CFMethod* cfnumMethod = dynamic_cast<ARM_CFMethod*>(&*GetNumMethod());
		if( !cfnumMethod)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : For the moment only the numerical method Closed Form is supported");
		//Integrals Parameters
		ARM_GP_Matrix IntegParameters = cfnumMethod->GetCFParameters();
		
		//strike
		double strike = strikePerState[0];
		double rhoFX1FX2 = matrix(1,0);//Correl (FX1,FX2)

		//fwds
		ARM_GP_VectorPtr fwd_payoff = payoffModel->Forward(marketPayoffName, evalTime, expiryTime, settlementTime, payTime, states );
		ARM_GP_VectorPtr fwd_quanto = quantoModel->Forward(marketQuantoName, evalTime, expiryTime, settlementTime, payTime, states );
	
		//pricing
		ARM_GaussReplic2D gaussReplic = ARM_FXVanillaFactory.Instance()->CreateFXVanillaAndGaussReplic(payCcy,
					modelName, callPut, strike,ARM_FXVanillaType::vanilla, quantoName, rhoFX1FX2);	
		double price = CallPrice(*payoffModel,*quantoModel,&gaussReplic,(*fwd_payoff)[0],(*fwd_quanto)[0],expiryTime,IntegParameters);

		//Normalisation	
		ARM_FXVanilla* fxvanilla = gaussReplic.GetFXVanilla();
		double norm = fxvanilla->GetIsInvG2() ? 1.0/(*fwd_quanto)[0] : (*fwd_quanto)[0];

		//payment ZC
		double zcT	= GetZeroCurve()->DiscountPrice(payTime/K_YEAR_LEN);
		return ARM_GP_VectorPtr(new ARM_GP_Vector(1,price* zcT/norm) );
	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Hybrid IRFX can be used just for its closed forms." );
		ARM_GP_VectorPtr dumyFwdVector;
		return dumyFwdVector;
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFX
///	Routine: CallSpreadVectorial
///	Returns: a vector of callspread vectorial
///	Action : computes fx callspread which can be paid
///		in one of the three currencies present 
///     the HybridIRFX model with 2FX as input
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HybridIRFX::CallSpreadVectorial(
		const string& model1Name,
		const string& model2Name,
		double evalTime,
		double expiryTime,
		double settlementTime1,
		double settlementTime2,
		double payTime,
		const std::vector<double>& strikePerState,
		double alpha,
		double beta,		
		const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context) const
{	
	if(	evalTime<=K_NEW_DOUBLE_TOL || states == ARM_PricingStatesPtr(NULL) )
	{
		/// init of call vectorial
		size_t payoffSize = states!=ARM_PricingStatesPtr(NULL)? states->size(): 1;
		ARM_GP_VectorPtr callspreadVect(new std::vector<double>(payoffSize,0.0));
	
		// Marketmodel corresponding to the FX of the payoff
		const ARM_ModelNameMap& modelMap = *GetModelMap();
		ARM_FXName  fxName(model1Name);
		ARM_ModelNameMap::const_iterator fxIter1 = modelMap[fxName.GetMktName()];
		ARM_EqFxBase* fxModel1 = dynamic_cast<ARM_EqFxBase*>(&*fxIter1->Model());
		if(!fxModel1)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The submodel corresponding to the first FX must be an FX model");

		fxName = ARM_FXName(model2Name);
		ARM_ModelNameMap::const_iterator fxIter2 = modelMap[fxName.GetMktName()];
		ARM_EqFxBase* fxModel2 = dynamic_cast<ARM_EqFxBase*>(&*fxIter2->Model());
		if(!fxModel2)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The submodel corresponding to the second FX must be an FX model");

		// pay Model
		string payCcy = GetPayModelName();
		if(modelMap.TestIfModelExisting(payCcy) ){
			ARM_ModelNameMap::const_iterator payIter = modelMap[payCcy];
			ARM_PricingModelIR* irModel = dynamic_cast<ARM_PricingModelIR*>(&*payIter->Model());
			if(!irModel)
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The submodel corresponding to the payment curve must be an IR model");

			// set de ir modelin the EqFXBase
			fxModel1->SetConvAdjustModel(irModel);
			fxModel2->SetConvAdjustModel(irModel);
		}
		
		//correls
		double rhoFX1FX2 = GetCorrelMatrix()->Interpolate(expiryTime)(0,1);//Correl (FX1,FX2)

		// forwards
		ARM_GP_VectorPtr fwd1 = fxModel1->Forward(fxName.GetMktName(), evalTime, expiryTime, settlementTime1, payTime, states );
		ARM_GP_VectorPtr fwd2 = fxModel2->Forward(fxName.GetMktName(), evalTime, expiryTime, settlementTime2, payTime, states );

		//NumMethod Validation
		ARM_CFMethod* cfnumMethod = dynamic_cast<ARM_CFMethod*>(&*GetNumMethod());
		if( !cfnumMethod)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : For the moment only the numerical method Closed Form is supported");
		//Integrals Parameters
		ARM_GP_Matrix IntegParameters = cfnumMethod->GetCFParameters();	

		// strike
		double strike = strikePerState[0];	
	
		//pricing
		ARM_GaussReplic2D gaussReplic = ARM_FXVanillaFactory.Instance()->CreateFXVanillaAndGaussReplic(payCcy,
					model1Name, 1.0, strike,ARM_FXVanillaType::spread, model2Name, rhoFX1FX2, alpha, beta);	
		double price = CallPrice(*fxModel1,*fxModel2,&gaussReplic,(*fwd1)[0],(*fwd2)[0],expiryTime,IntegParameters);

		//normalization
		ARM_FXVanilla* fxvanilla = gaussReplic.GetFXVanilla();	
		ARM_GaussReplic2D::QuantoType quantoType = gaussReplic.GetQuantoType();
		double norm = (quantoType == ARM_GaussReplic2D::Quanto1) ? (fxvanilla->GetIsInvG() ? 1.0/(*fwd1)[0]:(*fwd1)[0]) : 
						(quantoType == ARM_GaussReplic2D::Quanto2) ? (fxvanilla->GetIsInvG2() ? 1.0/(*fwd2)[0]:(*fwd2)[0]) : 1.0;

		//payment ZC
		double zcT	= GetZeroCurve()->DiscountPrice(payTime/K_YEAR_LEN);
		return ARM_GP_VectorPtr(new ARM_GP_Vector(1,price* zcT/norm) );
	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Hybrid IRFX can be used just for its closed forms." );
		ARM_GP_VectorPtr dumyFwdVector;
		return dumyFwdVector;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFX
///	Routine: CallSpreadVectorial
///	Returns: a vector of callspread vectorial
///	Action : computes fx callspread which can be paid
///		in one of the three currencies present 
///     the HybridIRFX model with 2FX as input
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HybridIRFX::CallQuotientVectorial(
		const string& model1Name,
		const string& model2Name,
		double evalTime,
		double expiryTime,
		double settlementTime1,
		double settlementTime2,
		double payTime,
		const std::vector<double>& strikePerState,	
		double alpha,
		double beta,
		double strike2,
		const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context) const
{	
	if(	evalTime<=K_NEW_DOUBLE_TOL || states == ARM_PricingStatesPtr(NULL) )
	{
		/// init of call vectorial
		size_t payoffSize = states!=ARM_PricingStatesPtr(NULL)? states->size(): 1;
		ARM_GP_VectorPtr callspreadVect(new std::vector<double>(payoffSize,0.0));
	
		// Marketmodel corresponding to the FX of the payoff
		const ARM_ModelNameMap& modelMap = *GetModelMap();
		ARM_FXName  fxName(model1Name);
		ARM_ModelNameMap::const_iterator fxIter1 = modelMap[fxName.GetMktName()];
		ARM_EqFxBase* fxModel1 = dynamic_cast<ARM_EqFxBase*>(&*fxIter1->Model());
		if(!fxModel1)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The submodel corresponding to the first FX must be an FX model");

		fxName = ARM_FXName(model2Name);
		ARM_ModelNameMap::const_iterator fxIter2 = modelMap[fxName.GetMktName()];
		ARM_EqFxBase* fxModel2 = dynamic_cast<ARM_EqFxBase*>(&*fxIter2->Model());
		if(!fxModel2)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The submodel corresponding to the second FX must be an FX model");

		// pay Model
		string payCcy = GetPayModelName();
		if(modelMap.TestIfModelExisting(payCcy) ){
			ARM_ModelNameMap::const_iterator payIter = modelMap[payCcy];
			ARM_PricingModelIR* irModel = dynamic_cast<ARM_PricingModelIR*>(&*payIter->Model());
			if(!irModel)
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The submodel corresponding to the payment curve must be an IR model");

			// set de ir modelin the EqFXBase
			fxModel1->SetConvAdjustModel(irModel);
			fxModel2->SetConvAdjustModel(irModel);
		}
		
		//correls
		double rhoFX1FX2 = GetCorrelMatrix()->Interpolate(expiryTime)(1,0);//Correl (FX1,FX2)
		// forwards
		ARM_GP_VectorPtr fwd1 = fxModel1->Forward(fxName.GetMktName(), evalTime, expiryTime, settlementTime1, payTime, states );
		ARM_GP_VectorPtr fwd2 = fxModel2->Forward(fxName.GetMktName(), evalTime, expiryTime, settlementTime2, payTime, states );

		//NumMethod Validation
		ARM_CFMethod* cfnumMethod = dynamic_cast<ARM_CFMethod*>(&*GetNumMethod());
		if( !cfnumMethod)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : For the moment only the numerical method Closed Form is supported");
		//Integrals Parameters
		ARM_GP_Matrix IntegParameters = cfnumMethod->GetCFParameters();	

		// strike
		double strike = strikePerState[0];	
	
		//pricing
		ARM_GaussReplic2D gaussReplic = ARM_FXVanillaFactory.Instance()->CreateFXVanillaAndGaussReplic(payCcy,
					model1Name, 1.0, strike,ARM_FXVanillaType::quotient, model2Name, rhoFX1FX2, alpha, beta, strike2);	
		double price = CallPrice(*fxModel1,*fxModel2,&gaussReplic,(*fwd1)[0],(*fwd2)[0],expiryTime,IntegParameters);

		//normalization
		ARM_FXVanilla* fxvanilla = gaussReplic.GetFXVanilla();	
		ARM_GaussReplic2D::QuantoType quantoType = gaussReplic.GetQuantoType();
		double norm = (quantoType == ARM_GaussReplic2D::Quanto1) ? (fxvanilla->GetIsInvG() ? 1.0/(*fwd1)[0]:(*fwd1)[0]) : 
						(quantoType == ARM_GaussReplic2D::Quanto2) ? (fxvanilla->GetIsInvG2() ? 1.0/(*fwd2)[0]:(*fwd2)[0]) : 1.0;

		//payment ZC
		double zcT	= GetZeroCurve()->DiscountPrice(payTime/K_YEAR_LEN);
		return ARM_GP_VectorPtr(new ARM_GP_Vector(1,price* zcT/norm) );
	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Hybrid IRFX can be used just for its closed forms." );
		ARM_GP_VectorPtr dumyFwdVector;
		return dumyFwdVector;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFX
///	Routine: RangeAccrualVectorial
///	Returns: a vector of rangeAccrual vectorial
///	Action : computes rangeAccrual ie 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HybridIRFX::RangeAccrualVectorial(
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
	//init of the Models
	const ARM_ModelNameMap& modelMap = *GetModelMap();
	ARM_FXName  fxName(fxModelName);
	ARM_ModelNameMap::const_iterator fxIter1 = modelMap[fxName.GetMktName()];
	ARM_EqFxBase* fxModel = dynamic_cast<ARM_EqFxBase*>(&*fxIter1->Model());
	if(!fxModel)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : RangeAccrualVectorial: FX model is required");
	ARM_ModelNameMap::const_iterator fxIter2 = modelMap[curveName];
	ARM_PricingModelIR* irModel = dynamic_cast<ARM_PricingModelIR*>(&*fxIter2->Model());
	if(!irModel)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : RangeAccrualVectorial: IR model is required");

	// set de ir model in the EqFXBasefor ConvexityAdjustment
	fxModel->SetConvAdjustModel(irModel);

	//notional and all the barriers are constant for the momemt
	double notional = notionals[0];
	double fxDownBarrier = fxDownBarriers[0];
	double fxUpBarrier = fxUpBarriers[0];
	double irDownBarrier = irDownBarriers[0];
	double irUpBarrier = irUpBarriers[0];

	//check the case
	//pure FX range accrual (no quanto here)
	if(irDownBarrier == irUpBarrier)
	{	
		ARM_CurveMatrix* correlMatrix = const_cast<ARM_HybridIRFX*>(this)->GetCorrelMatrix();
		fxModel->SetCorrelMatrix(*correlMatrix);
		return fxModel->RangeAccrualVectorial(fxModelName, evalTime, startTime, endTime, fixingTimes, payTime,
				fxDownBarriers, fxUpBarriers, notionals, states);
	}
	//pure IR range accrual
	else if(fxDownBarrier == fxUpBarrier)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : For pure rate range accrual, please use the keyword Corridor");
		ARM_GP_VectorPtr dumyFwdVector;
		return dumyFwdVector;
	}
	//hybrid case
	else
	{
		/// init of rangeAccrual vectorial
		size_t payoffSize = states!=ARM_PricingStatesPtr(NULL)? states->size(): 1;
		std::vector<double> rangeAccrualVect(payoffSize,0.0);

		//commom arguments
		double expiryTime,fwdPeriod, resetTime;  
		ARM_DensityFunctor* fxdensity = dynamic_cast<ARM_DensityFunctor*>(fxModel->GetDensityFunctor());
		ARM_DensityFunctor* irdensity = dynamic_cast<ARM_DensityFunctor*>(irModel->GetDensityFunctor());;
	
		//range accrual coupon calculation
		int nbFixing = fixingTimes.size();
		double price = 0;

		//xmin and xmax
		double x_inf = -6.0;
		double x_sup = +6.0;
		double fxXmin = -6.0;
		double irXmin = -6.0;
		double fxXmax = +6.0;
		double irXmax = +6.0;

		//Used for local model calibration for CRA double
		eachFixingPrices->push_back(nbFixing);

		for( int k = 0; k < nbFixing; ++k)
		{
			expiryTime = fixingTimes[k];
			resetTime  = irIndexResetTimes[k];
			startTime  = irIndexStartTimes[k];
			endTime    = irIndexEndTimes[k];
			fwdPeriod  = irIndexTerms[k];
			ARM_GP_VectorPtr fxfwdVect = fxModel->Forward(fxModelName, evalTime, expiryTime, expiryTime, expiryTime, states );
			ARM_GP_VectorPtr irfwdVect = irModel->Libor(curveName, evalTime, startTime, endTime, fwdPeriod, resetTime, endTime, states );
			double fxfwd = (*fxfwdVect)[0];
			double irfwd = (*irfwdVect)[0];

			//Update of the densities
			fxModel->UpdateDensityFunctor(fxfwd,expiryTime,0.0);
			irModel->UpdateDensityFunctor(irfwd,expiryTime,fwdPeriod);

			//get the correlation
			double rhoIRFX1 = fxModel->GetRho().Interpolate(expiryTime);//Correl (IR,FX1)			

			//change of variable to work in Gaussian world: FX = f(fxX) and IR = g(irX), with (fxX,irX)~bivariate(rho)
			double maturity = expiryTime/K_YEAR_LEN;
			// FX
			double fxXmin = ARM_GaussianAnalytics::cdfNormal_Inv( fxdensity->Proba(fxDownBarrier,fxfwd,maturity) ); 
			double fxXmax = ARM_GaussianAnalytics::cdfNormal_Inv( fxdensity->Proba(fxUpBarrier,fxfwd,maturity) ); 
			if(fxXmax < fxXmin) CC_NS(std,swap) (fxXmax, fxXmin);
			fxXmax = CC_Min (fxXmax, x_sup);
			fxXmin = CC_Max (fxXmin, x_inf);
			// IR
			double irXmin = ARM_GaussianAnalytics::cdfNormal_Inv( irdensity->Proba(irDownBarrier,irfwd,maturity) );//a revoir
			double irXmax = ARM_GaussianAnalytics::cdfNormal_Inv( irdensity->Proba(irUpBarrier,irfwd,maturity) ); 
			if(irXmax < irXmin) CC_NS(std,swap) (irXmax, irXmin);
			irXmax = CC_Min (irXmax, x_sup);
			irXmin = CC_Max (irXmin, x_inf);
		
			//computation of Proba( fxXmin<fxX<fxXmax ; irXmin<irX<irXmax ) = Phi(fxXmin,irXmin) - Phi(fxXmax,irXmin) - Phi(fxXmin,irXmax) + Phi(fxXmax,irXmax)
			double cpn1 = NormalCDF(fxXmin,irXmin,rhoIRFX1);//Phi(fxXmin,irXmin) (bivariate_cdfNormal is bugged)
			double cpn2 = NormalCDF(fxXmax,irXmin,rhoIRFX1);//Phi(fxXmax,irXmin) 
			double cpn3 = NormalCDF(fxXmin,irXmax,rhoIRFX1);//Phi(fxXmin,irXmax)
			double cpn4 = NormalCDF(fxXmax,irXmax,rhoIRFX1);//Phi(fxXmax,irXmax)
			double cpn = cpn1 - cpn2 - cpn3 + cpn4;

			//summation
			price += cpn;
			

			//store target volatilities to calibrate local model for CRA double condition
			eachFixingPrices->push_back(expiryTime);

			double price1 = fxdensity->Call_Option(fxDownBarrier, fxfwd, maturity);
			eachFixingPrices->push_back(price1);
			eachFixingPrices->push_back(fxXmin);

			double price2 = fxdensity->Call_Option(fxUpBarrier, fxfwd, maturity);
			eachFixingPrices->push_back(price2);
			eachFixingPrices->push_back(fxXmax);

			double price3 = irdensity->Call_Option(MAX(irDownBarrier,0.0), irfwd, maturity);
			eachFixingPrices->push_back(price3);
			eachFixingPrices->push_back(irXmin);

			double price4 = irdensity->Call_Option(irUpBarrier, irfwd, maturity);
			eachFixingPrices->push_back(price4);
			eachFixingPrices->push_back(irXmax);

			//end 

		}
		price /= nbFixing;
		//payment ZC
		double zcT	= GetZeroCurve()->DiscountPrice(payTime/K_YEAR_LEN);
		return ARM_GP_VectorPtr(new ARM_GP_Vector(1,price*zcT*notional) );
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFX
///	Routine: CallPriceTypeFx
///	Returns: double
///	Action : compute the call in the case "Fx"
///		
//////////////////////////////////////////////////// 
double ARM_HybridIRFX::CallPrice( const ARM_EqFxBase& payoffModel, 
		const ARM_EqFxBase& quantoModel,
		ARM_GaussReplic2D* gaussReplic,
		double  fwd_payoff,
		double  fwd_quanto,
		double	expiryTime,
		const ARM_GP_Matrix& IntegParameters) const
{
	ARM_FXVanilla* fxvanilla = gaussReplic->GetFXVanilla();	

	//Normal Mapping Fwd_payoff
	ARM_DensityFunctor* densityFunctor = payoffModel.GetDensityFunctor();
	densityFunctor->SetIsDirect(!fxvanilla->GetIsInvG());

	std::vector<double> glparams = IntegParameters.empty()? std::vector<double>() : IntegParameters.GetColumns(0);
	ARM_GP_Matrix payoffMatrix= payoffModel.Forward_Mapping(fwd_payoff,	expiryTime,	glparams);

	//Normal Mapping Fwd_Quanto
	densityFunctor = quantoModel.GetDensityFunctor();
	densityFunctor->SetIsDirect(!fxvanilla->GetIsInvG2());
	glparams = (!IntegParameters.empty() && IntegParameters.cols() == 2 )? IntegParameters.GetColumns(1) : std::vector<double>();
	ARM_GP_Matrix quantoMatrix = quantoModel.Forward_Mapping(fwd_quanto,expiryTime,glparams);

	//Set of the NumMethod
	gaussReplic->SetGLParams( payoffMatrix );
	gaussReplic->SetGLParams2( quantoMatrix );
	//Pricing
	double price = gaussReplic->Price();
	
	return price;
}

////////////////////////////////////////////////////
///	Class   : ARM_HybridIRFX
///	Routines: toString
///	Returns :
///	Action  : stringify the object
////////////////////////////////////////////////////
string ARM_HybridIRFX::toString(const string& indent, const string& nextIndent) const
{
 	const ARM_ModelNameMap& modelMap = *GetModelMap();

    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "HybridIRFX model (IRFX)\n";
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

    os << indent << "\n\n------> Stochastic FX1 Model <------\n";
    os << modelMap[0]->Model()->toString(indent,nextIndent);

    os << indent << "\n\n------> Stochastic FX2 Model  <------\n";
    os << modelMap[1]->Model()->toString(indent,nextIndent);

    return os.str();
}

#define ARM_CF_EPS 1.0e-13

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

