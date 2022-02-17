/*!
 *
 * Copyright (c) CM CIB January 2005 Paris
 *
 *	\file 2FXModel.cpp
 *
 *  \brief 2 FX model + 1 Payment model
 *
 *	\author  K.BELKHEIR
 *  \ reviewer E.Ezzine
 *	\version 1.0
 *	\date October 2006
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/2FXModel.h"

/// gpmodels
#include "gpmodels/EqFxBase.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/typedef.h"
#include "gpmodels/fxname.h"

/// gpbase
#include "gpbase/singleton.h"

/// gpinfra
#include "gpinfra/modelnamemap.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelparams.h"

/// gpnummethods
#include "gpnummethods/cfmethod.h"

//gpcalculator
#include "gpcalculators\forexvanilla.h"
#include "gpcalculators\fxvanillafactory.h"

#include <algorithm>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class   : ARM_2FXModel
///	Routines: Constructor
///	Returns :
///	Action  : builds the object
///				model types contains the integer on the various model
///	
////////////////////////////////////////////////////
ARM_2FXModel::ARM_2FXModel(	const ARM_ModelNameMap&	modelNameMap, 
	const ARM_CurveMatrix& correlationMatrix )
:
ARM_MultiAssetsModel(&modelNameMap,correlationMatrix.empty() ? NULL:&correlationMatrix)
{
}
///	Class   : ARM_2FXModel
///	Routines: Copy constructor
///	Returns :
///	Action  : constructor
ARM_2FXModel::ARM_2FXModel(const ARM_2FXModel& rhs)
:	ARM_MultiAssetsModel(rhs)
{
}

///////////////////////////////////////////////////
///	Class   : ARM_2FXModel
///	Routine : Init
///	Returns : ARM_PricingStatesPtr
///	Action  : intialise the model
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_2FXModel::Init(const string& payModelName,
			const ARM_TimeInfoPtrVector& timeInfos)
{
	///to complete
	/// Delegate to common code
	return ARM_MultiAssetsModel::Init( payModelName, timeInfos );
}

////////////////////////////////////////////////////
///	Class  : ARM_2FXModel
///	Routine: CallVectorial
///	Returns: a vector of call vectorial
///	Action : computes fx call which can be paid
///		in one of the three currencues present 
///     the 2FX model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_2FXModel::CallVectorial(
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
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+": 2IR FX needs either" + marketQuantoName + "or" + marketQuantoName2 +" to quanto adjustment, please advise");
		marketQuantoName = marketQuantoName2;
	}

	// check if we have only one FX
	if(marketPayoffName==marketQuantoName)
	{
		ARM_ModelNameMap::const_iterator payOffIter = modelMap[marketPayoffName];
		ARM_EqFxBase* payoffModel = dynamic_cast<ARM_EqFxBase*>(&*payOffIter->Model());
		return payoffModel->CallVectorial(modelName, evalTime, expiryTime, settlementTime, strikePerState, callPut, payTime, states, context);
	}

	/// init of call vectorial
	size_t payoffSize = states!=ARM_PricingStatesPtr(NULL)? states->size(): 1;
	ARM_GP_VectorPtr callVect(new ARM_GP_Vector(payoffSize,0.0));

	//NumMethod Validation
	ARM_CFMethod* cfnumMethod = dynamic_cast<ARM_CFMethod*>(&*GetNumMethod());
	if( !cfnumMethod)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : For the moment only the numerical method Closed Form is supported");
	//Integrals Parameters
	ARM_GP_Matrix IntegParameters = cfnumMethod->GetCFParameters();
	
	//strike
	double strike = strikePerState[0];
	//sub model names
	ARM_ModelNameMap::const_iterator payOffIter = modelMap[marketPayoffName];
	ARM_ModelNameMap::const_iterator quantoIter = modelMap[marketQuantoName];

	ARM_EqFxBase* payoffModel = dynamic_cast<ARM_EqFxBase*>(&*payOffIter->Model());
	ARM_EqFxBase* quantoModel = dynamic_cast<ARM_EqFxBase*>(&*quantoIter->Model());

	ARM_GP_VectorPtr fwd_payoff = payoffModel->Forward(marketPayoffName, evalTime, expiryTime, settlementTime, payTime, states );
	ARM_GP_VectorPtr fwd_quanto = quantoModel->Forward(marketQuantoName, evalTime, expiryTime, settlementTime, payTime, states );

	//pricing
	double rho = GetCorrelMatrix()->Interpolate(expiryTime)(0,1);
	ARM_GaussReplic2D* gaussReplic = ARM_FXVanillaFactory.Instance()->CreateFXVanillaAndGaussReplic(payCcy,
				modelName, callPut, strike,ARM_FXVanillaType::callquanto, quantoName, rho);	
	double price = CallPrice(*payoffModel,*quantoModel,gaussReplic,(*fwd_payoff)[0],(*fwd_quanto)[0],expiryTime,IntegParameters);

	//Normalisation	
	ARM_FXVanilla* fxvanilla = gaussReplic->GetFXVanilla();
	double norm = fxvanilla->GetIsInvG2() ? 1.0/(*fwd_quanto)[0] : (*fwd_quanto)[0];

	//payment ZC
	double zcT	= GetZeroCurve()->DiscountPrice(payTime/K_YEAR_LEN);
	return ARM_GP_VectorPtr(new ARM_GP_Vector(1,price* zcT/norm) );
}

////////////////////////////////////////////////////
///	Class  : ARM_2FXModel
///	Routine: CallSpreadVectorial
///	Returns: a vector of callspread vectorial
///	Action : computes fx callspread which can be paid
///		in one of the three currencies present 
///     the 2FX model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_2FXModel::CallSpreadVectorial(
		const string& Model1Name,
		const string& Model2Name,
		double evalTime,
		double expiryTime,
		double settlementTime1,
		double settlementTime2,
		double payTime,
		const ARM_GP_Vector& strikePerState,
		double alpha,
		double beta,		
		const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context) const
{	
	/// init of call vectorial
	size_t payoffSize = states!=ARM_PricingStatesPtr(NULL)? states->size(): 1;
	ARM_GP_VectorPtr callspreadVect(new ARM_GP_Vector(payoffSize,0.0));
	//NumMethod Validation
	ARM_CFMethod* cfnumMethod = dynamic_cast<ARM_CFMethod*>(&*GetNumMethod());
	if( !cfnumMethod)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : For the moment only the numerical method Closed Form is supported");
	//Integrals Parameters
	ARM_GP_Matrix IntegParameters = cfnumMethod->GetCFParameters();	
	//commom arguments

	// Marketmodel corresponding to the FX of the payoff
	const ARM_ModelNameMap& modelMap = *GetModelMap();
	ARM_FXName  fxName(Model1Name);
	ARM_ModelNameMap::const_iterator fxIter1 = modelMap[fxName.GetMktName()];
	ARM_EqFxBase* fxModel1 = dynamic_cast<ARM_EqFxBase*>(&*fxIter1->Model());
	ARM_GP_VectorPtr fwd1 = fxModel1->Forward(fxName.GetMktName(), evalTime, expiryTime, settlementTime1, payTime, states );
	fxName = ARM_FXName(Model2Name);
	ARM_ModelNameMap::const_iterator fxIter2 = modelMap[fxName.GetMktName()];
	ARM_EqFxBase* fxModel2 = dynamic_cast<ARM_EqFxBase*>(&*fxIter2->Model());
	ARM_GP_VectorPtr fwd2 = fxModel2->Forward(fxName.GetMktName(), evalTime, expiryTime, settlementTime2, payTime, states );

	double strike = strikePerState[0];	
	double rho = GetCorrelMatrix()->Interpolate(expiryTime)(0,1);

	string payCcy = GetPayModelName();
	ARM_GaussReplic2D* gaussReplic = ARM_FXVanillaFactory.Instance()->CreateFXVanillaAndGaussReplic(payCcy,
				Model1Name, 1.0, strike,ARM_FXVanillaType::spread, Model2Name, rho, alpha, beta);	
	double price = CallPrice(*fxModel1,*fxModel2,gaussReplic,(*fwd1)[0],(*fwd2)[0],expiryTime,IntegParameters);

	ARM_FXVanilla* fxvanilla = gaussReplic->GetFXVanilla();	
	ARM_GaussReplic2D::QuantoType quantoType = gaussReplic->GetQuantoType();
	double norm = (quantoType == ARM_GaussReplic2D::Quanto1) ? (fxvanilla->GetIsInvG() ? 1.0/(*fwd1)[0]:(*fwd1)[0]) : 
					(quantoType == ARM_GaussReplic2D::Quanto2) ? (fxvanilla->GetIsInvG2() ? 1.0/(*fwd2)[0]:(*fwd2)[0]) : 1.0;

	//payment ZC
	double zcT	= GetZeroCurve()->DiscountPrice(payTime/K_YEAR_LEN);
	return ARM_GP_VectorPtr(new ARM_GP_Vector(1,price* zcT/norm) );
}

////////////////////////////////////////////////////
///	Class  : ARM_2FXModel
///	Routine: RangeAccrualVectorial
///	Returns: a vector of rangeAccrual vectorial
///	Action : computes rangeAccrual ie 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_2FXModel::RangeAccrualVectorial(
		const string& model1Name,
		double evalTime,
		double startTime,
		double endTime,
		const ARM_GP_Vector& fixingTimes,
		double payTime,
		const ARM_GP_Vector& downBarrierVect,
		const ARM_GP_Vector& upBarrierVect,
		const ARM_GP_Vector& notionalVect,
		const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context) const
{	
	/// init of rangeAccrual vectorial
	size_t payoffSize = states!=ARM_PricingStatesPtr(NULL)? states->size(): 1;
	ARM_GP_Vector rangeAccrualVect(payoffSize,0.0);
	//NumMethod Validation
	ARM_CFMethod* cfnumMethod = dynamic_cast<ARM_CFMethod*>(&*GetNumMethod());
	if( !cfnumMethod)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : For the moment only the numerical method Closed Form is supported");

	//init of the Models
	const ARM_ModelNameMap& modelMap = *GetModelMap();
	ARM_FXName  fxName(model1Name);
	ARM_ModelNameMap::const_iterator fxIter1 = modelMap[fxName.GetMktName()];
	ARM_EqFxBase* fxModel1 = dynamic_cast<ARM_EqFxBase*>(&*fxIter1->Model());
	
	//commom arguments
	double expiryTime,settlementTime;  
	double cpnPayTime = 0.0;//to avoid the discount of the keyword DIGIT
	int callPut = 1;//to calculate Proba(ind1>=X) 
	ARM_GP_Vector downBarrierPerStates, upBarrierPerStates;
	double notional;

	//range accrual coupon calculation
	ARM_GP_VectorPtr cpnDown, cpnUp;
	ARM_GP_Vector valueVect;
	int nbFixing = fixingTimes.size();
	for( int k = 0; k < nbFixing; ++k)
	{
		expiryTime = fixingTimes[k];
		settlementTime = expiryTime;
		downBarrierPerStates = ARM_GP_Vector(payoffSize,downBarrierVect[0]);
		upBarrierPerStates = ARM_GP_Vector(payoffSize,upBarrierVect[0]);
		notional = notionalVect[0];
		cpnDown = fxModel1->DigitalVectorial(model1Name, evalTime, expiryTime, settlementTime, downBarrierPerStates, notional, callPut, cpnPayTime, states, context);
		cpnUp   = fxModel1->DigitalVectorial(model1Name, evalTime, expiryTime, settlementTime, upBarrierPerStates, notional, callPut, cpnPayTime, states, context);
		rangeAccrualVect += *cpnDown - *cpnUp;
	}
	rangeAccrualVect /= nbFixing;
	//payment ZC
	double zcT	= GetZeroCurve()->DiscountPrice(payTime/K_YEAR_LEN);
	rangeAccrualVect*=zcT;
	return ARM_GP_VectorPtr(new ARM_GP_Vector(rangeAccrualVect) );
}

////////////////////////////////////////////////////
///	Class  : ARM_2FXModel
///	Routine: CallPriceTypeFx
///	Returns: double
///	Action : compute the call in the case "Fx"
///		
//////////////////////////////////////////////////// 
double ARM_2FXModel::CallPrice( const ARM_EqFxBase& payoffModel, 
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

	ARM_GP_Vector glparams = IntegParameters.empty()? ARM_GP_Vector() : IntegParameters.GetColumns(0);
	ARM_GP_Matrix payoffMatrix= payoffModel.Forward_Mapping(fwd_payoff,	expiryTime,	glparams);

	//Normal Mapping Fwd_Quanto
	densityFunctor = quantoModel.GetDensityFunctor();
	densityFunctor->SetIsDirect(!fxvanilla->GetIsInvG2());
	glparams = (!IntegParameters.empty() && IntegParameters.cols() == 2 )? IntegParameters.GetColumns(1) : ARM_GP_Vector();
	ARM_GP_Matrix quantoMatrix = quantoModel.Forward_Mapping(fwd_quanto,expiryTime,glparams);

	//Set of the NumMethod
	gaussReplic->SetGLParams( payoffMatrix );
	gaussReplic->SetGLParams2( quantoMatrix );
	//Pricing
	double price = gaussReplic->Price();
	
	return price;
}

////////////////////////////////////////////////////
///	Class   : ARM_2FXModel
///	Routines: toString
///	Returns :
///	Action  : stringify the object
////////////////////////////////////////////////////
string ARM_2FXModel::toString(const string& indent, const string& nextIndent) const
{
 	const ARM_ModelNameMap& modelMap = *GetModelMap();

    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "2 FX model (2FX)\n";
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

