/*!
 *
 * Copyright (c) CDC IXIS CM July 2004 Paris
 *
 *	\file MultiAssets.cpp
 *
 *  \brief multi asset model
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2004
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/MultiAssets.h"
#include "gpmodels/forwardforex.h"
#include "gpmodels/ForwardMargin.h"


/// gpbase
#include "gpbase/gpmatrix.h"
#include "gpbase/vectormanip.h"
#include "gpbase/datestrip.h"
#include "gpinfra/zccrvfunctor.h"

/// gpinfra
#include "gpinfra/modelnamemap.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/timeinfo.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/discretisationscheme.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/pricingmodeltype.h"
#include "gpinfra/ModelParamsVec.h"

/// gpcalib
#include "gpcalib/modelfitter.h"
#include "gpcalib/calibmethod.h"

/// Kernel

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: Constructor
///	Returns :
///	Action  : builds the object with a correl
/// curve matrix
////////////////////////////////////////////////////

ARM_MultiAssetsModel::ARM_MultiAssetsModel(
		const ARM_ModelNameMap*	modelNameMap,
		const ARM_GP_Matrix* correlMatrix)
:	ARM_PricingModelIR(),
	ARM_PricingFunctionEquity(),
	itsModelMap( modelNameMap? static_cast<ARM_ModelNameMap*>(modelNameMap->Clone() ) : NULL )
{
	if(correlMatrix)
		correlMatrix->CheckCorrelMatrix();

	if (correlMatrix)
		itsCorrelMatrix = new ARM_CurveMatrix(*correlMatrix);
	else
		itsCorrelMatrix = NULL;

	InitMultiAsset();
	InitModelNb();
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: Constructor
///	Returns :
///	Action  : builds the object with a correl
/// curve matrix
////////////////////////////////////////////////////

ARM_MultiAssetsModel::ARM_MultiAssetsModel(
	const ARM_ModelNameMap*	modelNameMap, 
	const ARM_CurveMatrix* correlMatrix )
:	
	ARM_PricingModelIR(), 
	ARM_PricingFunctionEquity(),
	itsModelMap( modelNameMap? static_cast<ARM_ModelNameMap*>(modelNameMap->Clone() ) : NULL ),
	itsCorrelMatrix( correlMatrix? static_cast<ARM_CurveMatrix*>(correlMatrix->Clone() ) : NULL ),
	itsRefModel(NULL)
{
	/// check that the matrix is symmetric and with real correlation number
	if(correlMatrix)
		correlMatrix->CheckCorrelMatrix();

	InitMultiAsset();
	InitModelNb();
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: InitMultiAsset
///	Returns :
///	Action  : Initialise the multi asset for
/// the two constructors
////////////////////////////////////////////////////

void ARM_MultiAssetsModel::InitMultiAsset()
{
	if( itsModelMap )
	{
		if( !itsModelMap->size() )
			ARM_THROW( ERR_INVALID_ARGUMENT, "Multi needs at least one model!" );

		ARM_ModelNameMap::const_iterator 
		iter	= itsModelMap->begin(), 
		end		= itsModelMap->end();

		ARM_Date asOfRefDate = (*iter).Model()->GetAsOfDate();
		(*iter).Model()->SetFromMultiFactor(true);
		++iter;

		for( ; iter!= end; ++iter )
		{
			if(		(*iter).Model()->GetZeroCurve() != ARM_ZeroCurvePtr(NULL) 
				&&	(*iter).Model()->GetAsOfDate() != asOfRefDate )
			{
				CC_Ostringstream os;
				os  << "Model " << itsModelMap->begin()->ModelName() << " with asOfDate " << asOfRefDate.toString() 
					<< " while model " << (*iter).ModelName() << " with asOfDate " << (*iter).Model()->GetAsOfDate().toString();
				ARM_THROW( ERR_INVALID_ARGUMENT, os.str() );
			}
			(*iter).Model()->SetFromMultiFactor(true);
		}
	}
	UpdateSubModelLinks();
	SetModelParamsVec();
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: Copy Constructor
///	Returns :
///	Action  : builds the object
////////////////////////////////////////////////////
ARM_MultiAssetsModel::ARM_MultiAssetsModel(const ARM_MultiAssetsModel& rhs)
:	ARM_PricingModelIR( rhs ),
	ARM_PricingFunctionEquity( rhs ),
	itsModelMap( CreateClone(rhs.itsModelMap) ),
	itsCorrelMatrix( CreateClone(rhs.itsCorrelMatrix) ),
	itsRefModel(NULL)
{
	UpdateSubModelLinks();
	SetModelParamsVec();
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: UpdateSubModelLinks
///	Returns :
///	Action  : updates the ref model
////////////////////////////////////////////////////
void ARM_MultiAssetsModel::UpdateSubModelLinks()
{
	if( itsModelMap )
	{
		/// handle the case of fwd margin model
		ARM_ModelNameMap::iterator 
			iter	= itsModelMap->begin(), 
			end		= itsModelMap->end();

		for( ; iter!= end; ++iter )
			(*iter).Model()->UpdateLinks( *this );
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: SetPayModelName
///	Returns :
///	Action  : set PayName to subModels
////////////////////////////////////////////////////
void ARM_MultiAssetsModel::SetPayModelName(const string& modelName ) 
{ 
	ARM_PricingModel::SetPayModelName( modelName );
	if( itsModelMap )
	{
		/// handle the case of fwd margin model
		ARM_ModelNameMap::iterator 
			iter	= itsModelMap->begin(), 
			end		= itsModelMap->end();

		for( ; iter!= end; ++iter )
			(*iter).Model()->SetPayModelName( modelName );
	}
 }

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: ~ARM_MultiAssetsModel
///	Returns :
///	Action  : destructor
////////////////////////////////////////////////////

ARM_MultiAssetsModel::~ARM_MultiAssetsModel()
{	
	delete itsModelMap;  
	delete itsCorrelMatrix;
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: toString
///	Returns :
///	Action  : stringify the object
////////////////////////////////////////////////////
string ARM_MultiAssetsModel::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << " ----------- Hybrid multi-asset model\n";
	os << itsModelMap->toString(indent,nextIndent);
	
	if( itsRefModel )
		os << "Ref Model: " << itsRefModel->GetModelName() <<"\n";
	
	if( itsCorrelMatrix )
	{
		os << "Correlation matrix\n";
		os << itsCorrelMatrix->toString(indent,nextIndent);
	}
	
    return os.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : ValidateModelParams
///	Returns : throw an exception since this does not make any sense for this kind of model!
///	Action  : validate model params
////////////////////////////////////////////////////
bool ARM_MultiAssetsModel::ValidateModelParams(const ARM_ModelParams& params) const
{
#ifdef __GP_STRICT_VALIDATION
	if( !dynamic_cast<const ARM_ModelParamsVec*>( &params ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_MultiAssetsModel only support paremeters vector!" );
#endif
	return true;
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : CallVectorial
///	Returns : ARM_VectorPtr
///	Action  : call CallVectorial on  the sub model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MultiAssetsModel::CallVectorial(
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
	ARM_PricingModelPtr model = (*itsModelMap)[modelName]->Model();
	ARM_PricingFunctionEquity* equityFctor = dynamic_cast<ARM_PricingFunctionEquity*>(&*model);

	if( !equityFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + modelName + " does not support equity function CallVectorial!" );
	
	//return equityFctor->CallVectorial( modelName, evalTime, expiryTime, settlementTime, strikePerState, callPut, payTime, states, context );
	size_t statesSize = model->ModelStatesSize()?model->ModelStatesSize():1;
	return equityFctor->CallVectorial( 
										modelName, 
										evalTime, 
										std::vector<double>(statesSize,expiryTime), 
										std::vector<double>(statesSize,settlementTime), 
										strikePerState, 
										callPut, 
										std::vector<double>(statesSize,payTime), 
										states, 
										context );
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : Forward
///	Returns : ARM_VectorPtr
///	Action  : call Forward on the sub model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MultiAssetsModel::Forward(
	const string& modelName, 
    double evalTime,
	double expiryTime,
	double settlementTime,
	double payTime,
    const ARM_PricingStatesPtr& states) const
{
	ARM_PricingModelPtr model = (*itsModelMap)[modelName]->Model();
	ARM_PricingFunctionEquity* equityFctor = dynamic_cast<ARM_PricingFunctionEquity*>(&*model);

	if( !equityFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + modelName + " does not support equity function Forward!" );

	return equityFctor->Forward( modelName, evalTime, expiryTime, settlementTime, payTime, states );
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : Libor
///	Returns : ARM_VectorPtr
///	Action  : call Libor on the sub model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MultiAssetsModel::Libor( 
	const string& modelName, 
	double evalTime,
	double fwdStartTime,
	double fwdEndTime,
	double period,
	double fwdResetTime,
	double payTime,
	const ARM_PricingStatesPtr& states) const
{
	ARM_PricingModelPtr model = (*itsModelMap)[modelName]->Model();
	ARM_PricingFunctionIR* IRFctor = dynamic_cast<ARM_PricingFunctionIR*>(&*model);

	if( !IRFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + modelName + " does not support IR function Libor!" );

	return IRFctor->Libor( modelName, evalTime, fwdStartTime, fwdEndTime, period, fwdResetTime, payTime, states );
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : SwapRate
///	Returns : ARM_VectorPtr
///	Action  : call SwapRate on the sub model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MultiAssetsModel::SwapRate(
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
	const ARM_PricingStatesPtr& states) const
{
	ARM_PricingModelPtr model = (*itsModelMap)[modelName]->Model();
	ARM_PricingFunctionIR* IRFctor = dynamic_cast<ARM_PricingFunctionIR*>(&*model);
	if( !IRFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + modelName + " does not support IR function SwapRate!" );

	return IRFctor->SwapRate( modelName, evalTime, floatStartTime, floatEndTime, fixPayTimes, fixPayPeriods, fwdStartTimes,
	    fwdEndTimes, fwdPayPeriods, floatPayTimes, floatPayPeriods, margin, isDbleNotional, states );
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: Spread
///	Returns : a vector of spread values
///	Action  : Default Spread computation
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MultiAssetsModel::Spread(
		const string& modelName, 
		double evalTime,
		double coeff1,
		double floatStartTime1, 
		double floatEndTime1, 
		const std::vector<double>& fixPayTimes1,
		const std::vector<double>& fixPayPeriods1,
		const std::vector<double>& fwdStartTimes1,
        const std::vector<double>& fwdEndTimes1,
        const std::vector<double>& fwdPayPeriods1,
		const std::vector<double>& floatPayTimes1,
        const std::vector<double>& floatPayPeriods1,
        const std::vector<double>& margin1,
		double coeff2,
		double floatStartTime2, 
		double floatEndTime2, 
		const std::vector<double>& fixPayTimes2,
		const std::vector<double>& fixPayPeriods2,
		const std::vector<double>& fwdStartTimes2,
        const std::vector<double>& fwdEndTimes2,
        const std::vector<double>& fwdPayPeriods2,
		const std::vector<double>& floatPayTimes2,
        const std::vector<double>& floatPayPeriods2,
        const std::vector<double>& margin2,
		const ARM_PricingStatesPtr& states) const
{
	ARM_PricingModelPtr model = (*itsModelMap)[modelName]->Model();
	ARM_PricingFunctionIR* IRFctor = dynamic_cast<ARM_PricingFunctionIR*>(&*model);
	if( !IRFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + modelName + " does not support IR function SwapRate!" );

	return IRFctor->Spread(	modelName,
				evalTime,
				coeff1,
				floatStartTime1, 
				floatEndTime1, 
				fixPayTimes1,
				fixPayPeriods1,
				fwdStartTimes1,
				fwdEndTimes1,
				fwdPayPeriods1,
				floatPayTimes1,
				floatPayPeriods1,
				margin1,
				coeff2,
				floatStartTime2, 
				floatEndTime2, 
				fixPayTimes2,
				fixPayPeriods2,
				fwdStartTimes2,
				fwdEndTimes2,
				fwdPayPeriods2,
				floatPayTimes2,
				floatPayPeriods2,
				margin2,
				states);
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : NPVSwap
///	Returns : ARM_VectorPtr
///	Action  : call NPVSwap on the sub model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MultiAssetsModel::NPVSwap(
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
	ARM_PricingModelPtr model = (*itsModelMap)[modelName]->Model();

	ARM_PricingFunctionIR* IRFctor = dynamic_cast<ARM_PricingFunctionIR*>(&*model);
	if( !IRFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + modelName + " does not support IR function NPVSwap!" );

	return IRFctor->NPVSwap( modelName, evalTime, floatStartTime, floatEndTime, fixPayTimes, fixPayPeriods,
		fwdStartTimes, fwdEndTimes, fwdPayPeriods, floatPayTimes, floatPayPeriods, margin, isDbleNotional,
		FixNotional,FloatNotional, strikesPerState, payRec, states);
}

////////////////////////////////////////////////////
///	Class  : ARM_QModelBaseIR
///	Routine: DefaultNPVSwapLeg
///	Returns: a matrix of NPVSwap(t,F(R,Ti),K)
///	Action : 
/// Default: Default NPVSwalLeg computation
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_MultiAssetsModel::NPVSwapLeg(
		const string& curveName, 
		double evalTime,
		const std::vector<double>& fwdStartTimes, 
		const std::vector<double>& fwdEndTimes, 
		const std::vector<double>& fwdPayPeriods, 
		const std::vector<double>& payTimes, 
		const std::vector<double>& payPeriods, 
		const std::vector<double>& margin, 
		const std::vector<double>& notional, 
		const ARM_PricingStatesPtr& states) const
{
	ARM_PricingModelPtr model = (*itsModelMap)[curveName]->Model();
	ARM_PricingFunctionIR* IRFctor = dynamic_cast<ARM_PricingFunctionIR*>(&*model);
	if( !IRFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + curveName + " does not support IR function NPVSwapLeg!" );

	ARM_GP_MatrixPtr result= IRFctor->NPVSwapLeg(curveName, evalTime,fwdStartTimes,fwdEndTimes, fwdPayPeriods, 
		 payTimes, payPeriods, margin, notional,  states); 

	return result;
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: NPVFixLeg
///	Returns : vector ptr
///	Action  : To calculate a fix leg
ARM_GP_MatrixPtr ARM_MultiAssetsModel::NPVFixLeg(
		const string& curveName, 
		double evalTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const std::vector<double>& FixNotional,
		const ARM_GP_Matrix& strikesPerState,
		int   payRec,
		const ARM_PricingStatesPtr& states) const
{
	ARM_PricingModelPtr model = (*itsModelMap)[curveName]->Model();
	ARM_PricingFunctionIR* IRFctor = dynamic_cast<ARM_PricingFunctionIR*>(&*model);
	if( !IRFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + curveName + " does not support IR function NPVSwap!" );

	return IRFctor->NPVFixLeg(curveName,evalTime,fixPayTimes,fixPayPeriods,
		FixNotional, strikesPerState,payRec,states);
}

////////////////////////////////////////////////////
///	Class   : ARM_PricingModelIR
///	Routines: NPVBasisSwapWithNullStrikeAndNotionalExchange
///	Returns : a vector of swap rate values with null strike
///	Action  : Default Swap Rate computation with null strike and exchanging notionl flow by flow
///				WARNING: Only floating leg is computed with a notional exchanging
ARM_VectorPtr ARM_MultiAssetsModel::NPVBasisSwap( 
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
	/// Compute the fixed leg price
	int nbStates= states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;

	/// Check and validate 
	string DomCcyName(GetCurrency( domCurveName )->GetCcyName());
	string forCcyName(GetCurrency( forCurveName )->GetCcyName());

	string payCcyName(GetCurrency(GetPayModelName())->GetCcyName());
	if(DomCcyName != payCcyName )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Payment model and domestic model should have the same currency, please advise!");

	/// Forex forward computing
	ARM_VectorPtr forexVector = Forward(fxCurveName, evalTime,evalTime, evalTime, evalTime,states);

	///  foreign leg computing
	ARM_GP_MatrixPtr foreignLeg;
	if(isForFlottant){

		/// This leg is computed in own currency and discounted by payModel
		foreignLeg = NPVSwapLeg( forCurveName, evalTime,forFwdStartTimes, forFwdEndTimes, 
			forFwdPayPeriods, forPayTimes, forPayPeriods, forMarginVector,forNotionalVector,states);
	}else{
		/// fixed leg computing
		foreignLeg = NPVFixLeg( forCurveName, evalTime,forPayTimes, forPayPeriods,forNotionalVector, forStrikePerState, 1, states);
	}

	ARM_PricingFunctionIR* forIRFctor = dynamic_cast<ARM_PricingFunctionIR*>(&*(*itsModelMap)[forCurveName]->Model());
	
	/// computing notional exchange
	int i, size = forPayTimes.size();
	ARM_VectorPtr forZcStart = GetDiscountFunctor()->DiscountFactor(forCurveName,evalTime,forFlowStartTimes[0],states);
    ARM_VectorPtr forZcEnd   = GetDiscountFunctor()->DiscountFactor(forCurveName,evalTime,forFlowEndTimes[size-1],states);
	
	double alpha, beta;
	if(exNotionalType == "BOTH")
		alpha=1.0, beta = 1.0;
	else if(exNotionalType == "START")
		alpha=1.0, beta = 0.0;
	else if(exNotionalType == "END")
		alpha=0.0, beta = 1.0;
	else
		ARM_THROW(  ERR_INVALID_ARGUMENT, " EXNotional: Only BOTH,START or END are valid" );

	ARM_Vector forFlows(nbStates);
	for(i=0;i<nbStates;++i)
		forFlows[i]= ( beta*(*forZcEnd)[i]*(forNotionalVector)[size-1] - alpha*(*forZcStart)[i] *(forNotionalVector)[0] ) * (*forexVector)[i];

	for(i=0;i<nbStates;++i){
		std::vector<double>& column = foreignLeg->GetColumn(i);
		forFlows[i] += (column->sum()*(*forexVector)[i]) ;
		delete column;
	}

	ARM_GP_MatrixPtr domesticLeg;
	if(isDomFlottant){
		
		/// This leg is computed in own currency and discounted by payModel
		domesticLeg = NPVSwapLeg( domCurveName, evalTime,domFwdStartTimes, domFwdEndTimes, domFwdPayPeriods,
					domPayTimes, domPayPeriods, domMarginVector,domNotionalVector,states);
	}else{
		/// fixed leg computing
		domesticLeg = NPVFixLeg( domCurveName, evalTime,domPayTimes, domPayPeriods, domNotionalVector, domStrikePerState, 1, states);
	}

	ARM_PricingFunctionIR* domIRFctor = dynamic_cast<ARM_PricingFunctionIR*>(&*(*itsModelMap)[domCurveName]->Model());


	/// computing notional exchange
	size = domPayTimes.size();
	ARM_VectorPtr domZcStart = GetDiscountFunctor()->DiscountFactor(domCurveName,evalTime,domFlowStartTimes[0],states);
    ARM_VectorPtr domZcEnd = GetDiscountFunctor()->DiscountFactor(domCurveName,evalTime,domFlowEndTimes[size-1],states);
	
	std::vector<double> domFlows(nbStates);
	for(i=0;i<nbStates;++i)
		domFlows[i]= ( beta*(*domZcEnd)[i]*(domNotionalVector)[size-1] - alpha*(*domZcStart)[i] *(domNotionalVector)[0] );

	for(i=0;i<nbStates;++i){
		std::vector<double>& column = domesticLeg->GetColumn(i);
		domFlows[i] += column->sum() ;
		delete column;
	}

	ARM_GP_VectorPtr values(new std::vector<double>(nbStates));
	for(i=0;i<nbStates;++i)
		(*values)[i] = payRec*( domFlows[i] - forFlows [i]);

	return values;
}
////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : VanillaCorridorlet
///	Returns : ARM_VectorPtr
///	Action  : call VanillaCorridorlet on the sub model
////////////////////////////////////////////////////

ARM_VectorPtr ARM_MultiAssetsModel::VanillaCorridorlet(
	const   string& modelName, 
	double  evalTime,
    double  payTime,
    double  resetTime,
    double  startTime,
    double  endTime,
    int     indexPaymentType,
    double  fwdPaymentPeriod,
    const std::vector<double>& refIdxResetTimes,
    const std::vector<double>& refIdxStartTimes,
    const std::vector<double>& refIdxEndTimes,
    const std::vector<double>& refFwdPeriods,
    const std::vector<double>& refIndexWeight,
    double  couponMargin,
    const vector<const std::vector<double>*> downBarrierPerState,
    const vector<const std::vector<double>*> upBarrierPerState,
    double  payNotional,
    int     capFloor,
    const   ARM_PricingStatesPtr& states) const
{
	ARM_PricingModelPtr model = (*itsModelMap)[modelName]->Model();
	ARM_PricingFunctionIR* IRFctor = dynamic_cast<ARM_PricingFunctionIR*>(&*model);
	if( !IRFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + modelName + " does not support IR function VanillaCorridorlet!" );
	
	return IRFctor->VanillaCorridorlet(	modelName, evalTime, payTime, resetTime, startTime, endTime, indexPaymentType,
		fwdPaymentPeriod, refIdxResetTimes, refIdxStartTimes, refIdxEndTimes, refFwdPeriods, refIndexWeight, couponMargin,
		downBarrierPerState, upBarrierPerState, payNotional, capFloor, states);
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : VanillaDigital
///	Returns : ARM_VectorPtr
///	Action  : call VanillaDigital on the sub model
////////////////////////////////////////////////////

ARM_VectorPtr ARM_MultiAssetsModel::VanillaDigital(
	const string& modelName, 
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
	ARM_PricingModelPtr model = (*itsModelMap)[modelName]->Model();
	ARM_PricingFunctionIR* IRFctor = dynamic_cast<ARM_PricingFunctionIR*>(&*model);
	if( !IRFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + modelName + " does not support IR function VanillaDigital!" );
	
	return IRFctor->VanillaDigital( modelName, evalTime, payTime, period, payNotional, fwdResetTime, fwdStartTime,
		fwdEndTime, fwdPeriod, strikesPerState, capFloor, states);
}



////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : Vanilla	
///	Returns : ARM_VectorPtr
///	Action  : call VanillaCaplet on the sub model
////////////////////////////////////////////////////

ARM_VectorPtr ARM_MultiAssetsModel::VanillaCaplet(
	const string& modelName, 
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
	ARM_PricingModelPtr model = (*itsModelMap)[modelName]->Model();
	ARM_PricingFunctionIR* IRFctor = dynamic_cast<ARM_PricingFunctionIR*>(&*model);
	if( !IRFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + modelName + " does not support IR function VanillaDigital!" );
	
	return IRFctor->VanillaCaplet(	modelName, evalTime, payTime, period, payNotional, fwdResetTime,
		fwdStartTime, fwdEndTime, fwdPeriod, strikesPerState, capFloor, states );
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : VanillaCaplet
///	Returns : ARM_VectorPtr
///	Action  : call VanillaCaplet on the sub model
////////////////////////////////////////////////////

ARM_VectorPtr ARM_MultiAssetsModel::VanillaSwaption(
	const string& modelName,
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
	const std::vector<double>& fixTimes,
	const std::vector<double>& fixPayPeriods,
	const ARM_GP_Matrix& strikesPerState,
	int callPut,
	const ARM_PricingStatesPtr& states,
	bool isConstantNotional,
	bool isConstantSpread,
	bool isConstantStrike) const
{
	ARM_PricingModelPtr model = (*itsModelMap)[modelName]->Model();
	ARM_PricingFunctionIR* IRFctor = dynamic_cast<ARM_PricingFunctionIR*>(&*model);
	if( !IRFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + modelName + " does not support IR function VanillaDigital!" );

	return IRFctor->VanillaSwaption( modelName, evalTime, swapResetTime, fixNotional, floatNotional,floatStartTime, floatEndTime,        
		floatResetTimes,floatStartTimes,floatEndTimes,floatIntTerms,fixTimes, fixPayPeriods, strikesPerState, callPut, states,isConstantNotional, isConstantSpread,isConstantStrike);
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : CPISpot
///	Returns : ARM_VectorPtr
///	Action  : call CPISpot on the sub model
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_MultiAssetsModel::CPISpot( const string& InfcurveName, double evalTime, double CPITime, string DCFLag, long DailyInterp, string ResetLag, const ARM_PricingStatesPtr& states) const
{
	ARM_PricingModelPtr model = (*itsModelMap)[InfcurveName]->Model();
	ARM_PricingFuncInflation* InfFctor = dynamic_cast<ARM_PricingFuncInflation*>(&*model);

	if( !InfFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + InfcurveName + " does not support Infaltion function CPISpot!" );

	return InfFctor->CPISpot(InfcurveName, evalTime, CPITime, DCFLag, DailyInterp, ResetLag, states);
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : CPIForward
///	Returns : ARM_VectorPtr
///	Action  : call CPIForward on the sub model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MultiAssetsModel::CPIForward(
		const string& InfcurveName, 
		double evalTime, 
		double CPITime,
		double FixingTime, 
		const ARM_PricingStatesPtr& states) const
{
	ARM_PricingModelPtr model = (*itsModelMap)[InfcurveName]->Model();
	ARM_PricingFuncInflation* InfFctor = dynamic_cast<ARM_PricingFuncInflation*>(&*model);

	if( !InfFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + InfcurveName + " does not support Infaltion function CPISpot!" );

	return InfFctor->CPIForward(InfcurveName, evalTime, CPITime, FixingTime, states);
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : ConvexityAdjustment
///	Returns : ARM_VectorPtr
///	Action  : call ConvexityAdjustment on the sub model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MultiAssetsModel::ConvexityAdjustment(
		const string& InfcurveName,
        double evalTime,
		double tenor,
		double maturityTime,
        const ARM_PricingStatesPtr& states) const
{
	ARM_PricingModelPtr model = (*itsModelMap)[InfcurveName]->Model();
	ARM_PricingFuncInflation* InfFctor = dynamic_cast<ARM_PricingFuncInflation*>(&*model);

	if( !InfFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + InfcurveName + " does not support Infaltion function ConvexityAdjustment!" );

	return InfFctor->ConvexityAdjustment(InfcurveName,evalTime,tenor,maturityTime,states);
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : ForwardCPIRatio
///	Returns : ARM_VectorPtr
///	Action  : call ForwardCPIRatio on the sub model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MultiAssetsModel::ForwardCPIRatio(
		const string& InfcurveName,
        double evalTime,
		double tenor,
		double CPITime,
		double maturityTime,
        const ARM_PricingStatesPtr& states) const
{
	ARM_PricingModelPtr model = (*itsModelMap)[InfcurveName]->Model();
	ARM_PricingFuncInflation* InfFctor = dynamic_cast<ARM_PricingFuncInflation*>(&*model);

	if( !InfFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + InfcurveName + " does not support Infaltion function ForwardCPIRatio!" );

	return InfFctor->ForwardCPIRatio(InfcurveName,evalTime,tenor,CPITime, maturityTime,states);
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : YoYSwapRate
///	Returns : ARM_VectorPtr
///	Action  : call YoYSwapRate on the sub model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MultiAssetsModel::YoYSwapRate(
	const string& irCurveName,
	const string& InfCurveName, 
	double evalTime,
	const ARM_DateStripPtr& numDateStrip,
	const ARM_DateStripPtr& denomDateStrip,
	const ARM_DateStripPtr& fixedDateStrip,
	double itsSpread,
	const ARM_PricingStatesPtr& states) const
{
	ARM_PricingModelPtr model = (*itsModelMap)[InfCurveName]->Model();
	ARM_PricingFuncInflation* InfFctor = dynamic_cast<ARM_PricingFuncInflation*>(&*model);

	if( !InfFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + InfCurveName + " does not support Infaltion function YoYSwapRate!" );

	return InfFctor->YoYSwapRate(irCurveName, InfCurveName, evalTime, numDateStrip, denomDateStrip, fixedDateStrip, itsSpread, states);
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : YoYSwap
///	Returns : ARM_VectorPtr
///	Action  : call YoYSwap on the sub model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MultiAssetsModel::YoYSwap(
	const string& irCurveName,
	const string& InfCurveName, 
	double evalTime,
	double Strike,
	double FloatMargin, 
	const ARM_DateStripPtr& numDateStrip,
	const ARM_DateStripPtr& denomDateStrip,
	const ARM_DateStripPtr& fixedDateStrip,
	double itsSpread,
	const ARM_PricingStatesPtr& states) const
{
	ARM_PricingModelPtr model = (*itsModelMap)[InfCurveName]->Model();
	ARM_PricingFuncInflation* InfFctor = dynamic_cast<ARM_PricingFuncInflation*>(&*model);

	if( !InfFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + InfCurveName + " does not support Infaltion function YoYSwap!" );

	return InfFctor->YoYSwap(irCurveName, InfCurveName, evalTime, Strike, FloatMargin, numDateStrip,denomDateStrip, fixedDateStrip, itsSpread, states);
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : OATSwapRate
///	Returns : ARM_VectorPtr
///	Action  : call OATSwapRate on the sub model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MultiAssetsModel::OATSwapRate(
	const string& irCurveName,
	const string& InfCurveName, 
	double evalTime,
	const ARM_DateStripPtr& numDateStrip,
	const ARM_DateStripPtr& denomDateStrip,
	const ARM_DateStripPtr& fixedDateStrip,
	double itsCoupon,
	const ARM_PricingStatesPtr& states) const
{
	ARM_PricingModelPtr model = (*itsModelMap)[InfCurveName]->Model();
	ARM_PricingFuncInflation* InfFctor = dynamic_cast<ARM_PricingFuncInflation*>(&*model);

	if( !InfFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + InfCurveName + " does not support Infaltion function OATSwapRate!" );

	return InfFctor->OATSwapRate(irCurveName, InfCurveName, evalTime, numDateStrip, denomDateStrip, fixedDateStrip, itsCoupon, states);
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : OATSwap
///	Returns : ARM_VectorPtr
///	Action  : call OATSwap on the sub model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MultiAssetsModel::OATSwap(
	const string& irCurveName,
	const string& InfCurveName, 
	double evalTime,
	double Strike,
	double FloatMargin, 
	const ARM_DateStripPtr& numDateStrip,
	const ARM_DateStripPtr& denomDateStrip,
	const ARM_DateStripPtr& fixedDateStrip,
	double itsCoupon,
	const ARM_PricingStatesPtr& states) const
{
	ARM_PricingModelPtr model = (*itsModelMap)[InfCurveName]->Model();
	ARM_PricingFuncInflation* InfFctor = dynamic_cast<ARM_PricingFuncInflation*>(&*model);

	if( !InfFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + InfCurveName + " does not support Infaltion function OATSwap!" );

	return InfFctor->OATSwap(irCurveName, InfCurveName, evalTime, Strike, FloatMargin, numDateStrip,denomDateStrip, fixedDateStrip, itsCoupon,states);
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : YoYCap
///	Returns : ARM_VectorPtr
///	Action  : call YoYCap on the sub model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MultiAssetsModel::YoYCapFloor(const string& irCurveName, const string& infCurveName, 
	double evalTime,
	double Strike,
	double FloatMargin, 
	int CapFloor,
	const ARM_DateStripPtr& numDateStrip,
	const ARM_DateStripPtr& denomDateStrip,
	double itsSpread,
	const ARM_PricingStatesPtr& states) const
{
	ARM_PricingModelPtr model = (*itsModelMap)[infCurveName]->Model();
	ARM_PricingFuncInflation* InfFctor = dynamic_cast<ARM_PricingFuncInflation*>(&*model);

	if( !InfFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + infCurveName + " does not support Infaltion function YoYCapFloor!" );

	return InfFctor->YoYCapFloor(irCurveName, infCurveName, evalTime, Strike, FloatMargin, CapFloor, numDateStrip,denomDateStrip, itsSpread,states);
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : YoYCap
///	Returns : ARM_VectorPtr
///	Action  : call YoYCap on the sub model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MultiAssetsModel::OATCapFloor(const string& irCurveName, const string& infCurveName, 
	double evalTime,
	double Strike,
	double FloatMargin, 
	int CapFloor,
	const ARM_DateStripPtr& numDateStrip,
	const ARM_DateStripPtr& denomDateStrip,
	double itsSpread,
	const ARM_PricingStatesPtr& states) const
{
	ARM_PricingModelPtr model = (*itsModelMap)[infCurveName]->Model();
	ARM_PricingFuncInflation* InfFctor = dynamic_cast<ARM_PricingFuncInflation*>(&*model);

	if( !InfFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + infCurveName + " does not support Infaltion function OATCapFloor!" );

	return InfFctor->OATCapFloor(irCurveName, infCurveName, evalTime, Strike, FloatMargin, CapFloor, numDateStrip,denomDateStrip, itsSpread,states);
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : ZCCap
///	Returns : ARM_VectorPtr
///	Action  : call ZCCap on the sub model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MultiAssetsModel::ZCCap( ) const
{
	return ARM_GP_VectorPtr( new std::vector<double>( 0 ) );
}


////////////////////////////////////////////////////
///	Class  : ARM_MultiAssetsModel
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : call MCModelStatesFromToNextTime on the sub models
////////////////////////////////////////////////////

void ARM_MultiAssetsModel::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
	/// 1) does the rotation according to the rotation
	ARM_GP_MatrixPtr& numMethodStates = states->GetNumMethodStates();
	ARM_GP_MatrixPtr& modelStates = states->GetModelStates();
	size_t bucketSize =  states->size();

#ifdef __GP_STRICT_VALIDATION

	if( numMethodStates->rows() != FactorCount() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "numMethodStates->rows() != FactorCount()!" );
	/*if( modelStates->rows() != FactorCount() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelStates->rows() != FactorCount()!" );*/
	if( modelStates->cols() != bucketSize )
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelStates->cols() != bucketSize!" );
	if( numMethodStates->cols() != bucketSize )
		ARM_THROW( ERR_INVALID_ARGUMENT, "numMethodStates->cols() != bucketSize!" );

#endif

	/// delegates to individual models

	ARM_ModelNameMap::const_iterator 
		iter, 
		start	= itsModelMap->begin(), 
		end		= itsModelMap->end();

	for( iter=start; iter!=end; itsModelMap->getNextUsedIter(iter) )
	{
		(*iter).Model()->MCModelStatesFromToNextTime( states, timeIndex);
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : TreeStatesToModelStates
///	Returns : void
///	Action  : call TreeStatesToModelStates on the sub models
////////////////////////////////////////////////////

void ARM_MultiAssetsModel::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{
	/// delegates to individual models
	for( ARM_ModelNameMap::iterator iter=itsModelMap->begin(); iter!=itsModelMap->end(); itsModelMap->getNextUsedIter(iter) )
		(*iter).Model()->TreeStatesToModelStates( states, timeIndex);
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : FirstPricingStates
///	Returns : ARM_VectorPtr
///	Action  : computes the first pricing states
////////////////////////////////////////////////////

ARM_PricingStatesPtr ARM_MultiAssetsModel::FirstPricingStates( size_t bucketSize ) const
{
	ARM_PricingStatesPtr result( new ARM_PricingStates(bucketSize, ModelStatesSize(),0,FactorCount() ) );
	ARM_PricingStatesPtr tmpResult;
	size_t modelNb,modelFactorCount, modelStatesSize;
	size_t i,j,k;

	k=0;
	/// delegates to individual models
	for( ARM_ModelNameMap::iterator iter=itsModelMap->begin(); iter!=itsModelMap->end(); itsModelMap->getNextUsedIter(iter) )
	{
		modelNb			= (*iter).Model()->GetModelNb();
		modelFactorCount= (*iter).Model()->FactorCount();

		if(modelFactorCount == 0) continue;

		tmpResult= (*iter).Model()->FirstPricingStates( bucketSize );
		modelStatesSize = tmpResult->ModelStatesSize();

		result->AddPricingStatesContextVector( tmpResult->GetPricingStatesContextVector() );

#ifdef __GP_STRICT_VALIDATION
		if( tmpResult->NumMethodStatesSize() != modelFactorCount )
			ARM_THROW( ERR_INVALID_ARGUMENT, "tmpResult->NumMethodStatesSize() != modelFactorCount " );
#endif

		for( i=0; i<bucketSize; ++i )
		{
			for( j=0; j<modelFactorCount; ++j )
				result->SetNumMethodState( i, modelNb+j, tmpResult->GetNumMethodState( i, j ) );

			for( j=0; j<tmpResult->ModelStatesSize(); ++j )
				result->SetModelState( i, k+j, tmpResult->GetModelState( i, j ) );
		}
		k += modelStatesSize;
	}
	return result;
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : Induct
///	Returns : ARM_VectorPtr
///	Action  : induct the states from the current time to next time
////////////////////////////////////////////////////

ARM_PricingStatesPtr ARM_MultiAssetsModel::Induct(ARM_PricingStatesPtr& states,double toTime)
{
	/// standard induct!
	return ARM_PricingModel::Induct( states, toTime );
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : ComputeModelTimes
///	Returns : std::vector<double>&
///	Action  : computes the model times to stop at
////////////////////////////////////////////////////

std::vector<double>& ARM_MultiAssetsModel::ComputeModelTimes( const ARM_TimeInfoPtrVector& timeInfos )
{
	ARM_ModelNameMap::const_iterator 
		iter, 
		start	= itsModelMap->begin(), 
		end		= itsModelMap->end();

	/// 1) collect all the model times
	size_t usedModelsSize = itsModelMap->UsedModelsSize();
	std::vector<double>& tmpModelTimes=NULL;
	std::vector<double>& result = new std::vector<double>(0);

	for( iter=start; iter!= end; itsModelMap->getNextUsedIter(iter) )
	{
		tmpModelTimes = (*iter).Model()->ComputeModelTimes( timeInfos );
		if( tmpModelTimes )
		{
			size_t start = result->size();
			result->resize( result->size() + tmpModelTimes->size() );
			for( size_t i=0; i<tmpModelTimes->size(); ++i )
				(*result)[i+start] = tmpModelTimes->Elt(i);
			delete tmpModelTimes;
		}
	}

	/// 2) sort and remove duplicates
	result->sort();
	result->unique();

	return result;
}
////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : ComputeNumeraireTimes
///	Returns : std::vector<double>&
///	Action  : computes the model times to stop at
////////////////////////////////////////////////////

ARM_VectorPtr ARM_MultiAssetsModel::ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const
{
	ARM_ModelNameMap::const_iterator 
		iter, 
		start	= itsModelMap->begin(), 
		end		= itsModelMap->end();

	/// 1) collect all the model times
	size_t usedModelsSize = itsModelMap->UsedModelsSize();
	ARM_VectorPtr tmpModelTimes=ARM_VectorPtr(NULL);
	std::vector<double>& result = new std::vector<double>(0);

	for( iter=start; iter!= end; itsModelMap->getNextUsedIter(iter) )
	{
		tmpModelTimes = (*iter).Model()->ComputeNumeraireTimes( timeInfos );
		if( tmpModelTimes !=ARM_VectorPtr(NULL))
		{
			size_t start = result->size();
			result->resize( result->size() + tmpModelTimes->size() );
			for( size_t i=0; i<tmpModelTimes->size(); ++i )
				(*result)[i+start] = tmpModelTimes->Elt(i);
		}
	}

	/// 2) sort and remove duplicates
	if(result->size()>0)
	{
		result->sort();
		result->unique();
	}
	else
    {
        /// Free memory !
        delete result;
        result=NULL;
    }

    return ARM_VectorPtr(result);
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : SupportBackwardInduction
///	Returns : bool
///	Action  : tells whether it supports BackwardInduction
///				if any of the model does not support the backward induction
///				return false, else return true
////////////////////////////////////////////////////

bool ARM_MultiAssetsModel::SupportBackwardInduction() const
{
	ARM_ModelNameMap::const_iterator 
		iter, 
		start	= itsModelMap->begin(), 
		end		= itsModelMap->end();

	for( iter=start; iter!=end; itsModelMap->getNextUsedIter(iter) )
	{
		if( !(*iter).Model()->SupportBackwardInduction() )
			return false;
	}
	
	return true;
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : SupportForwardInduction
///	Returns : bool
///	Action  : tells whether it support ForwardInduction
///				if any of the model does not support the forward induction
///				return false, else return true
////////////////////////////////////////////////////

bool ARM_MultiAssetsModel::SupportForwardInduction() const
{
	ARM_ModelNameMap::const_iterator 
		iter, 
		start	= itsModelMap->begin(), 
		end		= itsModelMap->end();

	for( iter=start; iter!= end; itsModelMap->getNextUsedIter(iter) )
	{
		if( !(*iter).Model()->SupportForwardInduction() )
			return false;
	}
	
	return true;
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : SupportAnalyticMarginal
///	Returns : bool
///	Action  : tells whether it support AnalyticMarginal
///				if any of the model does not support AnalyticMarginal
///				return false, else return true
////////////////////////////////////////////////////

bool ARM_MultiAssetsModel::SupportAnalyticMarginal() const
{
	ARM_ModelNameMap::const_iterator 
		iter, 
		start	= itsModelMap->begin(), 
		end		= itsModelMap->end();

	for( iter=start; iter!= end; itsModelMap->getNextUsedIter(iter) )
	{
		if( !(*iter).Model()->SupportAnalyticMarginal() )
			return false;
	}
	
	return true;
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : BackwardLookingInit
///	Returns : ARM_PricingStatesPtr
///	Action  : Default initialisation for backward looking
///           numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_MultiAssetsModel::BackwardLookingInit(ARM_NumMethodPtr& numMethod, double firstInductTime)
{
    /// Delegate to numerical method
	return numMethod->Init(*this,firstInductTime);
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : ForwardLookingInit
///	Returns : ARM_PricingStatesPtr
///	Action  : Default initialisation for forward looking
///           numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_MultiAssetsModel::ForwardLookingInit(ARM_NumMethodPtr& numMethod, double firstInductTime)
{
	ARM_ModelNameMap::iterator 
		iter,
		start	= itsModelMap->begin(), 
		end		= itsModelMap->end();

	/// first notify the models about their indexes
	for( iter=start ; iter!=end; itsModelMap->getNextUsedIter(iter) )
	{
		numMethod->ComputeAndSetTimeSteps(*(*iter).Model());
	}

    /// Delegate to numerical method
	return numMethod->Init(*this,firstInductTime);
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: InitModelNb
///	Returns :
///	Action  : Initialize the model nb
////////////////////////////////////////////////////

void ARM_MultiAssetsModel::InitModelNb()
{
	ARM_ModelNameMap::iterator 
		iter,
		start	= itsModelMap->begin(), 
		end		= itsModelMap->end();
	size_t index = 0;
	size_t rank = 0;
	
	/// first notify the models about their indexes
	for( iter=start ; iter!=end; itsModelMap->getNextUsedIter(iter) )
	{
		(*iter).Model()->SetModelNb( index );
		index += (*iter).Model()->FactorCount();
		(*iter).Model()->SetModelRank( rank );
		rank++;
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : Init
///	Returns : ARM_PricingStatesPtr
///	Action  : intialise the model
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_MultiAssetsModel::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
	/// get the reference model and set the zero curve
	itsRefModel = &*(*itsModelMap)[payModelName]->Model();
	SetZeroCurve( itsRefModel->GetZeroCurve() );

	//itsRefModel->Init( payModelName, timeInfos );

	/// then std init
	int nbEvents=timeInfos.size();
    bool isSpotUse = nbEvents == 0 || (nbEvents==1 && timeInfos[0]->GetEventTime() <= K_NEW_DOUBLE_TOL);

	InitModelNb();

	ARM_ModelNameMap::iterator 
		iter,
		start	= itsModelMap->begin(), 
		end		= itsModelMap->end();

	/// second set time steps and delegate to the numerical method!
    if(!isSpotUse)
    {
		/// Numerical method and numeraire are needed to go on
        ARM_NumMethodPtr numMethod=GetNumMethod();
		if( numMethod == ARM_NumMethodPtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numerical method not set in multi asset model!");
		
		/// initialize the numeraire
		ARM_NumerairePtr numeraire = GetNumeraire();
		if( numeraire == ARM_NumerairePtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numeraire is not set in the multi-asset model!");
		
		/// creates the model schedule (smart pointor for exception safety!)
		std::vector<double>& ptimeSteps = new std::vector<double>(0), 
			*monoAssettimeSteps = NULL, *previousTimeSteps = NULL;
		
		/// Computes the merged timeSteps of all the models
		for( iter=start; iter!=end;iter++)
		{
            if(iter->UsedInPricing())
            {
			    previousTimeSteps = ptimeSteps;
			    monoAssettimeSteps = iter->Model()->PricingTimeSteps(timeInfos);
			    ptimeSteps = MergeSortedVectorNoDuplicates( *previousTimeSteps, *monoAssettimeSteps );
			    delete previousTimeSteps;
			    delete monoAssettimeSteps;
            }
		}

		/// Set the basic schedule in the numerical method
		numMethod->SetTimeSteps(*ptimeSteps);
		delete ptimeSteps;

		double firstInductTime = timeInfos[0]->GetEventTime();
		CC_NS(std,auto_ptr)<std::vector<double>> pModelTimes( (*this).ComputeModelTimes( timeInfos ) );
		numeraire->Init(*(GetDiscountFunctor()),payModelName,timeInfos,(*this).ComputeNumeraireTimes( timeInfos ));

		/// ...initialise it
		if( ARM_NumMethod::GP_BCKWDLOOKING == numMethod->GetPricingDirection() )
            return BackwardLookingInit(numMethod,firstInductTime);

		else
            return ForwardLookingInit(numMethod,firstInductTime);
	}
	else
	{
		/// analytic part
		ARM_PricingStatesPtr totalStates( new ARM_PricingStates(1,FactorCount(),0) );
		size_t j=0;
		for( iter=start ; iter!=end; itsModelMap->getNextUsedIter(iter) )
		{
			ARM_PricingStatesPtr modelStates = (*iter).Model()->Init( payModelName, timeInfos );
			for( size_t i=0; i<(*iter).Model()->FactorCount(); ++i, ++j )
			{
				double currentModelState = modelStates->GetModelState(0,i);
				totalStates->SetModelState(0,j,currentModelState);
			}
		}
		return totalStates;
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : ReInit
///	Returns : ARM_PricingStatesPtr
///	Action  : Re-initialise the model
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_MultiAssetsModel::ReInit()
{
	return ARM_PricingModel::ReInit();
}



////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : PostInit
///	Returns : ARM_PricingStatesPtr
///	Action  : PostInit the model
////////////////////////////////////////////////////
void ARM_MultiAssetsModel::PostInit()
{
	/// delegates to individual models
	for( ARM_ModelNameMap::iterator iter=itsModelMap->begin(); iter!=itsModelMap->end(); itsModelMap->getNextUsedIter(iter) )
		(*iter).Model()->PostInit();
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : double
///	Returns : ARM_PricingStatesPtr
///	Action  : Computes the VarianceToTime
////////////////////////////////////////////////////

double ARM_MultiAssetsModel::VarianceToTime(double var,double minTime,double maxTime) const
{
	double time = 0.0;
	for( ARM_ModelNameMap::iterator iter=itsModelMap->begin(); iter!=itsModelMap->end(); itsModelMap->getNextUsedIter(iter) )	
		time =(*iter).Model()->VarianceToTime(var,minTime,maxTime);
	return time;
	//ARM_THROW( ERR_INVALID_ARGUMENT, "unimplemented VarianceToTime!" );
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: UnderlyingCorrelation
///	Returns : void
///	Action  : No default implementation
////////////////////////////////////////////////////

double ARM_MultiAssetsModel::UnderlyingCorrelation(	   string	underlyingType,
													   double	fromTime,
													   double	toTime,
													   double	startTime1,
													   double   endTime1,
													   double	startTime2,
													   double   endTime2,
													   double	startTime3,
													   double   endTime3,
													   double	startTime4,
													   double   endTime4) const

{   
	if( underlyingType == "CPIFwd" )
	{
		ARM_PricingFuncInflation* InfFctor = NULL, *InfFctor2 = NULL;

		for( ARM_ModelNameMap::iterator iter=itsModelMap->begin(); iter!=itsModelMap->end(); itsModelMap->getNextUsedIter(iter) )
		{
			ARM_PricingFuncInflation* InfFctor2 = dynamic_cast<ARM_PricingFuncInflation*>(&*((*iter).Model()));
			if( InfFctor2 )
				if( !InfFctor )
					InfFctor = InfFctor2;
				else
					ARM_THROW( ERR_INVALID_ARGUMENT, "ambiguous call: 2 or more models are possible to compute this. " );
		}

		if( InfFctor )
			return UnderlyingCorrelation(underlyingType,fromTime,toTime,startTime1,endTime1,startTime2,endTime2,startTime3,endTime3,startTime4,endTime4);
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, "unimplemented UnderlyingCorrelation!" );

	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "unimplemented UnderlyingCorrelation!" );
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: UnderlyingCovariance
///	Returns : void
///	Action  : No default implementation
////////////////////////////////////////////////////

double ARM_MultiAssetsModel::UnderlyingCovariance(	   string	underlyingType,
													   double	fromTime,
													   double	toTime,
													   double	startTime1,
													   double   endTime1,
													   double	startTime2,
													   double   endTime2,
													   double	startTime3,
													   double   endTime3,
													   double	startTime4,
													   double   endTime4) const

{   
	if( underlyingType == "CPIFwd" )
	{
		ARM_PricingFuncInflation* InfFctor = NULL, *InfFctor2 = NULL;

		for( ARM_ModelNameMap::iterator iter=itsModelMap->begin(); iter!=itsModelMap->end(); itsModelMap->getNextUsedIter(iter) )
		{
			ARM_PricingFuncInflation* InfFctor2 = dynamic_cast<ARM_PricingFuncInflation*>(&*((*iter).Model()));
			if( InfFctor2 )
				if( !InfFctor )
					InfFctor = InfFctor2;
				else
					ARM_THROW( ERR_INVALID_ARGUMENT, "ambiguous call: 2 ore more models are possible to compute this. " );
		}

		if( InfFctor )
			return UnderlyingCovariance(underlyingType,fromTime,toTime,startTime1,endTime1,startTime2,endTime2,startTime3,endTime3,startTime4,endTime4);
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, "unimplemented UnderlyingCovariance!" );

	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "unimplemented UnderlyingCovariance!" );
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: VanillaSpreadOption
///	Returns : void
///	Action  : No default implementation
////////////////////////////////////////////////////
ARM_VectorPtr  ARM_MultiAssetsModel::VanillaSpreadOptionLet(const string& curveName,
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
	ARM_PricingModelPtr model = (*itsModelMap)[curveName]->Model();
	ARM_PricingFunctionIR* IRFctor = dynamic_cast<ARM_PricingFunctionIR*>(&*model);
	if( !IRFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + curveName + " does not support IR function VanillaDigital!" );
	
	return IRFctor->VanillaSpreadOptionLet(	curveName, 
											evalTime, 
											callPut, 
											startTime, 
											endTime, 
											resetTime, 
											payTime, 
											payPeriod, 
											notional, 
											coeffLong, 
											coeffShort, 
											strikes, 
											swapLongFloatStartTime,
											swapLongFloatEndTime,
											swapLongFixPayTimes,
											swapLongFixPayPeriods,
											swapShortFloatStartTime,
											swapShortFloatEndTime,
											swapShortFixPayTimes,
											swapShortFixPayPeriods,
											states) ;
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: VanillaSpreadOption
///	Returns : void
///	Action  : No default implementation
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MultiAssetsModel::VanillaCMSCorridorlet(
		const string& curveName,
		double evalTime,
		double payTime,
		double resetTime,
		double startTime,
		double endTime,
		const std::vector<double>& refIdxResettimes,
		const std::vector<double>& refIndexWeights,
		const std::vector<double>& coeff1,
		const ARM_SwapRatePtrVector& firstIndex,
		const std::vector<double>& coeff2,
		const ARM_SwapRatePtrVector& secondIndex,
		int		payIndexType,			/// K_FIXED, K_LIBOR or K_CMS
		double	coupon,					/// in case of a fixed payment (K_FIXED)
		const	ARM_SwapRate& payRate,	/// rate description (K_LIBOR or K_CMS)	
		double  payIndexLeverage,
		const std::vector<double>& downBarriers,
        const std::vector<double>& upBarriers,
        double  payNotional,
        int     rcvPay,
		const ARM_SwapRatePtrVector& thirdIndex, // Single rate index for double condition
		const std::vector<double>& downBarriers3,
		const std::vector<double>& upBarriers3,
        const   ARM_PricingStatesPtr& states) const
{
	ARM_PricingModelPtr model = (*itsModelMap)[curveName]->Model();
	ARM_PricingFunctionIR* IRFctor = dynamic_cast<ARM_PricingFunctionIR*>(&*model);
	if( !IRFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + curveName + " does not support IR function VanillaDigital!" );
	
	return IRFctor->VanillaCMSCorridorlet(
				curveName,
				evalTime,
				payTime,
				resetTime,
				startTime,
				endTime,
				refIdxResettimes,
				refIndexWeights,
				coeff1,
				firstIndex,
				coeff2,
				secondIndex,
				payIndexType,
				coupon,
				payRate,
				payIndexLeverage,
				downBarriers,
				upBarriers,
				payNotional,
				rcvPay,
				thirdIndex,
				downBarriers3,
				upBarriers3,
				states);
}

////////////////////////////////////////////////////
///	Class  : ARM_MultiAssetsModel
///	Routine: RangeAccrualVectorial
///	Returns: void (go to the local normal model)
///	Action : compute the corridor double condition: one condition
///				on FX, one condition on Libor, both fixing 
///				at the same date
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MultiAssetsModel::RangeAccrualVectorial(
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
	ARM_PricingModelPtr model = (*itsModelMap)[curveName]->Model();
	ARM_PricingFunctionEquity* EQFctor = dynamic_cast<ARM_PricingFunctionEquity*>(&*model);
	if( !EQFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + curveName + " does not support IR function Range Accrual Vectorial!" );
	
	return EQFctor->RangeAccrualVectorial(curveName,
										evalTime,
										startTime,
										endTime,
										payTime,
										fixingTimes,
										payIndexType, 
										payIndexTerm,
										fxModelName,
										irIndexType, 
										irIndexResetTimes,
										irIndexStartTimes,
										irIndexEndTimes,
										irIndexTerms,
										fxDownBarriers,
										fxUpBarriers,
										irDownBarriers,
										irUpBarriers,
										notionals,
										states,
										eachFixingPrices,
										context);
}



////////////////////////////////////////////////////
///	Class  : ARM_Local_Normal_Model
///	Routine: DoubleDigital
///	Returns: double digital condition values
///	Action : Implicitly option expiry = eval date
///			 Computes the double digital price on
///			 both rates given their values
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MultiAssetsModel::DoubleDigital(
		const string& modelName, 
		double evalTime,
		const ARM_VectorPtr& firstRate,
        const std::vector<double>& firstStrikeDown,
        const std::vector<double>& firstStrikeUp,
		double firstStrikeSpread,
		const ARM_VectorPtr& secondRate,
        const std::vector<double>& secondStrikeDown,
        const std::vector<double>& secondStrikeUp,
		double secondStrikeSpread,
        const ARM_PricingStatesPtr& states) const
{
	ARM_PricingModelPtr model = (*itsModelMap)[modelName]->Model();
	ARM_PricingFunctionIR* IRFctor = dynamic_cast<ARM_PricingFunctionIR*>(&*model);
	if( !IRFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + modelName + " does not support IR function VanillaDigital!" );
	
	return IRFctor->DoubleDigital(modelName,evalTime,
		firstRate,firstStrikeDown,firstStrikeUp,firstStrikeSpread,
		secondRate,secondStrikeDown,secondStrikeUp,secondStrikeSpread,states);
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : GetCurrency
///	Returns : ARM_Currency*
///	Action  : Gets the currency corresponding to the modelName
////////////////////////////////////////////////////

ARM_Currency* ARM_MultiAssetsModel::GetCurrency( const string& modelName  ) const
{
	return (*itsModelMap)[modelName]->Model()->GetZeroCurve()->GetCurrencyUnit(); 	
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : ProcessPaidPayoffs
///	Returns : void
///	Action  : ProcessPaidPayoffs
////////////////////////////////////////////////////
void ARM_MultiAssetsModel::ProcessPaidPayoffs(const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, const ARM_PricingStatesPtr& states ) const
{
	itsRefModel->ProcessPaidPayoffs( payModelName, payoffs, evalTime, states );
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : ProcessUnPaidPayoffs
///	Returns : void
///	Action  : ProcessUnPaidPayoffs
////////////////////////////////////////////////////
void ARM_MultiAssetsModel::ProcessUnPaidPayoffs(const string& payModelName, ARM_VectorPtr& payoffs, double evalTime, const ARM_PricingStatesPtr& states ) const
{
	itsRefModel->ProcessUnPaidPayoffs( payModelName, payoffs, evalTime, states );
}



////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : NeedArrowDebreuPrices
///	Returns : bool
///	Action  : teels whether it requires arrow debreu prices!
////////////////////////////////////////////////////

bool ARM_MultiAssetsModel::NeedArrowDebreuPrices() const
{
	ARM_ModelNameMap::const_iterator 
		iter, 
		start	= itsModelMap->begin(), 
		end		= itsModelMap->end();

	for( iter=start; iter!= end; itsModelMap->getNextUsedIter(iter) )
	{
		if( !(*iter).Model()->NeedArrowDebreuPrices() )
			return false;
	}
	
	return true;
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : NeedsToCholeskyDecomposeFactors
///	Returns : bool
///	Action  : teels whether it Needs To CholeskyDecompose Factors
////////////////////////////////////////////////////

bool ARM_MultiAssetsModel::NeedsToCholeskyDecomposeFactors( ) const
{
	ARM_ModelNameMap::const_iterator 
		iter, 
		start	= itsModelMap->begin(), 
		end		= itsModelMap->end();

	for( iter=start; iter!= end; itsModelMap->getNextUsedIter(iter) )
	{
		if( !(*iter).Model()->NeedsToCholeskyDecomposeFactors() )
			return false;
	}
	return true;
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : NeedMCIntegProcess
///	Returns : bool
///	Action  : to say how MC method give simulated
///           processes
////////////////////////////////////////////////////

ARM_BoolVector ARM_MultiAssetsModel::NeedMCIntegProcess() const
{
	ARM_BoolVector ret(FactorCount());

	ARM_ModelNameMap::const_iterator 
		iter, 
		start	= itsModelMap->begin(), 
		end		= itsModelMap->end();

	ARM_BoolVector tmpFlags;
	int tmpCount,tmpModelNb;
	size_t i;

	for( iter=start; iter!= end; itsModelMap->getNextUsedIter(iter) )
	{
		tmpCount = (*iter).Model()->FactorCount();
		tmpFlags = (*iter).Model()->NeedMCIntegProcess();
		if(tmpFlags.empty())
			tmpFlags = ARM_BoolVector(tmpCount, false);
		
		tmpModelNb = (*iter).Model()->GetModelNb();
		for (i = 0; i < tmpCount; ++i)
			ret[tmpModelNb+i] = tmpFlags[i];
	}
	return ret;
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : NumMethodStateLocalVariances
///	Returns : void
///	Action  : compute NumMethodStateLocalVariances
////////////////////////////////////////////////////
void ARM_MultiAssetsModel::NumMethodStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const
{
	size_t nbSteps	= timeSteps.size();
	size_t factorNb = FactorCount();

	/// aggregates the local variances with the use of the correlation!
#if defined(__GP_STRICT_VALIDATION)
	if( localVariances.size()!= 0 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "localVariances.size()!= 0" );
#endif

	localVariances.resize(nbSteps-1);

	size_t i,j,k;
	for( i=0; i<nbSteps-1; ++i )
		localVariances[i] = new ARM_GP_TriangularMatrix(factorNb,0.0);

	ARM_MatrixVector tmpLocalVariances;
	tmpLocalVariances.reserve((nbSteps-1)*factorNb);
	size_t modelNb, modelFactor;
	

	/// delegates to individual models and get the corresponding variances and local variances
	for( ARM_ModelNameMap::iterator iter=itsModelMap->begin(); iter!=itsModelMap->end(); itsModelMap->getNextUsedIter(iter) )
	{
		(*iter).Model()->NumMethodStateLocalVariances( timeSteps, tmpLocalVariances );
		modelNb		= (*iter).Model()->GetModelNb();
		modelFactor = (*iter).Model()->FactorCount();
		
		for( i=0; i<nbSteps-1; ++i )
		{
			for( j=0; j<modelFactor; ++j )
				for( k=0; k<modelFactor; ++k )
					localVariances[i]->Elt(modelNb+j,modelNb+k) = tmpLocalVariances[(nbSteps-1)*modelNb+i]->Elt(j,k);
		}
	}

	/// destroy the tmpLocalVariances
	DeletePointorVector<ARM_GP_Matrix>( tmpLocalVariances );
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : NumMethodStateGlobalVariances
///	Returns : void
///	Action  : compute the nummethod states global variances
///				To be improved!!!
////////////////////////////////////////////////////
void ARM_MultiAssetsModel::NumMethodStateGlobalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& variances ) const
{
	size_t nbSteps	= timeSteps.size();
	size_t factorNb = FactorCount();

	/// aggregates the local and variances with the use of the correlation!
#if defined(__GP_STRICT_VALIDATION)
	if( variances.size()!= 0 ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "variances.size()!= 0" );
#endif

	variances.resize(nbSteps);

	size_t i,j,k;
	for( i=0; i<nbSteps-1; ++i )
	{
		variances[i] = new ARM_GP_TriangularMatrix(factorNb,0.0);
	}
	variances[i] = new ARM_GP_TriangularMatrix(factorNb,0.0);

	ARM_MatrixVector tmpVariances;
	tmpVariances.reserve(nbSteps*factorNb);
	size_t modelNb, modelFactor;
	ARM_GP_T_Matrix<bool> hasBeenSetted (factorNb, factorNb, false);

	/// delegates to individual models and get the corresponding variances and local variances
	for( ARM_ModelNameMap::iterator iter=itsModelMap->begin(); iter!=itsModelMap->end(); itsModelMap->getNextUsedIter(iter) )
	{
		(*iter).Model()->NumMethodStateGlobalVariances( timeSteps, tmpVariances );
		modelNb		= (*iter).Model()->GetModelNb();
		modelFactor = (*iter).Model()->FactorCount();
		
		for( i=0; i<nbSteps-1; ++i )
		{
			for( j=0; j<modelFactor; ++j )
			{
				for( k=0; k<modelFactor; ++k )
				{
					variances[i]->Elt(modelNb+j,modelNb+k)		= tmpVariances[nbSteps*modelNb+i]->Elt(j,k);
				}
			}
		}

		/// last but not least, take the element of variances
		for( j=0; j<modelFactor; ++j )
			for( k=0; k<modelFactor; ++k )
				variances[i]->Elt(modelNb+j,modelNb+k) = tmpVariances[nbSteps*modelNb+i]->Elt(j,k);
	}

	/// destroy the tmpVariances and tmpLocalVariances
	DeletePointorVector<ARM_GP_Matrix>( tmpVariances );
}



////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : VolatilitiesAndCorrelations
///	Returns : void
///	Action  : compute the volatilities and correlation
////////////////////////////////////////////////////

void ARM_MultiAssetsModel::VolatilitiesAndCorrelations( const std::vector<double>& timeSteps, 
	ARM_GP_MatrixPtr& vols,
	ARM_GP_MatrixPtr& d1Vols,
	ARM_GP_MatrixPtr& correls,
	bool linearVol) const
{
	ARM_GP_MatrixPtr tmpVols;
	ARM_GP_MatrixPtr tmpD1Vols;
	ARM_GP_MatrixPtr tmpCorrels;

	size_t factorNb = FactorCount();
	vols = 	ARM_GP_MatrixPtr( new ARM_GP_Matrix(factorNb,timeSteps.size(),0.0) );
	d1Vols	= ARM_GP_MatrixPtr( new ARM_GP_Matrix(factorNb,timeSteps.size(),0.0) );
	correls	= ARM_GP_MatrixPtr( new ARM_GP_Matrix(factorNb*(factorNb-1)/2,timeSteps.size(),0.0) );
	size_t i,j,k,
		modelNb,
		modelFactorCount,
		timeStepsNb = timeSteps.size();;
	
	/// delegates to individual models
	for( ARM_ModelNameMap::iterator iter=itsModelMap->begin(); iter!=itsModelMap->end(); itsModelMap->getNextUsedIter(iter) )
	{
		(*iter).Model()->VolatilitiesAndCorrelations( timeSteps, tmpVols, tmpD1Vols, tmpCorrels, linearVol );

		modelNb			= (*iter).Model()->GetModelNb();
		modelFactorCount= (*iter).Model()->FactorCount();

		/// copy to the 
		for( i=0; i<modelFactorCount; ++i )
		{
			for( j=0; j<timeStepsNb; ++j )
			{
				(*vols)(modelNb+i,j)  = (*tmpVols)(i,j);
				
				if ( !tmpD1Vols.IsNull() ) 
					(*d1Vols)(modelNb+i,j)= (*tmpD1Vols)(i,j);
			}
		}
	}

	/// creates the correlation matrix from the used correl matrix
	size_t offset=0,offseti=0;

	ARM_GP_Matrix correlMatrix;

	for( i=0; i<factorNb; ++i )
	{
		for( j=0; j<timeStepsNb; ++j )
		{
			correlMatrix = GetCorrelMatrix()->Interpolate(timeSteps[i]);
			for( k=i+1,offset=offseti; k<factorNb; ++k,++offset )
				(*correls)(offset,j) = correlMatrix(i,k);
		}

		offseti += factorNb-1-i;
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : ModelStateLocalVariances
///	Returns : void
///	Action  : computes ModelStateLocalVariances
////////////////////////////////////////////////////

void ARM_MultiAssetsModel::ModelStateLocalVariances( const std::vector<double>& timeSteps,
	ARM_MatrixVector& localVariances ) const
{
	size_t nbSteps = timeSteps.size();
	localVariances.reserve( (nbSteps-1)*itsModelMap->UsedModelsSize());	

	ARM_ModelNameMap::iterator 
		iter,
		start	= itsModelMap->begin(), 
		end		= itsModelMap->end();

	/// first notify the models about their indexes
	for( iter=start ; iter!=end; itsModelMap->getNextUsedIter(iter) )
		(*iter).Model()->ModelStateLocalVariances(timeSteps,localVariances);
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : SetModelStateLocalStdDevs
///	Returns : void
///	Action  : Sets Model State Local Std Devs
///////////////////////////////////////////////////

void ARM_MultiAssetsModel::SetModelStateLocalStdDevs( const ARM_MatrixVector& stateLocalStdDevs, bool shareStateLocalStdDevs )
{
	/// delegates to individual models
	for( ARM_ModelNameMap::iterator iter=itsModelMap->begin(); iter!=itsModelMap->end(); itsModelMap->getNextUsedIter(iter) )
		(*iter).Model()->SetModelStateLocalStdDevs( stateLocalStdDevs, true );

	/// finally set the local std devto itself
	ARM_PricingModel::SetModelStateLocalStdDevs(stateLocalStdDevs, false );
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : SetModelStateLocalVars
///	Returns : void
///	Action  : Sets Model State Local Vars
///////////////////////////////////////////////////

void ARM_MultiAssetsModel::SetModelStateLocalVars( const ARM_MatrixVector& stateLocalVars, bool shareStateLocalVars )
{
	/// delegates to individual models
	for( ARM_ModelNameMap::iterator iter=itsModelMap->begin(); iter!=itsModelMap->end(); itsModelMap->getNextUsedIter(iter) )
		(*iter).Model()->SetModelStateLocalVars( stateLocalVars, true );

	/// finally set the local var to itself
	ARM_PricingModel::SetModelStateLocalVars(stateLocalVars, false );
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : SetModelStateLocalVars
///	Returns : void
///	Action  : Sets Model State Local Vars
///////////////////////////////////////////////////

void ARM_MultiAssetsModel::SetNumMethodStateLocalVars( const ARM_MatrixVector& stateLocalVars, bool isShared )
{
	/// delegates to individual models
	for( ARM_ModelNameMap::iterator iter=itsModelMap->begin(); iter!=itsModelMap->end(); itsModelMap->getNextUsedIter(iter) )
		(*iter).Model()->SetNumMethodStateLocalVars( stateLocalVars, !isShared);

	/// finally set the local var to itself
	ARM_PricingModel::SetNumMethodStateLocalVars(stateLocalVars, isShared);
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : SetModelStateLocalVars
///	Returns : void
///	Action  : Sets Model State Local Vars
///////////////////////////////////////////////////

void ARM_MultiAssetsModel::SetNumMethodStateLocalStdDevs( const ARM_MatrixVector& stateLocalVars, bool isShared )
{
	/// delegates to individual models
	for( ARM_ModelNameMap::iterator iter=itsModelMap->begin(); iter!=itsModelMap->end(); itsModelMap->getNextUsedIter(iter) )
		(*iter).Model()->SetNumMethodStateLocalStdDevs( stateLocalVars,!isShared);

	/// finally set the local var to itself
	ARM_PricingModel::SetNumMethodStateLocalStdDevs(stateLocalVars,isShared);
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : LocalDiscounts
///	Returns : ARM_VectorPtr
///	Action  : computes LocalDiscounts
////////////////////////////////////////////////////


ARM_VectorPtr ARM_MultiAssetsModel::LocalDiscounts( size_t timeIdx, double dt, 
	const ARM_PricingStatesPtr& states) const
{
	return itsRefModel->LocalDiscounts(timeIdx,dt,states);
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : DiscountFactor
///	Returns : ARM_VectorPtr
///	Action  : computes discount factor for a given model
////////////////////////////////////////////////////

ARM_VectorPtr ARM_MultiAssetsModel::DiscountFactor( 
	const string& modelName,
    double evalTime, 
	double maturityTime,
    const ARM_PricingStatesPtr& states) const
{
	return (*itsModelMap)[modelName]->Model()->DiscountFactor(	modelName, evalTime, maturityTime, states); 	
}



////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : Re_InitialiseCalibParams
///	Returns : void
///	Action  : 
////////////////////////////////////////////////////

void ARM_MultiAssetsModel::Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter)
{
	if( modelFitter.GetFactorNb() >= itsModelMap->size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": modelFitter.GetFactorNb() >= itsModelMap->size()!" );
	
	(*itsModelMap)[ modelFitter.GetFactorNb() ]->Model()->Re_InitialiseCalibParams( modelFitter );
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : PreProcessing
///	Returns : void
///	Action  : 
////////////////////////////////////////////////////

void ARM_MultiAssetsModel::PreProcessing(ARM_ModelFitter& modelFitter)
{
	if( modelFitter.GetFactorNb() >= itsModelMap->size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": modelFitter.GetFactorNb() >= itsModelMap->size()!" );
	
	(*itsModelMap)[ modelFitter.GetFactorNb() ]->Model()->PreProcessing( modelFitter );
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : PostProcessing
///	Returns : void
///	Action  : 
////////////////////////////////////////////////////

void ARM_MultiAssetsModel::PostProcessing(const ARM_ModelFitter& modelFitter)
{
	if( modelFitter.GetFactorNb() >= itsModelMap->size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": modelFitter.GetFactorNb() >= itsModelMap->size()!" );
	
	(*itsModelMap)[ modelFitter.GetFactorNb() ]->Model()->PostProcessing( modelFitter );
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : AdviseCurrentCalibSecIndex
///	Returns : void
///	Action  : 
////////////////////////////////////////////////////

void ARM_MultiAssetsModel::AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter)
{
	if( modelFitter.GetFactorNb() >= itsModelMap->size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": modelFitter.GetFactorNb() >= itsModelMap->size()!" );
	
	(*itsModelMap)[ modelFitter.GetFactorNb() ]->Model()->AdviseCurrentCalibSecIndex( index, modelFitter );
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : AdviseCurrentCalib
///	Returns : void
///	Action  : 
////////////////////////////////////////////////////

void ARM_MultiAssetsModel::AdviseCurrentCalib(ARM_ModelFitter& modelFitter)
{
	if( modelFitter.GetFactorNb() >= itsModelMap->size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": modelFitter.GetFactorNb() >= itsModelMap->size()!" );
	
	(*itsModelMap)[ modelFitter.GetFactorNb() ]->Model()->AdviseCurrentCalib( modelFitter );
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : AdviseBreakPointTimes
///	Returns : void
///	Action  : 
////////////////////////////////////////////////////
void ARM_MultiAssetsModel::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, 
							  ARM_ModelParam* modelParam, 
							  size_t factorNb )
{
	if( factorNb >= itsModelMap->size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": factorNb  >= itsModelMap->size()!" );
	
	(*itsModelMap)[ factorNb ]->Model()->AdviseBreakPointTimes( portfolio, modelParam,factorNb );
}




////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : ValidateCalibMethod
///	Returns : void
///	Action  : 
////////////////////////////////////////////////////
void ARM_MultiAssetsModel::ValidateCalibMethod(ARM_CalibMethod& calibMethod)
{
	if( calibMethod.GetFactorNb() >= itsModelMap->size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": calibMethod.GetFactorNb() >= itsModelMap->size()!" );
	
	(*itsModelMap)[ calibMethod.GetFactorNb() ]->Model()->ValidateCalibMethod( calibMethod );
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : GetAsOfDate
///	Returns : ARM_Date
///	Action  : returns the asOfDate of the first model
////////////////////////////////////////////////////
ARM_Date ARM_MultiAssetsModel::GetAsOfDate() const
{
	ARM_ModelNameMap::const_iterator iter	= itsModelMap->begin();
	return (*iter).Model()->GetAsOfDate();
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : FactorCount
///	Returns : size_t
///	Action  : returns the number of factors
////////////////////////////////////////////////////
size_t ARM_MultiAssetsModel::FactorCount() const
{
	ARM_ModelNameMap::const_iterator 
		iter, 
		start	= itsModelMap->begin(), 
		end		= itsModelMap->end();
	size_t nbFactors=0;

	for( iter=start; iter!= end; itsModelMap->getNextUsedIter(iter) )
	{
		nbFactors += (*iter).Model()->FactorCount();
	}
	return nbFactors;
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : ModelStatesCount
///	Returns : size_t
///	Action  : returns the number of factors
////////////////////////////////////////////////////
size_t ARM_MultiAssetsModel::ModelStatesSize() const
{
	ARM_ModelNameMap::const_iterator 
		iter, 
		start	= itsModelMap->begin(), 
		end		= itsModelMap->end();
	size_t nbModelStates=0;

	for( iter=start; iter!= end; itsModelMap->getNextUsedIter(iter) )
	{
		nbModelStates += (*iter).Model()->ModelStatesSize();
	}
	return nbModelStates;
}



////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: GetSettlementCalendar
///	Returns : string
///	Action  : get the settlement calendar
////////////////////////////////////////////////////
string ARM_MultiAssetsModel::GetSettlementCalendar(const string& modelName) const
{
	ARM_PricingModelPtr model = (*itsModelMap)[modelName]->Model();
	ARM_PricingFunctionEquity* equityFctor = dynamic_cast<ARM_PricingFunctionEquity*>(&*model);

	if( !equityFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + modelName + " does not support equity function CallVectorial!" );
	
	return equityFctor->GetSettlementCalendar( modelName );
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: GetSettlementGap
///	Returns : double
///	Action  : get the settlement gap
////////////////////////////////////////////////////
double ARM_MultiAssetsModel::GetSettlementGap(const string& modelName) const
{
	ARM_PricingModelPtr model = (*itsModelMap)[modelName]->Model();
	ARM_PricingFunctionEquity* equityFctor = dynamic_cast<ARM_PricingFunctionEquity*>(&*model);

	if( !equityFctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "the sub model " + modelName + " does not support equity function CallVectorial!" );
	
	return equityFctor->GetSettlementGap( modelName );
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: SetNumeraire
///	Returns : void
///	Action  : sets the numeraire to the pricing model
////////////////////////////////////////////////////

void ARM_MultiAssetsModel::SetNumeraire(const ARM_NumerairePtr& numerairePtr)
{
	ARM_PricingModel::SetNumeraire( numerairePtr);

	/// propagates the numeraire to the base models
	/// Initialises the various numeraires of the models
	ARM_ModelNameMap::iterator iter;
	for( iter=itsModelMap->begin(); iter!=itsModelMap->end(); ++iter )
		(*iter).Model()->SetNumeraire(numerairePtr);
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: SetNumMethod
///	Returns : void
///	Action  : sets the numerical method to the pricing model
////////////////////////////////////////////////////
void ARM_MultiAssetsModel::SetNumMethod(const ARM_NumMethodPtr& numMethodPtr )
{
	ARM_PricingModel::SetNumMethod( numMethodPtr );

	/// propagates the numeraire to the base models
	ARM_ModelNameMap::iterator iter;
	for( iter=itsModelMap->begin(); iter!=itsModelMap->end(); ++iter )
		(*iter).Model()->SetNumMethod(numMethodPtr);
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: GetRefModel
///	Returns : ARM_PricingModel*
///	Action  : returns the reference model
////////////////////////////////////////////////////

const ARM_PricingModel* ARM_MultiAssetsModel::GetRefModel() const
{
	return itsRefModel;
}

ARM_PricingModel* ARM_MultiAssetsModel::GetRefModel()
{
	return itsRefModel;
}



////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: ComputeDriftCommon
///	Returns : void
///	Action  : computes the relative and absolute drifts
////////////////////////////////////////////////////

void ARM_MultiAssetsModel::ComputeDriftCommon( 
	const std::vector<double>& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts,
	size_t rowSize,
	const ComputeDriftFunc& func ) const
{
	size_t factorNb = FactorCount();

	relativeDrifts	= ARM_GP_MatrixPtr( new ARM_GP_Matrix( rowSize, factorNb, 0.0 ) );
	absoluteDrifts	= ARM_GP_MatrixPtr( new ARM_GP_Matrix( rowSize, factorNb, 0.0 ) );

	/// delegates to individual models
	ARM_GP_MatrixPtr  tmpRelativeDrifts, tmpAbsoluteDrifts;

	size_t i,j,modelFactor, modelNb;
	for( ARM_ModelNameMap::iterator iter=itsModelMap->begin(); iter!=itsModelMap->end(); itsModelMap->getNextUsedIter(iter) )
	{
		(*((*iter).Model()).*func)( timeSteps, tmpRelativeDrifts, tmpAbsoluteDrifts );
		modelFactor = (*iter).Model()->FactorCount();
		modelNb		= (*iter).Model()->GetModelNb();

		if( tmpRelativeDrifts != ARM_GP_MatrixPtr(NULL) && tmpAbsoluteDrifts != ARM_GP_MatrixPtr(NULL) )
		{
			for( i=0; i<rowSize; ++i)
			{
				for( j=0; j<modelFactor; ++j )
				{
					(*relativeDrifts)(i,modelNb+j)=(*tmpRelativeDrifts)(i,j);
					(*absoluteDrifts)(i,modelNb+j)=(*tmpAbsoluteDrifts)(i,j);
				}
			}
		}
		else if( tmpRelativeDrifts != ARM_GP_MatrixPtr(NULL) )
		{
			for( i=0; i<rowSize; ++i)
			{
				for( j=0; j<modelFactor; ++j )
				{
					(*relativeDrifts)(i,modelNb+j)=(*tmpRelativeDrifts)(i,j);
				}
			}

		}
		else if( tmpAbsoluteDrifts != ARM_GP_MatrixPtr(NULL) )
		{
			for( i=0; i<rowSize; ++i)
			{
				for( j=0; j<modelFactor; ++j )
				{
					(*absoluteDrifts)(i,modelNb+j)=(*tmpAbsoluteDrifts)(i,j);
				}
			}
		}

		tmpRelativeDrifts=ARM_GP_MatrixPtr(NULL);
		tmpAbsoluteDrifts=ARM_GP_MatrixPtr(NULL);
	}
}




////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : IntegratedLocalDrifts
///	Returns : void
///	Action  : compute integrated relative and absolute drifts
////////////////////////////////////////////////////
void ARM_MultiAssetsModel::IntegratedLocalDrifts(
	const std::vector<double>& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts) const
{
	size_t nbSteps	= timeSteps.size()-1;
	ComputeDriftCommon( timeSteps, relativeDrifts, absoluteDrifts, nbSteps, &ARM_PricingModel::IntegratedLocalDrifts );
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: EulerLocalDrifts
///	Returns : void
///	Action  : computes the Euler relative and absolute drifts
////////////////////////////////////////////////////

void ARM_MultiAssetsModel::EulerLocalDrifts( const std::vector<double>& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const
{
	size_t nbSteps	= timeSteps.size();
	ComputeDriftCommon( timeSteps, relativeDrifts, absoluteDrifts, nbSteps, &ARM_PricingModel::EulerLocalDrifts );
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: MarkovianDrift
///	Returns : ARM_GP_MatrixPtr
///	Action  : computes the Markovian drift
////////////////////////////////////////////////////

ARM_GP_MatrixPtr ARM_MultiAssetsModel::MarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates ) const
{
	/// FIX FIX
	return ARM_GP_MatrixPtr( new ARM_GP_Matrix( numMethodStates->rows() , numMethodStates->cols(), 0.0 ) );
}


////////////////////////////////////////////////////
///	Class  : ARM_MultiAssetsModel
///	Routine: GetType
///	Returns: int
///	Action : A multi assets model links stochastic models
///			 but is not stcohastic itself !
////////////////////////////////////////////////////

int ARM_MultiAssetsModel::GetType() const
{
	return MT_MULTIASSET_MODEL | MT_NON_STOCHASTIC_MODEL;
}


////////////////////////////////////////////////////
///	Class  : ARM_MultiAssetsModel
///	Routine: SetModelMap, SetModelMapNoClone
///	Returns: void
///	Action : function to set a model map to a model
////////////////////////////////////////////////////

void ARM_MultiAssetsModel::SetModelMap(ARM_ModelNameMap* RefModelMap)
{
	delete itsModelMap;  
	itsModelMap = RefModelMap? static_cast<ARM_ModelNameMap*>( RefModelMap->Clone() ): NULL;
	SetModelParamsVec();
	SetMultiFactorFlagOnModel();
}

void ARM_MultiAssetsModel::SetModelMapNoClone(ARM_ModelNameMap* RefModelMap)
{
	delete itsModelMap;  
	itsModelMap = RefModelMap;
	SetModelParamsVec();
	SetMultiFactorFlagOnModel();
}

////////////////////////////////////////////////////
///	Class  : ARM_MultiAssetsModel
///	Routine: SetMultiFactorFlagOnModel
///	Returns: void
///	Action : tells the basic model that they are part of multi-factor ones
////////////////////////////////////////////////////
void ARM_MultiAssetsModel::SetMultiFactorFlagOnModel()
{
	if( itsModelMap )
	{
		ARM_ModelNameMap::iterator 
			iter	= itsModelMap->begin(), 
			end		= itsModelMap->end();

		for( ; iter != end; ++iter )
			(*iter).Model()->SetFromMultiFactor(true);
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_MultiAssetsModel
///	Routine: SetModelParamsVec
///	Returns: void
///	Action : function to set the corresponding model params vec according to the mode name map
////////////////////////////////////////////////////
void ARM_MultiAssetsModel::SetModelParamsVec()
{
	if( itsModelMap )
	{
		vector<ARM_ModelParams*> paramsVec( itsModelMap->size() );
		for( size_t i=0; i<itsModelMap->size(); ++i )
			paramsVec[i] = (*itsModelMap)[i]->Model()->GetModelParams();
		ARM_ModelParamsVec modelParamsVec( paramsVec );
		SetModelParams( modelParamsVec );
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_MultiAssetsModel
///	Routine: GetModelNamePerFactor
///	Returns: string
///	Action : to get model name per factor
////////////////////////////////////////////////////
string ARM_MultiAssetsModel::GetModelNamePerFactor(size_t factorNb ) const
{
	if( itsModelMap )
	{
		if( factorNb >= itsModelMap->size() )
			ARM_THROW( ERR_INVALID_ARGUMENT, ": out or range for factor nb!" );
		return (*itsModelMap)[factorNb]->Model()->GetModelName();
	}
	else
		return GetModelName();
}


////////////////////////////////////////////////////
///	Class  : ARM_MultiAssetsModel
///	Routine: GetCorrelSubMatrix
///	Returns: void
///	Action : returns the submatrix of correlation between
///  model named modelName1 and modelName2
////////////////////////////////////////////////////
const ARM_GP_MatrixPtr ARM_MultiAssetsModel::GetCorrelSubMatrix( const string& ModelName1, const string& ModelName2 ) const
{
	if( itsModelMap )
	{
		size_t factorNb1, factorNb2, totalFactorNb;
		size_t model1Idx = 0, model2Idx = 0;
		ARM_ModelNameMap::const_iterator 
		iter	= itsModelMap->begin(), 
		end		= itsModelMap->end();

		/// Get Index and FactorNb of the first Model
		for( ; iter!= end; ++iter )
		{
				factorNb1 = (*iter).Model()->FactorCount();
				if ((*iter).ModelName() == ModelName1 )
					break;
				else
					model1Idx += factorNb1;
		}

		/// Get Index and FactorNb of the second Model
		iter	= itsModelMap->begin();

		for( ; iter!= end; ++iter )
		{
				factorNb2 = (*iter).Model()->FactorCount();
				if ((*iter).ModelName() == ModelName2 )
					break;
				else
					model2Idx += factorNb2;
		}

		if( model1Idx == FactorCount() || model2Idx == FactorCount() )
			ARM_THROW( ERR_INVALID_ARGUMENT, "Model Not found in modelName Map" );

		if(!itsCorrelMatrix)
			ARM_THROW( ERR_INVALID_ARGUMENT, " Invalid correlation matrix" );

		ARM_GP_Matrix correlMatrix = itsCorrelMatrix->Interpolate(0);

		/// Build and fill the correlSubMatrix
		totalFactorNb = factorNb1 + factorNb2;
		ARM_GP_Matrix * returnedMatrix = new ARM_GP_Matrix( totalFactorNb, totalFactorNb );

		size_t i,j;
		i=0;j=0;

		/// Diagonal terms of the matrix (correl between model1's own factors and model2's own factors)
		for(i=0;i<factorNb1;i++)
			for(j=0;j<factorNb1;j++)
				returnedMatrix->Elt(i,j) = correlMatrix(i+model1Idx,j+model1Idx);

		for(i=0;i<factorNb2;i++)
			for(j=0;j<factorNb2;j++)
				returnedMatrix->Elt(i+factorNb1,j+factorNb1) = correlMatrix(i+model2Idx,j+model2Idx);

		/// Correl between model1 and model2
		for(i=0;i<factorNb2;i++)
			for(j=0;j<factorNb1;j++)
			{
				returnedMatrix->Elt(i+factorNb1,j) = correlMatrix(i+model2Idx,j+model1Idx);
				returnedMatrix->Elt(i,j+factorNb1) = correlMatrix(i+model2Idx,j+model1Idx);
			}

		return returnedMatrix;

	}

	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "Multi needs a modelMap!" );
}

////////////////////////////////////////////////////
///	Class  : ARM_MultiAssetsModel
///	Routine: VolatilitiesAndCorrelationTimesSteps
///	Returns: void
///	Action : Volatilities and correlation time steps for PDEs
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_MultiAssetsModel::VolatilitiesAndCorrelationTimesSteps() const
{
	std::vector<double> *result = new std::vector<double>(0);
	std::vector<double> *intermed = NULL;

	ARM_ModelNameMap::const_iterator 
	iter	= itsModelMap->begin(), 
	end		= itsModelMap->end();

	for( ; iter!=end ; itsModelMap->getNextUsedIter(iter))
	{
		intermed = MergeSortedVectorNoDuplicates( *((*iter).Model()->VolatilitiesAndCorrelationTimesSteps()) , *result );
		delete result;
		result = intermed;
	}

	return ARM_GP_VectorPtr(result);
}

////////////////////////////////////////////////////
///	Class  : ARM_MultiAssetsModel
///	Routine: ImpliedVol
///	Returns: int
///	Action : tells the type of the model
////////////////////////////////////////////////////
double ARM_MultiAssetsModel::ImpliedVol(const ARM_VanillaArg& arg) const
{
    return DefaultImpliedVol(arg);
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

